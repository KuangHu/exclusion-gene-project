"""Extract short-protein centers and ordered neighborhoods from PLSDB/IMGPR.

No EEX calls, functional labels, anchor selection, or downloads are performed.
All CDS remain in the ordering, including invalid and partial translations.
"""
import argparse
from contextlib import contextmanager, closing
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import sqlite3
import sys

from Bio import SeqIO, __version__ as BIO_VERSION

AA = set('ACDEFGHIKLMNPQRSTVWY')
MODULE = Path(__file__).resolve().parent
REPO = MODULE.parents[1]

def digest(text):
    return hashlib.sha256(text.encode()).hexdigest()

def file_sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

@contextmanager
def open_text(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf8') as f:
        yield f

def protein(sequence):
    s = ''.join(str(sequence).split()).upper()
    if s.endswith('*'):
        s = s[:-1]
    if not s:
        return s, '', 'missing_translation'
    if not set(s) <= AA:
        return s, '', 'nonstandard_translation'
    return s, digest(s), 'valid'

def gene(start0, end0, strand, sequence, source_id, partial='00', product='unknown', name=''):
    if start0 < 0 or end0 <= start0 or strand not in (-1, 1):
        raise ValueError('Invalid CDS coordinates or strand: ' + str(source_id))
    s, h, status = protein(sequence)
    if partial not in ('00', '', 'unknown'):
        status = 'partial_CDS'
    return dict(start0=start0, end0=end0, strand=strand, sequence=s,
                sequence_sha256=h if status == 'valid' else '', sequence_status=status,
                aa_length=len(s), source_gene_id=source_id, partial=partial,
                gene=name, product=product or 'unknown')

def parse_protein_header(record, fmt):
    if fmt == 'project-cds-faa':
        # scripts/109_full_cds_cache.py: accession|ordinal|start1|end1|strand|aa_length
        fields = record.id.rsplit('|', 5)
        if len(fields) != 6:
            raise ValueError('Expected complete project CDS cache header: ' + record.id)
        accession, ordinal, start, end, strand, length = fields
        g = gene(int(start)-1, int(end), int(strand), record.seq, record.id, 'unknown')
        if len(protein(record.seq)[0]) != int(length):
            raise ValueError('Header/sequence length disagreement: ' + record.id)
    elif fmt == 'prodigal-faa':
        fields = record.description.split(' # ')
        if len(fields) != 5 or '_' not in record.id:
            raise ValueError('Expected Prodigal coordinate header: ' + record.id)
        accession, ordinal = record.id.rsplit('_', 1)
        match = re.search(r'(?:^|;)partial=([01]{2})(?:;|$)', fields[4])
        if not match:
            raise ValueError('Prodigal header lacks partial flag: ' + record.id)
        g = gene(int(fields[1])-1, int(fields[2]), int(fields[3]), record.seq,
                 record.id, match.group(1))
    else:
        raise ValueError(fmt)
    g['source_ordinal'] = int(ordinal)
    return accession, g

def ordered_genes(genes, check_ordinals=False):
    genes = sorted(genes, key=lambda g: (g['start0'], g['end0'], g.get('source_ordinal', 0)))
    if check_ordinals:
        ords = sorted(g['source_ordinal'] for g in genes)
        if not ords or ords[0] not in (0, 1) or ords != list(range(ords[0], ords[0]+len(ords))):
            raise ValueError('Non-contiguous or duplicated CDS ordinals: input may be a filtered protein set')
    for i, g in enumerate(genes):
        g['cds_index'] = i
    return genes

def protein_records(path, fmt):
    seen = set()
    current, genes = None, []
    with open_text(path) as f:
        for record in SeqIO.parse(f, 'fasta'):
            accession, g = parse_protein_header(record, fmt)
            if current is not None and accession != current:
                seen.add(current)
                ordered = ordered_genes(genes, True)
                yield current, digest(json.dumps(ordered, sort_keys=True)), ordered
                genes = []
            if accession in seen:
                raise ValueError('Contig reappears non-contiguously in protein FASTA: ' + accession)
            current = accession
            genes.append(g)
    if current is not None:
        ordered = ordered_genes(genes, True)
        yield current, digest(json.dumps(ordered, sort_keys=True)), ordered

def records(path, fmt):
    if fmt in ('project-cds-faa', 'prodigal-faa'):
        yield from protein_records(path, fmt)
        return
    with open_text(path) as f:
        for rec in SeqIO.parse(f, 'fasta' if fmt == 'genome-fasta' else 'genbank'):
            payload = str(rec.seq).upper() if fmt == 'genome-fasta' else rec.format('genbank')
            yield rec.id, digest(payload), rec

def decode_genes(payload, fmt):
    if fmt in ('project-cds-faa', 'prodigal-faa'):
        return payload
    if fmt == 'genome-fasta':
        import orf_caller
        calls = orf_caller.call(str(payload.seq))
        return ordered_genes([gene(c['start']-1, c['end'], c['strand'], c['aa'],
                                  payload.id + '|' + str(i), '00') for i, c in enumerate(calls)])
    genes = []
    for i, feat in enumerate(payload.features):
        if feat.type != 'CDS':
            continue
        if feat.location is None:
            raise ValueError('CDS without location; cannot safely preserve ordering')
        q = feat.qualifiers
        strand = feat.location.strand
        if strand not in (-1, 1):
            raise ValueError('CDS without unambiguous strand')
        partial = '00' if type(feat.location.start).__name__ == type(feat.location.end).__name__ == 'ExactPosition' else 'partial'
        g = gene(int(feat.location.start), int(feat.location.end), strand,
                 q.get('translation', [''])[0], q.get('protein_id', [payload.id+'|'+str(i)])[0],
                 partial, q.get('product', ['unknown'])[0], q.get('gene', [''])[0])
        if len(feat.location.parts) != 1:
            g['sequence_status'] = 'compound_CDS'
            g['sequence_sha256'] = ''
        genes.append(g)
    return ordered_genes(genes)

def windows_for_record(dataset, accession, key, genes, min_aa=40, max_aa=300):
    for i, center in enumerate(genes):
        if center['sequence_status'] != 'valid' or not min_aa <= center['aa_length'] <= max_aa:
            continue
        slots = []
        for rel in range(-7, 8):
            j = i + rel * center['strand']
            if 0 <= j < len(genes):
                g = dict(genes[j])
                g.pop('sequence')
            else:
                g = dict(start0='', end0='', strand=0, sequence_sha256='',
                         sequence_status='outside_record', aa_length=0, source_gene_id='',
                         partial='', gene='', product='unknown', cds_index=j)
            g['relative_cds'] = rel
            g['relative_strand'] = g['strand'] * center['strand']
            slots.append(g)
        needed = [g for g in slots if abs(g['relative_cds']) <= 3]
        valid = all(g['sequence_status'] == 'valid' for g in needed)
        duplicate = any(g['relative_cds'] and g['sequence_sha256'] == center['sequence_sha256'] for g in needed)
        case_id = digest(key + '\0' + str(i))
        row = dict(case_id=case_id, dataset=dataset, record_id=accession, record_key=key,
                   center_cds=i, center_strand=center['strand'], center_sha256=center['sequence_sha256'],
                   center_aa_length=center['aa_length'], center_start0=center['start0'],
                   center_end0=center['end0'], center_gene=center['gene'], center_product=center['product'],
                   full_pm3=valid, full_pm7=all(g['sequence_status']=='valid' for g in slots),
                   center_repeated_pm3=duplicate, context_eligible=valid and not duplicate,
                   functional_label='unassigned', source_annotation_mode='source_preserved',
                   window_sequence_id=digest('|'.join(g['sequence_sha256'] for g in slots)))
        yield row, slots

def connect(path):
    db = sqlite3.connect(path, uri=True)
    db.executescript('''PRAGMA journal_mode=WAL;
    CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY, value TEXT);
    CREATE TABLE IF NOT EXISTS records(key TEXT PRIMARY KEY, dataset TEXT, accession TEXT, digest TEXT, n_cds INTEGER);
    CREATE TABLE IF NOT EXISTS occurrences(file_index INTEGER, ordinal INTEGER, record_key TEXT, PRIMARY KEY(file_index,ordinal));
    CREATE TABLE IF NOT EXISTS cases(id TEXT PRIMARY KEY, record_key TEXT, payload TEXT, slots TEXT);
    CREATE TABLE IF NOT EXISTS sequences(hash TEXT PRIMARY KEY, sequence TEXT NOT NULL);
    CREATE TABLE IF NOT EXISTS needed(hash TEXT PRIMARY KEY);
    ''')
    return db

def export(db, out, embedding_radius):
    case_cursor = db.execute('SELECT payload,slots,record_key FROM cases ORDER BY id')
    case_fields = None
    slot_fields = ['case_id','dataset','record_id','relative_cds','cds_index','sequence_sha256',
                   'sequence_status','aa_length','start0','end0','strand','relative_strand','partial',
                   'source_gene_id','gene','product']
    copies = dict(db.execute('SELECT record_key,COUNT(*) FROM occurrences GROUP BY record_key'))
    counts = dict(records=db.execute('SELECT COUNT(*) FROM records').fetchone()[0], centers=0, context_eligible=0)
    with (out/'centers.tsv').open('w',encoding='utf8',newline='') as cf, (out/'slots.tsv').open('w',encoding='utf8',newline='') as sf, (out/'embedding_slots.tsv').open('w',encoding='utf8',newline='') as ef, (out/'windows.jsonl').open('w',encoding='utf8') as jf:
        sw=csv.DictWriter(sf,fieldnames=slot_fields,delimiter='\t',lineterminator='\n');sw.writeheader()
        ew=csv.DictWriter(ef,fieldnames=slot_fields,delimiter='\t',lineterminator='\n');ew.writeheader()
        for data, raw_slots, key in case_cursor:
            c=json.loads(data);slots=json.loads(raw_slots);c['source_occurrence_count']=copies[key]
            if case_fields is None:
                case_fields=list(c);cw=csv.DictWriter(cf,fieldnames=case_fields,delimiter='\t',lineterminator='\n');cw.writeheader()
            cw.writerow(c);counts['centers']+=1;counts['context_eligible']+=c['context_eligible']
            jf.write(json.dumps(dict(center=c,slots=slots),ensure_ascii=False)+'\n')
            for g in slots:
                row={k:g.get(k,'') for k in slot_fields};row.update(case_id=c['case_id'],dataset=c['dataset'],record_id=c['record_id']);sw.writerow(row)
                if abs(g['relative_cds']) <= embedding_radius:ew.writerow(row)
    counts['unique_embedding_sequences']=db.execute('SELECT COUNT(*) FROM needed').fetchone()[0]
    with (out/'all_proteins.faa').open('w',encoding='utf8') as f:
        for h,s in db.execute('SELECT s.hash,s.sequence FROM sequences s JOIN needed n ON n.hash=s.hash ORDER BY s.hash'):
            if digest(s)!=h:raise ValueError('Stored sequence hash mismatch')
            f.write('>SHA256_'+h+'\n'+s+'\n')
    counts['distinct_window_sequences']=db.execute("SELECT COUNT(DISTINCT json_extract(payload,'$.window_sequence_id')) FROM cases").fetchone()[0]
    return counts

def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--inputs',type=Path,required=True,help='TSV: dataset, path, format; paths relative to this TSV')
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--min-aa',type=int,default=40);ap.add_argument('--max-aa',type=int,default=300)
    ap.add_argument('--embedding-radius',type=int,choices=(3,7),default=3)
    ap.add_argument('--shard-index',type=int,default=0);ap.add_argument('--num-shards',type=int,default=1)
    args=ap.parse_args(argv)
    if not 0 <= args.shard_index < args.num_shards or not 1 <= args.min_aa <= args.max_aa:ap.error('Invalid length range or shard')
    with args.inputs.open(encoding='utf8') as f:inputs=list(csv.DictReader(f,delimiter='\t'))
    if not inputs:raise ValueError('Empty input manifest')
    provenance=[]
    for r in inputs:
        if r.get('format') not in ('genome-fasta','genbank','project-cds-faa','prodigal-faa'):raise ValueError('Unrecognized explicit input format')
        if not r.get('dataset','').strip():raise ValueError('Dataset namespace required')
        path=Path(r['path']);path=(args.inputs.parent/path).resolve() if not path.is_absolute() else path.resolve()
        provenance.append(dict(dataset=r['dataset'],format=r['format'],path=str(path),sha256=file_sha(path)))
    caller_note={}
    if any(r['format']=='genome-fasta' for r in inputs):
        sys.path.insert(0,str(REPO/'scripts/lib'))
        import adapters,orf_caller,pyrodigal
        adapters.verify_all(verbose=True)
        caller_note=dict(meta=orf_caller.ANON,min_gene_nt=orf_caller.MIN_GENE,closed=orf_caller.CLOSED,translation_table=orf_caller.TRANSL_TABLE,pyrodigal_version=pyrodigal.__version__,source_sha256=file_sha(REPO/'scripts/lib/orf_caller.py'))
    config=dict(schema=1,extractor_sha256=file_sha(Path(__file__)),biopython_version=BIO_VERSION,inputs=provenance,min_aa=args.min_aa,max_aa=args.max_aa,embedding_radius=args.embedding_radius,shard=args.shard_index,num_shards=args.num_shards,caller=caller_note,wrap_around=False)
    args.out.mkdir(parents=True,exist_ok=True)
    # Reuse the same cross-platform OS-owned lock used by embedding inference.
    from neighborhood_io import output_lock
    with output_lock(args.out), closing(connect(args.out/'extraction.sqlite')) as db:
        encoded=json.dumps(config,sort_keys=True);old=db.execute("SELECT value FROM meta WHERE key='config'").fetchone()
        if old and old[0]!=encoded:raise ValueError('Inputs/configuration changed; use a new output directory')
        db.execute("INSERT OR IGNORE INTO meta VALUES ('config',?)",(encoded,));db.commit()
        for file_index,r in enumerate(provenance):
            for ordinal,(accession,raw_digest,payload) in enumerate(records(r['path'],r['format'])):
                key=digest(r['dataset']+'\0'+accession)
                if int(key,16)%args.num_shards!=args.shard_index:continue
                old=db.execute('SELECT digest FROM records WHERE key=?',(key,)).fetchone()
                if old:
                    if old[0]!=raw_digest:raise ValueError('Conflicting content under dataset/record ID: '+accession)
                    with db:db.execute('INSERT OR IGNORE INTO occurrences VALUES (?,?,?)',(file_index,ordinal,key))
                    continue
                genes=decode_genes(payload,r['format'])
                with db:
                    db.execute('INSERT INTO records VALUES (?,?,?,?,?)',(key,r['dataset'],accession,raw_digest,len(genes)))
                    db.execute('INSERT INTO occurrences VALUES (?,?,?)',(file_index,ordinal,key))
                    for c,slots in windows_for_record(r['dataset'],accession,key,genes,args.min_aa,args.max_aa):
                        db.execute('INSERT INTO cases VALUES (?,?,?,?)',(c['case_id'],key,json.dumps(c),json.dumps(slots)))
                        for slot in slots:
                            h=slot['sequence_sha256']
                            if not h:continue
                            seq=genes[slot['cds_index']]['sequence'];db.execute('INSERT OR IGNORE INTO sequences VALUES (?,?)',(h,seq))
                            if abs(slot['relative_cds'])<=args.embedding_radius:db.execute('INSERT OR IGNORE INTO needed VALUES (?)',(h,))
                if ordinal%100==0:print(json.dumps(dict(dataset=r['dataset'],record=accession,ordinal=ordinal,shard=args.shard_index)),flush=True)
        counts=export(db,args.out,args.embedding_radius);db.execute('PRAGMA wal_checkpoint(TRUNCATE)')
        if counts['records']==0 and args.num_shards==1:
            raise ValueError('No source records parsed; check input format before analysis')
        manifest=dict(status='complete',config=config,counts=counts,copy_definition='One distinct dataset/record ID and center CDS; repeated source records counted separately, never treated as independent biological samples.',ordering='All CDS sorted by coordinate, including invalid/partial translations; offsets oriented to center strand.',annotation_status='Unassigned; source gene/product preserved. No EEX or novelty classification.',empty_shard=counts['centers']==0)
        (args.out/'extraction_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf8');print(json.dumps(manifest,indent=2))

if __name__=='__main__':main()
