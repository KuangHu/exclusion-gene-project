"""Merge a complete extraction shard set, deduplicating required ESM sequences globally."""
import argparse,json
from pathlib import Path
from extract import connect,export
from neighborhood_io import output_lock

def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--shards-root',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args(argv)
    files=sorted(a.shards_root.glob('shard_*/extraction_manifest.json'))
    if not files:raise ValueError('No shard_*/extraction_manifest.json found')
    manifests=[json.loads(p.read_text()) for p in files]
    if any(m['status']!='complete' for m in manifests):raise ValueError('Incomplete extraction shard')
    canonical={k:v for k,v in manifests[0]['config'].items() if k!='shard'}
    if any({k:v for k,v in m['config'].items() if k!='shard'}!=canonical for m in manifests):raise ValueError('Shard configurations differ')
    ids=[m['config']['shard'] for m in manifests]
    if sorted(ids)!=list(range(canonical['num_shards'])):raise ValueError('Missing or duplicate shard indices')
    a.out.mkdir(parents=True,exist_ok=True)
    if (a.out/'extraction.sqlite').exists():raise ValueError('Use a new merge output directory')
    with output_lock(a.out):
        db=connect(a.out/'extraction.sqlite')
        try:
            for p in files:
                # ATTACH is parameterized; source is opened read-only.
                uri=(p.parent/'extraction.sqlite').resolve().as_uri()+'?mode=ro'
                db.execute('ATTACH DATABASE ? AS shard',(uri,))
                with db:
                    for table in ('records','occurrences','cases'):db.execute(f'INSERT INTO {table} SELECT * FROM shard.{table}')
                    conflict=db.execute('SELECT s.hash FROM sequences s JOIN shard.sequences t ON t.hash=s.hash WHERE s.sequence != t.sequence LIMIT 1').fetchone()
                    if conflict:raise ValueError('Sequence hash/content conflict')
                    for table in ('sequences','needed'):db.execute(f'INSERT OR IGNORE INTO {table} SELECT * FROM shard.{table}')
                db.execute('DETACH DATABASE shard')
            counts=export(db,a.out,canonical['embedding_radius']);db.execute('PRAGMA wal_checkpoint(TRUNCATE)')
        finally:db.close()
        report=dict(status='complete',config=canonical,merged_shards=ids,counts=counts,annotation_status='unassigned; no EEX calls',source_manifests=[str(p.resolve()) for p in files],empty_shard=counts['centers']==0)
        (a.out/'extraction_manifest.json').write_text(json.dumps(report,indent=2),encoding='utf8');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
