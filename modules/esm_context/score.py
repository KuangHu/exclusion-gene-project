"""Score extracted windows with two separate frozen heads and late fusion.

Computational decision flags are exploratory model calls, never source annotations.
This script does not create a novelty shortlist or infer conjugation anchor classes.
"""
import argparse,csv,json,hashlib,sqlite3
from contextlib import closing
from pathlib import Path
import numpy as np
from extract import file_sha
from neighborhood_io import output_lock

class Scorer:
    def __init__(self,folder):
        folder=Path(folder);self.card=json.loads((folder/'model_card.json').read_text());self.heads={}
        for row in self.card['heads']:
            path=folder/row['path']
            if file_sha(path)!=row['sha256']:raise ValueError('Head checksum mismatch')
            with np.load(path,allow_pickle=False) as z:
                self.heads[row['name']]={k:z[k].copy() for k in z.files}
        for name,positions in [('center',[0]),('masked_pm3',[-3,-2,-1,1,2,3])]:
            m=self.heads[name]
            if m['positions'].tolist()!=positions or int(m['dimensions'])!=2560*len(positions):raise ValueError('Unexpected model positions/dimensions')
        path=folder/self.card['training_hash_file']
        if file_sha(path)!=self.card['training_hash_file_sha256']:raise ValueError('Training hash manifest mismatch')
        self.training=set(path.read_text().splitlines())
    def head(self,name,vectors):
        m=self.heads[name];a=np.asarray(vectors,dtype='float32');norm=np.linalg.norm(a,axis=1,keepdims=True)
        if not np.isfinite(a).all() or (norm<=0).any():raise ValueError('Invalid input vector')
        x=(a/norm).ravel().astype(float)
        if len(x)!=int(m['dimensions']):raise ValueError('Head dimension mismatch')
        z=(x-m['mean'])/m['scale']/np.sqrt(len(x));return float(z@m['coef'][0]+m['intercept'][0])
    def fusion(self,c,n):
        f=self.card['fusion'];z=(np.array([c,n])-f['mean'])/f['scale'];return float(f['alpha']*z[0]+(1-f['alpha'])*z[1])

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--extraction',type=Path,required=True);p.add_argument('--embeddings',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--model',type=Path,default=Path(__file__).parent/'model');a=p.parse_args(argv)
    em=json.loads((a.extraction/'extraction_manifest.json').read_text())
    if em['status']!='complete':raise ValueError('Extraction is incomplete')
    validation=json.loads((a.embeddings/'embedding_validation.json').read_text());progress=json.loads((a.embeddings/'progress.json').read_text());scorer=Scorer(a.model)
    if validation['status']!='passed' or progress['status']!='complete' or progress['signature']!=scorer.card['signature']:raise ValueError('Embedding validation/protocol mismatch')
    reference=json.loads((a.embeddings/'reference_validation.json').read_text())
    if reference.get('passed') is not True or reference.get('reference_sha256')!=scorer.card['embedding_reference_sha256']:raise ValueError('Classifier numerical embedding reference was not validated')
    if file_sha(a.extraction/'all_proteins.faa')!=validation['full_input_sha256']:raise ValueError('Embedding input does not match extraction')
    if file_sha(a.embeddings/'embeddings.npy')!=validation['embeddings_file_sha256']:raise ValueError('Embedding matrix checksum mismatch')
    x=np.load(a.embeddings/'embeddings.npy',mmap_mode='r');index={}
    with (a.embeddings/'sequence_manifest.tsv').open(encoding='utf8') as f:
        for i,r in enumerate(csv.DictReader(f,delimiter='\t')):
            if int(r['row_0based'])!=i or r['sequence_sha256'] in index:raise ValueError('Embedding manifest order/duplicates invalid')
            index[r['sequence_sha256']]=i
    if x.shape!=(len(index),2560):raise ValueError('Embedding matrix shape mismatch')
    a.out.mkdir(parents=True,exist_ok=True)
    if (a.out/'window_scores.tsv').exists():raise ValueError('Use a new scoring output directory')
    with output_lock(a.out), closing(sqlite3.connect(a.out/'score_aggregation.sqlite')) as agg:
        fields=['case_id','dataset','record_id','center_cds','center_sha256','center_aa_length','center_gene','center_product','source_occurrence_count','window_sequence_id','center_score','center_margin','center_model_positive','context_score','context_margin','context_model_positive','fusion_score','fusion_margin','fusion_model_positive','context_status','exact_training_center','functional_label','novelty_shortlist_status']
        agg.execute('CREATE TABLE IF NOT EXISTS variants(hash TEXT,window_id TEXT,c REAL,n REAL,f REAL,PRIMARY KEY(hash,window_id))');n=0;both=0
        with (a.extraction/'windows.jsonl').open(encoding='utf8') as src,(a.out/'window_scores.tsv').open('w',encoding='utf8',newline='') as dst:
            writer=csv.DictWriter(dst,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
            for line in src:
                w=json.loads(line);c=w['center'];slots={g['relative_cds']:g for g in w['slots']};h=c['center_sha256']
                if h not in index:raise ValueError('Center vector missing from completed input')
                cs=scorer.head('center',[x[index[h]]]);cm=cs-float(scorer.heads['center']['threshold']);ns=fs=None;status='eligible' if c['context_eligible'] else 'incomplete_or_repeated_center_pm3'
                if c['context_eligible']:
                    hs=[slots[i]['sequence_sha256'] for i in (-3,-2,-1,1,2,3)]
                    if any(v not in index for v in hs):raise ValueError('Required context vector missing')
                    ns=scorer.head('masked_pm3',[x[index[v]] for v in hs]);fs=scorer.fusion(cs,ns);both+=1
                row={k:c.get(k,'') for k in fields};row.update(center_score=cs,center_margin=cm,center_model_positive=cm>=0,context_score='' if ns is None else ns,context_margin='' if ns is None else ns-float(scorer.heads['masked_pm3']['threshold']),context_model_positive='' if ns is None else ns>=float(scorer.heads['masked_pm3']['threshold']),fusion_score='' if fs is None else fs,fusion_margin='' if fs is None else fs-scorer.card['fusion']['threshold'],fusion_model_positive='' if fs is None else fs>=scorer.card['fusion']['threshold'],context_status=status,exact_training_center=h in scorer.training,functional_label='unassigned',novelty_shortlist_status='not_evaluated_requires_existing_class_matched_anchor_pm7')
                writer.writerow(row);agg.execute('INSERT OR IGNORE INTO variants VALUES (?,?,?,?,?)',(h,c['window_sequence_id'],cs,ns,fs));n+=1
        agg.commit()
        with (a.out/'center_scores.tsv').open('w',encoding='utf8',newline='') as dst:
            writer=csv.writer(dst,delimiter='\t',lineterminator='\n');writer.writerow(['sequence_sha256','distinct_windows','scorable_context_windows','mean_center_score','mean_context_score','mean_fusion_score','center_model_positive','context_model_positive','fusion_model_positive','exact_training_center','functional_label'])
            for h,nn,nc,c,ns,fs in agg.execute('SELECT hash,COUNT(*),COUNT(n),AVG(c),AVG(n),AVG(f) FROM variants GROUP BY hash ORDER BY hash'):
                writer.writerow([h,nn,nc,c,'' if ns is None else ns,'' if fs is None else fs,c>=float(scorer.heads['center']['threshold']),'' if ns is None else ns>=float(scorer.heads['masked_pm3']['threshold']),'' if fs is None else fs>=scorer.card['fusion']['threshold'],h in scorer.training,'unassigned'])
        unique=agg.execute('SELECT COUNT(DISTINCT hash) FROM variants').fetchone()[0]
        summary=dict(status='complete',model_id=scorer.card['model_id'],model_card_sha256=file_sha(a.model/'model_card.json'),case_windows=n,context_scored_cases=both,unique_center_sequences=unique,annotation_status='No functional annotations assigned',aggregation='Equal weight per distinct full +/-7 sequence window; context/fusion means use only context-scorable windows.',limitations=scorer.card['limitations']);(a.out/'scoring_manifest.json').write_text(json.dumps(summary,indent=2),encoding='utf8');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
