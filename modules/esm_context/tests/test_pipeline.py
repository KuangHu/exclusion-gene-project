import contextlib
import csv
import io
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest
from types import SimpleNamespace
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import extract
import merge
import embed
from score import Scorer
import score

class PipelineTests(unittest.TestCase):
    def genes(self,n=17):
        aa='ACDEFGHIKLMNPQRSTVWY'
        return extract.ordered_genes([extract.gene(i*200,i*200+153,-1 if i==8 else 1,'M'+aa[i%20]*49,str(i)) for i in range(n)])
    def manifest(self,root,nrecords=2,invalid=False):
        path=root/'proteins.faa'
        with path.open('w') as f:
            for rec in range(nrecords):
                for i,g in enumerate(self.genes()):
                    s=g['sequence'] if not invalid or i!=6 else 'M'+'X'*49
                    f.write(f'>ACC{rec}|{i}|{g["start0"]+1}|{g["end0"]}|{g["strand"]}|50\n{s}\n')
        manifest=root/'inputs.tsv';manifest.write_text('dataset\tpath\tformat\nPLSDB\tproteins.faa\tproject-cds-faa\n')
        return manifest
    def run_extract(self,manifest,out,*args):
        with contextlib.redirect_stdout(io.StringIO()):extract.main(['--inputs',str(manifest),'--out',str(out),*args])
    def test_reverse_orientation_is_not_reindexed(self):
        gs=self.genes();r,s=next((r,s) for r,s in extract.windows_for_record('D','A','key',gs) if r['center_cds']==8)
        self.assertEqual([g['cds_index'] for g in s],list(range(15,0,-1)))
        self.assertTrue(r['context_eligible']);self.assertEqual(s[7]['sequence_sha256'],gs[8]['sequence_sha256'])
    def test_invalid_translation_keeps_its_slot(self):
        gs=self.genes();gs[6]=extract.gene(1200,1353,1,'MXX','bad');gs[6]['cds_index']=6
        r,s=next((r,s) for r,s in extract.windows_for_record('D','A','key',gs) if r['center_cds']==8)
        self.assertEqual(s[9]['sequence_status'],'nonstandard_translation');self.assertEqual(s[10]['cds_index'],5);self.assertFalse(r['context_eligible'])
    def test_no_wrap_and_no_center_copy_in_context(self):
        gs=self.genes();r,s=next(extract.windows_for_record('D','A','k',gs));self.assertEqual(s[0]['sequence_status'],'outside_record');self.assertFalse(r['context_eligible'])
        gs[7].update(sequence=gs[8]['sequence'],sequence_sha256=gs[8]['sequence_sha256'])
        r,_=next((r,s) for r,s in extract.windows_for_record('D','A','k',gs) if r['center_cds']==8);self.assertTrue(r['center_repeated_pm3']);self.assertFalse(r['context_eligible'])
    def test_prodigal_header_partial_and_coordinate(self):
        r=SeqRecord(Seq('MAA*'),id='CP_001.1_1',description='CP_001.1_1 # 4 # 15 # -1 # ID=1_1;partial=10;start_type=ATG;')
        acc,g=extract.parse_protein_header(r,'prodigal-faa');self.assertEqual(acc,'CP_001.1');self.assertEqual((g['start0'],g['end0'],g['strand']),(3,15,-1));self.assertEqual(g['sequence_status'],'partial_CDS');self.assertEqual(g['sequence'],'MAA')
    def test_filtered_faa_ordinals_rejected(self):
        gs=self.genes(2);gs[0]['source_ordinal']=1;gs[1]['source_ordinal']=3
        with self.assertRaisesRegex(ValueError,'ordinals'):extract.ordered_genes(gs,True)
    def test_resume_identical_and_changed_input_refused(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);m=self.manifest(root);out=root/'out';self.run_extract(m,out);before=(out/'centers.tsv').read_bytes();self.run_extract(m,out);self.assertEqual(before,(out/'centers.tsv').read_bytes())
            with contextlib.closing(sqlite3.connect(out/'extraction.sqlite')) as db:self.assertEqual(db.execute('SELECT COUNT(*) FROM occurrences').fetchone()[0],2)
            (root/'proteins.faa').write_text((root/'proteins.faa').read_text().replace('ACC0','NEW0'))
            with self.assertRaisesRegex(ValueError,'changed'):self.run_extract(m,out)
    def test_shards_merge_equal_unsharded_and_missing_shard_refused(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);m=self.manifest(root,4);self.run_extract(m,root/'full');shards=root/'shards'
            self.run_extract(m,shards/'shard_0','--num-shards','2','--shard-index','0')
            with self.assertRaisesRegex(ValueError,'Missing'):merge.main(['--shards-root',str(shards),'--out',str(root/'bad')])
            self.run_extract(m,shards/'shard_1','--num-shards','2','--shard-index','1')
            with contextlib.redirect_stdout(io.StringIO()):merge.main(['--shards-root',str(shards),'--out',str(root/'merged')])
            for name in ('centers.tsv','slots.tsv','all_proteins.faa'):self.assertEqual((root/'full'/name).read_bytes(),(root/'merged'/name).read_bytes())
    def test_coordinate_golden_existing_genbank(self):
        rec=SeqIO.read(extract.REPO/'data/seed/genbank/RP4__BN000925.1.gb','genbank');genes=extract.decode_genes(rec,'genbank');hit=[g for g in genes if (g['start0'],g['end0'],g['strand'])==(27449,27659,1)];self.assertEqual(len(hit),1);self.assertEqual(hit[0]['aa_length'],69)
    def test_frozen_caller_sentinel_and_raw_fasta_path(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);rec=SeqIO.read(extract.REPO/'data/seed/genbank/RP4__BN000925.1.gb','genbank');SeqIO.write(rec,root/'rp4.fna','fasta');m=root/'inputs.tsv';m.write_text('dataset\tpath\tformat\nPLSDB\trp4.fna\tgenome-fasta\n');out=root/'out';self.run_extract(m,out)
            with (out/'centers.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
            hit=[r for r in rows if r['center_start0']=='27449' and r['center_end0']=='27659'];self.assertEqual(len(hit),1);self.assertEqual(hit[0]['center_aa_length'],'69')
    def test_vocab_none_falls_back_and_invalid_vocab_fails(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);tokens=['<cls>','<pad>','<eos>','<unk>','<mask>']+sorted(embed.AA);(root/'vocab.txt').write_text('\n'.join(tokens)+'\n');c=SimpleNamespace(vocab_list=None,vocab_size=len(tokens));self.assertEqual(embed.load_vocab(c,root)['<pad>'],1)
            c.vocab_size+=1
            with self.assertRaises(ValueError):embed.load_vocab(c,root)
    def test_packaged_head_dimension_mask_and_checksum(self):
        s=Scorer(extract.MODULE/'model');self.assertEqual(s.heads['center']['positions'].tolist(),[0]);self.assertEqual(s.heads['masked_pm3']['positions'].tolist(),[-3,-2,-1,1,2,3]);v=np.ones((1,2560),dtype='float32');self.assertTrue(np.isfinite(s.head('center',v)));self.assertTrue(np.isfinite(s.fusion(0,0)))
    def test_scoring_plumbing_with_explicit_synthetic_vector_fixture(self):
        # This is an I/O/arithmetic test, not a claim of biological embedding compatibility.
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);m=self.manifest(root);ex=root/'extract';self.run_extract(m,ex);emb=root/'embed';emb.mkdir();seqs=list(SeqIO.parse(ex/'all_proteins.faa','fasta'));rng=np.random.default_rng(42);vectors=rng.normal(size=(len(seqs),2560)).astype('float32');np.save(emb/'embeddings.npy',vectors)
            with (emb/'sequence_manifest.tsv').open('w',newline='') as f:
                w=csv.writer(f,delimiter='\t');w.writerow(['row_0based','sequence_sha256']);w.writerows((i,extract.digest(str(r.seq))) for i,r in enumerate(seqs))
            card=json.loads((extract.MODULE/'model/model_card.json').read_text());(emb/'embedding_validation.json').write_text(json.dumps(dict(status='passed',full_input_sha256=extract.file_sha(ex/'all_proteins.faa'),embeddings_file_sha256=extract.file_sha(emb/'embeddings.npy'))));(emb/'progress.json').write_text(json.dumps(dict(status='complete',signature=card['signature'])));(emb/'reference_validation.json').write_text(json.dumps(dict(passed=True,reference_sha256=card['embedding_reference_sha256'],test_fixture_only=True)))
            with contextlib.redirect_stdout(io.StringIO()):score.main(['--extraction',str(ex),'--embeddings',str(emb),'--out',str(root/'scores')])
            with (root/'scores/window_scores.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
            self.assertEqual(len(rows),34);self.assertTrue(all(r['functional_label']=='unassigned' for r in rows));self.assertTrue(any(r['fusion_score']=='' for r in rows));self.assertTrue(any(r['fusion_score']!='' for r in rows))
            with (root/'scores/center_scores.tsv').open() as f:centers=list(csv.DictReader(f,delimiter='\t'))
            self.assertEqual(len(centers),17);self.assertTrue(all(r['distinct_windows']=='1' for r in centers))
            for r in rows:
                if r['fusion_score']!='':self.assertAlmostEqual(float(r['fusion_score']),Scorer(extract.MODULE/'model').fusion(float(r['center_score']),float(r['context_score'])),places=10)

if __name__=='__main__':unittest.main()
