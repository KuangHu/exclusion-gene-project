#!/usr/bin/env python3
"""Six-frame every ADJACENT anchor gap. Three deliverables from one pass.

Why now: the probe settled it. The 206 plasmids with a callable VirB7 have median
length 56 aa (Q1 43, Q3 56) -- the same short range as pKM101's uncalled 48 aa
traN -- so length is NOT what distinguishes them. But 94% of them overlap no
neighbouring CDS, while pKM101's overlaps by 11 bp. Pyrodigal's floor is 30 aa, so
it can call proteins this short; OVERLAP is the discriminator, and the overlap
fraction <=50% rule is the matching remedy.

Deliverables:
  1. VirB7 recovery -- the position fails on 97.2% of plasmids and recall
     recovered exactly zero, because the proteins are not in the cache at all.
  2. unresolved_gaps at ADJACENT-PAIR granularity. The catalogue currently counts
     genes between the FIRST and LAST anchor, so a plasmid with scattered anchors
     counts nearly everything (max 2,815; tier A Q3 = 32). That number overstates
     the candidate pool and must not feed anything downstream.
  3. the per-gap material the per-class slot criteria will consume.

Controls, both mandatory:
  POSITIVE  pKM101 traN 7899..8045 recovered from U09868.1 -- a real gene the
            caller really misses.
  BASELINE  six-frame must re-find every Pyrodigal gene inside the window,
            matched on (strand, stop_coord). Stop-to-stop enumeration is a strict
            superset of any one-ORF-per-stop caller.
"""
import argparse
import collections
import csv
import os
import re
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
import decontaminate as DC
import orf_caller
import sixframe
from assertions import require_compute_node
from Bio.Seq import Seq

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
PLSDB = "/global/scratch/users/kh36969/plsdb/sequences.fasta"
EXTRA = "/global/scratch/users/kh36969/exclusion_gene/mash/extra_controls.fa"
MIN_AA, MAX_OVL = 40, 0.50
LIPO = re.compile(r"[LVI][ASTVIG][GASN]C")


def novel(seq, lo, hi, calls, py_stops):
    if hi - lo + 1 < 3 * MIN_AA + 3:
        return []
    out = []
    for o in sixframe.call(seq[lo - 1:hi], min_aa=MIN_AA, mode="first_start"):
        o["start"] += lo - 1; o["end"] += lo - 1
        o["stop_coord"] = o["end"] if o["strand"] == 1 else o["start"]
        if (o["strand"], o["stop_coord"]) in py_stops:
            continue
        span = o["end"] - o["start"] + 1
        ov = sum(max(0, min(o["end"], c["end"]) - max(o["start"], c["start"]) + 1)
                 for c in calls)
        if ov > MAX_OVL * span:
            continue
        out.append(o)
    return out


def aa_of(seq, o):
    s = Seq(seq[o["start"] - 1:o["end"]])
    if o["strand"] == -1:
        s = s.reverse_complement()
    return str(s.translate(table=11)).rstrip("*")


def positive_control():
    from Bio import SeqIO
    rec = next(SeqIO.parse(os.path.join(PROJ, "data", "seed", "genbank",
                                        "pKM101__U09868.1.gb"), "genbank"))
    seq = str(rec.seq).upper()
    calls = orf_caller.call(seq)
    py = {(c["strand"], c["stop_coord"]) for c in calls}
    got = novel(seq, 7500, 8400, calls, py)
    ok = any((o["strand"], o["stop_coord"]) == (1, 8045) for o in got)
    print("POSITIVE CONTROL pKM101 traN 7899..8045: %s"
          % ("RECOVERED" if ok else "*** NOT RECOVERED ***"))
    if not ok:
        raise SystemExit("positive control failed; nothing else may be read")


def main():
    require_compute_node()
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--cpus", type=int, default=8)
    a = ap.parse_args()
    import hashlib
    import pyhmmer
    from Bio import SeqIO
    positive_control()
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    seqs, owner, idx = [], [], []
    aas = []
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"):
            continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">":
                nm = line[1:].rstrip("\n")
            else:
                f = nm.split("|")
                owner.append(f[0]); idx.append(int(f[1])); aas.append(line.rstrip("\n"))
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(seqs)).encode(),
                    sequence=line.rstrip("\n")).digitize(alpha))
    anch = collections.defaultdict(dict)
    for name, fam, ev, minaa in DC.ANCHOR_FAMILIES:
        p = os.path.join(HMM, fam + ".hmm")
        if not os.path.exists(p):
            continue
        with pyhmmer.plan7.HMMFile(p) as fh:
            model = next(iter(fh))
        kw = {"E": ev} if ev else {"bit_cutoffs": "gathering"}
        for top in pyhmmer.hmmsearch([model], seqs, cpus=a.cpus, **kw):
            for h in top:
                i = int(dec(h.name))
                if minaa and len(aas[i]) <= minaa:
                    continue
                cur = anch[owner[i]].get(idx[i])
                if cur is None or h.score > cur[1]:
                    anch[owner[i]][idx[i]] = (name, h.score)
    print("anchored plasmids: %d" % len(anch), flush=True)

    _h = lambda x: int(hashlib.md5(x.encode()).hexdigest(), 16)
    want = {k for k in anch if a.nshards == 1 or _h(k) % a.nshards == a.shard}
    rows = []
    stat = collections.Counter()
    for path in (PLSDB, EXTRA):
        if not os.path.exists(path):
            continue
        for rec in SeqIO.parse(path, "fasta"):
            if rec.id not in want:
                continue
            seq = str(rec.seq).upper()
            calls = orf_caller.call(seq)
            py = {(c["strand"], c["stop_coord"]) for c in calls}
            A = sorted(anch[rec.id])
            for u, v in zip(A, A[1:]):
                if v - u < 2:
                    continue                       # adjacent genes, no gap
                lo = calls[u]["end"] + 1
                hi = calls[v]["start"] - 1
                if hi < lo:
                    continue
                between = calls[u + 1:v]
                # BASELINE inside the anchor-to-anchor span
                span = sixframe.call(seq[calls[u]["start"] - 1:calls[v]["end"]],
                                     min_aa=MIN_AA, mode="first_start")
                sst = set()
                for o in span:
                    s2 = o["start"] + calls[u]["start"] - 1
                    e2 = o["end"] + calls[u]["start"] - 1
                    sst.add((o["strand"], e2 if o["strand"] == 1 else s2))
                for c in between:
                    if c["aa_len"] < MIN_AA:
                        continue
                    if (c["strand"], c["stop_coord"]) in sst:
                        stat["baseline_ok"] += 1
                    else:
                        stat["BASELINE_MISS"] += 1
                nv = novel(seq, lo, hi, calls, py)
                stat["gaps"] += 1
                if nv:
                    stat["gaps_with_novel"] += 1
                rows.append({
                    "accession": rec.id,
                    "anchor_left": anch[rec.id][u][0], "anchor_right": anch[rec.id][v][0],
                    "gap_bp": hi - lo + 1, "n_called_between": len(between),
                    "called_lengths": ";".join(str(c["aa_len"]) for c in between),
                    "n_sixframe_novel": len(nv),
                    "novel_lengths": ";".join(str(o["aa_len"]) for o in nv),
                    "novel_lipobox": ";".join(
                        str(int(bool(LIPO.search(aa_of(seq, o)[:40])))) for o in nv),
                    "novel_coords": ";".join("%d..%d%s" % (o["start"], o["end"],
                                             "+" if o["strand"] == 1 else "-") for o in nv),
                    "novel_seqs": ";".join(aa_of(seq, o) for o in nv)})
    out = os.path.join(PROJ, "data", "positional")
    os.makedirs(out, exist_ok=True)
    sfx = "" if a.nshards == 1 else "_%03d" % a.shard
    dest = os.path.join(out, "interanchor_gaps%s.tsv" % sfx)
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    for k, v in sorted(stat.items()):
        print("%-20s %d" % (k, v))
    if stat["BASELINE_MISS"]:
        print("*** BASELINE FAILED (%d) -- do not read any other number ***"
              % stat["BASELINE_MISS"])
    print("wrote %s (%d gaps)" % (dest, len(rows)))


if __name__ == "__main__":
    main()
