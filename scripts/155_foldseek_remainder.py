#!/usr/bin/env python3
"""Structural comparison of the unnamed remainder against the named families.

Closes the one gap sequence search cannot: category 2 members whose fold is
conserved but whose sequence has diverged past HMM detection. The hold-out showed
this is real -- TrbK, the family known to sit at 41% identity between RP4 and
R751, is recovered at only 50-75%, while the large well-sampled families reach
94-100%. A genuinely new family is by definition sparsely sampled, i.e.
TrbK-shaped, so the sensitivity that applies to the candidate pool is the low one.

NO STRUCTURE PREDICTION IS NEEDED. Foldseek's ProstT5 mode predicts the 3Di
structural alphabet directly from sequence. ESMFold/AF2 are unnecessary, and AFDB
lookup was never possible anyway: these are ab-initio Pyrodigal calls with no
UniProt identifiers.

DESIGN. Query = the unnamed remainder. Target = the 801 unique NAMED slot
occupants, which carry family labels. This is self-contained -- it asks exactly
"is this unnamed protein structurally a member of a family we already have?"
without needing AFDB or PDB.

THE CONTROL, and it gates the result: the named set is also searched against
itself with self-hits excluded. A named protein must recover its OWN family
structurally. If that fails, the 3Di representation is not informative for
proteins of this size (40-260 aa) and no negative result about the remainder can
be read from it.

TM-score >= 0.5 is the conventional same-fold threshold and is what the
operational definition of "new" uses.
"""
import collections, csv, os, subprocess, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
SOF = os.path.join(PROJ, "data", "release", "v1.1", "slot_occupant_families.tsv")
FS = os.environ.get("FOLDSEEK_BIN", "/global/scratch/users/kh36969/bin/foldseek/bin/foldseek")
W = "/global/scratch/users/kh36969/foldseek_db/prostt5_weights"
WORK = "/global/scratch/users/kh36969/foldseek_db/work"
EEX = ["TIGR04359", "NF033894", "NF041429", "NF033891", "PF10624", "PF14729"]
# ProstT5 databases contain 3Di SEQUENCES, not coordinates -- there is no _ca
# file, so alntmscore cannot be computed (job 26102692). TM-score would require
# actually folding all 3,098 proteins (ESMFold), a far larger job.
#
# Instead the 3Di-space E-value is used, and the threshold is CALIBRATED on our
# own labelled data: the 801 named proteins carry family labels, so the
# same-family vs different-family score distributions say where to cut. That is
# stronger than importing the TM>=0.5 convention, which was never measured on
# proteins of this size anyway.
TM = None


def main():
    require_compute_node()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    os.makedirs(WORK, exist_ok=True)

    rows = [r for r in csv.DictReader(open(SOF), delimiter="\t")
            if r["slot_status"] in ("candidate", "eex_occupied") and r["slot_occupant_coords"]]
    want = collections.defaultdict(set)
    for r in rows: want[r["accession"]].add(r["slot_occupant_coords"])
    seq_of = {}
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"): continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">": nm = line[1:].rstrip("\n")
            else:
                f = nm.split("|"); acc = f[0]
                if acc in want:
                    k = "%s..%s" % (f[2], f[3])
                    if k in want[acc]: seq_of[(acc, k)] = line.rstrip("\n")
    named, cand = {}, set()
    for r in rows:
        s = seq_of.get((r["accession"], r["slot_occupant_coords"]))
        if not s: continue
        if r["slot_status"] == "eex_occupied":
            named[s] = r["slot_occupant_family"].split(";")[0]
        else: cand.add(s)
    cand = sorted(cand)

    dig = lambda L: [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=s).digitize(alpha)
                     for i, s in enumerate(L)]
    C = dig(cand)
    drop = set()
    for f in EEX:
        p = os.path.join(HMM, f + ".hmm")
        if not os.path.exists(p): continue
        with pyhmmer.plan7.HMMFile(p) as fh: m = next(iter(fh))
        for top in pyhmmer.hmmsearch([m], C, cpus=8, E=1e-3):
            for h in top: drop.add(int(dec(h.name)))
    rem = [cand[i] for i in range(len(cand)) if i not in drop]
    nm_list = sorted(named)
    print("named %d unique | candidates %d unique | remainder after HMM %d"
          % (len(nm_list), len(cand), len(rem)), flush=True)

    def write(path, seqs, pref):
        with open(path, "w") as fh:
            for i, s in enumerate(seqs): fh.write(">%s%d\n%s\n" % (pref, i, s))
    qf, tf = os.path.join(WORK, "rem.fa"), os.path.join(WORK, "named.fa")
    write(qf, rem, "r"); write(tf, nm_list, "n")

    def run(cmd):
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode:
            print("  FAILED: %s\n  %s" % (" ".join(cmd[:4]), r.stderr[-400:])); sys.exit(1)

    # ProstT5 3Di prediction is the expensive step -- it timed out at 1 h on CPU
    # (job 26100374). --gpu 1 offloads it; the flag is ignored if no CUDA device
    # is present, so the same script runs either way.
    gpu = ["--gpu", "1"] if os.environ.get("FS_GPU", "1") == "1" else []
    for fa, db in ((qf, "remdb"), (tf, "nameddb")):
        run([FS, "createdb", fa, os.path.join(WORK, db), "--prostt5-model", W] + gpu)
    print("  3Di databases built (ProstT5, no structure prediction)", flush=True)

    # --- CONTROL FIRST: can a named protein find its own family structurally? ---
    # easy-search expects STRUCTURE files. We have prebuilt 3Di databases from
    # createdb --prostt5-model, so the correct path is search + convertalis.
    ctl = os.path.join(WORK, "ctl.tsv")
    ndb = os.path.join(WORK, "nameddb"); rdb = os.path.join(WORK, "remdb")
    run([FS, "search", ndb, ndb, os.path.join(WORK, "ctlres"),
         os.path.join(WORK, "tmp1"), "-e", "10", "--exhaustive-search", "1"])
    run([FS, "convertalis", ndb, ndb, os.path.join(WORK, "ctlres"), ctl,
         "--format-output", "query,target,evalue,bits,fident,alnlen"])
    same, diff = [], []
    for line in open(ctl):
        f = line.rstrip("\n").split("\t")
        q, t, ev = f[0], f[1], float(f[2])
        if q == t: continue
        fq = named[nm_list[int(q[1:])]]; ft = named[nm_list[int(t[1:])]]
        (same if fq == ft else diff).append(ev)
    print("\n=== CONTROL / CALIBRATION: 3Di scores within vs across known families ===")
    print("  same-family pairs      %7d" % len(same))
    print("  different-family pairs %7d" % len(diff))
    if not same:
        print("  CONTROL FAILED -- no same-family structural pair found at all.")
        print("  3Di is not informative here; no negative result can be read.")
        return 1
    same.sort(); diff.sort()
    def pct(v, p): return v[min(len(v)-1, int(len(v)*p))]
    print("  same-family E:  median %.2e  90th %.2e" % (pct(same,0.5), pct(same,0.9)))
    if diff:
        print("  diff-family E:  median %.2e  10th %.2e" % (pct(diff,0.5), pct(diff,0.1)))
    # threshold = the E that keeps 90% of same-family pairs
    CUT = pct(same, 0.9)
    fp = sum(1 for e in diff if e <= CUT)
    print("  chosen cut E<=%.2e keeps 90%% of same-family pairs" % CUT)
    print("  it also admits %d of %d different-family pairs (%.1f%% false rate)"
          % (fp, len(diff), 100.0*fp/max(1,len(diff))))
    if diff and fp/max(1,len(diff)) > 0.5:
        print("  CONTROL FAILED -- the cut cannot separate families; over half of")
        print("  different-family pairs pass. No negative result can be read.")
        return 1
    print("  CONTROL PASSES -- proceeding with E<=%.2e" % CUT)

    res = os.path.join(WORK, "rem_vs_named.tsv")
    run([FS, "search", rdb, ndb, os.path.join(WORK, "remres"),
         os.path.join(WORK, "tmp2"), "-e", "10", "--exhaustive-search", "1"])
    run([FS, "convertalis", rdb, ndb, os.path.join(WORK, "remres"), res,
         "--format-output", "query,target,evalue,bits,fident,alnlen"])
    best = {}
    for line in open(res):
        f = line.rstrip("\n").split("\t")
        q, t, ev = f[0], f[1], float(f[2])
        if q not in best or ev < best[q][1]: best[q] = (t, ev)
    struct = {q: v for q, v in best.items() if v[1] <= CUT}
    print("\n=== REMAINDER vs NAMED, structural (3Di) ===")
    print("  remainder sequences               %6d" % len(rem))
    print("  with a structural hit E<=%.1e  %6d  (%.1f%%)  <- CATEGORY 2, structural"
          % (CUT, len(struct), 100.0*len(struct)/len(rem)))
    print("  no structural hit              %6d  (%.1f%%)  <- survives to candidate pool"
          % (len(rem)-len(struct), 100.0*(len(rem)-len(struct))/len(rem)))
    fam = collections.Counter(named[nm_list[int(t[1:])]] for t, _ in struct.values())
    print("  assigned to: %s" % ", ".join("%s %d" % kv for kv in fam.most_common()))
    dest = os.path.join(PROJ, "data", "anchors", "foldseek_remainder.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["query_idx", "aa_len", "best_named_family", "evalue_3di", "assigned"])
        for i, s in enumerate(rem):
            q = "r%d" % i
            b = best.get(q)
            w.writerow([i, len(s), named[nm_list[int(b[0][1:])]] if b else "",
                        ("%.3e" % b[1]) if b else "", int(bool(b and b[1] <= CUT))])
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    sys.exit(main() or 0)
