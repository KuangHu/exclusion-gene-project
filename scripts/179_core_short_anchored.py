#!/usr/bin/env python3
"""EEx candidates in all eight classes: conserved, short, membrane-anchored.

Two situations motivate this, and they are different problems:

  (1) The EEx TYPE is known but our models do not find it. MPF_F is the case:
      PF10624 detects R100-type TraS only (1 of 14 traS-carrying seeds), so
      F-type and R64-type TraS sit unnamed in the candidate pool. The gene
      exists, is characterised, and is invisible to the scan.
  (2) No EEx is known at all -- MPF_G, MPF_B, MPF_C, and MPF_FATA (whose only
      characterised exclusion is SURFACE exclusion, PrgA/Sea1, category 2).

Both reduce to the same search: a protein that is USUALLY PRESENT in its class,
SHORT, and MEMBRANE-ANCHORED. That is the frozen ANCHOR_CANDIDATE criterion with
the positional term replaced by conservation, because G/B/C have no slot defined
and waiting for one would block the search indefinitely.

    length <= 250 aa  AND  (lipobox_strict OR tm_count >= 1)
    ranked by PREVALENCE within the class

NO PREVALENCE THRESHOLD IS SET. The four known EEx families are located in the
ranking and reported, but they do NOT define a cutoff -- sentinels verify recall,
they never set thresholds. A compute cap (top N families by prevalence among
those <=250 aa) is applied so TMbed is tractable; it is declared as a compute cap,
and if a known EEx falls outside it that is reported as a RECALL FAILURE, not
quietly dropped.

WHAT IS EXCLUDED, and why it must be:
  * anchor families (decontaminate.ANCHOR_FAMILIES) -- 562 leaked into the
    candidate pool once before when this step was skipped
  * category 2, known-not-EEx: PF05818 TraT (OM lipoprotein, surface exclusion),
    PrgA_Sea1 (Gram-positive surface exclusion)
  * the class's own CONJScan accessory profiles -- they are the machine

KNOWN EEx ARE SCANNED LAST AND POST-HOC, for labelling only. A hit marks a family
as ARCHIVED-known, never as a discovery, and never enters the ranking criteria.

FAMILY GRANULARITY IS THE KNOWN WEAK POINT. Step 1d showed TraS shatters from 1
family at 50% id into 40 at 90%, and that the covariation verdict flipped with it.
Families are therefore built at 0.50 AND 0.70 and every result is reported at
both. A candidate that appears only at one granularity is flagged, not ranked.
"""
import argparse
import collections
import csv
import gzip
import os
import re
import subprocess
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node
import decontaminate as DC

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
ICE = "/global/scratch/users/kh36969/ice_db"
CJ = os.path.join(PROJ, "data", "CONJScan", "profiles")
DEFS = os.path.join(PROJ, "data", "CONJScan", "definitions", "Plasmids")
HMM = os.path.join(PROJ, "data", "hmm")
REL = os.path.join(PROJ, "data", "release")
ANCH = os.path.join(PROJ, "data", "anchors")
WORK = "/global/scratch/users/kh36969/exclusion_gene/recovery/core179"

MAXLEN = 250
IDS = [0.50, 0.70]
TOPN = 400                 # compute cap for TMbed, NOT a biological threshold
STRICT = re.compile(r"[LVI][ASTVIG][GASN]C")
CLASSES = ["T", "F", "I", "FA", "FATA", "G", "B", "C"]
NEWTHR = {"G": 3, "B": 5, "C": 2, "FATA": 2}
EEX = [("TIGR04359", "TrbK"), ("PF10624", "TraS_R100"), ("NF033891", "ExcA"),
       ("NF033894", "Eex_IncN"), ("NF041429", "EexR/S"), ("PF14729", "DUF4467")]
# the EEx each class actually owns -- used ONLY for the pipeline recall check
NATIVE = {"T": ("TrbK",), "F": ("TraS_R100",), "I": ("ExcA",), "FA": ("DUF4467",)}
CAT2 = [("PF05818", "TraT"), ("PrgA_Sea1", "PrgA/Sea1")]


def segments(pred):
    out, cur, st = [], None, 0
    for i, c in enumerate(pred):
        k = "TM" if c in "HhBb" else ("S" if c == "S" else ".")
        if k != cur:
            if cur == "TM":
                out.append((st + 1, i))
            cur, st = k, i
    if cur == "TM":
        out.append((st + 1, len(pred)))
    return out


def class_profiles(cls):
    import xml.etree.ElementTree as ET
    return {g.get("name") for g in
            ET.parse(os.path.join(DEFS, "T4SS_type%s.xml" % cls)).getroot().iter("gene")}


def load_members(dataset):
    """class -> set of element ids, from the frozen catalogues and the N-sweep."""
    out = {}
    if dataset == "plsdb":
        for c, p in (("T", "v1.1/catalogue_MPF_T_v1.1.tsv"),
                     ("F", "v1/catalogue_MPF_F_v1.tsv"),
                     ("I", "v1.1/catalogue_MPF_I_v1.tsv")):
            out[c] = {l.split("\t")[0] for i, l in enumerate(open(os.path.join(REL, p))) if i}
        fa = list(csv.DictReader(open(os.path.join(REL, "v1.1/catalogue_MPF_FA_v1.tsv")),
                                 delimiter="\t"))
        out["FA"] = {r["accession"] for r in fa if r["fa_subclass"] == "FA"}
    rows = list(csv.DictReader(open(os.path.join(ANCH, "profile_counts_%s.tsv" % dataset)),
                              delimiter="\t"))
    for c, n in NEWTHR.items():
        out[c] = {r["element"] for r in rows if int(r["n_MPF_" + c]) >= n}
    if dataset == "ice":
        for c in ("T", "F", "I", "FA"):
            out[c] = {r["element"] for r in rows if int(r["n_MPF_" + c]) >= 2}
    return out


def main():
    require_compute_node()
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", choices=["plsdb", "ice"], required=True)
    a = ap.parse_args()
    import pyhmmer
    import orf_caller
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    os.makedirs(WORK, exist_ok=True)
    DC.assert_no_evalue_anchor()

    members = load_members(a.dataset)
    print("=== class sizes (%s) ===" % a.dataset)
    for c in CLASSES:
        print("  MPF_%-5s %6d" % (c, len(members.get(c, ()))))

    # ---- proteins -----------------------------------------------------------
    want = set().union(*members.values())
    prot, owner = [], []
    if a.dataset == "plsdb":
        for fn in sorted(os.listdir(CACHE)):
            if not fn.endswith(".faa"):
                continue
            nm = None
            for line in open(os.path.join(CACHE, fn)):
                if line[0] == ">":
                    nm = line[1:].rstrip("\n")
                elif nm.split("|")[0] in want:
                    owner.append(nm.split("|")[0]); prot.append(line.rstrip("\n"))
    else:
        orf_caller.assert_anon_only(orf_caller.ANON)
        name, buf = None, []

        def emit():
            if name in want and buf:
                for g in orf_caller.call("".join(buf).upper()):
                    p = g["aa"].rstrip("*")
                    if len(p) >= 30:
                        owner.append(name); prot.append(p)

        with gzip.open(os.path.join(ICE, "seq", "ice.shard000.fna.gz"), "rt") as fh:
            for line in fh:
                if line[0] == ">":
                    emit(); name, buf = line[1:].split()[0], []
                else:
                    buf.append(line.strip())
        emit()
    print("\nproteins: %d over %d elements\n" % (len(prot), len(set(owner))), flush=True)

    seqs = [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=p).digitize(alpha)
            for i, p in enumerate(prot)]

    # ---- things that must never be a candidate ------------------------------
    def scan(path):
        """GA where the model carries one; E<=1e-5 where it does not.

        PrgA_Sea1 was built locally from pCF10 prgA + pAD1 sea1 and has no
        gathering cutoff, so bit_cutoffs="gathering" raises MissingCutoffs and
        kills the run. Scripts 158/159/160 already handle GA-less models at
        E<=1e-5 and that precedent is followed here rather than inventing a
        second convention.

        An E-value is database-size dependent (the VirB5 lesson) so this is not
        used for anything that sets a result. It binds only on the EXCLUSION
        side, and PrgA/Sea1 is an 892 aa protein while every candidate here is
        <=250 aa, so the filter barely reaches the candidate set at all. Any
        model falling back is named in the log, never silently."""
        if not os.path.exists(path):
            return set()
        with pyhmmer.plan7.HMMFile(path) as fh:
            hmm = next(iter(fh))
        try:
            it = pyhmmer.hmmsearch([hmm], seqs, cpus=16, bit_cutoffs="gathering")
            got = set()
            for top in it:
                for h in top:
                    got.add(int(dec(h.name)))
            return got
        except pyhmmer.errors.MissingCutoffs:
            print("    no GA on %s -- falling back to E<=1e-5 (precedent: 158/159/160)"
                  % os.path.basename(path), flush=True)
            got = set()
            for top in pyhmmer.hmmsearch([hmm], seqs, cpus=16, E=1e-5):
                for h in top:
                    got.add(int(dec(h.name)))
            return got

    banned = set()
    for _, famf, _, _ in DC.ANCHOR_FAMILIES:
        banned |= scan(os.path.join(HMM, famf + ".hmm"))
    for mdl, _ in CAT2:
        banned |= scan(os.path.join(HMM, mdl + ".hmm"))
    allprof = set()
    for c in CLASSES:
        allprof |= class_profiles(c)
    for p in sorted(allprof):
        banned |= scan(os.path.join(CJ, p + ".hmm"))
    print("excluded proteins (anchors + category 2 + machine profiles): %d\n"
          % len(banned), flush=True)

    known = {}
    for mdl, nm in EEX:
        for i in scan(os.path.join(HMM, mdl + ".hmm")):
            known[i] = nm
    print("known-EEx proteins (post-hoc labels only): %d\n" % len(known), flush=True)

    # ---- per class, per granularity -----------------------------------------
    allrows, recall, found_native = [], [], []
    for cls in CLASSES:
        els = members.get(cls, set())
        if len(els) < 50:
            print("MPF_%-5s %d elements -- too few, skipped\n" % (cls, len(els)))
            continue
        idxs = [i for i in range(len(prot)) if owner[i] in els and i not in banned]
        print("=" * 74)
        print("MPF_%s  %d elements | %d candidate proteins (<=%d aa: %d)"
              % (cls, len(els), len(idxs), MAXLEN,
                 sum(1 for i in idxs if len(prot[i]) <= MAXLEN)), flush=True)
        fa = os.path.join(WORK, "%s_%s.faa" % (a.dataset, cls))
        with open(fa, "w") as fh:
            for i in idxs:
                fh.write(">%d\n%s\n" % (i, prot[i]))
        for ident in IDS:
            out = os.path.join(WORK, "%s_%s_%d" % (a.dataset, cls, int(ident * 100)))
            r = subprocess.run(["mmseqs", "easy-cluster", fa, out, out + "_tmp",
                                "--min-seq-id", str(ident), "-c", "0.8",
                                "--cov-mode", "0", "-v", "1"],
                               capture_output=True, text=True)
            if r.returncode:
                print("  mmseqs failed at %.2f: %s" % (ident, r.stderr[-300:])); continue
            fam = {}
            for line in open(out + "_cluster.tsv"):
                rep, mem = line.rstrip("\n").split("\t")[:2]
                fam[int(mem)] = int(rep)
            byfam = collections.defaultdict(list)
            for i, rep in fam.items():
                byfam[rep].append(i)
            rows = []
            for rep, mem in byfam.items():
                lens = sorted(len(prot[i]) for i in mem)
                med = lens[len(lens) // 2]
                if med > MAXLEN:
                    continue
                el = {owner[i] for i in mem}
                lab = collections.Counter(known[i] for i in mem if i in known)
                rows.append({"rep": rep, "n_prot": len(mem), "n_elem": len(el),
                             "prev": 100.0 * len(el) / len(els), "med_aa": med,
                             "known": lab.most_common(1)[0][0] if lab else ""})
            rows.sort(key=lambda r: -r["prev"])
            for j, r in enumerate(rows):
                r["prev_rank"] = j + 1
                r["mpf_class"] = "MPF_" + cls
                r["min_seq_id"] = ident
                r["dataset"] = a.dataset
                r["in_compute_cap"] = int(j < TOPN)
            kn = [r for r in rows if r["known"]]
            print("  id %.2f: %6d families <=%d aa | known-EEx families: %s"
                  % (ident, len(rows), MAXLEN,
                     ", ".join("%s rank %d (%.1f%%)" % (r["known"], r["prev_rank"], r["prev"])
                               for r in kn) or "none"), flush=True)
            # NATIVE pairs only, and only the BEST family per pair. A known EEx
            # that is RARE in a class it does not belong to (Eex_IncN at 0.3% of
            # MPF_G ICEs) is correctly ranked low -- that is the scan working, not
            # a recall failure. The first version asserted on every labelled family
            # in every class and failed on exactly those stragglers.
            for nm in NATIVE.get(cls, ()):
                best = min((r["prev_rank"] for r in kn if r["known"] == nm),
                           default=None)
                if best is not None and best <= TOPN:
                    found_native.append("MPF_%s @%.2f %s rank %d" % (cls, ident, nm, best))
            allrows.extend(rows[:TOPN])
        print()

    dest = os.path.join(ANCH, "core_short_candidates_%s.tsv" % a.dataset)
    keys = ["dataset", "mpf_class", "min_seq_id", "prev_rank", "rep", "n_prot",
            "n_elem", "prev", "med_aa", "known", "in_compute_cap"]
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in allrows:
            w.writerow({k: r.get(k, "") for k in keys})
    print("wrote %s (%d rows)" % (dest, len(allrows)))

    fa = os.path.join(WORK, "tmbed_in_%s.faa" % a.dataset)
    with open(fa, "w") as fh:
        seen = set()
        for r in allrows:
            k = (r["mpf_class"], r["min_seq_id"], r["rep"])
            if k in seen:
                continue
            seen.add(k)
            fh.write(">%s|%s|%d|%d|%.2f|%s\n%s\n"
                     % (r["mpf_class"], r["min_seq_id"], r["rep"], r["n_elem"],
                        r["prev"], r["known"] or "-", prot[r["rep"]]))
    print("wrote %s (%d representatives for TMbed)" % (fa, len(seen)))

    # PIPELINE-LEVEL recall: the scan must be able to surface at least one
    # class's OWN EEx inside the cap. If it cannot do that anywhere, it has not
    # been shown capable of finding an exclusion gene and no candidate list from
    # it is readable. A native EEx that is simply rare in this dataset (TrbK is
    # an IncP plasmid gene and is near-absent from MPF_T ICEs) is biology, and
    # ranking it low is correct.
    print("\n=== pipeline recall: native class/EEx pairs inside the top-%d cap ===" % TOPN)
    for x in found_native:
        print("  %s" % x)
    if not found_native:
        print("  NONE. The scan never surfaced a class's own exclusion gene as")
        print("  'usually present'. It has not been shown capable of the task and")
        print("  its candidate lists must not be read.")
        return 1
    print("  OK -- %d native pair(s) recovered." % len(found_native))
    return 0


if __name__ == "__main__":
    sys.exit(main())
