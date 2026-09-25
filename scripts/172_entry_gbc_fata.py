#!/usr/bin/env python3
"""MPF_G / B / C / FATA §2 -- entry criterion, measured. Plasmids and ICEs.

The pilot (script 171, job 26151746) showed the four prototypes are separated by
their accessory profiles: ICEHin1056 16/20 G, CTnDOT 12/12 B, pCC7120alpha 8/8 C,
pCF10 8/27 FATA, with two cross-class hits in the entire 4x8 matrix. That is
RECALL on known material. It says nothing about what a threshold admits from a
database, which is what this measures.

WHY BOTH DATASETS. Guglielmini's G and B are ICE-dominant -- ICEHin1056 is an ICE,
CTnDOT is a conjugative transposon. Running PLSDB alone would return near-zero for
both and that zero would read as "not detectable" when it actually means "not on
plasmids". The ICE set is scanned in the same run, with the same profiles, so the
two numbers are comparable.

ALL EIGHT CLASSES ARE SCANNED, not just the four new ones. Each class is then the
others' background, which is the only way a per-class threshold is readable. The
MPF_FA precedent is the warning: at `>=1` accessory profile the FA residual looked
strong at 73.9% until the MPF_F background was printed next to it at 35.2%.

SHARED PROFILES ARE REMOVED. Three profiles are claimed by more than one class
definition -- T4SS_G_tfc7 (F,G), T4SS_I_traE (I,G), T4SS_T_virB1 (T,F,G). They
cannot discriminate, so they are dropped from every class's count. MPF_G therefore
scores out of 17, not 20; the prototype's 16/20 includes tfc7.

  BASELINE, PLSDB: the 1,170 plasmids already called `fa_subclass == FATA` in
  catalogue_MPF_FA_v1.tsv must be re-admitted by the independent FATA count here.
  That call was made by a different script from a different profile tally, so it
  is a genuine outside check. If recall on it is not near-total the scan is wrong
  and no threshold table may be read.

  PRE-REGISTERED, ICE: host phylum is recorded in the ICEberg table and played no
  part in any profile. Written down BEFORE the run, from the published class
  descriptions -- B -> Bacteroidota, C -> Cyanobacteria, G -> Pseudomonadota,
  FATA -> Bacillota/Actinomycetota. These are CHECKS, never filters: a class whose
  hits scatter across unrelated phyla is picking up noise, but a class that lands
  where it should has an orthogonal line of support.

No threshold is chosen here. The full table is printed for every N and the choice
is a separate, documented step -- sentinels verify recall, they never set cutoffs.
"""
import argparse
import collections
import csv
import gzip
import os
import sys
import xml.etree.ElementTree as ET

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node
import orf_caller

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
ICE = "/global/scratch/users/kh36969/ice_db"
CJ = os.path.join(PROJ, "data", "CONJScan", "profiles")
DEFS = os.path.join(PROJ, "data", "CONJScan", "definitions", "Plasmids")
REL = os.path.join(PROJ, "data", "release")
OUT = os.path.join(PROJ, "data", "anchors")

CLASSES = ["T", "F", "I", "FA", "FATA", "G", "B", "C"]
NEW = ["G", "B", "C", "FATA"]
GENERIC = ("T4SS_virb4", "T4SS_I_traU", "T4SS_t4cp1", "T4SS_t4cp2", "T4SS_tcpA")
CATS = [("MPF_T", os.path.join(REL, "v1.1", "catalogue_MPF_T_v1.1.tsv")),
        ("MPF_F", os.path.join(REL, "v1", "catalogue_MPF_F_v1.tsv")),
        ("MPF_I", os.path.join(REL, "v1.1", "catalogue_MPF_I_v1.tsv"))]
FACAT = os.path.join(REL, "v1.1", "catalogue_MPF_FA_v1.tsv")
# The PREDICTION is the taxon, not the string. ICEberg labels phyla with the older
# NCBI names while these were written in GTDB names, so the first run scored three
# of four as "NOT as predicted" purely on nomenclature -- Pseudomonadota IS
# Proteobacteria, Cyanobacteriota IS Cyanobacteria, Bacillota IS Firmicutes (which
# GTDB further splits into Firmicutes, Firmicutes_A, ...). The taxa asserted before
# the run are unchanged; only the synonyms they are matched against are recorded here.
PREDICT = {"B": ("Bacteroidota", "Bacteroidetes"),
           "C": ("Cyanobacteriota", "Cyanobacteria"),
           "G": ("Pseudomonadota", "Proteobacteria"),
           "FATA": ("Bacillota", "Firmicutes", "Firmicutes_A", "Firmicutes_B",
                    "Firmicutes_C", "Actinomycetota", "Actinobacteriota")}
NMAX = 6


def profile_sets():
    """Per-class accessory profiles, with any profile claimed by >1 class removed."""
    acc = {}
    for c in CLASSES:
        names = {g.get("name") for g in
                 ET.parse(os.path.join(DEFS, "T4SS_type%s.xml" % c)).getroot().iter("gene")}
        acc[c] = sorted(n for n in names
                        if n not in GENERIC and not n.startswith("T4SS_MOB"))
    owner = collections.Counter(p for c in CLASSES for p in acc[c])
    shared = sorted(p for p, n in owner.items() if n > 1)
    disc = {c: [p for p in acc[c] if owner[p] == 1] for c in CLASSES}
    return disc, shared


def load_plsdb():
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    seqs, owner = [], []
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"):
            continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">":
                nm = line[1:].rstrip("\n")
            else:
                owner.append(nm.split("|")[0])
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(seqs)).encode(),
                    sequence=line.rstrip("\n")).digitize(alpha))
    return seqs, owner, {}


def load_ice():
    """ORF-call the ICE shard with the frozen caller and keep host taxonomy."""
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    orf_caller.assert_anon_only(orf_caller.ANON)
    meta = {r["element_id"]: r for r in
            csv.DictReader(open(os.path.join(ICE, "subset", "ice_elements.tsv")),
                           delimiter="\t")}
    seqs, owner = [], []

    def emit(name, buf):
        if not (name and buf):
            return
        for g in orf_caller.call("".join(buf).upper()):
            p = (g["aa"] if isinstance(g, dict) else g).rstrip("*")
            if len(p) < 30:
                continue
            owner.append(name)
            seqs.append(pyhmmer.easel.TextSequence(
                name=str(len(seqs)).encode(), sequence=p).digitize(alpha))

    name, buf = None, []
    with gzip.open(os.path.join(ICE, "seq", "ice.shard000.fna.gz"), "rt") as fh:
        for line in fh:
            if line[0] == ">":
                emit(name, buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    emit(name, buf)
    return seqs, owner, meta


def main():
    require_compute_node()
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", choices=["plsdb", "ice"], required=True)
    a = ap.parse_args()
    import pyhmmer
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    disc, shared = profile_sets()
    print("=== discriminating profile sets ===")
    print("  shared, dropped from every class: %s" % ", ".join(shared))
    for c in CLASSES:
        print("  MPF_%-5s %2d" % (c, len(disc[c])))

    seqs, owner, meta = (load_plsdb if a.dataset == "plsdb" else load_ice)()
    elements = sorted(set(owner))
    NDB = len(elements)
    print("\n%s: %d proteins over %d elements\n" % (a.dataset, len(seqs), NDB), flush=True)
    if NDB == 0:
        print("no elements loaded -- refusing to report"); return 1

    # ---- scan every discriminating profile once ----------------------------
    count = collections.defaultdict(collections.Counter)   # element -> class -> n
    allp = sorted({p for c in CLASSES for p in disc[c]})
    missing = [p for p in allp if not os.path.exists(os.path.join(CJ, p + ".hmm"))]
    if missing:
        print("profiles missing from disk: %s" % missing)
        print("A silent skip reports a class as absent. Stopping.")
        return 1
    cls_of = {p: c for c in CLASSES for p in disc[c]}
    for i, prof in enumerate(allp):
        with pyhmmer.plan7.HMMFile(os.path.join(CJ, prof + ".hmm")) as fh:
            hmm = next(iter(fh))
        hits = set()
        for top in pyhmmer.hmmsearch([hmm], seqs, cpus=16, bit_cutoffs="gathering"):
            for h in top:
                hits.add(owner[int(dec(h.name))])
        for e in hits:
            count[e][cls_of[prof]] += 1
        if (i + 1) % 10 == 0:
            print("  scanned %d/%d profiles" % (i + 1, len(allp)), flush=True)

    # ---- baseline, PLSDB only ----------------------------------------------
    if a.dataset == "plsdb":
        fata_called = {r["accession"] for r in
                       csv.DictReader(open(FACAT), delimiter="\t")
                       if r["fa_subclass"] == "FATA"}
        print("\n=== BASELINE: re-admitting the %d plasmids already called FATA ==="
              % len(fata_called))
        for n in (1, 2, 3):
            got = sum(1 for e in fata_called if count[e]["FATA"] >= n)
            print("  FATA count >= %d : %5d / %d  (%.1f%%)"
                  % (n, got, len(fata_called), 100.0 * got / len(fata_called)))
        rec2 = sum(1 for e in fata_called if count[e]["FATA"] >= 2) / len(fata_called)
        if rec2 < 0.95:
            print("\n  BASELINE FAILED: %.1f%% recall at >=2 on a set already called"
                  " FATA by an independent tally." % (100 * rec2))
            print("  STOPPING -- the threshold table below would not be readable.")
            return 1
        print("  BASELINE PASSED (%.1f%% at >=2).\n" % (100 * rec2))

    # ---- admission by threshold --------------------------------------------
    admitted = {}
    for name, path in CATS:
        admitted[name] = {l.split("\t")[0] for i, l in enumerate(open(path)) if i}
    prev = set().union(*admitted.values()) if a.dataset == "plsdb" else set()

    print("=== admission by threshold, %s ===" % a.dataset)
    print("  %-6s %2s %9s %8s %10s %s"
          % ("class", "N", "admitted", "%db", "contam%",
             "  ".join("%5s" % ("+" + c) for c in CLASSES)))
    rows = []
    for c in CLASSES:
        for n in range(1, NMAX + 1):
            sel = [e for e in elements if count[e][c] >= n]
            if not sel:
                continue
            contam = (100.0 * len(set(sel) & prev) / len(sel)) if prev else float("nan")
            # how many of these ALSO clear the same bar in another class
            co = {o: sum(1 for e in sel if o != c and count[e][o] >= n) for o in CLASSES}
            print("  MPF_%-3s %2d %9d %7.2f%% %9s  %s"
                  % (c, n, len(sel), 100.0 * len(sel) / NDB,
                     ("%.1f%%" % contam) if prev else "n/a",
                     "  ".join("%5d" % co[o] for o in CLASSES)), flush=True)
            r = {"dataset": a.dataset, "class": "MPF_" + c, "threshold": n,
                 "admitted": len(sel), "pct_db": round(100.0 * len(sel) / NDB, 3),
                 "contamination_pct": (round(contam, 1) if prev else "")}
            r.update({"also_MPF_" + o: co[o] for o in CLASSES})
            rows.append(r)
        print()

    # ---- pre-registered taxonomy check, ICE only ---------------------------
    if a.dataset == "ice" and meta:
        print("=== PRE-REGISTERED host-phylum check (recorded before the run) ===")
        print("  profiles carry no taxonomic information; phylum comes from ICEberg.")
        for c in NEW:
            sel = [e for e in elements if count[e][c] >= 2]
            if not sel:
                print("  MPF_%-5s expected %-16s no elements at >=2" % (c, PREDICT[c]))
                continue
            ph = collections.Counter(meta.get(e, {}).get("phylum", "?") for e in sel)
            top, ntop = ph.most_common(1)[0]
            # score the TAXON, summing every synonym/GTDB split of it
            hit = 100.0 * sum(k for p, k in ph.items() if p in PREDICT[c]) / len(sel)
            print("  MPF_%-5s expected %-16s n=%4d  top=%-18s in-taxon %.1f%%  %s"
                  % (c, PREDICT[c][0], len(sel), top, hit,
                     "AS PREDICTED" if top in PREDICT[c] else "NOT as predicted"))
            for p, k in ph.most_common(4):
                print("        %-24s %4d  %5.1f%%" % (p, k, 100.0 * k / len(sel)))

    dest = os.path.join(OUT, "entry_gbc_fata_%s.tsv" % a.dataset)
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("\nwrote %s (%d rows)" % (dest, len(rows)))

    per = os.path.join(OUT, "profile_counts_%s.tsv" % a.dataset)
    with open(per, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["element"] + ["n_MPF_" + c for c in CLASSES] +
                   (["phylum"] if a.dataset == "ice" else []))
        for e in elements:
            if any(count[e][c] for c in CLASSES):
                row = [e] + [count[e][c] for c in CLASSES]
                if a.dataset == "ice":
                    row.append(meta.get(e, {}).get("phylum", ""))
                w.writerow(row)
    print("wrote %s" % per)
    return 0


if __name__ == "__main__":
    sys.exit(main())
