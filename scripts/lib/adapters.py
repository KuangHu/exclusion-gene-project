#!/usr/bin/env python3
"""A2 -- every caller/parser adapter must reproduce a known CDS before use.

The bug this prevents already happened: orfipy's BED is 0-based half-open AND
excludes the stop codon, while GenBank CDS coordinates include it. Verified on
ColE1 mob2 -- orfipy 1842..2187 (len 345) vs GenBank 1843..2190 (348 nt). The
off-by-3 made orfipy score 0/10 on a recall benchmark, and the first run was
measuring that bug rather than the caller.

Rule: an adapter is USABLE only if it is registered here with a golden test that
reproduces one known CDS at exact coordinates. `verify_all()` fails on any
adapter that lacks a golden test as well as on any that fails one -- an
unregistered adapter is not "untested", it is unusable.

This matters more going forward, not less: extending to MPF_I and MPF_F will add
adapters, and each new one is another chance to repeat the orfipy error.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
GB = os.path.join(PROJ, "data", "seed", "genbank")


class AdapterUnverified(RuntimeError):
    pass


_REGISTRY = {}


def register(name, fn, golden, note=""):
    """golden = (record_filename, start, end, strand, aa_len) -- 1-based inclusive,
    GenBank CDS convention, stop codon INCLUDED."""
    if golden is None:
        raise AdapterUnverified(
            "adapter %r registered without a golden test. An adapter with no "
            "golden test is unusable, not merely untested -- this is the check "
            "that catches coordinate-convention errors like orfipy's "
            "stop-codon-excluding BED." % name)
    _REGISTRY[name] = {"fn": fn, "golden": golden, "note": note}


def verify(name):
    """Run one adapter's golden test. Returns (ok, detail)."""
    from Bio import SeqIO
    a = _REGISTRY[name]
    fn, (recfile, s, e, st, aa) = a["fn"], a["golden"]
    path = os.path.join(GB, recfile)
    if not os.path.exists(path):
        return False, "golden record missing: %s" % recfile
    rec = next(SeqIO.parse(path, "genbank"))
    try:
        calls = fn(str(rec.seq))
    except Exception as exc:
        return False, "adapter raised: %s" % exc
    hit = [c for c in calls
           if c.get("start") == s and c.get("end") == e and c.get("strand") == st]
    if not hit:
        near = sorted(calls, key=lambda c: abs(c.get("start", 0) - s))[:1]
        n = near[0] if near else {}
        return False, ("no call at %d..%d(%s); nearest %s..%s(%s) -- a constant "
                       "offset here means a coordinate-convention mismatch"
                       % (s, e, "+" if st == 1 else "-", n.get("start"),
                          n.get("end"), n.get("strand")))
    if hit[0].get("aa_len") != aa:
        return False, "coords match but aa_len %s != expected %s" % (
            hit[0].get("aa_len"), aa)
    return True, "%d..%d exact, %d aa" % (s, e, aa)


def verify_all(verbose=True):
    if not _REGISTRY:
        raise AdapterUnverified("no adapters registered")
    bad = []
    for name in sorted(_REGISTRY):
        ok, detail = verify(name)
        if verbose:
            print("  %-20s %-6s %s" % (name, "OK" if ok else "FAIL", detail))
        if not ok:
            bad.append((name, detail))
    if bad:
        raise AdapterUnverified(
            "%d adapter(s) failed their golden test: %s"
            % (len(bad), "; ".join("%s (%s)" % b for b in bad)))
    return True


# ---------------------------------------------------------------- registrations
def _pyrodigal_anon(seq):
    import orf_caller
    return orf_caller.call(seq)


register(
    "pyrodigal_anon", _pyrodigal_anon,
    # RP4 trbK: the shortest slot gene in the set, 69 aa, and the one whose
    # exact start the whole slot hypothesis rests on
    golden=("RP4__BN000925.1.gb", 27450, 27659, 1, 69),
    note="THE production caller. docs/orf_caller_decision.md",
)


if __name__ == "__main__":
    print("A2 -- adapter golden tests\n")
    try:
        verify_all()
        print("\nPASS -- %d adapter(s) verified" % len(_REGISTRY))
    except AdapterUnverified as e:
        print("\nFAIL: %s" % e)
        sys.exit(1)

    print("\n=== A2 self-test: an adapter with a coordinate offset must fail ===")
    def _offset_by_three(seq):
        # simulates exactly the orfipy bug: stop codon excluded from the 3' end
        return [{**c, "end": c["end"] - 3, "aa_len": c["aa_len"]}
                for c in _pyrodigal_anon(seq)]
    register("buggy_offset3", _offset_by_three,
             golden=("RP4__BN000925.1.gb", 27450, 27659, 1, 69))
    ok, detail = verify("buggy_offset3")
    print("  buggy_offset3        %-6s %s" % ("OK" if ok else "CAUGHT", detail))
    if ok:
        print("\n  A2 IS NOT WORKING -- it passed a known-broken adapter")
        sys.exit(1)

    print("\n=== and registering without a golden test must fail ===")
    try:
        register("no_golden", _pyrodigal_anon, golden=None)
        print("  NOT CAUGHT -- A2 is not working")
        sys.exit(1)
    except AdapterUnverified as e:
        print("  CAUGHT: %s" % str(e).split(".")[0])
