"""Hard-fail guards against the silent-failure pattern.

Four silent failures across two levels, every one of the form "it looked like it
ran, the result was empty":

  L1  an undefined gene was skipped, reported as "no-baseline"
  L2  a string-replace edit did not match; provenance columns silently absent
  L2  name_regex was tested only against /product, never /gene
  L2  a config edit wrote a note into restrict_inc_group, silently restricting
      a gene to a nonexistent Inc group

None were caught by a test. All three guards below turn the empty result into an
exception at the point it happens.
"""

import os
import sys


class AssertionFailure(RuntimeError): pass


# ---------------------------------------------------------------- guard 7
# A7 is the most valuable of these, because the bug it prevents produced a FALSE
# SIGNAL rather than a missing one.
#
# csv.writer defaults to \r\n. The last field of every row therefore held a
# stray "\r", so an EMPTY anchor field was truthy and 146/146 plasmids were
# scored as "has an anchor". A1-A6 all guard against "matched nothing, returned
# empty" -- an absence, which at least has a chance of being noticed. This one
# manufactured a 100% occupancy rate.
#
# 100% is exactly the value at which pipeline v2 6.1 says to ABANDON THE
# DIRECTION ("near 100% -> slot occupancy is not a variable -> change the
# criterion"). A line-terminator character came within one smoke test of killing
# the project by way of its own decision rule.

_WHITESPACE_ONLY = {"", "\r", "\n", "\r\n", "\t", " "}


def is_blank(value):
    """A7: the ONLY sanctioned emptiness test for a parsed field.

    Never use bare truthiness on a field read from a file. "\r" is truthy.
    """
    if value is None:
        return True
    if isinstance(value, str):
        return value.strip() == ""
    return False


def assert_clean_field(value, where=""):
    """A7: a field carrying stray line-terminator characters is a hard failure.

    Catches the writer bug at the READ side too, so a bad file produced by some
    other tool cannot slip through.
    """
    if isinstance(value, str) and value != value.strip() and value.strip() != "":
        raise AssertionFailure(
            "%s: field %r carries leading/trailing whitespace or line-terminator "
            "characters. This is how an empty field became truthy and produced a "
            "100%% occupancy rate. Strip on read, and write with "
            "lineterminator='\\n'." % (where, value))
    if isinstance(value, str) and value in _WHITESPACE_ONLY and value != "":
        raise AssertionFailure(
            "%s: field is whitespace-only (%r), not empty. Under bare truthiness "
            "this reads as PRESENT." % (where, value))
    return value


# ---------------------------------------------------------------- guard 8
# A8 exists because the same mistake was made TWICE:
#   146 -> "~8,900"   : a shard taken in file order, treated as a random sample
#   sorted()[:800]    : 694 CP + 95 AP + ZERO NZ_, against a set that is 80% NZ_
# Both produced a plausible number from an unrepresentative slice. A7 blocks a
# MANUFACTURED PRESENCE; A8 blocks MANUFACTURED REPRESENTATIVENESS. Both are
# false signals rather than missing ones, and both propagate straight into a
# decision -- the 80.7% landed inside the 60-90% band that decides the project's
# direction.

SAMPLING_METHODS = {"random", "hash", "contiguous", "full"}


def declare_sampling(method, n_taken, n_total, produces_statistic=True, where=""):
    """A8: any --limit / truncation / subsample must declare how it sampled.

    `contiguous` is a hard failure whenever the script produces a statistic.
    It remains legal for debugging, where produces_statistic=False.
    """
    if method not in SAMPLING_METHODS:
        raise AssertionFailure(
            "%s: sampling method %r not declared. One of %s is required."
            % (where, method, sorted(SAMPLING_METHODS)))
    if method == "contiguous" and produces_statistic:
        raise AssertionFailure(
            "%s: CONTIGUOUS sampling (%d of %d) may not produce a statistic. "
            "PLSDB is ordered by submission, so a contiguous block is a "
            "stratum, not a sample -- sorted()[:800] returned 694 CP + 95 AP "
            "and zero NZ_ from a set that is 80%% NZ_. Use hash or random "
            "sampling, or pass produces_statistic=False for a debug run."
            % (where, n_taken, n_total))
    if method != "full" and n_total and n_taken / n_total < 0.5:
        sys.stderr.write(
            "[A8] %s: %s sample of %d/%d (%.1f%%) -- report the sampling "
            "method alongside any figure derived from it.\n"
            % (where, method, n_taken, n_total, 100.0 * n_taken / n_total))
    return method


def assert_no_cr(path):
    """A7: refuse to read a TSV/CSV containing carriage returns."""
    with open(path, "rb") as fh:
        head = fh.read(1 << 20)
    if b"\r" in head:
        raise AssertionFailure(
            "%s contains carriage returns. Rewrite it with "
            "lineterminator='\\n'; a CR in the final field makes an empty field "
            "truthy." % path)
    return path

# ---------------------------------------------------------------- guard 1
ENUMS={
 "side":{"exclusion","partner"},
 "role_exclusion":{"positive","negative","ambiguous","candidate"},
 "role_adjacency":{"positive","negative","unmeasured"},
 "record_scope":{"chromosome","plasmid","ICE","fragment"},
 "status":{"PINNED","PROVISIONAL","ASSEMBLE","CHROMOSOMAL","UNRESOLVED"},
 "evidence_tier":{"E1","E2","E3","E4"},
 "hmm_hit":{"HIT","no_hit","no_seq","unknown"},
 "partner_status":{"verified","predicted","none","self"},
 "type":{"entry","surface","both","undetermined"},
}
INC_GROUPS={"IncFI","IncFII","IncFV","IncI1","IncI-gamma","IncIgamma","IncN",
 "IncP-alpha","IncW","IncHI1","IncHI2","IncC","IncA","IncQ","IncP/IncN-like",
 "IncC (also IncA)","not_applicable","NA",""}

def check_enum(field, value, where=""):
    """Guard 1: a config value outside its enum is a hard failure, not a warning."""
    if field not in ENUMS: return value
    if value in (None,""): return value
    if value not in ENUMS[field]:
        raise AssertionFailure(
            f"{where}: field '{field}' has value {value!r} which is not in its "
            f"enum {sorted(ENUMS[field])}. This is how a note written into the "
            f"wrong column silently disabled a filter.")
    return value

def check_inc_group(value, where=""):
    for v in str(value).split(","):
        if v.strip() and v.strip() not in INC_GROUPS:
            raise AssertionFailure(
                f"{where}: inc_group {v.strip()!r} unknown. If this is a real new "
                f"Inc group, add it to INC_GROUPS; if it is prose, it is in the "
                f"wrong column.")
    return value

# ---------------------------------------------------------------- guard 2
def require_nonempty(result, what, where=""):
    """Guard 2: a query/regex that matches nothing fails loudly."""
    if result is None or (hasattr(result,"__len__") and len(result)==0):
        raise AssertionFailure(
            f"{where}: {what} matched 0 records. An empty result is a failure "
            f"here, not a valid answer -- pass allow_empty=True at the call site "
            f"if emptiness is genuinely expected.")
    return result

# ---------------------------------------------------------------- guard 3
def replace_or_fail(text, old, new, where="", count=None):
    """Guard 3: a str.replace that changes nothing fails loudly."""
    n=text.count(old)
    if n==0:
        raise AssertionFailure(
            f"{where}: replacement target not found, 0 substitutions made.\n"
            f"  looked for: {old[:110]!r}")
    if count is not None and n!=count:
        raise AssertionFailure(f"{where}: expected {count} occurrences, found {n}")
    return text.replace(old,new)

# ---------------------------------------------------------------- guard 9
# A9: a curated file and a generated file must never share a path.
#
# 80_determine_anchors.py read Pfam accessions off the control plasmids and wrote
# them to data/anchors/anchor_set.tsv -- which was the hand-curated anchor
# definition, including the VirB5 rule (E<=1e-5 AND aa>150) that GA cannot
# express and the VirB7/PF20898 correction. The curation was gone, and nothing
# failed: the file still parsed, still had the right columns, and the pipeline
# ran on a machine read-off of three plasmids. Silent, and upstream of every
# occupancy number.
#
# Two mechanisms, because either alone leaks:
#   sealed hashes   -- detects the clobber before the next stage consumes it
#   write guard     -- a script that opens a curated path for writing dies there
import hashlib as _hashlib

CURATED_REGISTRY = "config/curated_files.tsv"


def _sha256(path):
    h = _hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def _registry(proj):
    reg = os.path.join(proj, CURATED_REGISTRY)
    out = {}
    if not os.path.exists(reg):
        return out
    for line in open(reg):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if f[0] == "path":
            continue
        out[f[0]] = (f[1], f[2] if len(f) > 2 else "")
    return out


def assert_curated_unchanged(proj):
    """Every registered curated file must still hash to its sealed value."""
    problems = []
    for rel, (want, note) in sorted(_registry(proj).items()):
        p = os.path.join(proj, rel)
        if not os.path.exists(p):
            problems.append(f"A9: curated file {rel} is MISSING (sealed {want[:12]}). "
                            f"{note}")
            continue
        got = _sha256(p)
        if got != want:
            problems.append(
                f"A9 VIOLATION: {rel} no longer matches its sealed hash.\n"
                f"      sealed {want}\n      actual {got}\n"
                f"    This is the anchor_set.tsv clobber signature: a generator "
                f"wrote over curation and nothing else would have noticed. If the "
                f"edit was deliberate, re-seal with:\n"
                f"      python3 scripts/lib/assertions.py --reseal {rel}\n"
                f"    {note}")
    return problems


def curated_guard(proj, path):
    """Call before opening `path` for writing. Refuses if it is curated."""
    rel = os.path.relpath(os.path.abspath(path), proj)
    if rel in _registry(proj):
        raise AssertionFailure(
            f"A9: refusing to write {rel} -- it is a CURATED file, not a "
            f"generated one. Write to a *_RAW_readoff.tsv sibling and merge by "
            f"hand, which is what 80_determine_anchors.py now does.")
    return path


def assert_no_generated_collision(proj):
    """Static check: no script may name a curated path on an open-for-write line."""
    problems = []
    curated = _registry(proj)
    bases = {os.path.basename(r): r for r in curated}
    for root, _dirs, files in os.walk(os.path.join(proj, "scripts")):
        for fn in files:
            if not fn.endswith(".py"):
                continue
            p = os.path.join(root, fn)
            for i, line in enumerate(open(p, errors="replace"), 1):
                if '"w"' not in line and "'w'" not in line:
                    continue
                for b, rel in bases.items():
                    if b in line and "curated_guard" not in line:
                        problems.append(
                            f"A9: {os.path.relpath(p, proj)}:{i} opens curated "
                            f"file {rel} for writing.")
    return problems


# ---------------------------------------------------------------------------
# A14: a file must have CONTENT before it can be sealed.
#
# A9 answers "has this file changed since sealing". It has no notion of whether
# the file ever had content. A failed heredoc left a 0-BYTE
# docs/SUMMARY_unnamed_stratification.md and _reseal hashed it happily -- the act
# of sealing then made it look verified.
#
# Same shape as the VirB5 threshold: sealed with a rationale that had never been
# reproducible. A seal is an integrity guarantee, never a correctness one, and it
# cannot even guarantee non-emptiness unless asked to.
#
# This cannot check that content is MEANINGFUL. It closes the narrower class:
# "a process failed and left a shell behind".
# ---------------------------------------------------------------------------
_A14_MIN_BYTES = 32
_A14_MIN_LINES = 2


def assert_sealable(proj, rels):
    """A14: refuse to seal empty or structurally broken files."""
    bad = []
    for rel in rels:
        p = os.path.join(proj, rel)
        if not os.path.exists(p):
            bad.append((rel, "does not exist")); continue
        n = os.path.getsize(p)
        if n < _A14_MIN_BYTES:
            bad.append((rel, "only %d bytes" % n)); continue
        with open(p, "rb") as fh:
            lines = fh.read().decode("utf-8", "replace").splitlines()
        body = [l for l in lines if l.strip()]
        if len(body) < _A14_MIN_LINES:
            bad.append((rel, "only %d non-blank lines" % len(body))); continue
        if rel.endswith(".md") and not any(l.lstrip().startswith("#") for l in body):
            bad.append((rel, "markdown with no heading"))
        elif rel.endswith(".tsv"):
            hdr = next((l for l in body if not l.startswith("#")), "")
            data = [l for l in body if not l.startswith("#")][1:]
            if "\t" not in hdr:
                bad.append((rel, "tsv header has no tab"))
            elif not data:
                bad.append((rel, "tsv has a header but no data row"))
    if bad:
        raise AssertionError(
            "A14: %d file(s) are not sealable -- a seal on an empty or broken\n"
            "  file makes it LOOK verified. Fix the file, do not seal it.\n%s"
            % (len(bad), "\n".join("      %-56s %s" % b for b in bad)))
    return True


def _reseal(proj, rels):
    assert_sealable(proj, rels)          # A14 -- never seal a shell
    reg = os.path.join(proj, CURATED_REGISTRY)
    cur = _registry(proj)
    for rel in rels:
        cur[rel] = (_sha256(os.path.join(proj, rel)),
                    cur.get(rel, ("", ""))[1])
    with open(reg, "w", newline="") as fh:
        fh.write("# A9: curated files, sealed. A generated file must never take "
                 "one of these paths.\npath\tsha256\tnote\n")
        for rel in sorted(cur):
            fh.write("%s\t%s\t%s\n" % (rel, cur[rel][0], cur[rel][1]))
    return reg


if __name__ == "__main__":
    import sys as _sys
    _proj = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    if len(_sys.argv) > 2 and _sys.argv[1] == "--reseal":
        print("re-sealed: " + _reseal(_proj, _sys.argv[2:]))
    else:
        for _p in assert_curated_unchanged(_proj) + assert_no_generated_collision(_proj):
            print("  * " + _p)


# ---------------------------------------------------------------- guard 10
# A10: an accession that is cited must resolve to a model we actually hold.
#
# `TIGR04239` was written for `TIGR04359` (TrbK_RP4). It is well-formed, plausible,
# one digit off, and no existing guard covers it: A1-A9 check enums, emptiness,
# replacements, CR bytes, sampling, metric names and curated-file hashes -- none
# looks at accession STRINGS. This project has already catalogued six
# attribute-mismatch errors in the literature and its own notes; a silent digit
# swap in an accession is the same failure class and propagates just as quietly.
#
# Rule: every Pfam/TIGRFAM/NCBIfam accession appearing in docs/ or in the curated
# anchor set must be resolvable -- PF* against the local Pfam-A catalogue, NF*/
# TIGR* against the models actually present in the HMM directory. Citing a model
# we do not hold means the claim about it was never checked.
import re as _re

PFAM_DAT = "/global/scratch/users/kh36969/funcannot_dbs/pfam/Pfam-A.hmm.dat"
HMM_DIR = "/global/scratch/users/kh36969/exclusion_gene/hmm"
ACC_RE = _re.compile(r"\b(PF\d{5}|TIGR\d{5}|NF\d{6})\b")
A10_SCOPE = ("docs", os.path.join("data", "anchors"))


def _known_accessions():
    pf, other = set(), set()
    if os.path.exists(PFAM_DAT):
        for line in open(PFAM_DAT):
            if line.startswith("#=GF AC"):
                pf.add(line.split()[2].split(".")[0])
    if os.path.isdir(HMM_DIR):
        for fn in os.listdir(HMM_DIR):
            if not fn.endswith(".hmm"):
                continue
            # NB: every HMMER3 file opens with `HMMER3/f [...]`, so breaking on
            # a "HMM" prefix exits before ACC is ever seen. Read the header block
            # by line count instead; ACC sits within the first ~25 lines.
            with open(os.path.join(HMM_DIR, fn)) as fh:
                for n, line in enumerate(fh):
                    if line.startswith("ACC "):
                        other.add(line.split()[1].split(".")[0])
                    if n > 30:
                        break
    return pf, other


def assert_accessions_resolve(proj):
    pf, held = _known_accessions()
    if not pf and not held:
        return ["A10: no local model catalogue found; cannot verify accessions"]
    problems = []
    targets = []
    for sc in A10_SCOPE:
        p = os.path.join(proj, sc)
        if os.path.isdir(p):
            for root, _d, files in os.walk(p):
                targets += [os.path.join(root, f) for f in files
                            if f.endswith((".md", ".tsv"))]
        elif os.path.exists(p):
            targets.append(p)
    for path in targets:
        for i, line in enumerate(open(path, errors="replace"), 1):
            for acc in set(ACC_RE.findall(line)):
                ok = acc in pf if acc.startswith("PF") else acc in held
                if not ok:
                    problems.append(
                        "A10: %s:%d cites %s, which resolves to no model we hold. "
                        "A one-digit accession slip is well-formed and silent -- "
                        "TIGR04239 for TIGR04359 is the case this guard exists for. "
                        "Fetch the model or fix the accession."
                        % (os.path.relpath(path, proj), i, acc))
    return sorted(set(problems))


# ---------------------------------------------------------------- guard 11
# A11: heavy work must not run on a login node.
#
# Scripts 83 and 91-95 were submitted through sbatch. From script 96 onward the
# analyses were launched with `nohup python3` directly because it iterated faster,
# and nothing re-checked. The result was multi-hour 16-20 CPU jobs running on
# n0001.scs00 next to fifteen other users, the worst being a phmmer sweep at
# 2 days 4 hours and 5 GB RSS.
#
# A convention note would not have prevented this -- the same session had already
# recorded "a safeguard written as prose and not as code" as a failure mode, in
# the context of running recall without its precision control. So this is a check.
#
# Call require_compute_node() at the top of anything that scans the CDS cache or
# runs hmmsearch/phmmer over it.
def require_compute_node(allow_env="ALLOW_LOGIN_NODE"):
    """Refuse to proceed unless inside a SLURM allocation."""
    if os.environ.get("SLURM_JOB_ID"):
        return True
    if os.environ.get(allow_env):
        sys.stderr.write("A11: running outside SLURM because %s is set. "
                         "Only acceptable for a small smoke test.\n" % allow_env)
        return False
    raise AssertionFailure(
        "A11: no SLURM_JOB_ID -- this looks like a login node (%s).\n"
        "  Heavy scans must go through sbatch. Write a wrapper in slurm/ and "
        "submit it, or set %s=1 for a deliberately small smoke test.\n"
        "  This guard exists because scripts 96-105 were run with `nohup python3` "
        "on the login node for two sessions, including a 2-day phmmer sweep."
        % (os.uname().nodename, allow_env))


# ---------------------------------------------------------------- guard 12
# A12: a lookup key must exist in the database it is being looked up in.
#
# Third variant of the same failure. A9 hashes curated files, A10 checks that a
# cited accession resolves to a model we hold -- neither asks whether the
# accession used as a KEY is present in the target.
#
#   1. `accession_related` carries a different coordinate frame, so coordinates
#      validated against it were meaningless.
#   2. The MPF_F entry criterion reported 0/5 seeds for all 13 families. PLSDB is
#      RefSeq-based and the seed records are INSDC: AP001918.1, AP000342.1,
#      CP033514.1, AF250878.1, KJ817376.1 are none of them in PLSDB. Four have
#      RefSeq equivalents; SXT is an ICE and legitimately absent, so it had to be
#      injected.
#
# The next collision is close to certain rather than hypothetical: ICEBs1 is a
# chromosomal element absent from PLSDB entirely; pLS20/pCF10/pAD1 are Gram+
# elements of unknown PLSDB coverage; and R64/R621a seeds are AP005147.1 /
# AP011954.1, INSDC again. MPF_I and MPF_FA hit all three at once.
def assert_keys_present(keys, universe, what="lookup keys", where="",
                        allow_missing=()):
    """Fail loudly when a lookup key is absent from the set it indexes.

    `keys` may be a dict {label: accession} or an iterable of accessions.
    `allow_missing` names keys that are known-absent for a stated reason (an ICE
    in a plasmid database, say) and must be injected rather than silently dropped.
    """
    items = keys.items() if hasattr(keys, "items") else [(k, k) for k in keys]
    u = set(universe)
    missing = [(lab, k) for lab, k in items if k not in u and lab not in allow_missing]
    if missing:
        raise AssertionFailure(
            "A12: %d of %d %s are absent from the target set%s:\n%s\n"
            "  A key that is not in the database returns nothing, and 'nothing' "
            "reads as a biological result. The MPF_F entry criterion reported "
            "0/5 seeds this way while the scan itself was working perfectly.\n"
            "  Check for a RefSeq/INSDC mismatch, or inject the record if it is "
            "legitimately absent (an ICE in a plasmid database)."
            % (len(missing), len(items), what, (" in " + where) if where else "",
               "\n".join("      %-10s %s" % (lab, k) for lab, k in missing)))
    return True


# ---------------------------------------------------------------------------
# A13: an anchor threshold may never be an E-value.
#
# hmmsearch reports E = P * N, so an E-value threshold is a property of the
# protein AND the size of the database searched. Measured on this project's own
# data (job 25913252): PF07996 at `E<=1e-5 AND aa>150`, over IDENTICAL sequences,
# admitted 5,368 plasmids in a 1,099,911-protein search and 5,341 in a 7,725,419-
# protein one. Bitscore cutoffs (GA) are invariant -- VirB1 (5,984), VirB2 (5,884)
# and VirB3 (6,769) reproduced to the digit across both.
#
# This is not a theoretical concern. `anchor_set.tsv` carried VirB5 at
# `E<=1e-5 AND aa>150` into release v1. Its stated rationale was that the rule
# rescues IncP TrbJ, which GA loses. Under the pipeline's own search those
# proteins score E = 0.22 (RP4) and E = 0.0093 (R751): the rule FAILED both
# controls it was written for. The rescuing E-values had come from a
# protein-vs-Pfam scan (N ~ 30,134 models) rather than a model-vs-proteome
# search. The rule was also STRICTER than the GA it replaced (E=1e-5 falls at
# ~31.7 bits vs GA 24.7), costing ~177 plasmids while documented as a loosening.
#
# A9 sealed that file and could not catch it: a seal proves a file has not
# CHANGED, never that a number in it was ever derivable. A13 is the check A9
# cannot be.
#
# See docs/virb5_evalue_threshold_invalid.md.
# ---------------------------------------------------------------------------
A13_ANCHOR_FILES = (os.path.join("data", "anchors", "anchor_set.tsv"),
                    os.path.join("data", "anchors", "anchor_set_MPF_F.tsv"))
_A13_EVALUE = _re.compile(r"(?:^|[^A-Za-z])[Ee]\s*(?:<=|<|=)\s*[0-9]|[Ee]-?value", _re.I)


def assert_no_evalue_thresholds(proj, files=None):
    """A13: reject any E-value in an anchor set's `threshold` column.

    GA, or an explicit bitscore. Never an E-value: it is database-size dependent
    and therefore not reproducible across cache scopes.
    """
    bad = []
    for rel in (files or A13_ANCHOR_FILES):
        path = os.path.join(proj, rel)
        if not os.path.exists(path):
            continue
        with open(path) as fh:
            head, cols = None, None
            for ln in fh:
                if ln.startswith("#") or not ln.strip():
                    continue
                f = ln.rstrip("\n").split("\t")
                if head is None:
                    head = f
                    cols = [i for i, c in enumerate(f) if c.strip() == "threshold"]
                    continue
                for i in cols:
                    if i < len(f) and _A13_EVALUE.search(f[i]):
                        bad.append((rel, f[0], f[i]))
    if bad:
        raise AssertionError(
            "A13: %d anchor threshold(s) are E-value based. hmmsearch reports\n"
            "  E = P * N, so these are database-size dependent and do NOT\n"
            "  reproduce across cache scopes (measured: 5,368 vs 5,341 hits for\n"
            "  the same rule on the same sequences at 1.1M vs 7.7M proteins).\n"
            "  Use GA, or an explicit bitscore derived from measurement.\n%s\n"
            "  See docs/virb5_evalue_threshold_invalid.md"
            % (len(bad), "\n".join("      %s  %s -> %r" % b for b in bad)))
    return True
