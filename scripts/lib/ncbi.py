"""E-utilities access for the exclusion-gene seed set.

Design notes
------------
* Dual accession by default. For these historically-studied elements the ORIGINAL
  GenBank submission is the better annotation source and RefSeq/PGAP is the
  cross-check -- the reverse of the usual default. Demonstrated on R100: PGAP
  (NC_002134.1, re-annotated 2025-06-04) dropped the /gene="traS" qualifier that
  the original submission (AP000342.1) carries. Coordinates and protein sequence
  are identical; only the annotation degraded.
* RefSeq -> GenBank backlink lives in the COMMENT block, not DBSOURCE/DBLINK:
      "REFSEQ INFORMATION: The reference sequence is identical to AP000342.1."
* Protein identity: GenBank protein_id (e.g. BAA78881.1) is ELEMENT-SPECIFIC.
  RefSeq WP_ accessions are non-redundant across species and identify a SEQUENCE,
  not a locus -- never use a WP_ alone to say which element a protein came from.
"""
import hashlib, re, time, urllib.parse, urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_LAST = [0.0]
_MIN_INTERVAL = 0.36          # <3 req/s, the un-keyed E-utilities limit

def _throttle():
    dt = time.time() - _LAST[0]
    if dt < _MIN_INTERVAL:
        time.sleep(_MIN_INTERVAL - dt)
    _LAST[0] = time.time()

def _get(endpoint, params, tries=4, timeout=120):
    url = f"{EUTILS}/{endpoint}?" + urllib.parse.urlencode(params)
    last = None
    for attempt in range(tries):
        _throttle()
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:                      # transient 429/5xx/timeout
            last = e
            time.sleep(2 ** attempt)
    raise RuntimeError(f"E-utilities failed after {tries} tries: {endpoint} {params}: {last}")

def esearch(term, db="nuccore", retmax=20):
    """Return list of UIDs. Caller must still verify -- a hit is a candidate."""
    xml = _get("esearch.fcgi", {"db": db, "term": term, "retmax": retmax})
    return re.findall(r"<Id>(\d+)</Id>", xml)

def esummary(uids, db="nuccore"):
    """Return [{uid, accver, title, slen, organism}] for candidate review."""
    if not uids:
        return []
    xml = _get("esummary.fcgi", {"db": db, "id": ",".join(uids), "version": "2.0"})
    out = []
    for blk in re.findall(r"<DocumentSummary[ >].*?</DocumentSummary>", xml, re.S):
        def field(tag):
            m = re.search(rf"<{tag}>(.*?)</{tag}>", blk, re.S)
            return m.group(1).strip() if m else ""
        out.append({
            "source":   "ncbi",
            "uid":      re.search(r'uid="(\d+)"', blk).group(1) if re.search(r'uid="(\d+)"', blk) else "",
            "accver":   field("AccessionVersion"),
            "title":    field("Title"),
            "slen":     field("Slen"),
            "organism": field("Organism"),
        })
    return out

def efetch_gb(accession, rettype="gbwithparts"):
    """Fetch a GenBank flatfile. Accession SHOULD carry .version."""
    return _get("efetch.fcgi", {"db": "nuccore", "id": accession,
                                "rettype": rettype, "retmode": "text"})

def refseq_to_genbank(gb_text):
    """Extract the original GenBank accession.version a RefSeq record derives from.

    The backlink is prose inside COMMENT and wraps across lines, so normalise
    whitespace before matching. Returns None for primary (non-RefSeq) records.
    """
    m = re.search(r"^COMMENT\s+(.*?)(?=^\S)", gb_text, re.S | re.M)
    if not m:
        return None
    comment = " ".join(m.group(1).split())
    m2 = re.search(r"reference sequence is identical to ([A-Z]{1,2}\d{5,8}\.\d+)", comment)
    if m2:
        return m2.group(1)
    m3 = re.search(r"derived from ([A-Z]{1,2}\d{5,8}\.\d+)", comment)
    return m3.group(1) if m3 else None

def aa_sha256(seq):
    """Freeze-value for a protein. Uppercased, stripped -- the regression anchor."""
    return hashlib.sha256(seq.strip().upper().encode()).hexdigest()

def record_md5(text):
    return hashlib.md5(text.encode()).hexdigest()
