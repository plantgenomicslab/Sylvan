#!/usr/bin/env python3
"""Bound miniprot evidence depth before EVM consumes it.

EVM 2.1.0 accumulates evidence weight with no normalisation:

    add_match_coverage   $CODING_SCORES[$i] = $current + $weight   (per base)
    add_introns          $INTRONS_TO_SCORE{$key} += $intron_score  (per intron)

So N alignments over one base contribute N x weight. With the full OrthoDB
Viridiplantae set as protein evidence a 1 Mb Ath_ODB partition reaches 297x
mean depth over covered bases and 2,316x at the peak, against 1,075 GeneWise
rows -- database redundancy becomes vote count, and a weight scalar cannot
correct it.

This module caps both accumulation paths:

    per genomic base                 retained coding segments <= k_base
    per exact intron key             supporting alignments    <= k_intron

Selection is atomic per alignment ID. `ID` is one genomic hit; `Target` is the
query protein and may have several hits at different loci, so capping by Target
would remove a protein's evidence everywhere at once. Identity and score are
per-segment values -- they vary between segments of the same alignment in 99%
of multi-segment alignments -- so both the ranking and the identity floor use
length-weighted alignment-level aggregates. Applying a floor per row would
amputate weak exons and emit chimeric partial alignments.

Input is NOT coordinate-sorted (real files jump Chr5:13203494 -> 13201546 ->
Chr2:...), so callers streaming huge files should sort by ID first; this module
groups by ID itself and does not rely on input order.

Design: docs/superpowers/specs/2026-08-25-evm-selftune-design.md
"""
import argparse
import math
import os
import shlex
import subprocess
import tempfile
import sys
import re
from collections import Counter, defaultdict


class Alignment:
    """One genomic hit: every segment sharing an ID."""

    __slots__ = ("aln_id", "seqid", "strand", "segments", "introns",
                 "identity_w", "aligned_aa", "score_total", "rank",
                 "family")

    def __init__(self, aln_id, seqid, strand, segments, introns, identity_w,
                 aligned_aa, score_total, rank, family=None):
        self.aln_id = aln_id
        self.seqid = seqid
        self.strand = strand
        self.segments = segments
        self.introns = introns
        self.identity_w = identity_w
        self.aligned_aa = aligned_aa
        self.score_total = score_total
        self.rank = rank
        self.family = family

    def priority(self):
        """Ordering key -- lower sorts first (i.e. is admitted first).

        Length-weighted identity leads: raw score scales with alignment length,
        so a long mediocre hit would otherwise evict a short excellent one.
        aligned_aa breaks ties so a short conserved domain does not monopolise
        capacity. Rank is only meaningful within one Target. The ID tail makes
        the result independent of input order.
        """
        return (-self.identity_w, -self.aligned_aa, self.rank,
                -self.score_total, self.aln_id)


def _attr(attrs, key):
    for field in attrs.split(";"):
        field = field.strip()
        if field.startswith(key + "="):
            return field[len(key) + 1:]
    return None


def _target_span(attrs):
    """Query interval from `Target=<id> <start> <end>`, normalised.

    Ranking must count query residues: summing genomic segment lengths instead
    lets long introns inflate an alignment past a compact one that actually
    aligned more of the protein.
    """
    raw = _attr(attrs, "Target")
    if not raw:
        return None
    parts = raw.split()
    if len(parts) < 3:
        return None
    try:
        t_start, t_end = int(parts[1]), int(parts[2])
    except ValueError:
        return None
    return (t_start, t_end) if t_start <= t_end else (t_end, t_start)


def family_key(target, aln_id):
    """Conservatively group obvious isoforms of one Target protein gene.

    OrthoDB-derived Targets such as ``AL1006U10010.t1`` encode the isoform in
    a terminal ``.tN`` suffix.  We also recognize the equally explicit
    ``.pN`` and ``-R[A-Z]+`` conventions.  Anything else is not safely
    interpretable as a family identifier, so it receives a singleton key
    based on the alignment ID.  Thus enabling the cap cannot accidentally
    merge unrelated proteins merely because they share a broad text prefix.
    """
    if target:
        base = re.sub(r"(?:\.(?:t|p)\d+|-R[A-Z]+)$", "", target,
                      flags=re.IGNORECASE)
        if base != target and base:
            return "isoform:" + base
    return "singleton:" + aln_id


def _parse(line, lineno=None):
    """Parse one segment. None for comments/blanks; ValueError for garbage.

    Passing an unparseable row through would hand EVM evidence that never went
    through the cap, so malformed input fails fast rather than being copied to
    the output.
    """
    if not line or line.startswith("#"):
        return None
    where = f"line {lineno}: " if lineno else ""
    f = line.rstrip("\n").split("\t")
    if len(f) < 9:
        raise ValueError(f"{where}expected 9 tab-separated columns, got {len(f)}")
    aln_id = _attr(f[8], "ID")
    if not aln_id:
        raise ValueError(f"{where}no ID attribute")
    try:
        start, end = int(f[3]), int(f[4])
    except ValueError:
        raise ValueError(f"{where}non-integer coordinates {f[3]!r} {f[4]!r}")
    if start < 1 or end < 1:
        raise ValueError(f"{where}coordinates must be >= 1")
    if start > end:
        start, end = end, start
    if f[6] not in ("+", "-"):
        raise ValueError(f"{where}strand must be + or -, got {f[6]!r}")
    try:
        score = float(f[5])
    except ValueError:
        score = 0.0
    if not math.isfinite(score):
        raise ValueError(f"{where}score is not finite")
    raw_identity = _attr(f[8], "Identity")
    try:
        identity = float(raw_identity)
    except (TypeError, ValueError):
        identity = 0.0
    if not math.isfinite(identity):
        raise ValueError(f"{where}Identity is not finite: {raw_identity!r}")
    try:
        rank = int(_attr(f[8], "Rank"))
    except (TypeError, ValueError):
        rank = 1
    target = _attr(f[8], "Target")
    query = target.split()[0] if target else None
    return (aln_id, f[0], f[6], start, end, identity, score, rank,
            _target_span(f[8]), query)


def _union_length(intervals):
    """Total length covered by a set of inclusive intervals."""
    total, cur_s, cur_e = 0, None, None
    for s, e in sorted(intervals):
        if cur_e is None or s > cur_e + 1:
            if cur_e is not None:
                total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        total += cur_e - cur_s + 1
    return total


def _build_alignment(aln_id, meta, rows):
    """Turn one ID's segments into an Alignment, or raise if they disagree."""
    seqid, strand, rank, query = meta
    rows.sort()
    spans = [(s, e) for s, e, _, _, _ in rows]

    # Query residues, not genomic span. All-or-nothing: mixing rows that carry
    # a Target with rows that do not would silently drop the latter from the
    # identity average.
    tspans = [t for _, _, _, _, t in rows if t]
    if tspans and len(tspans) != len(rows):
        raise ValueError(
            f"alignment {aln_id}: Target interval present on "
            f"{len(tspans)} of {len(rows)} segments -- refusing to mix "
            "query-coordinate and genomic weighting")
    if tspans:
        aligned = _union_length(tspans)
        weight_of = [(t[1] - t[0] + 1) for t in tspans]
    else:
        aligned = _union_length(spans)
        weight_of = [e - s + 1 for s, e in spans]
    idents = [idn for _, _, idn, _, _ in rows]
    denom = sum(weight_of)
    weighted = sum(idn * w for idn, w in zip(idents, weight_of))

    introns = [(e1 + 1, s2 - 1) for (_, e1), (s2, _) in zip(spans, spans[1:])
               if s2 - 1 >= e1 + 1]
    return Alignment(
        aln_id=aln_id, seqid=seqid, strand=strand, segments=spans,
        introns=introns,
        identity_w=(weighted / denom) if denom else 0.0,
        aligned_aa=aligned,
        score_total=sum(sc for _, _, _, sc, _ in rows),
        rank=rank,
        family=family_key(query, aln_id),
    )


def _check_meta(aln_id, seen, incoming):
    """An ID is one genomic hit; disagreement means the file is not what we
    assume. `_target_span` drops the query name, so a mixed ID would union
    coordinates belonging to two different proteins."""
    if seen == incoming:
        return
    fields = ("seqid", "strand", "Rank", "Target")
    diff = [f"{name}={a!r} vs {b!r}"
            for name, a, b in zip(fields, seen, incoming) if a != b]
    raise ValueError(f"alignment {aln_id}: inconsistent {', '.join(diff)}")


def summarize_alignments(lines):
    """Group segments by ID into alignment-level records."""
    segs = defaultdict(list)
    meta = {}
    for lineno, line in enumerate(lines, 1):
        parsed = _parse(line, lineno)
        if parsed is None:
            continue
        (aln_id, seqid, strand, start, end, identity, score, rank,
         tspan, query) = parsed
        incoming = (seqid, strand, rank, query)
        if aln_id in meta:
            _check_meta(aln_id, meta[aln_id], incoming)
        else:
            meta[aln_id] = incoming
        segs[aln_id].append((start, end, identity, score, tspan))

    return [_build_alignment(aln_id, meta[aln_id], rows)
            for aln_id, rows in segs.items()]


def _components(alignments):
    """Split one (seqid, strand) group into connected overlap components.

    Capacity is only contended between alignments that overlap, so admission
    can be decided independently per component. This also keeps the coordinate
    arrays small.
    """
    ordered = sorted(alignments,
                     key=lambda a: (min(s for s, _ in a.segments),
                                    max(e for _, e in a.segments)))
    comp, reach = [], None
    for aln in ordered:
        lo = min(s for s, _ in aln.segments)
        hi = max(e for _, e in aln.segments)
        if comp and lo > reach:
            yield comp
            comp = []
        comp.append(aln)
        reach = hi if reach is None or lo > reach else max(reach, hi)
    if comp:
        yield comp


def _admit_component(comp, k_base, k_intron, intron_used, family_cap=None):
    """Greedy admission: best quality first, whole alignment or nothing.

    Deciding per alignment (rather than evicting per crowded window) keeps the
    result independent of traversal order and guarantees ID atomicity. It is
    not a global optimum -- this is a weighted packing problem -- but a high
    quality hit always reserves capacity before a weaker one.
    """
    bounds = sorted({b for aln in comp for s, e in aln.segments
                     for b in (s, e + 1)})
    slot = {b: i for i, b in enumerate(bounds)}
    depth = [0] * max(len(bounds), 1)

    accepted, reasons, conflict, rescuable = set(), {}, {}, {}
    family_used = Counter()
    for aln in sorted(comp, key=Alignment.priority):
        if family_cap is not None and family_used[aln.family] >= family_cap:
            reasons[aln.aln_id] = "family_cap"
            continue
        # Count the alignment's own contribution: its segments may overlap each
        # other, and checking only the current depth would let a single ID push
        # a base past the cap by itself.
        need = Counter()
        for s, e in aln.segments:
            for cell in range(slot[s], slot[e + 1]):
                need[cell] += 1
        blocked = [c for c, n in need.items() if depth[c] + n > k_base]
        if blocked:
            # How much of the alignment actually collided decides whether
            # atomicity is cheap or expensive here: losing a whole alignment
            # over a sliver of overlap is the expensive case.
            span = sum(bounds[c + 1] - bounds[c] for c in need
                       if c + 1 < len(bounds))
            hit = sum(bounds[c + 1] - bounds[c] for c in blocked
                      if c + 1 < len(bounds))
            reasons[aln.aln_id] = "base_cap"
            conflict[aln.aln_id] = (hit, span)
            # The part that had room: an upper bound on what any non-atomic
            # policy could have admitted here.
            free = set(need) - set(blocked)
            rescuable[aln.aln_id] = [
                (bounds[c], bounds[c + 1] - 1) for c in sorted(free)
                if c + 1 < len(bounds)]
            continue

        want = Counter((aln.seqid, aln.strand, s, e) for s, e in aln.introns)
        if any(intron_used[k] + n > k_intron for k, n in want.items()):
            reasons[aln.aln_id] = "intron_cap"
            continue

        for cell, n in need.items():
            depth[cell] += n
        for key, n in want.items():
            intron_used[key] += n
        accepted.add(aln.aln_id)
        family_used[aln.family] += 1
    return accepted, reasons, conflict, rescuable


class SelectionReport:
    """Which alignments survived, why the rest did not, and what it cost.

    `cap_lost_coverage_bases` is raw coverage minus retained coverage. It is a
    useful operating statistic but NOT a measure of atomicity cost: with no
    identity floor every non-accepted alignment is a cap reject, so the value
    equals the total loss by construction. The atomicity question is answered
    by `conflict` -- how much of a rejected alignment was actually saturated.
    """

    __slots__ = ("accepted", "reasons", "cap_lost_coverage_bases",
                 "covered_bases", "conflict", "raw_covered_bases",
                 "n_alignments", "rescuable_bp_sum",
                 "rescuable_new_coverage_bases")

    def __init__(self, accepted, reasons, cap_lost_coverage_bases,
                 covered_bases, conflict=None, raw_covered_bases=0,
                 n_alignments=0, rescuable_bp_sum=0,
                 rescuable_new_coverage_bases=0):
        self.accepted = accepted
        self.reasons = reasons
        # Named for what it is: coverage the cap gave up. It is NOT an
        # atomicity cost -- with min_identity=0 every non-accepted alignment is
        # a cap reject, so this equals raw minus accepted by construction.
        self.cap_lost_coverage_bases = cap_lost_coverage_bases
        self.covered_bases = covered_bases
        self.conflict = conflict or {}
        self.raw_covered_bases = raw_covered_bases
        self.n_alignments = n_alignments
        # Upper bound on what relaxing ID atomicity could add: the unsaturated
        # part of rejected alignments (summed, then unioned and net of what
        # accepted alignments already cover).
        self.rescuable_bp_sum = rescuable_bp_sum
        self.rescuable_new_coverage_bases = rescuable_new_coverage_bases


def _run_selection(alignments, k_base, k_intron, min_identity,
                   family_cap=None):
    """Shared admission pass -> (accepted ids, reason per rejected id)."""
    if k_base < 1 or k_intron < 1:
        raise ValueError(
            f"caps must be >= 1 (got k_base={k_base}, k_intron={k_intron})")
    if family_cap is not None and family_cap < 1:
        raise ValueError(f"family_cap must be >= 1 (got {family_cap})")
    groups = defaultdict(list)
    reasons, conflict = {}, {}
    rescuable_by_key, rescuable_bp = defaultdict(list), [0]
    for aln in alignments:
        if aln.identity_w < min_identity:
            reasons[aln.aln_id] = "identity_floor"
            continue
        groups[(aln.seqid, aln.strand)].append(aln)

    accepted = set()
    for _, group in sorted(groups.items()):
        intron_used = defaultdict(int)
        for comp in _components(group):
            ok, why, conf, resc = _admit_component(
                comp, k_base, k_intron, intron_used, family_cap)
            accepted |= ok
            reasons.update(why)
            conflict.update(conf)
            rescuable_by_key[(comp[0].seqid, comp[0].strand)].extend(
                iv for ivs in resc.values() for iv in ivs)
            rescuable_bp[0] += sum(e - s2 + 1 for ivs in resc.values()
                                   for s2, e in ivs)
    return accepted, reasons, conflict, rescuable_by_key, rescuable_bp[0]


def select_alignments(alignments, k_base, k_intron, min_identity=0.0,
                      family_cap=None):
    """Return the set of alignment IDs that fit under both caps."""
    accepted = _run_selection(
        alignments, k_base, k_intron, min_identity, family_cap)[0]
    return accepted


def audit_selection(alignments, k_base, k_intron, min_identity=0.0,
                    family_cap=None):
    """Run selection and attribute the coverage it gave up."""
    (accepted, reasons, conflict, rescuable_by_key,
     rescuable_bp) = _run_selection(
        alignments, k_base, k_intron, min_identity, family_cap)

    kept_by_key = defaultdict(list)
    lost_by_key = defaultdict(list)
    for aln in alignments:
        key = (aln.seqid, aln.strand)
        if aln.aln_id in accepted:
            kept_by_key[key].extend(aln.segments)
        elif reasons.get(aln.aln_id) in ("base_cap", "intron_cap"):
            lost_by_key[key].extend(aln.segments)

    collateral = 0
    for key, lost in lost_by_key.items():
        kept = kept_by_key.get(key, [])
        collateral += (_union_length(lost + kept) - _union_length(kept))

    covered = sum(_union_length(segs) for segs in kept_by_key.values())
    raw_by_key = defaultdict(list)
    for aln in alignments:
        raw_by_key[(aln.seqid, aln.strand)].extend(aln.segments)
    raw_covered = sum(_union_length(v) for v in raw_by_key.values())

    rescuable_new = 0
    for key, ivs in rescuable_by_key.items():
        kept = kept_by_key.get(key, [])
        rescuable_new += _union_length(ivs + kept) - _union_length(kept)

    return SelectionReport(accepted, reasons, collateral, covered, conflict,
                           raw_covered, len(alignments), rescuable_bp,
                           rescuable_new)


def cap_miniprot_lines(lines, k_base, k_intron, min_identity=0.0,
                       family_cap=None):
    """Filter GFF lines down to the alignments that fit under both caps."""
    lines = list(lines)
    accepted = select_alignments(
        summarize_alignments(lines), k_base, k_intron, min_identity,
        family_cap)

    kept = []
    for line in lines:
        if not line or line.startswith("#"):
            kept.append(line)
            continue
        parsed = _parse(line)
        if parsed is None or parsed[0] in accepted:
            kept.append(line)
    return kept


def _summarize_stream(path, tmpdir):
    """Stream alignment summaries without holding the file in memory.

    Real input is neither coordinate- nor ID-sorted (`Chr5:13203494 ->
    13201546 -> Chr2:...`), and at 7 GB / 60M rows the whole-file approach
    needs many times the input size in Python objects. Sorting by ID puts a
    single alignment's segments next to each other, so each one can be folded
    and released as it goes by.
    """
    prefixed = (
        "awk -F'\\t' 'BEGIN{OFS=\"\\t\"} /^#/{next} NF>=9 "
        "{ if (match($9, /ID=[^;]+/)) "
        "{ print substr($9, RSTART+3, RLENGTH-3), $0 } "
        "else { print \"\", $0 } }' " + shlex.quote(path)
        + " | LC_ALL=C sort -t'\t' -k1,1 -S 512M -T " + shlex.quote(tmpdir)
    )
    proc = subprocess.Popen(prefixed, shell=True, stdout=subprocess.PIPE,
                            text=True)
    current, meta, rows = None, None, []
    try:
        for line in proc.stdout:
            key, _, raw = line.rstrip("\n").partition("\t")
            parsed = _parse(raw)
            if parsed is None:
                continue
            (aln_id, seqid, strand, start, end, identity, score, rank,
             tspan, query) = parsed
            incoming = (seqid, strand, rank, query)
            if aln_id != current:
                if current is not None:
                    yield _build_alignment(current, meta, rows)
                current, meta, rows = aln_id, incoming, []
            else:
                _check_meta(aln_id, meta, incoming)
            rows.append((start, end, identity, score, tspan))
        if current is not None:
            yield _build_alignment(current, meta, rows)
    finally:
        proc.stdout.close()
        if proc.wait() != 0:
            raise RuntimeError(f"sort pipeline failed for {path}")


def cap_miniprot_file(infile, outfile, k_base, k_intron, min_identity=0.0,
                      tmpdir=None, family_cap=None):
    """Cap a GFF on disk. Returns the SelectionReport.

    Three passes: fold segments into alignments over an ID-sorted stream,
    decide admission, then copy through the rows of accepted IDs. The output
    is written to a sibling temporary file and renamed, so a failure part-way
    leaves any previous output untouched rather than half-written.
    """
    tmpdir = tmpdir or os.path.dirname(os.path.abspath(outfile)) or "."
    alignments = list(_summarize_stream(infile, tmpdir))
    report = audit_selection(alignments, k_base, k_intron, min_identity,
                             family_cap)

    fd, tmp = tempfile.mkstemp(dir=tmpdir, suffix=".capped.tmp")
    try:
        with os.fdopen(fd, "w") as out, open(infile) as src:
            for lineno, line in enumerate(src, 1):
                line = line.rstrip("\n")
                parsed = _parse(line, lineno)
                if parsed is None or parsed[0] in report.accepted:
                    out.write(line + "\n")
        os.replace(tmp, outfile)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise
    return report


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("infile")
    ap.add_argument("outfile")
    ap.add_argument("--k-base", type=int, required=True,
                    help="max retained coding segments over any base")
    ap.add_argument("--k-intron", type=int, required=True,
                    help="max alignments supporting any exact intron")
    ap.add_argument("--min-identity", type=float, default=0.0,
                    help="alignment-level length-weighted identity floor")
    ap.add_argument("--tmpdir", default=None,
                    help="scratch for external sort (GPFS, not compute /tmp "
                         "which is tmpfs and eats node RAM)")
    ap.add_argument("--family-cap", type=int, default=None, metavar="N",
                    help="opt in: admit at most N alignments from an obvious "
                         "Target isoform family per overlap component")
    args = ap.parse_args(argv)

    report = cap_miniprot_file(args.infile, args.outfile, args.k_base,
                               args.k_intron, args.min_identity, args.tmpdir,
                               args.family_cap)

    by_reason = Counter(report.reasons.values())
    lost = report.raw_covered_bases - report.covered_bases

    print(f"alignments {report.n_alignments} -> {len(report.accepted)}",
          file=sys.stderr)
    print(f"bases      {report.raw_covered_bases} -> {report.covered_bases} "
          f"(lost {lost}, "
          f"{100.0 * lost / report.raw_covered_bases if report.raw_covered_bases else 0:.1f}%)",
          file=sys.stderr)
    for reason, n in sorted(by_reason.items()):
        print(f"  rejected {reason}: {n}", file=sys.stderr)

    bound = report.rescuable_new_coverage_bases
    print(f"  atomicity upper bound: {bound} bp "
          f"({100.0 * bound / lost if lost else 0:.1f}% of loss); "
          f"summed unsaturated {report.rescuable_bp_sum} bp",
          file=sys.stderr)

    if report.conflict:
        ratios = sorted(hit / span for hit, span in report.conflict.values()
                        if span)
        n = len(ratios)
        def pct(q):
            return 100.0 * ratios[min(int(q * n), n - 1)]
        print("  conflicting bp share of rejected alignments: "
              f"median {pct(0.5):.1f}%  p90 {pct(0.9):.1f}%  "
              f"<10%: {100.0 * sum(1 for r in ratios if r < 0.1) / n:.1f}% "
              "of rejects", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
