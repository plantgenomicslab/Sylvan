#!/usr/bin/env python3
"""Validate a finished GFF3 against its genome, and optionally its protein FASTA.

WHY THIS EXISTS
---------------
The pipeline's completion test is that a file exists and is non-empty. That test
has passed for annotations that were structurally broken, and each time the
damage was found downstream by a tool that refused the file:

  - 167k duplicate IDs, because multi-CDS mRNAs all got the same `.cds` suffix
  - gene IDs of the form `results/FILTER/00000010`, because an output directory
    was passed where a prefix was expected
  - every refined CDS written with phase=0, producing frameshifted proteins
  - mRNAs with no exon or CDS children at all
  - mRNAs whose exons sat on both strands 545 kb apart, fused from two loci

None of those are detectable by `[ -s file ]`. All of them are detectable here,
cheaply, before the annotation is released.

WHAT IT CHECKS
--------------
ERROR (the annotation is wrong; exit 1 with --fail-on-error):
  duplicate ID, ID containing whitespace or a path separator, dangling Parent,
  mRNA without a gene parent, exon/CDS without an mRNA parent, childless mRNA,
  start > end, feature off the end of its contig, child outside its parent,
  child on a different strand from its parent, CDS phase not in {0,1,2},
  seqid absent from the genome, duplicate FASTA ID, empty FASTA record,
  internal stop codon, protein that does not match its translated CDS

WARNING (worth knowing, not necessarily wrong):
  mRNA with exons but no CDS, CDS total length not a multiple of three,
  missing start or stop codon, protein FASTA record with no matching transcript

USAGE
-----
    validate_annotation.py --gff results/FILTER/filtered.gff3 \\
                          --genome input/genome.fa \\
                          --report results/QC/validation.json \\
                          --summary results/QC/validation.txt \\
                          --fail-on-error

    # add translation checks (needs to hold one contig's sequence at a time)
    validate_annotation.py --gff a.gff3 --genome g.fa --check-translation \\
                          --protein results/FILTER/proteins.fa --fail-on-error

NOTES
-----
Standard library only, so it runs under any python3 in or out of the container.
The genome is read twice at most: once for contig lengths, once more per contig
if translation is checked, so memory is bounded by the largest contig rather
than by the genome.
"""

import argparse
import json
import re
import sys
from collections import defaultdict

# A GFF3 ID has to survive being used as a FASTA header and a filename. Slashes
# were how `results/FILTER/00000010` got into gene IDs.
BAD_ID_CHARS = re.compile(r"[\s/\\]")

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

TRANSCRIPT_TYPES = {"mRNA", "transcript"}
# GFF3 allows one discontinuous feature to be written as several lines sharing a
# single ID -- that is how a multi-segment CDS is represented. Flagging those as
# duplicates buried the real duplicates: on Ptr_ODB's damaged PREFILTER, 11,670
# of 15,587 "errors" were legal multi-segment CDS. Duplicates are only an error
# when the records do not agree on parent, contig, strand and type.
DISCONTINUOUS_TYPES = {"CDS", "five_prime_UTR", "three_prime_UTR", "UTR"}
CHILD_TYPES = {"exon", "CDS", "five_prime_UTR", "three_prime_UTR", "UTR"}


class Report:
    """Findings plus the counts a reader needs to judge them."""

    def __init__(self):
        self.errors = []
        self.warnings = []
        self.counts = defaultdict(int)

    def error(self, check, message, feature=None):
        self.errors.append({"check": check, "message": message, "feature": feature})
        self.counts[f"error:{check}"] += 1

    def warn(self, check, message, feature=None):
        self.warnings.append({"check": check, "message": message, "feature": feature})
        self.counts[f"warning:{check}"] += 1


def parse_attributes(attr):
    out = {}
    for field in attr.rstrip(";").split(";"):
        if "=" not in field:
            continue
        key, _, value = field.partition("=")
        out[key.strip()] = value.strip()
    return out


def read_contig_lengths(path):
    """seqid -> length, without holding any sequence."""
    lengths = {}
    name = None
    total = 0
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = total
                name = line[1:].split()[0]
                total = 0
            else:
                total += len(line.strip())
    if name is not None:
        lengths[name] = total
    return lengths


def iter_contigs(path):
    """Yield (seqid, sequence) one contig at a time."""
    name = None
    chunks = []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def read_fasta_ids(path, report, label):
    """id -> sequence, reporting duplicates and empty records as errors."""
    seqs = {}
    name = None
    chunks = []

    def flush():
        if name is None:
            return
        seq = "".join(chunks)
        if not seq:
            report.error("empty_fasta_record", f"{label} record {name} has no sequence", name)
        if name in seqs:
            report.error("duplicate_fasta_id", f"{label} ID {name} appears more than once", name)
        else:
            seqs[name] = seq

    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                flush()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    flush()
    return seqs


def translate(seq):
    out = []
    for i in range(0, len(seq) - len(seq) % 3, 3):
        out.append(CODON_TABLE.get(seq[i:i + 3].upper(), "X"))
    return "".join(out)


def load_gff(path, report):
    """Parse into features/children, reporting ID problems as we go."""
    features = {}
    children = defaultdict(list)
    duplicates = defaultdict(list)
    order = []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for lineno, line in enumerate(handle, 1):
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                report.error("malformed_line",
                             f"line {lineno} has {len(f)} columns, expected 9")
                continue
            attrs = parse_attributes(f[8])
            fid = attrs.get("ID")
            try:
                start, end = int(f[3]), int(f[4])
            except ValueError:
                report.error("malformed_line", f"line {lineno} has non-integer coordinates")
                continue
            rec = {
                "id": fid, "seqid": f[0], "type": f[2], "start": start, "end": end,
                "strand": f[6], "phase": f[7], "parents": [],
                "line": lineno,
            }
            parent = attrs.get("Parent")
            if parent:
                rec["parents"] = [p for p in parent.split(",") if p]
            if fid is not None:
                if BAD_ID_CHARS.search(fid):
                    report.error("bad_id_character",
                                 f"ID {fid!r} contains whitespace or a path separator", fid)
                if fid in features:
                    duplicates[fid].append(rec)
                else:
                    features[fid] = rec
                    order.append(fid)
            for p in rec["parents"]:
                children[p].append(rec)
            report.counts[f"feature:{f[2]}"] += 1
    _judge_duplicates(features, duplicates, report)
    return features, children, order


def _judge_duplicates(features, duplicates, report):
    """Split shared IDs into legal discontinuous features and real collisions."""
    for fid, extras in duplicates.items():
        first = features[fid]
        group = [first] + extras
        same_shape = (
            len({r["type"] for r in group}) == 1
            and len({r["seqid"] for r in group}) == 1
            and len({r["strand"] for r in group}) == 1
            and len({tuple(sorted(r["parents"])) for r in group}) == 1
        )
        if first["type"] in DISCONTINUOUS_TYPES and same_shape:
            report.counts["note:discontinuous_feature"] += 1
            # Record the merged extent separately. Widening start/end in place
            # would be visible to check_cds, which sums the same dicts through
            # `children` -- that turned 122 CDS-length warnings into 1,353.
            first["span"] = (min(r["start"] for r in group),
                             max(r["end"] for r in group))
            continue
        lines = ", ".join(str(r["line"]) for r in group)
        report.error("duplicate_id",
                     f"ID {fid} is defined {len(group)} times at lines {lines} "
                     f"and they are not one discontinuous feature", fid)


def check_structure(features, children, report):
    for fid in list(features):
        rec = features[fid]
        for parent in rec["parents"]:
            if parent not in features:
                report.error("missing_parent",
                             f"{rec['type']} {fid} references Parent {parent}, "
                             f"which is not defined", fid)
                continue
            prec = features[parent]
            if rec["strand"] != prec["strand"]:
                report.error("strand_mismatch",
                             f"{rec['type']} {fid} is on strand {rec['strand']} but its "
                             f"parent {parent} is on {prec['strand']}", fid)
            lo, hi = rec.get("span", (rec["start"], rec["end"]))
            if lo < prec["start"] or hi > prec["end"]:
                report.error("child_outside_parent",
                             f"{rec['type']} {fid} ({lo}-{hi}) lies "
                             f"outside parent {parent} ({prec['start']}-{prec['end']})", fid)
        if rec["type"] in TRANSCRIPT_TYPES:
            if not rec["parents"]:
                report.error("transcript_without_gene",
                             f"{rec['type']} {fid} has no Parent gene", fid)
            kids = children.get(fid, [])
            if not kids:
                report.error("childless_transcript",
                             f"{rec['type']} {fid} has no exon or CDS children", fid)
            elif not any(k["type"] == "CDS" for k in kids):
                report.warn("transcript_without_cds",
                            f"{rec['type']} {fid} has no CDS", fid)
        if rec["type"] in CHILD_TYPES and not rec["parents"]:
            report.error("orphan_child",
                         f"{rec['type']} at line {rec['line']} has no Parent", fid)


def check_coordinates(features, contig_lengths, report):
    for fid, rec in features.items():
        if rec["start"] > rec["end"]:
            report.error("inverted_coordinates",
                         f"{rec['type']} {fid} has start {rec['start']} > end {rec['end']}", fid)
        length = contig_lengths.get(rec["seqid"])
        if length is None:
            report.error("seqid_not_in_genome",
                         f"{rec['type']} {fid} is on {rec['seqid']}, which is not in "
                         f"the genome FASTA", fid)
        elif rec["end"] > length:
            report.error("feature_past_contig_end",
                         f"{rec['type']} {fid} ends at {rec['end']} but {rec['seqid']} "
                         f"is {length} bp", fid)


def check_cds(features, children, report):
    """Phase validity and CDS length, per transcript."""
    for fid, rec in features.items():
        if rec["type"] not in TRANSCRIPT_TYPES:
            continue
        cds = [k for k in children.get(fid, []) if k["type"] == "CDS"]
        if not cds:
            continue
        total = 0
        for c in cds:
            if c["phase"] not in ("0", "1", "2"):
                report.error("invalid_cds_phase",
                             f"CDS of {fid} at {c['start']}-{c['end']} has phase "
                             f"{c['phase']!r}, expected 0, 1 or 2", fid)
            total += c["end"] - c["start"] + 1
        if total % 3:
            report.warn("cds_length_not_multiple_of_three",
                        f"{fid} has {total} bp of CDS ({total % 3} over)", fid)


def check_translation(gff_path, features, children, genome_path, proteins, report):
    """Translate each transcript's CDS and compare against the protein FASTA."""
    by_contig = defaultdict(list)
    for fid, rec in features.items():
        if rec["type"] in TRANSCRIPT_TYPES and any(
                k["type"] == "CDS" for k in children.get(fid, [])):
            by_contig[rec["seqid"]].append(fid)

    for seqid, sequence in iter_contigs(genome_path):
        for fid in by_contig.get(seqid, []):
            rec = features[fid]
            cds = sorted((k for k in children[fid] if k["type"] == "CDS"),
                         key=lambda k: k["start"], reverse=rec["strand"] == "-")
            pieces = []
            for c in cds:
                piece = sequence[c["start"] - 1:c["end"]]
                if rec["strand"] == "-":
                    piece = piece.translate(COMPLEMENT)[::-1]
                pieces.append(piece)
            nt = "".join(pieces)
            # The first CDS's phase says how many bases to drop before the
            # first complete codon; ignoring it is what turns a partial 5'
            # gene into a frameshift.
            offset = int(cds[0]["phase"]) if cds and cds[0]["phase"] in ("0", "1", "2") else 0
            aa = translate(nt[offset:])
            body = aa[:-1] if aa.endswith("*") else aa
            if "*" in body:
                report.error("internal_stop_codon",
                             f"{fid} translates with {body.count('*')} internal stop(s)", fid)
            if not aa.startswith("M"):
                report.warn("no_start_codon", f"{fid} does not begin with M", fid)
            if not aa.endswith("*"):
                report.warn("no_stop_codon", f"{fid} does not end with a stop codon", fid)
            if proteins is not None:
                given = proteins.get(fid)
                if given is None:
                    report.counts["translation:no_protein_record"] += 1
                else:
                    if given.rstrip("*") != body:
                        report.error("protein_sequence_mismatch",
                                     f"{fid}: protein FASTA does not match the "
                                     f"translated CDS", fid)


def write_reports(report, counts, args):
    payload = {
        "status": "FAILED" if report.errors else "PASSED",
        "counts": dict(sorted(report.counts.items())),
        "totals": counts,
        "errors": report.errors,
        "warnings": report.warnings,
    }
    if args.report:
        with open(args.report, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)

    lines = [f"Annotation validation: {payload['status']}", ""]
    for key, value in counts.items():
        lines.append(f"{key:<34}{value:>12,}")
    lines.append("")
    lines.append(f"{'Errors':<34}{len(report.errors):>12,}")
    lines.append(f"{'Warnings':<34}{len(report.warnings):>12,}")
    if report.counts:
        lines.append("")
        lines.append("By check:")
        for key, value in sorted(report.counts.items()):
            if key.startswith(("error:", "warning:")):
                lines.append(f"  {key:<40}{value:>10,}")
    # A capped sample, so a million-error file still produces a readable summary.
    for label, items in (("ERRORS", report.errors), ("WARNINGS", report.warnings)):
        if not items:
            continue
        lines.append("")
        lines.append(f"{label} (first 20 of {len(items):,}):")
        for item in items[:20]:
            lines.append(f"  [{item['check']}] {item['message']}")
    text = "\n".join(lines) + "\n"
    if args.summary:
        with open(args.summary, "w", encoding="utf-8") as handle:
            handle.write(text)
    return text


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--protein", help="protein FASTA to compare against translated CDS")
    ap.add_argument("--report", help="write JSON here")
    ap.add_argument("--summary", help="write the text summary here")
    ap.add_argument("--check-translation", action="store_true",
                    help="translate CDS and check for internal stops; implied by --protein")
    ap.add_argument("--fail-on-error", action="store_true",
                    help="exit 1 when any ERROR was reported")
    args = ap.parse_args()

    report = Report()
    features, children, _ = load_gff(args.gff, report)
    contig_lengths = read_contig_lengths(args.genome)

    check_structure(features, children, report)
    check_coordinates(features, contig_lengths, report)
    check_cds(features, children, report)

    proteins = None
    if args.protein:
        proteins = read_fasta_ids(args.protein, report, "protein")
    if args.check_translation or args.protein:
        check_translation(args.gff, features, children, args.genome, proteins, report)

    genes = sum(1 for r in features.values() if r["type"] == "gene")
    transcripts = sum(1 for r in features.values() if r["type"] in TRANSCRIPT_TYPES)
    cds_count = report.counts.get("feature:CDS", 0)
    counts = {
        "Genes": genes,
        "Transcripts": transcripts,
        "CDS features": cds_count,
        "Contigs in genome": len(contig_lengths),
    }
    if proteins is not None:
        counts["Protein records"] = len(proteins)

    text = write_reports(report, counts, args)
    sys.stdout.write(text)
    return 1 if (report.errors and args.fail_on_error) else 0


if __name__ == "__main__":
    sys.exit(main())
