#!/usr/bin/env python3
"""Unit tests for validate_annotation.py.

Each case is a synthetic GFF3 plus a tiny genome, written to a temp dir and run
through the validator. The point of every case is that it FAILS on the defect it
names -- a checker that only ever passes is worse than none, and the first draft
of this validator did exactly that for multi-segment CDS (it called 11,670 legal
discontinuous features "duplicate IDs" on real data, burying 3,383 real ones).
"""

import json
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
VALIDATOR = HERE / "validate_annotation.py"

# Two designed reading frames, so the translation cases are not accidents:
#   1..15   ATG GGG CCC AAA TAA   -- clean, starts M, ends *
#   21..35  ATG GGG TAA CCC TAA   -- stop at codon 3, i.e. an internal stop
# The rest is N padding so out-of-bounds cases have room to be out of bounds.
GENOME = ">chr1\n" + "ATGGGGCCCAAATAANNNNNATGGGGTAACCCTAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" + "\n"

HEADER = "##gff-version 3\n"


def gff(*rows):
    return HEADER + "".join("\t".join(map(str, r)) + "\n" for r in rows)


def gene(start, end, strand="+", gid="g1"):
    return ("chr1", "test", "gene", start, end, ".", strand, ".", f"ID={gid}")


def mrna(start, end, strand="+", mid="g1.t1", gid="g1"):
    return ("chr1", "test", "mRNA", start, end, ".", strand, ".",
            f"ID={mid};Parent={gid}")


def cds(start, end, strand="+", cid="cds1", mid="g1.t1", phase="0"):
    return ("chr1", "test", "CDS", start, end, ".", strand, phase,
            f"ID={cid};Parent={mid}")


def exon(start, end, strand="+", eid="e1", mid="g1.t1"):
    return ("chr1", "test", "exon", start, end, ".", strand, ".",
            f"ID={eid};Parent={mid}")


def run(gff_text, extra=()):
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        (tmp / "a.gff3").write_text(gff_text)
        (tmp / "g.fa").write_text(GENOME)
        report = tmp / "r.json"
        proc = subprocess.run(
            [sys.executable, str(VALIDATOR), "--gff", str(tmp / "a.gff3"),
             "--genome", str(tmp / "g.fa"), "--report", str(report),
             "--fail-on-error", *extra],
            capture_output=True, text=True)
        payload = json.loads(report.read_text()) if report.exists() else {}
        return proc.returncode, payload


def checks(payload):
    return {e["check"] for e in payload.get("errors", [])}


def case(name, gff_text, expect_check, extra=()):
    code, payload = run(gff_text, extra)
    got = checks(payload)
    if expect_check is None:
        assert code == 0, f"{name}: expected pass, got errors {got}"
        assert not got, f"{name}: expected no errors, got {got}"
    else:
        assert code == 1, f"{name}: expected exit 1, got {code} (errors {got})"
        assert expect_check in got, f"{name}: expected {expect_check}, got {got}"
    print(f"  PASS {name}")


def main():
    print("validate_annotation.py")

    # A clean single-CDS gene must produce no errors at all.
    case("clean gene", gff(gene(1, 30), mrna(1, 30), exon(1, 30), cds(1, 30)), None)

    # Multi-segment CDS sharing one ID is legal GFF3, not a duplicate.
    case("multi-segment CDS shares an ID",
         gff(gene(1, 60), mrna(1, 60),
             exon(1, 30, eid="e1"), exon(40, 60, eid="e2"),
             cds(1, 30, cid="c"), cds(40, 60, cid="c", phase="0")), None)

    # Two unrelated features with one ID is a real collision.
    case("real duplicate ID",
         gff(gene(1, 30), mrna(1, 30), exon(1, 30),
             cds(1, 30, cid="c"),
             ("chr1", "test", "gene", 40, 60, ".", "+", ".", "ID=c")),
         "duplicate_id")

    # The Ptr_ODB chimera shape: a child on the opposite strand.
    case("child on the opposite strand",
         gff(gene(1, 60), mrna(1, 60),
             exon(1, 30), cds(1, 30),
             cds(40, 60, strand="-", cid="c2")),
         "strand_mismatch")

    case("childless mRNA", gff(gene(1, 30), mrna(1, 30)), "childless_transcript")

    case("dangling Parent",
         gff(gene(1, 30), mrna(1, 30), exon(1, 30, mid="nope")),
         "missing_parent")

    case("child outside its parent",
         gff(gene(1, 30), mrna(1, 30), exon(1, 30), cds(1, 45, cid="c")),
         "child_outside_parent")

    case("feature past the end of the contig",
         gff(gene(1, 500), mrna(1, 500), exon(1, 500), cds(1, 500)),
         "feature_past_contig_end")

    case("phase outside 0/1/2",
         gff(gene(1, 30), mrna(1, 30), exon(1, 30), cds(1, 30, phase="7")),
         "invalid_cds_phase")

    case("ID carrying a path separator",
         gff(("chr1", "test", "gene", 1, 30, ".", "+", ".",
              "ID=results/FILTER/00000010"),
             ("chr1", "test", "mRNA", 1, 30, ".", "+", ".",
              "ID=m1;Parent=results/FILTER/00000010"),
             exon(1, 30, mid="m1"), cds(1, 30, mid="m1")),
         "bad_id_character")

    case("seqid absent from the genome",
         gff(("chrX", "test", "gene", 1, 30, ".", "+", ".", "ID=g1"),
             ("chrX", "test", "mRNA", 1, 30, ".", "+", ".", "ID=g1.t1;Parent=g1"),
             ("chrX", "test", "exon", 1, 30, ".", "+", ".", "ID=e1;Parent=g1.t1"),
             ("chrX", "test", "CDS", 1, 30, ".", "+", "0", "ID=c1;Parent=g1.t1")),
         "seqid_not_in_genome")

    # 1..15 is the clean frame: M G P K *
    case("clean translation",
         gff(gene(1, 15), mrna(1, 15), exon(1, 15), cds(1, 15)),
         None, extra=("--check-translation",))

    # 21..35 stops at codon 3, so the protein carries a stop before its end.
    case("internal stop codon",
         gff(gene(21, 35), mrna(21, 35), exon(21, 35), cds(21, 35)),
         "internal_stop_codon", extra=("--check-translation",))

    print("\nall cases passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
