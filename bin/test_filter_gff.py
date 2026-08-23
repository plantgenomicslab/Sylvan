#!/usr/bin/env python3
"""Unit tests for Filter.py's filter_gff.

The case that matters is the last one. filter_gff used a single flat name set
for gene rows, mRNA rows and exon/CDS rows. Gene IDs were added to it (by
stripping `.tN`) so gene rows would survive a partial discard, but the parent
test cannot tell a gene row from an mRNA row, so a *discarded* mRNA whose Parent
is that gene matched it too. Its children matched nothing and went to the
discard file -- leaving an mRNA with no exon or CDS in the released annotation
(issue #37: 138 across the nine head-to-head runs; in 9 of Osa's 10 cases the
empty record had taken `.t1`, which primary_fasta.py picks as the
representative isoform).
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from Filter import filter_gff  # noqa: E402  (needs the path insert above)


def row(ftype, start, end, attrs):
    return ["chr1", "test", ftype, start, end, ".", "+", ".", attrs]


def frame(*rows):
    return pd.DataFrame(list(rows))


def keep_frame(*new_ids):
    return pd.DataFrame({"New_ID": list(new_ids)})


def ids(df):
    import re
    out = []
    for attr in df[8]:
        m = re.search(r"ID=([^;]+)", attr)
        if m:
            out.append(m.group(1))
    return out


def childless(df):
    """mRNA IDs in df with no exon/CDS child in df -- the #37 symptom."""
    import re
    from collections import defaultdict
    kids = defaultdict(int)
    mrna = []
    for ftype, attr in zip(df[2], df[8]):
        if ftype == "mRNA":
            m = re.search(r"ID=([^;]+)", attr)
            if m:
                mrna.append(m.group(1))
        elif ftype in ("exon", "CDS"):
            m = re.search(r"Parent=([^;]+)", attr)
            if m:
                for p in m.group(1).split(","):
                    kids[p] += 1
    return [m for m in mrna if kids[m] == 0]


def two_isoform_gene():
    """g1 with t1 and t2, each with its own exon and CDS."""
    return frame(
        row("gene", 100, 900, "ID=g1"),
        row("mRNA", 100, 500, "ID=g1.t1;Parent=g1"),
        row("exon", 100, 500, "ID=g1.t1.exon1;Parent=g1.t1"),
        row("CDS", 150, 450, "ID=g1.t1.cds1;Parent=g1.t1"),
        row("mRNA", 200, 900, "ID=g1.t2;Parent=g1"),
        row("exon", 200, 900, "ID=g1.t2.exon1;Parent=g1.t2"),
        row("CDS", 250, 850, "ID=g1.t2.cds1;Parent=g1.t2"),
    )


def test_all_kept():
    kept, discarded = filter_gff(two_isoform_gene(), keep_frame("g1.t1", "g1.t2"))
    assert len(discarded) == 0, f"nothing should be discarded, got {ids(discarded)}"
    assert not childless(kept)
    print("  PASS both isoforms kept")


def test_all_discarded():
    kept, discarded = filter_gff(two_isoform_gene(), keep_frame())
    assert len(kept) == 0, f"nothing should survive, got {ids(kept)}"
    assert "g1" in ids(discarded), "the gene row must follow its transcripts out"
    print("  PASS whole gene discarded")


def test_shared_exon_multi_parent():
    """#20.5: an exon listing two parents survives if either is kept."""
    df = frame(
        row("gene", 100, 900, "ID=g1"),
        row("mRNA", 100, 500, "ID=g1.t1;Parent=g1"),
        row("mRNA", 200, 900, "ID=g1.t2;Parent=g1"),
        row("exon", 200, 500, "ID=shared.exon1;Parent=g1.t1,g1.t2"),
        row("CDS", 250, 450, "ID=g1.t2.cds1;Parent=g1.t2"),
    )
    kept, _ = filter_gff(df, keep_frame("g1.t2"))
    assert "shared.exon1" in ids(kept), "shared exon dropped though g1.t2 is kept"
    print("  PASS shared exon with two parents")


def test_partial_discard_leaves_no_orphan_mrna():
    """Issue #37. t1 is discarded, t2 is kept; t1 must leave whole."""
    kept, discarded = filter_gff(two_isoform_gene(), keep_frame("g1.t2"))

    orphans = childless(kept)
    assert not orphans, (
        f"a discarded transcript was resurrected without its children: {orphans}")

    assert "g1" in ids(kept), "the gene row must survive with t2"
    assert "g1.t2" in ids(kept)
    assert "g1.t1" not in ids(kept), "t1 was discarded and must not appear in kept"
    assert "g1.t1" in ids(discarded)
    assert "g1.t1.exon1" in ids(discarded) and "g1.t1.cds1" in ids(discarded), (
        "t1's children must travel with it")
    print("  PASS partial discard leaves no orphan mRNA (#37)")


def main():
    print("Filter.py filter_gff")
    test_all_kept()
    test_all_discarded()
    test_shared_exon_multi_parent()
    test_partial_discard_leaves_no_orphan_mrna()
    print("\nall cases passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
