#!/usr/bin/env python3
"""Drop childless mRNA shells from the annotation GFF3 before sequence extraction.

The PREFILTER GFF3 carries mRNA rows that have **no exon and no CDS children**
(tracker wyim-pgl/Sylvan-EGAPx#51). They are all ``agat-mrna-*`` records left
behind when the EVM->PREFILTER merge re-parented a model's features onto a
neighbouring gene: the features moved, the mRNA row did not. Every run has
some (Ath 79, Ptr 169, Spe_A 175, Dca 543, Osa 2,521).

They are not merely dead weight -- they break the two consumers of this file
in different ways:

**EVM's gff3_file_to_proteins.pl (#50).** A childless shell is often placed on
the *opposite strand* from the real mRNA under the same gene, so the gene ends
up with isoforms of conflicting orientation -- invalid GFF3. Gene_obj then
assembles a sequence far larger than the locus: Osa produced an 82,338,402 bp
"transcript" for ``evm.model.Osa_Chr1.326`` (correct answer 1,746 bp) on a
43.3 Mb chromosome, and a 63,683,521 aa "protein" for
``evm.model.Osa_Chr1.2051``. The rule then runs until it exhausts memory --
raising the memory limit does not help, because the output is unbounded.
Strand-mixed genes track the shells exactly: Osa 1,261, every other run 0,
which is why only Osa hit the wall.

**gffread -w (#49).** gffread does *not* skip a childless mRNA; it falls back
to the mRNA row's span and emits raw genomic sequence, introns and all. Those
enter the RSEM reference as fake transcripts. That is the source of the
``SAM/BAM file declares less reference sequences (44506) than RSEM knows
(45049)`` warning -- the 543 extra are exactly these shells, which STAR's
transcriptome index omits. Removing them makes the two agree.

Dropping a shell cannot change a real model's sequence: gff3_file_to_proteins
already skipped these (it emits ``-warning, no cDNA sequence for ...``), so
``.cdna``/``.pep`` are unchanged, and only ``.mrna`` loses the fake entries.

An mRNA is dropped together with any child feature that points at it (UTR rows
survive re-parenting even when exon/CDS do not). A gene row is dropped only
once no mRNA remains under it.

Usage:
    python sanitize_filter_gff.py input.gff3 output.gff3
"""
import re
import sys
from collections import defaultdict

ID_RE = re.compile(r"ID=([^;\s]+)")
PARENT_RE = re.compile(r"Parent=([^;\s]+)")

# Feature types that make an mRNA usable for sequence extraction. UTR rows are
# deliberately excluded: a shell that kept only its UTRs still yields no
# transcript, and gffread would fall back to the genomic span for it.
STRUCTURAL = ("exon", "CDS")


def _attr(regex, attrs):
    match = regex.search(attrs)
    return match.group(1) if match else None


def _parents(attrs):
    """Return every Parent of a row. GFF3 allows a comma-separated list."""
    raw = _attr(PARENT_RE, attrs)
    return raw.split(",") if raw else []


def find_childless_mrnas(lines):
    """Return the set of mRNA IDs that have no exon and no CDS child."""
    all_mrnas = set()
    with_structure = set()
    for line in lines:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        if fields[2] == "mRNA":
            mrna_id = _attr(ID_RE, fields[8])
            if mrna_id is not None:
                all_mrnas.add(mrna_id)
        elif fields[2] in STRUCTURAL:
            with_structure.update(_parents(fields[8]))
    return all_mrnas - with_structure


def count_mixed_strand_genes(lines, dropped=frozenset()):
    """Count genes whose surviving mRNAs disagree on strand.

    Invalid GFF3, and the condition that makes EVM's Gene_obj emit sequences
    longer than the chromosome. Reported before and after so the fix is
    auditable in the rule log.
    """
    strands = defaultdict(set)
    for line in lines:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "mRNA":
            continue
        mrna_id = _attr(ID_RE, fields[8])
        if mrna_id in dropped:
            continue
        for parent in _parents(fields[8]):
            strands[parent].add(fields[6])
    return sum(1 for seen in strands.values() if len(seen) > 1)


def sanitize(lines, dropped):
    """Yield the GFF3 with dropped mRNAs, their children, and emptied genes gone."""
    surviving_genes = set()
    for line in lines:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "mRNA":
            continue
        if _attr(ID_RE, fields[8]) in dropped:
            continue
        surviving_genes.update(_parents(fields[8]))

    for line in lines:
        if line.startswith("#"):
            yield line
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            yield line
            continue
        feature, attrs = fields[2], fields[8]
        if feature == "gene":
            if _attr(ID_RE, attrs) in surviving_genes:
                yield line
            continue
        if feature == "mRNA":
            if _attr(ID_RE, attrs) not in dropped:
                yield line
            continue
        # Any other row is a child: drop it when every Parent was dropped.
        parents = _parents(attrs)
        if parents and all(parent in dropped for parent in parents):
            continue
        yield line


def find_strand_inconsistent_mrnas(lines):
    """Return mRNA IDs whose exon/CDS children disagree on strand.

    This is the condition that makes EVM's extractor emit sequences longer than
    the chromosome, and it cannot be repaired here: such an mRNA is two real
    gene models fused under one ID by the upstream merge (issue #51), so any
    local "fix" would be a guess at which half to keep. Abort instead -- the
    nine-species comparison assumes one protocol, and silently discarding 2.3%
    of one species' models would break that quietly.
    """
    mrna_strand = {}
    child_strands = defaultdict(set)
    for line in lines:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        if fields[2] == "mRNA":
            mrna_id = _attr(ID_RE, fields[8])
            if mrna_id is not None:
                mrna_strand[mrna_id] = fields[6]
        elif fields[2] in STRUCTURAL:
            for parent in _parents(fields[8]):
                child_strands[parent].add(fields[6])
    bad = set()
    for mrna_id, strands in child_strands.items():
        if len(strands) > 1:
            bad.add(mrna_id)
        elif mrna_id in mrna_strand and strands - {mrna_strand[mrna_id]}:
            bad.add(mrna_id)
    return bad


def main(argv):
    if len(argv) != 3:
        sys.exit(__doc__)
    in_path, out_path = argv[1], argv[2]

    with open(in_path) as handle:
        lines = handle.readlines()

    chimeric = find_strand_inconsistent_mrnas(lines)
    if chimeric:
        sample = ", ".join(sorted(chimeric)[:3])
        sys.exit(
            f"sanitize_filter_gff: ERROR {len(chimeric)} mRNA(s) in {in_path} have "
            f"exon/CDS children on conflicting strands (e.g. {sample}).\n"
            "These are two gene models fused under one ID by the PASA/EVM merge "
            "(wyim-pgl/Sylvan-EGAPx#51). Extracting sequence from them produces "
            "output larger than the chromosome and exhausts memory (#50).\n"
            "Regenerate PREFILTER with the fixed merge_pasa_evm.py rather than "
            "filtering around this."
        )

    dropped = find_childless_mrnas(lines)
    mixed_before = count_mixed_strand_genes(lines)
    mixed_after = count_mixed_strand_genes(lines, dropped)

    with open(out_path, "w") as handle:
        handle.writelines(sanitize(lines, dropped))

    print(
        f"sanitize_filter_gff: dropped {len(dropped)} childless mRNA(s); "
        f"strand-mixed genes {mixed_before} -> {mixed_after}",
        file=sys.stderr,
    )
    if mixed_after:
        # A gene may legitimately hold same-ID-free isoforms on both strands
        # once every mRNA is internally consistent; that is odd but not the
        # unbounded-output trigger, so report rather than abort.
        print(
            f"sanitize_filter_gff: WARNING {mixed_after} gene(s) still carry "
            "mRNAs on both strands",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
