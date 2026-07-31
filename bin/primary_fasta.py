#!/usr/bin/env python3
"""Extract the primary isoform from a FASTA, optionally renaming it to the gene ID.

WHY THIS EXISTS
---------------
Downstream tools that score an annotation -- OMArk, BUSCO, OrthoFinder, InterProScan
-- expect one representative protein (or CDS) per gene. Feeding them every isoform
inflates duplication metrics and wastes runtime. `TidyGFF.py` numbers isoforms per
gene as `<gene>.t1`, `<gene>.t2`, ... and always emits `.t1`, so selecting records
whose ID ends in `.t1` yields exactly one record per gene.

Verified on B. tournefortii (`results/FILTER/filtered.gff3`, 43,981 genes):
43,981 `.t1` records against 4,102 `.t2`, 84 `.t3` and 2 `.t4` -- one `.t1` per gene,
no gene without one.

This complements `Pick_Primaries.py`, which selects primaries by TransDecoder score
and therefore needs a TransDecoder run first. Use this one when the annotation has
already been through `TidyGFF.py` and the `.tN` numbering is authoritative.

USAGE
-----
    # protein FASTA -> one protein per gene, IDs left as <gene>.t1
    primary_fasta.py filtered.pep primary.pep

    # same, but IDs renamed to <gene> -- what OMArk/BUSCO/OrthoFinder want
    primary_fasta.py filtered.pep primary.pep --rename

    # CDS works identically; the ID convention is shared
    primary_fasta.py filtered.cds primary.cds --rename

    # a different isoform-numbering scheme
    primary_fasta.py in.pep out.pep --suffix .t01

    # stdin/stdout, so it composes in a pipeline
    gffread filtered.gff3 -g genome.fa -y - | primary_fasta.py - primary.pep --rename

    # keep the header text after the first whitespace (dropped by default)
    primary_fasta.py filtered.pep primary.pep --rename --keep-description

TYPICAL PIPELINE
----------------
    gffread results/FILTER/filtered.gff3 -g input/genome.fa -y all.pep
    primary_fasta.py all.pep primary.pep --rename
    busco -i primary.pep -l brassicales_odb10 -m proteins -o busco_primary

OUTPUT
------
Filtered FASTA on `output`; a one-line count and any warnings on stderr, so stdout
stays clean when writing to `-`. Warnings are raised when nothing matched the
suffix, when the suffix has no leading dot (`t1` would also match `sampleXt1`), and
when renaming produces duplicate gene IDs -- that last one means the suffix does not
identify a unique isoform per gene and the selection is wrong.

NOTES
-----
Streams the input and depends only on the standard library, so it runs under any
python3 and on files larger than memory. `fasta_utils.readFasta` is not used here
because it materialises the whole file as a dict.

The legacy positional call `primary_fasta.py in out [suffix]` still works, but the
default suffix is `.t1` rather than the older dotless `t1`.
"""
import argparse
import sys


def stream_fasta(path):
    """Yield (header_without_'>', sequence_lines) without holding the file in memory."""
    handle = sys.stdin if path == "-" else open(path)
    try:
        header, seq = None, []
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, seq
                header, seq = line[1:].rstrip("\n"), []
            elif header is not None:
                seq.append(line.rstrip("\n"))
        if header is not None:
            yield header, seq
    finally:
        if handle is not sys.stdin:
            handle.close()


def main():
    ap = argparse.ArgumentParser(
        description="Keep only primary-isoform records (default: IDs ending in .t1).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="examples:\n"
               "  primary_fasta.py filtered.pep primary.pep --rename\n"
               "  primary_fasta.py filtered.cds primary.cds --rename\n"
               "  gffread ann.gff3 -g gen.fa -y - | primary_fasta.py - out.pep --rename")
    ap.add_argument("input", help="input FASTA ('-' for stdin)")
    ap.add_argument("output", help="output FASTA ('-' for stdout)")
    ap.add_argument("suffix", nargs="?", default=".t1",
                    help="isoform suffix to keep (default: .t1)")
    ap.add_argument("--suffix", dest="suffix_opt", default=None,
                    help="same as the positional suffix; takes precedence")
    ap.add_argument("--rename", action="store_true",
                    help="strip the suffix so the ID becomes the gene ID")
    ap.add_argument("--keep-description", action="store_true",
                    help="keep header text after the first whitespace "
                         "(dropped by default -- most downstream tools split on it)")
    args = ap.parse_args()

    suffix = args.suffix_opt if args.suffix_opt is not None else args.suffix
    if not suffix:
        ap.error("suffix must not be empty")
    # A bare "t1" also matches IDs like "sampleXt1". That was the old default, so
    # warn rather than silently widening the selection.
    if not suffix.startswith("."):
        print(f"WARNING: suffix {suffix!r} has no leading dot -- it will also match "
              f"IDs that merely end in those characters. Did you mean '.{suffix}'?",
              file=sys.stderr)

    out = sys.stdout if args.output == "-" else open(args.output, "w")
    total = kept = dup = 0
    seen = set()
    try:
        for header, seq in stream_fasta(args.input):
            total += 1
            seq_id = header.split(None, 1)[0]
            if not seq_id.endswith(suffix):
                continue
            kept += 1
            desc = header.split(None, 1)[1:]
            if args.rename:
                seq_id = seq_id[: -len(suffix)]
                if seq_id in seen:
                    dup += 1
                seen.add(seq_id)
            header = " ".join([seq_id] + desc) if args.keep_description else seq_id
            out.write(">" + header + "\n")
            out.write("\n".join(seq) + "\n")
    finally:
        if out is not sys.stdout:
            out.close()

    print(f"{kept} / {total} records kept (suffix {suffix!r}"
          f"{', renamed to gene ID' if args.rename else ''})", file=sys.stderr)
    if kept == 0:
        print(f"WARNING: nothing matched {suffix!r} -- check the ID format.",
              file=sys.stderr)
    if dup:
        print(f"WARNING: {dup} duplicate gene IDs after renaming -- the suffix does "
              f"not uniquely identify one isoform per gene.", file=sys.stderr)


if __name__ == "__main__":
    main()
