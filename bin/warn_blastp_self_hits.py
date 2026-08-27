#!/usr/bin/env python3
"""Report a likely circular-homology signal without changing BLAST results."""

import argparse
import sys


def summarize(path, identity_cutoff):
	total = 0
	high_identity = 0
	malformed = 0
	with open(path, encoding="utf-8") as handle:
		for line in handle:
			if not line.strip() or line.startswith("#"):
				continue
			fields = line.rstrip("\n").split("\t")
			try:
				pident = float(fields[2])
			except (IndexError, ValueError):
				malformed += 1
				continue
			total += 1
			high_identity += pident >= identity_cutoff
	return total, high_identity, malformed


def main():
	parser = argparse.ArgumentParser(
		description="Warn when near-identical BLASTP best hits suggest self-species proteins"
	)
	parser.add_argument("blast_table", help="BLAST outfmt 6 table (pident in column 3)")
	parser.add_argument("--log", required=True, help="diagnostic log file to write")
	parser.add_argument("--identity-cutoff", type=float, default=99.5)
	parser.add_argument("--fraction-threshold", type=float, default=0.20)
	args = parser.parse_args()

	total, high_identity, malformed = summarize(args.blast_table, args.identity_cutoff)
	fraction = high_identity / total if total else 0.0
	summary = (
		"BLASTP circular-homology diagnostic: "
		f"{high_identity}/{total} hits ({fraction:.1%}) have "
		f"pident >= {args.identity_cutoff:g}; warning threshold "
		f"{args.fraction_threshold:.1%}; malformed rows skipped: {malformed}."
	)
	warning = ""
	if total and fraction >= args.fraction_threshold:
		warning = (
		"WARNING: A high fraction of near-identical BLASTP best hits suggests that "
		"the protein database includes the target species itself. Homology evidence "
		"may therefore be circular and can allow overpredicted gene models to pass. "
		"This diagnostic does not exclude hits or change filtering. For the separate "
		"opt-in strict homology path, see feat/filter-evidence-quality (PR #27; "
		"--strict-single-exon)."
		)
		print(f"{summary}\n{warning}", file=sys.stderr)

	with open(args.log, "w", encoding="utf-8") as handle:
		handle.write(summary + "\n")
		if warning:
			handle.write(warning + "\n")
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
