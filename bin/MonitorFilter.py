#!/usr/bin/env python3
"""Plot semi-supervised random-forest training diagnostics across iterations.

Reads the ``process`` dictionary emitted by ``Filter.semiSupRandomForest``
(keys: ``kept``, ``discarded``, ``kept_buscos``, ``discarded_buscos``, ``OOB``
-- each a list with one entry per training iteration) from a JSON file and
writes a multi-panel figure.

Previously this was a notebook leftover: it plotted a hard-coded demo dict at
import time and called ``plt.show()``, which hangs on a headless cluster/CI
node (no ``$DISPLAY``) and rendered nothing useful.
"""
import argparse
import json
import sys

import matplotlib
matplotlib.use("Agg")  # headless-safe: no DISPLAY on the cluster or in CI
import matplotlib.pyplot as plt  # noqa: E402  (must follow matplotlib.use)
from matplotlib.ticker import FuncFormatter  # noqa: E402

# (process key -> human-readable panel title), in display order.
PANELS = [
	("OOB", "OOB Error"),
	("kept", "Kept Genes"),
	("discarded", "Discarded Genes"),
	("kept_buscos", "Kept BUSCOs"),
	("discarded_buscos", "Discarded BUSCOs"),
]


def format_with_commas(value, _pos):
	"""Thousands separators for counts; trimmed decimals for the OOB fraction."""
	if value >= 1:
		return "{:,.0f}".format(value)
	return "{:.10f}".format(value).rstrip("0")


def select_panels(process):
	"""Return the (key, title) panels that actually have data in ``process``."""
	return [(key, title) for key, title in PANELS if process.get(key)]


def plot_process(process, out_path):
	"""Render one line panel per populated series to ``out_path``."""
	panels = select_panels(process)
	if not panels:
		raise SystemExit("No plottable series found in the process data.")

	fig, axs = plt.subplots(3, 2, figsize=(8, 7))
	axs = axs.flatten()
	for i, (key, title) in enumerate(panels):
		values = process[key]
		ax = axs[i]
		ax.plot(range(1, len(values) + 1), values, marker="o")
		ax.set_xlabel("Iteration")
		ax.yaxis.set_major_formatter(FuncFormatter(format_with_commas))
		ax.set_title(title)

	# Drop any unused panels so the grid does not show empty axes.
	for j in range(len(panels), len(axs)):
		fig.delaxes(axs[j])

	fig.tight_layout()
	fig.savefig(out_path)
	print(f"Wrote {out_path}", file=sys.stderr)


def main():
	ap = argparse.ArgumentParser(description=__doc__)
	ap.add_argument("process_json", help="JSON file with the semiSupRandomForest process dict")
	ap.add_argument("-o", "--out", default="filter_iterations.png", help="Output image path")
	args = ap.parse_args()
	with open(args.process_json) as fh:
		process = json.load(fh)
	plot_process(process, args.out)


if __name__ == "__main__":
	main()
