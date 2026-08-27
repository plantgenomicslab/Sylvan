#!/usr/bin/env python3
"""Merge per-partition EVM outputs at whole-mRNA granularity.

Stock `recombine_EVM_partial_outputs.pl` assigns every gene/mRNA/CDS ROW to a
partition by its row midpoint, so a model straddling the 20 kb partition
overlap is torn apart: of 87 reference multi-exon models crossing an Ath_ODB
boundary, 10 survived intact (chain Sn 11.5%) and 68 that were correct in some
partition vanished from the merge. This module assigns each mRNA -- with its
gene row and all child rows -- to ONE owner partition instead:

    1. core interval per partition = [lo, next partition's lo - 1]
       (the last / an unsegmented partition owns through its end)
    2. a model belongs to the partition whose core contains the midpoint of
       its child-row span
    3. where the overlap holds two versions of the same locus (same
       chrom/strand, >= 50% span overlap), the copy with more child rows --
       the more complete gene structure -- wins, wherever it was computed

Input partition GFF3s are in partition-local coordinates, exactly what
`EVM_to_GFF3.pl <dir>/evm.out <chrom>` emits; coordinates are shifted by
lo - 1 here. Output IDs are prefixed with the partition basename because
EVM numbering restarts per partition (`evm.TU.<chrom>.1` exists in every
partition of that chromosome).

Validated against the Ath_ODB whole-genome experiment (remerge.py v2):
row-midpoint merge 73.30 IC F1, this merge 73.54.

Usage:
    merge_evm_partitions.py <partitions_list.out> <out.gff3>

partitions_list.out lines (from partition_EVM_inputs.pl):
    chrom<TAB>chrom_dir<TAB>Y<TAB>partition_dir   segmented chromosome
    chrom<TAB>chrom_dir<TAB>N                     unsegmented chromosome
"""
import argparse
import bisect
import os
import re
import sys
from collections import defaultdict

PART_NAME_RE = re.compile(r"(.+)_(\d+)-(\d+)$")
PARTITION_GFF = "evm.partition.gff3"
MIN_OVERLAP_FRACTION = 0.5


def _attr(attrs, key):
    match = re.search(r"(?:^|;)" + re.escape(key) + r"=([^;]+)", attrs)
    return match.group(1) if match else None


def _prefix_ids(attrs, tag):
    attrs = re.sub(r"(^|;)ID=", r"\1ID=" + tag + "_", attrs)
    return re.sub(r"(^|;)Parent=", r"\1Parent=" + tag + "_", attrs)


def read_partitions(listing_path):
    """Parse partitions_list.out into partition records."""
    parts = []
    with open(listing_path) as fh:
        for lineno, line in enumerate(fh, 1):
            fields = line.rstrip("\n").split("\t")
            if not line.strip():
                continue
            if len(fields) >= 4 and fields[2] == "Y":
                pdir = fields[3]
                name = os.path.basename(pdir.rstrip("/"))
                match = PART_NAME_RE.match(name)
                if not match:
                    raise ValueError(
                        f"{listing_path}:{lineno}: segmented partition dir "
                        f"{name!r} does not end in _<lo>-<hi>")
                lo, hi = int(match.group(2)), int(match.group(3))
            elif len(fields) >= 3:
                pdir = fields[1]
                name = os.path.basename(pdir.rstrip("/"))
                lo, hi = 1, float("inf")
            else:
                raise ValueError(
                    f"{listing_path}:{lineno}: expected 3 or 4 tab-separated "
                    f"fields, got {len(fields)}")
            parts.append({"name": name, "dir": pdir, "chrom": fields[0],
                          "lo": lo, "hi": hi, "off": lo - 1})
    return parts


def assign_cores(parts):
    """core = [lo, next partition's lo - 1]; the last one owns through hi."""
    by_chrom = defaultdict(list)
    for part in parts:
        by_chrom[part["chrom"]].append(part)
    for chrom_parts in by_chrom.values():
        chrom_parts.sort(key=lambda p: p["lo"])
        for i, part in enumerate(chrom_parts):
            part["core_lo"] = part["lo"]
            if i + 1 < len(chrom_parts):
                part["core_hi"] = chrom_parts[i + 1]["lo"] - 1
            else:
                part["core_hi"] = part["hi"]


def load_models(part, warn_counter):
    """Read one partition's GFF3 into whole-mRNA model records."""
    path = os.path.join(part["dir"], PARTITION_GFF)
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return []
    off = part["off"]
    genes, mrnas, children = {}, {}, defaultdict(list)
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                warn_counter[0] += 1
                continue
            try:
                fields[3] = str(int(fields[3]) + off)
                fields[4] = str(int(fields[4]) + off)
            except ValueError:
                warn_counter[0] += 1
                continue
            if fields[2] == "gene":
                genes[_attr(fields[8], "ID")] = fields
            elif fields[2] == "mRNA":
                mrnas[_attr(fields[8], "ID")] = fields
            else:
                parent = _attr(fields[8], "Parent")
                if parent:
                    children[parent].append(fields)

    models = []
    for mrna_id, mrna_row in mrnas.items():
        segs = children.get(mrna_id, [])
        if not segs:
            continue
        # Deterministic child order regardless of row order in the input.
        segs.sort(key=lambda s: (int(s[3]), int(s[4]), s[2], s[8]))
        lo = min(int(s[3]) for s in segs)
        hi = max(int(s[4]) for s in segs)
        gene_id = _attr(mrna_row[8], "Parent")
        models.append({
            "chrom": mrna_row[0], "strand": mrna_row[6],
            "lo": lo, "hi": hi, "mid": (lo + hi) / 2,
            "ncds": len(segs), "part": part["name"],
            "gene": genes.get(gene_id), "gene_id": gene_id,
            "mrna": mrna_row, "mrna_id": mrna_id, "children": segs,
        })
    return models


def _model_sort_key(model):
    return (model["chrom"], model["strand"], model["lo"], model["hi"],
            model["ncds"], model["part"], model["mrna_id"])


def pick_owned(models, core_of):
    kept = []
    for model in models:
        chrom, core_lo, core_hi = core_of[model["part"]]
        if model["chrom"] == chrom and core_lo <= model["mid"] <= core_hi:
            kept.append(model)
    return kept


def prefer_complete(kept, models):
    """Swap in a more complete overlapping copy from a neighbour partition.

    Candidates are windowed with bisect over lo-sorted lists: any overlapping
    candidate satisfies cand.lo >= m.lo - max_span and cand.lo <= m.hi, so the
    window sees exactly the candidates a full scan would, in the same order.
    """
    by_pos = defaultdict(list)
    for model in models:
        by_pos[(model["chrom"], model["strand"])].append(model)
    max_span = {}
    for key, group in by_pos.items():
        group.sort(key=_model_sort_key)
        group.sort(key=lambda m: m["lo"])
        max_span[key] = max(m["hi"] - m["lo"] for m in group)

    replaced = 0
    final = []
    for model in kept:
        key = (model["chrom"], model["strand"])
        group = by_pos[key]
        los = [m["lo"] for m in group]
        start = bisect.bisect_left(los, model["lo"] - max_span[key])
        best = model
        for idx in range(start, len(group)):
            cand = group[idx]
            if cand["lo"] > model["hi"]:
                break
            if cand["part"] == model["part"]:
                continue
            overlap = (min(model["hi"], cand["hi"])
                       - max(model["lo"], cand["lo"]))
            if overlap <= 0:
                continue
            if overlap / max(1, model["hi"] - model["lo"]) < MIN_OVERLAP_FRACTION:
                continue
            if cand["ncds"] > best["ncds"]:
                best = cand
        if best is not model:
            replaced += 1
        final.append(best)

    seen, deduped = set(), []
    for model in final:
        dedup_key = (model["chrom"], model["strand"], model["lo"],
                     model["hi"], model["ncds"])
        if dedup_key in seen:
            continue
        seen.add(dedup_key)
        deduped.append(model)
    return deduped, replaced, len(final) - len(deduped)


def emit(kept, out_path):
    kept = sorted(kept, key=lambda m: (m["chrom"], m["lo"], m["hi"],
                                       m["strand"], m["part"], m["mrna_id"]))
    emitted_genes = set()
    with open(out_path, "w") as out:
        for model in kept:
            tag = model["part"]

            def write(row, tag=tag):
                fields = list(row)
                fields[8] = _prefix_ids(fields[8], tag)
                out.write("\t".join(fields) + "\n")

            gene_key = (tag, model["gene_id"])
            if model["gene"] is not None and gene_key not in emitted_genes:
                emitted_genes.add(gene_key)
                write(model["gene"])
            write(model["mrna"])
            for child in model["children"]:
                write(child)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n", maxsplit=1)[0])
    ap.add_argument("partitions_list")
    ap.add_argument("outfile")
    args = ap.parse_args(argv)

    parts = read_partitions(args.partitions_list)
    assign_cores(parts)
    core_of = {p["name"]: (p["chrom"], p["core_lo"], p["core_hi"])
               for p in parts}

    warn_counter = [0]
    models = []
    for part in parts:
        models.extend(load_models(part, warn_counter))
    # Deterministic regardless of row order inside partition files.
    models.sort(key=_model_sort_key)

    kept = pick_owned(models, core_of)
    kept, replaced, deduped = prefer_complete(kept, models)
    emit(kept, args.outfile)

    if warn_counter[0]:
        print(f"merge_evm_partitions: skipped {warn_counter[0]} malformed "
              "rows", file=sys.stderr)
    print(f"merge_evm_partitions: partitions {len(parts)}, models "
          f"{len(models)}, kept {len(kept)} (replaced {replaced}, "
          f"deduped {deduped}) -> {args.outfile}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
