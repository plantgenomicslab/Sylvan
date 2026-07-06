#!/usr/bin/env python3
"""Curate filtered.gff3 -> keep RF-Keep genes + evidenced RF-None genes.

Reproduces (approximately) the historical Sylvan.curated.gff3 reduction. The exact
historical procedure is not recorded; this is a documented label-based approximation.
"""
import argparse, re, sys
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--filtered", required=True)
    ap.add_argument("--features", required=True)   # keep_data.tsv
    ap.add_argument("--out", required=True)
    ap.add_argument("--cov-thresh", type=float, default=0.0)
    ap.add_argument("--evidence", default="HELIXER,MINIPROT")  # cols that, if present, rescue a None gene
    ap.add_argument("--report", default="")
    a = ap.parse_args()

    df = pd.read_csv(a.features, sep="\t", dtype=str).fillna("")
    ev_cols = [c for c in a.evidence.split(",") if c]
    def truthy(v):  # "1.0"/"1"/"True" -> True
        return str(v).strip() in ("1", "1.0", "True", "true")
    is_keep = df["label"] == "Keep"
    none_mask = df["label"] == "None"
    ev = pd.Series(False, index=df.index)
    for c in ev_cols:
        ev = ev | df[c].map(truthy)
    if a.cov_thresh > 0:
        ev = ev | (pd.to_numeric(df["COVERAGE"], errors="coerce").fillna(0) > a.cov_thresh)
    keep_tx = set(df.loc[is_keep | (none_mask & ev), "transcript_id"])
    # mRNA evm.model.X.N -> gene evm.TU.X.N
    keep_genes = {t.replace("evm.model.", "evm.TU.", 1) for t in keep_tx}

    def gid_of(attr, ftype):
        if ftype == "gene":
            m = re.search(r'ID=([^;]+)', attr);  return m.group(1) if m else None
        m = re.search(r'Parent=([^;]+)', attr)   # mRNA->gene ; lower features handled by mRNA gate below
        return m.group(1) if m else None

    # stream: keep a gene block iff gene id in keep_genes. Track current gene via Parent chain.
    id_parent = {}
    lines = []
    with open(a.filtered) as fh:
        for ln in fh:
            if ln.startswith("#"):
                lines.append(("hdr", ln, None)); continue
            c = ln.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            fid = (re.search(r'ID=([^;]+)', c[8]) or [None, None])[1] if "ID=" in c[8] else None
            par = (re.search(r'Parent=([^;]+)', c[8]) or [None, None])[1] if "Parent=" in c[8] else None
            if fid: id_parent[fid] = par
            lines.append(("feat", c, (fid, par)))

    def root_gene(fid, par, ftype):
        if ftype == "gene":
            return fid
        cur = par
        for _ in range(20):
            if cur is None: return None
            nxt = id_parent.get(cur)
            if nxt is None: return cur
            cur = nxt
        return cur

    kept_genes = 0
    with open(a.out, "w") as out:
        for kind, payload, ids in lines:
            if kind == "hdr":
                out.write(payload); continue
            c = payload; fid, par = ids
            g = root_gene(fid, par, c[2])
            if g in keep_genes:
                out.write("\t".join(c) + "\n")
                if c[2] == "gene": kept_genes += 1
    sys.stderr.write(f"[curate_genes] kept genes: {kept_genes}\n")
    if a.report:
        with open(a.report, "w") as r:
            r.write(f"kept_genes\t{kept_genes}\n")

if __name__ == "__main__":
    main()
