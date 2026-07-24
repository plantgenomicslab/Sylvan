#!/usr/bin/env python3
"""
Score gene models using simple logistic regression and random forest.

Features:
- pfam_hit: 1 if PfamScan reports a domain
- homolog_hit: 1 if BLASTP to homolog DB hits
- rex_hit: 1 if BLASTP to RexDB hits
- tpm: TPM from RSEM (isoforms.results)
- cov_aug/helixer/repeat: fraction coverage from bedtools coverage outputs

Pseudo-labels (weak supervision):
- positive if pfam_hit or homolog_hit
- negative if rex_hit
Only entries with a pseudo-label are used for training.

Threshold selection:
- For each model, pick threshold maximizing F1 on pseudo-labels (computed from PR curve).

Outputs:
- scores.csv: per-gene scores and features
- scores.metrics.txt: AUC/PR/F1 and chosen thresholds
"""
import argparse
import os
import sys
import pandas as pd
import numpy as np

from sklearn.base import clone
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, precision_recall_curve, average_precision_score
from sklearn.model_selection import cross_val_predict, StratifiedKFold

# Features that deterministically DEFINE the pseudo-labels (label_pos = pfam|homolog,
# label_neg = rex). Feeding them to the model lets it memorise the labelling rule
# (train AUC ~= 1, useless generalisation), so they are excluded from X.
LEAKAGE_FEATURES = ["pfam_hit", "homolog_hit", "rex_hit"]
# Evidence the model is actually allowed to learn from.
MODEL_FEATURES = [
    "tpm",
    "cov_aug",
    "cov_helixer",
    "cov_repeat",
    "lnc_prob",
    "aa_len",
    "internal_stop",
    "short_orf",
]


def load_pfam(path):
    if not os.path.exists(path):
        return pd.DataFrame(columns=["gene_id", "pfam_hit"])
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    # PfamScan output: query id is column 0
    df["gene_id"] = df.iloc[:, 0].astype(str)
    return df[["gene_id"]].assign(pfam_hit=1).drop_duplicates()


def load_blast(path, col=0, flag_name="hit"):
    if not os.path.exists(path):
        return pd.DataFrame(columns=["gene_id", flag_name])
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    df["gene_id"] = df.iloc[:, col].astype(str)
    return df[["gene_id"]].assign(**{flag_name: 1}).drop_duplicates()


def load_rsem(path):
    if not os.path.exists(path):
        return pd.DataFrame(columns=["gene_id", "tpm"])
    df = pd.read_csv(path, sep="\t")
    # isoforms.results has "transcript_id" or "transcript_id(s)"
    if "transcript_id" in df.columns:
        id_col = "transcript_id"
    else:
        id_col = df.columns[0]
    tpm_col = "TPM" if "TPM" in df.columns else df.columns[-1]
    df["gene_id"] = df[id_col].astype(str)
    return df[["gene_id", tpm_col]].rename(columns={tpm_col: "tpm"})


def load_cov(path):
    if not os.path.exists(path):
        return pd.DataFrame(columns=["gene_id", "cov"])
    # The prefilter BED handed to `bedtools coverage -a` has 5 columns
    # (chrom, start, end, strand, transcript_id); bedtools then appends
    # (num_overlaps, bp_covered, feature_len, frac_covered). So the gene ID is
    # column index 4 and the coverage fraction is column index 8 -- matching
    # Filter.py's read_abInitio_cov (usecols=[4, 8]). The previous code read
    # column index 3 (strand) as the gene ID, so the join key never matched and
    # cov_aug/cov_helixer/cov_repeat were all NaN -> 0.
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    df["gene_id"] = df.iloc[:, 4].astype(str)
    frac = df.iloc[:, 8]
    return df[["gene_id"]].assign(cov=frac.values)


def select_threshold(y_true, scores):
    precision, recall, thresholds = precision_recall_curve(y_true, scores)
    f1 = 2 * precision * recall / (precision + recall + 1e-9)
    best_idx = f1.argmax()
    thresh = thresholds[best_idx] if best_idx < len(thresholds) else 0.5
    return thresh, precision[best_idx], recall[best_idx], f1[best_idx]


def load_orf_flags(pep_path, len_cutoff=100):
    """
    Simple ORF quality checks from peptide FASTA:
    - aa_len
    - internal_stop flag (any '*' before the end)
    - short_orf flag (length < len_cutoff)
    """
    if not pep_path or not os.path.exists(pep_path):
        return pd.DataFrame(columns=["gene_id", "aa_len", "internal_stop", "short_orf"])
    records = []
    gene_id = None
    seq_chunks = []
    with open(pep_path) as f:
        for line in f:
            if line.startswith(">"):
                if gene_id is not None:
                    seq = "".join(seq_chunks).strip()
                    aa_len = len(seq.replace("*", ""))
                    internal_stop = 1 if "*" in seq[:-1] else 0
                    short_orf = 1 if aa_len < len_cutoff else 0
                    records.append((gene_id, aa_len, internal_stop, short_orf))
                gene_id = line[1:].strip().split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if gene_id is not None:
            seq = "".join(seq_chunks).strip()
            aa_len = len(seq.replace("*", ""))
            internal_stop = 1 if "*" in seq[:-1] else 0
            short_orf = 1 if aa_len < len_cutoff else 0
            records.append((gene_id, aa_len, internal_stop, short_orf))
    return pd.DataFrame(records, columns=["gene_id", "aa_len", "internal_stop", "short_orf"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pfam", required=True)
    ap.add_argument("--blast", required=True)
    ap.add_argument("--blast_rex", required=True)
    ap.add_argument("--rsem", required=True)
    ap.add_argument("--cov_aug", required=True)
    ap.add_argument("--cov_helixer", required=True)
    ap.add_argument("--cov_repeat", required=True)
    ap.add_argument("--lnc_pred", required=False, default=None)
    ap.add_argument("--pep", required=False, default=None)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--out_metrics", required=True)
    args = ap.parse_args()

    pfam = load_pfam(args.pfam)
    blast = load_blast(args.blast, flag_name="homolog_hit")
    blast_rex = load_blast(args.blast_rex, flag_name="rex_hit")
    rsem = load_rsem(args.rsem)
    cov_aug = load_cov(args.cov_aug).rename(columns={"cov": "cov_aug"})
    cov_helixer = load_cov(args.cov_helixer).rename(columns={"cov": "cov_helixer"})
    cov_repeat = load_cov(args.cov_repeat).rename(columns={"cov": "cov_repeat"})
    # LncDC predictions: CSV with gene IDs and probabilities/labels
    if args.lnc_pred and os.path.exists(args.lnc_pred):
        try:
            lnc = pd.read_csv(args.lnc_pred)
            # Accept columns: gene_id, prob or label
            if "gene_id" not in lnc.columns:
                lnc.rename(columns={lnc.columns[0]: "gene_id"}, inplace=True)
            if "prob" not in lnc.columns:
                # best-effort: take last column as probability if numeric
                for col in lnc.columns[::-1]:
                    if np.issubdtype(lnc[col].dtype, np.number):
                        lnc["prob"] = lnc[col]
                        break
            lnc = lnc[["gene_id"] + ([c for c in ["prob"] if "prob" in lnc.columns])].copy()
            lnc.rename(columns={"prob": "lnc_prob"}, inplace=True)
        except Exception:
            lnc = pd.DataFrame(columns=["gene_id", "lnc_prob"])
    else:
        lnc = pd.DataFrame(columns=["gene_id", "lnc_prob"])

    orf_flags = load_orf_flags(args.pep)

    # Merge features
    dfs = [pfam, blast, blast_rex, rsem, cov_aug, cov_helixer, cov_repeat, lnc, orf_flags]
    feats = dfs[0]
    for df in dfs[1:]:
        feats = feats.merge(df, on="gene_id", how="outer")

    feats = feats.fillna(
        {
            "pfam_hit": 0,
            "homolog_hit": 0,
            "rex_hit": 0,
            "tpm": 0,
            "cov_aug": 0,
            "cov_helixer": 0,
            "cov_repeat": 0,
            "lnc_prob": 0,
            "aa_len": 0,
            "internal_stop": 0,
            "short_orf": 0,
        }
    )

    # Quantile-based adaptive cutoffs
    q_tpm_keep = feats["tpm"].quantile(0.05) if len(feats) else 0
    q_cov_aug_keep = feats["cov_aug"].quantile(0.05) if len(feats) else 0
    q_cov_helixer_keep = feats["cov_helixer"].quantile(0.05) if len(feats) else 0
    q_cov_repeat_high = feats["cov_repeat"].quantile(0.95) if len(feats) else 1
    # Apply lower bounds so very low-coverage datasets don't set near-zero thresholds
    q_tpm_keep = max(q_tpm_keep, 1.0)
    q_cov_aug_keep = max(q_cov_aug_keep, 0.05)
    q_cov_helixer_keep = max(q_cov_helixer_keep, 0.05)

    exp_pass = feats["tpm"] >= q_tpm_keep
    cov_pass = (feats["cov_aug"] >= q_cov_aug_keep) | (feats["cov_helixer"] >= q_cov_helixer_keep)
    lnc_low = feats["lnc_prob"] >= 0.6  # more aggressive: drop if lncRNA prob is moderate-high
    pseudo_flag = (feats["internal_stop"] == 1) | (feats["short_orf"] == 1)
    # If expression is low, require Pfam OR homolog; expression+coverage only helps when Pfam/homolog evidence exists.
    soft_keep = (feats["pfam_hit"] == 1) | (feats["homolog_hit"] == 1) | ((exp_pass & cov_pass) & ((feats["pfam_hit"] == 1) | (feats["homolog_hit"] == 1)))
    soft_drop = ((feats["rex_hit"] == 1) & (feats["cov_repeat"] >= q_cov_repeat_high)) | lnc_low | pseudo_flag
    feats["soft_vote"] = np.where(soft_drop, "drop", np.where(soft_keep, "keep", "unsure"))

    # Pseudo-labels
    feats["label_pos"] = (feats["pfam_hit"] + feats["homolog_hit"] > 0).astype(int)
    feats["label_neg"] = (feats["rex_hit"] > 0).astype(int)
    # Use only rows where label is defined (pos or neg)
    has_label = (feats["label_pos"] == 1) | (feats["label_neg"] == 1)
    train = feats[has_label].copy()
    train["label"] = train["label_pos"]

    metrics_lines = []
    if len(train) < 10 or train["label"].nunique() < 2:
        feats["logit_score"] = 0.0
        feats["rf_score"] = 0.0
        metrics_lines.append("Insufficient labels for training; scores set to 0.")
    else:
        # Exclude the label-defining features (pfam/homolog/rex) so the model
        # learns from independent evidence instead of memorising the rule.
        X_train = train[MODEL_FEATURES]
        y_train = train["label"]

        models = {
            "logit": LogisticRegression(max_iter=1000),
            "rf": RandomForestClassifier(
                n_estimators=200,
                max_depth=None,
                n_jobs=-1,
                class_weight="balanced",
                random_state=42,
            ),
        }

        # Metrics on the same rows a model trained on are optimistic (AUC ~= 1).
        # Estimate generalisation with out-of-fold cross-validated predictions;
        # fall back to in-sample only when there are too few minority labels to
        # split into folds.
        min_class = int(y_train.value_counts().min())
        n_splits = min(5, min_class)
        can_cv = n_splits >= 2

        for name, est in models.items():
            # Fit on the full labelled set to score EVERY gene for the CSV.
            est.fit(X_train, y_train)
            feats[f"{name}_score"] = est.predict_proba(feats[MODEL_FEATURES])[:, 1]

            if can_cv:
                oof = cross_val_predict(
                    clone(est),
                    X_train,
                    y_train,
                    cv=StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42),
                    method="predict_proba",
                )[:, 1]
                eval_scores = pd.Series(oof, index=y_train.index)
                note = f"(out-of-fold {n_splits}-fold CV)"
            else:
                eval_scores = feats.loc[has_label, f"{name}_score"]
                note = "(in-sample; too few minority labels for CV)"

            auc_roc = roc_auc_score(y_train, eval_scores)
            # average_precision_score is the correct PR-AUC. The previous
            # auc(precision, recall) integrated precision (non-monotonic) as the
            # x-axis, which auc() is not defined for -> garbage/negative values.
            pr_auc = average_precision_score(y_train, eval_scores)
            thresh, p, r, f1 = select_threshold(y_train, eval_scores)
            metrics_lines.append(
                f"{name}: ROC_AUC={auc_roc:.4f} PR_AUC={pr_auc:.4f} "
                f"threshold={thresh:.4f} precision={p:.4f} recall={r:.4f} F1={f1:.4f} {note}"
            )
            feats[f"{name}_threshold"] = thresh

    metrics_lines.append(
        f"quantiles: tpm_keep(p05)={q_tpm_keep:.4f} "
        f"cov_aug_keep(p05)={q_cov_aug_keep:.4f} cov_helixer_keep(p05)={q_cov_helixer_keep:.4f} "
        f"cov_repeat_high(p95)={q_cov_repeat_high:.4f} "
        f"lnc_drop_threshold=0.6 "
        f"pseudo: short_orf<len 100 aa or internal_stop"
    )

    # Save outputs
    feats[
        [
            "gene_id",
            "pfam_hit",
            "homolog_hit",
            "rex_hit",
            "tpm",
            "cov_aug",
            "cov_helixer",
            "cov_repeat",
            "lnc_prob",
            "aa_len",
            "internal_stop",
            "short_orf",
            "soft_vote",
            "logit_score",
            "rf_score",
        ]
    ].to_csv(args.out_csv, index=False)
    with open(args.out_metrics, "w") as f:
        for line in metrics_lines:
            f.write(line + "\n")


if __name__ == "__main__":
    main()
