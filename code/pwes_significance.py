"""
Permutation-significance test for PWES clusters (the `df_wap` /
`df_pwes_shuffled_sums` / `df_summary` / `df_summary_revised` /
`WAP_histograms.png` / `aas_dict_FullPDB.pml` outputs of the old
Figure4_PWES_* folders).

Why this module exists
-----------------------
251122_bescan_Plots_Github.ipynb ("NB1") does not contain the code that
produced these files -- they exist only as outputs in
previous_results/Figure4_PWES_*/. The definitions below were
reverse-engineered from those outputs and reproduce them exactly for the
LOF+Stability run (see errors.md, "Output-parity pass"):

  obs_wap(cluster)  = sum over residue pairs i<j inside the cluster of
                      |PWES_ij|, where PWES = df_pairwise * df_gauss
                      (i.e. df_pwes_unsorted from be_scan.pwes_clustering).
  null iteration i  = np.random.seed(i); p = np.random.permutation(n);
                      PWES* = df_pairwise[p][:, p] * df_gauss, then the same
                      cluster sum. This is exactly be_scan's shuffle_pwes()
                      permutation scheme, run for n_perm iterations.
  mean / std        = mean and sample std (ddof=1) of the null.
  95ci_min / max    = 2.5th / 97.5th percentiles of the null.
  obs_gt            = number of null iterations strictly greater than obs.
  1t_pval           = obs_gt / n_perm.
  df_summary_revised adds a Bonferroni test across clusters:
  "Signficiant at 0.05" (sic, the old column name) = p * n_clusters < 0.05,
  and "Signficiance" (sic) stars from the Bonferroni-adjusted p:
  **** < 1e-4, *** < 1e-3, ** < 1e-2, * < 0.05, else ns. (Only "****",
  "*" and "ns" occur in the old files, so the ** / *** cut-offs are the
  conventional ones, not verified.)
  df_pwes_shuffled_sums = element-wise mean of the n_perm shuffled PWES
                      matrices, in df_pwes_unsorted row/column order.
"""

import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def pwes_permutation_test(df_pairwise, df_gauss, df_clus, label_col="label",
                          cluster_col="cl_new", n_perm=10000):
    """Run the cluster-WAP permutation test.

    df_pairwise : the `_df_pairwise.csv` frame from pwes_clustering (rows and
                  columns in df_clus[label_col] order; its own index may be
                  the internal sgRNA_* ids).
    df_gauss    : the `_df_gauss.csv` frame (indexed by residue label).
    df_clus     : the `_df_clus.csv` frame (one row per residue, in the same
                  order as df_pairwise).

    Returns (df_wap, df_shuffled_mean, df_summary, df_summary_revised).
    """
    labels = df_clus[label_col].tolist()
    pair = np.asarray(df_pairwise.values, dtype=float)
    gauss = df_gauss.loc[labels, labels].values.astype(float)
    n = len(labels)
    assert pair.shape == (n, n), "df_pairwise must be in df_clus residue order"

    clusters = sorted(df_clus[cluster_col].unique())
    lab = df_clus[cluster_col].values
    onehot = np.stack([(lab == c).astype(float) for c in clusters], axis=1)  # n x k
    def cluster_wap(pwes):
        # PWES is symmetric, so the i<j sum is (full block sum - diagonal) / 2.
        a = np.abs(pwes)
        full = ((a @ onehot) * onehot).sum(axis=0)
        diag = (np.diagonal(a)[:, None] * onehot).sum(axis=0)
        return (full - diag) / 2.0

    obs = cluster_wap(pair * gauss)
    null = np.empty((len(clusters), n_perm))
    shuffled_sum = np.zeros((n, n))
    for i in range(n_perm):
        np.random.seed(i)
        p = np.random.permutation(n)
        s = pair[np.ix_(p, p)] * gauss
        shuffled_sum += s
        null[:, i] = cluster_wap(s)

    df_wap = pd.DataFrame(null, columns=[f"iter{i}" for i in range(n_perm)])
    df_wap.insert(0, "obs_wap", obs)
    df_wap.insert(0, "cluster", clusters)

    df_shuffled = pd.DataFrame(shuffled_sum / n_perm, index=labels, columns=labels)

    obs_gt = (null > obs[:, None]).sum(axis=1)
    df_summary = pd.DataFrame({
        "cluster": clusters,
        "obs_wap": obs,
        "mean": null.mean(axis=1),
        "std": null.std(axis=1, ddof=1),
        "95ci_min": np.percentile(null, 2.5, axis=1),
        "95ci_max": np.percentile(null, 97.5, axis=1),
        "obs_gt": obs_gt,
        "1t_pval": obs_gt / n_perm,
    })
    df_rev = df_summary.copy()
    adj = df_rev["1t_pval"] * len(clusters)
    df_rev["Signficiant at 0.05"] = adj < 0.05
    df_rev["Signficiance"] = np.select(
        [adj < 1e-4, adj < 1e-3, adj < 1e-2, adj < 0.05],
        ["****", "***", "**", "*"], default="ns")
    return df_wap, df_shuffled, df_summary, df_rev


def plot_wap_histograms(df_wap, out_png, ncols=3, bins=30):
    """Grid of null-WAP histograms with the observed WAP as a red dashed line
    (layout of the old *_WAP_histograms.png)."""
    iters = df_wap.filter(like="iter").values
    k = len(df_wap)
    nrows = math.ceil(k / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(20, 20 * nrows / 12))
    axes = np.atleast_1d(axes).ravel()
    for i, (cl, obs) in enumerate(zip(df_wap["cluster"], df_wap["obs_wap"])):
        ax = axes[i]
        ax.hist(iters[i], bins=bins, color="blue", alpha=0.7, edgecolor="black")
        ax.axvline(obs, color="red", ls="--", lw=1, label="Observed WAP")
        ax.set_title(f"Cluster: {cl}"); ax.set_xlabel("WAP Value"); ax.set_ylabel("Frequency")
        lo, hi = min(iters[i].min(), obs), max(iters[i].max(), obs)
        pad = 0.03 * (hi - lo)
        ax.set_xlim(lo - pad, hi + pad)
        ax.legend()
    for ax in axes[k:]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)
