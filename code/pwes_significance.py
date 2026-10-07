import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def pwes_permutation_test(df_pairwise, df_gauss, df_clus, label_col="label",
                          cluster_col="cl_new", n_perm=10000):
    """
    Run the cluster-WAP permutation test.

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
    """
    Grid of null-WAP histograms with the observed WAP as a red dashed line
    """
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
