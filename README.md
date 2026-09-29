# PRC2 base-editor screen — analysis notebooks

Consolidates the analysis previously split across
`251122_bescan_Plots_Github.ipynb` and `figures_ycm_xyh_v2.ipynb` into one
pipeline, built on [`be_scan`](https://github.com/liaulab/be-scan).

**Status:** both notebooks are built and execute end-to-end with zero
errors. `master_analysis.ipynb`'s static outputs were validated against
the old figures during migration (see `errors.md` for the full record and
current open items — most validation history has been condensed out of
that file; ask if you want it back). Both notebooks currently have their
cell outputs cleared (no baked-in plots) — run them fresh to regenerate
figures.

## Layout

- `master_analysis.ipynb` — the publication-facing notebook. Static
  (matplotlib/seaborn, via `be_scan`) figures only. Loess plots are
  intentionally excluded.
- `interactive_analysis.ipynb` — the same analyses as interactive Plotly
  HTML instead of static SVGs, via
  `be_scan.figure_plot.figure_interactive`. Not for publication — feeds
  `docs/`. Produces 75 HTML figures under
  `new_outputs_interactive/<Category>/` (34 Scatterplots, 12 Boxplots, 20
  Lollipops, 2 ClusteredHeatmaps, 7 PWES). Several plot categories were
  trimmed down from earlier, larger sets per Calvin's requests (predictor-
  based FoldX/ThermoMPNN lollipops and scatterplots removed, several
  redundant boxplot/heatmap/scatterplot variants dropped) — see `errors.md`
  for the full trim history. A real `be_scan` bug (interactive lollipops'
  negative-threshold hit filter never actually filtered) was found and
  monkey-patched here; only affected the interactive notebook.
- `code/` — shared Python modules, styled to match `be_scan.figure_plot`:
  - `config.py` — every shared constant: gene lists, colors, **input file
    paths (all under `inputs/`, see below)**, the timepoint-selection
    table, and k-means/PWES parameters (including `KMEANS_RANDOM_STATE`,
    chosen by grid search — see `errors.md` for the achieved match and its
    limits). Both notebooks import from here so they can't drift apart.
  - `data_loading.py` — reshapes `PRC2_Screens_Counts_Scores_ZScores.xlsx`
    and the companion mapping/cluster workbooks into the long-format
    tables the plotting code expects. Column mappings were validated
    against the reference workbooks; two real bugs (a Q575R column swap,
    an MTF2 position-column error) were found this way — see `errors.md`.
  - `pymol_export.py` — residue-score → PyMOL b-factor `.txt` export.
  - `structural_enrichment.py` — odds-ratio/forest-plot/histogram/
    substitution-matrix/quantile-boxplot/violin-plot functions, ported
    from `figures_ycm_xyh_v2.ipynb` and restyled to `be_scan` conventions.
  - `pwes_significance.py` — permutation-test significance calculation for
    PWES clusters, rebuilt from the old output files and verified to
    reproduce them.
- `inputs/` — every raw data file either notebook reads, consolidated into
  one flat folder (2026-09-29 cleanup; everything unrelated to this
  analysis was deleted, and the old `FoldX/`/`YoochanStuff/` subfolders
  were flattened into it): `PRC2_Screens_Counts_Scores_ZScores.xlsx`,
  `PRC2_Screens_AA_Mapping.xlsx`, `PRC2_Screens_KMeans_Clusters.xlsx`,
  `PRC2_Screens_PWES_Clusters.xlsx`, `mutation_level_data.tsv`, `6wkr.pdb`,
  `sgRNA_level_{ABE,CBE}_251112.tsv`, and `EZH2.tsv`/`EED.tsv`/`SUZ12.tsv`
  (FoldX predictions, three genes only). `code/config.py` is the only file
  with a literal path in it — if this folder ever moves, only
  `INPUTS_DIR` there needs to change.
- `environment.yml` — dependencies; the first cell of both notebooks
  installs from this on every run (Colab has no persistent environment).
- `errors.md` — every bug, inconsistency, and validation result found
  while migrating the two original notebooks, and how each was resolved
  (or why it's still open). Read it before trusting any specific figure.
- `new_outputs/` — figures from `master_analysis.ipynb` (Scatterplots,
  Boxplots, Lollipops, ClusteredHeatmaps, Stability, PWES [combined
  LOF+Stability run only], PyMOL `.txt` exports).
- `docs/` — the interactive data page built from
  `new_outputs_interactive/` (named `docs/`, not `webpage/`, so it can be
  served directly as the repo's GitHub Pages source — main branch,
  `/docs`): a static `index.html` (light/dark toggle, sidebar of
  figure-type categories, one shared legend bar per group of plots
  instead of a legend on every plot) plus `manifest.json` and `data/` (a
  processed copy of every interactive figure — autosized to fill its
  card, with domain-highlight background bands fixed to span the full
  plot height). `Scatterplots` splits into 5 sidebar entries (WT, Q575R,
  EZH2/EED Stability, Correlation) all backed by one `data/Scatterplots/`
  folder. No build step needed to view it — see `docs/README.md` for
  regenerating `data/`/`manifest.json` after a fresh notebook run, or
  deploying via GitHub Pages.
- No longer here (moved or deleted in the 2026-09-29 cleanup):
  `previous_results/` (the old-figure comparison target used during
  validation — no longer needed to run the pipeline; validation results
  are recorded in `errors.md`), `load_bfact.py`, `COSMIC/`, and the two
  original pre-migration notebooks.

## Figure categories (in `master_analysis.ipynb`'s order)

1. Scatterplots + correlation scatterplots — WT/Q575R/MTF2/Stability
   dropout & abundance Z-scores, per gene and pairwise correlations.
2. Boxplots / strip plots — editor × gene and editor × mutation-type
   distributions.
3. Lollipop plots — per-residue LOF hits and LOF-vs-stability comparisons
   (screen data only). Ordered EZH2, EED, SUZ12, AEBP2, JARID2.
4. Clustered heatmaps — k-means clustering of multi-screen hit patterns
   (ABE and CBE, k=6).
5. Stability-related analysis — clustered heatmaps, forest plots, violin
   plots, strip plots of structural predictors vs. abundance.
6. PWES section — 3D-proximity clustering against the 6WKR structure
   (scatterplots, heatmap, boxplots). Only the combined LOF+Stability run
   is kept.
7. PyMOL export — residue-level scores → per-residue `.txt` files for
   b-factor coloring, run separately inside a PyMOL session (not called
   from the notebook).

## Setup

Both notebooks run fresh in Google Colab every time (no persistent
environment). The first cell installs everything in `environment.yml`,
including `be_scan` from GitHub, and (for `master_analysis.ipynb`)
fetches `Arial.ttf`. Both expect to run with this folder as the working
directory, with `inputs/` present alongside them.

## Known open issues

See `errors.md` for the full list. Still needs a decision:
- Two Q575R-CBE39 timepoints (`PRC2i-D0`, `DMSO-D0`) are unconfirmed —
  every other ambiguous condition is now confirmed against an old figure;
  these two have no vector reference to check against (`errors.md` E-01).
- K-means cluster labels are the best achievable match (ARI 0.93 ABE /
  0.91 CBE against the reference workbook), not exact — a scikit-learn
  version-level reproducibility limit. Use
  `PRC2_Screens_KMeans_Clusters.xlsx` directly if exact published labels
  matter more than a fresh clustering.
- One ABE cluster subset heatmap and a handful of old CBE lollipop plots
  can't be reproduced exactly — both are explained (not bugs here) in
  `errors.md`.
