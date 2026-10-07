# Base editing charts a sequence-function atlas of Polycomb Repressive Complex 2 at amino acid resolution. 

## Layout

- `master_analysis.ipynb` — Static `be_scan` figures only. 
- `interactive_analysis.ipynb` — The same analyses as interactive Plotly HTML instead of static SVGs, via `be_scan.figure_plot.figure_interactive`. Not for publication. 
- `code/` shared Python modules: 
  - `config.py` shared constants. 
  - `data_loading.py` formats `PRC2_Screens_Counts_Scores_ZScores.xlsx`
  - `pymol_export.py` maps residue-score PyMOL to b-factor `.txt`.
  - `structural_enrichment.py` functions for odds-ratio/forest-plot/histogram/substitution-matrix/quantile-boxplot/violin-plot.
  - `pwes_significance.py` permutation significance test for PWES clusters. 
- `inputs/` raw data files consolidated into one folder. `code/config.py` is the only file with paths, if this folder ever moves only `INPUTS_DIR` needs to change. This folder is not included in this repo until the paper is published. 
- `environment.yml` dependencies. Both notebooks installs from this on every run (Colab has no persistent environment).
- `errors.md` bugs and inconsistency. 
- `docs/` — the interactive data page deploying via GitHub Pages.

## Figure categories (in `master_analysis.ipynb`'s order)

1. Scatterplots and correlation plots. 
2. Boxplots and strip plots. 
3. Lollipop plots. 
4. K-means clustered heatmaps. 
5. Stability-related analysis (clustered heatmaps, forest plots, violin
   plots, strip plots of structural predictors vs. abundance). 
6. PWES against the 6WKR structure. 
7. PyMOL mapping scores to per-residue `.txt` files for b-factor coloring. Run separately inside a PyMOL session. 

## Setup

Both notebooks run fresh in Google Colab every time. The first cell installs everything in `environment.yml` including `be_scan` from GitHub. 
