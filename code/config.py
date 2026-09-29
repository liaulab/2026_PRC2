"""
Shared configuration for the PRC2 base-editor screen analysis notebooks
(master_analysis.ipynb and interactive_analysis.ipynb).

Every notebook cell that needs a gene list, a color, a file path, or a
condition name should import from here rather than re-declaring it, so
that both notebooks (and any script in code/) stay in sync.

NOTE ON TIMEPOINT SELECTION (read this before trusting any WT/Q575R output)
----------------------------------------------------------------------------
PRC2_Screens_Counts_Scores_ZScores.xlsx stores more than one timepoint's
Z-score for some conditions. Calvin was not sure, as of the initial
migration pass, which timepoint the manuscript figures used for:
  - WT screen, ABE editor      (Wk3 vs Wk4 available)
  - WT screen, CBE6b editor    (Wk3 vs Wk6 available)
  - Q575R screen, CBE39 editor (Wk3 vs Wk4 available)
This file defaults all three to Wk3 as a PLACEHOLDER. This default is
almost certainly wrong for at least one of them and MUST be confirmed
against the manuscript / the person who generated the original Table_S1
export before these results are treated as final. See errors.md, item
E-01, for full detail. To change it, edit TIMEPOINTS below -- every
loader function in data_loading.py reads from this dict, so one edit
fixes every downstream plot.
"""

from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
# All paths are relative to the notebook's working directory, which should
# be the RevisionAnalysis/ folder itself (both on a local machine and after
# the Colab setup cell clones/mounts this folder). Do not hardcode absolute
# paths here -- that was the single biggest source of breakage in the two
# original notebooks (see errors.md).
DATA_DIR = Path(".")
CODE_DIR = Path("code")
OUTPUT_DIR = Path("new_outputs")
INTERACTIVE_OUTPUT_DIR = Path("new_outputs_interactive")
PREVIOUS_RESULTS_DIR = Path("previous_results")

# All raw inputs were consolidated into a single flat inputs/ folder
# (2026-09-29 cleanup) -- everything unrelated to this analysis was deleted
# and the old FoldX/ and YoochanStuff/ subfolders were flattened into it.
# If this ever needs to move again, this is the only line that has to change.
INPUTS_DIR = Path("inputs")

SCREEN_DATA_XLSX = INPUTS_DIR / "PRC2_Screens_Counts_Scores_ZScores.xlsx"
AA_MAPPING_XLSX = INPUTS_DIR / "PRC2_Screens_AA_Mapping.xlsx"
KMEANS_CLUSTERS_XLSX = INPUTS_DIR / "PRC2_Screens_KMeans_Clusters.xlsx"
PWES_CLUSTERS_XLSX = INPUTS_DIR / "PRC2_Screens_PWES_Clusters.xlsx"
PDB_FILE = INPUTS_DIR / "6wkr.pdb"
# FoldX/{gene}.tsv used to live in its own FoldX/ subfolder; those three
# per-gene files (EZH2.tsv, EED.tsv, SUZ12.tsv) now sit directly in inputs/,
# so FOLDX_DIR is just inputs/ itself -- code/data_loading.py's
# `config.FOLDX_DIR / f"{g}.tsv"` pattern is unchanged.
FOLDX_DIR = INPUTS_DIR

# Supplied 2026-09-28 (previously missing -- see errors.md E-02 history).
# Loaders that need it still check for its existence and raise a clear,
# actionable error if it's ever absent again, rather than failing deep
# inside a plotting function.
MUTATION_LEVEL_DATA_TSV = INPUTS_DIR / "mutation_level_data.tsv"

# These already have very close to the schema the structural-enrichment
# section expects (Gene Symbol, sgRNA ID, FoldX, FoldX_aa, ThermoMPNN,
# GFP-MAX, ...) -- see errors.md E-03 for the exact column comparison.
# Used to live under a YoochanStuff/ subfolder; now flat in inputs/.
SGRNA_LEVEL_ABE_TSV = INPUTS_DIR / "sgRNA_level_ABE_251112.tsv"
SGRNA_LEVEL_CBE_TSV = INPUTS_DIR / "sgRNA_level_CBE_251112.tsv"

# ---------------------------------------------------------------------------
# Genes / subunits
# ---------------------------------------------------------------------------
PRC2_GENES = ["EZH2", "EED", "SUZ12", "JARID2", "AEBP2"]
MTF2_GENE = ["MTF2"]
ALL_GENES = PRC2_GENES + MTF2_GENE

# Protein lengths used to build the full 1..N residue scaffold for the
# PyMOL b-factor export (see code/pymol_export.py). Confirm these against
# the isoforms actually used for the screen library before relying on them.
PROTEIN_LENGTHS = {
    # EZH2: CORRECTED 2026-09-28 (see errors.md V-03) -- the screen uses the
    # 751-residue isoform (NB1's xmax_dict agreed), not 746. 746 truncated
    # the last 5 residues (positions 742-744 have non-zero scores in
    # PRC2_Screens_AA_Mapping.xlsx that were being dropped).
    "EZH2": 751,
    "EED": 441,
    "SUZ12": 739,
    "JARID2": 1246,
    "AEBP2": 517,  # NB1 truncated this to <=503 for PDB purposes -- see errors.md E-06
    "MTF2": 593,
}

# ---------------------------------------------------------------------------
# Colors (carried over from the original notebooks' custom_palette / mut_pal
# / subunit_colors -- consolidated here so both notebooks agree)
# ---------------------------------------------------------------------------
SUBUNIT_COLORS = {
    "EZH2": "#4C72B0",
    "EED": "#DD8452",
    "SUZ12": "#55A868",
    "JARID2": "#C44E52",
    "AEBP2": "#8172B2",
    "MTF2": "#937860",
}

MUTATION_COLORS = {
    "Missense": "#4C72B0",
    "Silent": "#B0B0B0",
    "Nonsense": "#C44E52",
    "Splice site": "#DD8452",
    "Intron": "#DDDDDD",
    "UTR": "#EEEEEE",
    "No mutation": "#F0F0F0",
}

EDITOR_COLORS = {
    "ABE": "#4C72B0",
    "CBE": "#C44E52",
}

# ---------------------------------------------------------------------------
# Timepoint selection -- see the module docstring above before changing
# ---------------------------------------------------------------------------
TIMEPOINTS = {
    ("WT", "ABE"): "Wk3",      # PLACEHOLDER -- unconfirmed, see errors.md E-01a
    ("WT", "CBE6b"): "Wk3",    # PLACEHOLDER -- unconfirmed, see errors.md E-01b
    ("Q575R", "CBE39"): "Wk3",  # PLACEHOLDER -- unconfirmed, see errors.md E-01c
    ("MTF2", "ABE"): "Wk3",    # Confirmed from previous_results/Figure3_* filenames
}

# Per-comparison overrides of TIMEPOINTS, keyed (screen, editor, comparison).
# Checked by data_loading._z() before TIMEPOINTS.
#   WT / CBE6b / PRC2i-D0 -> Wk6: ADDED 2026-09-28 (errors.md, "Output-parity
#   pass"). The old Figure1_Boxplot/WT_Resistance_Boxplot.svg
#   (a vector plot) has CBE6b box quartiles/whiskers that match
#   WT-CBE6b-PRC2i-D0-Wk6_Z to the drawn coordinates, and do not match Wk3.
#   Its ABE box matches WT-ABE-PRC2i-D0-Wk3_Z, not Wk4. This settles E-01b for
#   the PRC2i-D0 comparison only. WT CBE6b DMSO-D0 stays Wk3, which the k-means
#   reference reproduces exactly. The only consumer of this column is the
#   WT_Resistance boxplot.
TIMEPOINT_OVERRIDES = {
    ("WT", "CBE6b", "PRC2i-D0"): "Wk6",
}

# ---------------------------------------------------------------------------
# K-means clustering / PWES defaults carried over from NB1
# ---------------------------------------------------------------------------
KMEANS_K = 6
LOF_STABILITY_CUTOFF = 3
LOESS_N_REPEATS = 1000  # NB1 referenced an undefined `n_repeats` -- see errors.md E-04

# K-means random_state for the clustered-heatmap section. The notebook passes
# it to be_scan's kmeans_pca_scatterplot(random_state=...), which runs
# sklearn's KMeans. kmeans_clustered_heatmap only draws precomputed labels and
# takes no random_state.
#
# CHOSEN BY GRID SEARCH, 2026-09-28, redone from scratch in the output-parity
# pass (errors.md "K-means seed lock-in"). Setup: sklearn 1.8.0,
# KMeans(n_clusters=6, random_state=s), n_init="auto" (one k-means++ init).
# Input: the notebook's cutoff-3 hit-guide matrix, sorted by sgRNA ID. Each
# seed s in 0-199 was scored by adjusted Rand index (ARI) against
# PRC2_Screens_KMeans_Clusters.xlsx (ABE_KMeans 473 guides, CBE_KMeans 237).
# One seed serves both editors; the criterion was max of min(ARI_ABE, ARI_CBE).
# Seed 60 also has the highest mean ARI.
#   Seed 60:  ABE ARI = 0.9275, CBE ARI = 0.9082   (seed 0: 0.479 / 0.740)
#   Seeds 0-199, ABE min / median / max = 0.418 / 0.651 / 0.963
#   Seeds 0-199, CBE min / median / max = 0.444 / 0.716 / 0.965
#   Seeds 0-99 only: ABE 0.418 / 0.577 / 0.928, CBE 0.493 / 0.742 / 0.948,
#                    and the best seed is still 60.
#   Seeds with both ARIs >= 0.8: 8%. Seeds with both >= 0.9: 1.5% (3 seeds).
# The best seed per editor differs: 119 for ABE (0.963) and 151 for CBE
# (0.965), and each does worse on the other editor (0.868 / 0.837).
# HONEST READING: no seed reproduces the reference exactly (ARI 1.0). Seed 60
# is the best of 200 tries, not a reproduction. The typical seed gives about
# 0.65 / 0.72. With seed 60, 16 of 473 ABE guides and 14 of 237 CBE
# guides differ from the reference in cluster assignment. If exact manuscript
# labels are needed, use PRC2_Screens_KMeans_Clusters.xlsx directly. The
# cutoff-2 *_K6 / *_K6_Raw heatmaps use the same seed, but it was not tuned
# for them.
KMEANS_RANDOM_STATE = 60

# Number of permutations for the PWES cluster-WAP significance test
# (code/pwes_significance.py). 10000 matches the old df_wap outputs.
PWES_N_PERMUTATIONS = 10000
