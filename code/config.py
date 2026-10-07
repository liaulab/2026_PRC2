"""
Shared configuration for the PRC2 base-editor screen analysis notebooks
"""

from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

# All paths are relative to the notebook's working directory
DATA_DIR = Path(".")
CODE_DIR = Path("code")
OUTPUT_DIR = Path("new_outputs")
INTERACTIVE_OUTPUT_DIR = Path("new_outputs_interactive")
PREVIOUS_RESULTS_DIR = Path("previous_results")

# All raw inputs were consolidated into a single flat inputs/ 
INPUTS_DIR = Path("inputs")

SCREEN_DATA_XLSX = INPUTS_DIR / "PRC2_Screens_Counts_Scores_ZScores.xlsx"
AA_MAPPING_XLSX = INPUTS_DIR / "PRC2_Screens_AA_Mapping.xlsx"
KMEANS_CLUSTERS_XLSX = INPUTS_DIR / "PRC2_Screens_KMeans_Clusters.xlsx"
PWES_CLUSTERS_XLSX = INPUTS_DIR / "PRC2_Screens_PWES_Clusters.xlsx"
PDB_FILE = INPUTS_DIR / "6wkr.pdb"
FOLDX_DIR = INPUTS_DIR

MUTATION_LEVEL_DATA_TSV = INPUTS_DIR / "mutation_level_data.tsv"
SGRNA_LEVEL_ABE_TSV = INPUTS_DIR / "sgRNA_level_ABE_251112.tsv"
SGRNA_LEVEL_CBE_TSV = INPUTS_DIR / "sgRNA_level_CBE_251112.tsv"

# ---------------------------------------------------------------------------
# Genes / subunits
# ---------------------------------------------------------------------------

PRC2_GENES = ["EZH2", "EED", "SUZ12", "JARID2", "AEBP2"]
MTF2_GENE = ["MTF2"]
ALL_GENES = PRC2_GENES + MTF2_GENE

PROTEIN_LENGTHS = {
    "EZH2": 751, # truncated this to 746 for PDB purposes
    "EED": 441,
    "SUZ12": 739,
    "JARID2": 1246,
    "AEBP2": 517,  # truncated this to 503 for PDB purposes
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
    ("WT", "ABE"): "Wk3", 
    ("WT", "CBE6b"): "Wk3", 
    ("Q575R", "CBE39"): "Wk3",  
    ("MTF2", "ABE"): "Wk3", 
}

TIMEPOINT_OVERRIDES = {
    ("WT", "CBE6b", "PRC2i-D0"): "Wk6",
}

# ---------------------------------------------------------------------------
# K-means clustering / PWES defaults carried over from NB1
# ---------------------------------------------------------------------------

KMEANS_K = 6
LOF_STABILITY_CUTOFF = 3
KMEANS_RANDOM_STATE = 60
PWES_N_PERMUTATIONS = 10000
