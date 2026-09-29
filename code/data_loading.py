"""
Reshapes PRC2_Screens_Counts_Scores_ZScores.xlsx (the current single source
of screen data) into the long-format, per-condition tables that the
plotting code (both the be_scan calls carried over from
251122_bescan_Plots_Github.ipynb and the ports of figures_ycm_xyh_v2.ipynb)
expects.

Why this module exists
-----------------------
Both original notebooks read a collection of small, per-condition files
(Table_S1_ScreenData_251214.csv, Table_S3_MTF2_ScreenData_251215.csv,
stability_*_controlZscore.csv, WT_*_controlZscore.csv,
Q575R_*_controlZscore.csv, sgRNA_level_{GENE}_{EDITOR}.tsv) that no longer
exist anywhere in this folder. Per Calvin's instructions, the single
source of truth is now PRC2_Screens_Counts_Scores_ZScores.xlsx, which
already contains raw counts, log-fold-changes, and per-condition Z-scores
for every sgRNA in four sheets (WT_Screens, Q575R_Screens,
Stability_Screens, MTF2_Screens). This module selects and renames the
columns each analysis needs so the rest of the pipeline (which was
written against the old file names/column names) needs minimal changes.

This module never reads the full workbook into a notebook cell's output --
it opens the file once per process (cached via functools.lru_cache) and
hands back only the columns asked for.
"""

from functools import lru_cache
from pathlib import Path
from typing import Optional

import pandas as pd

from . import config


# ---------------------------------------------------------------------------
# Low-level sheet access (cached so the ~30MB workbook is only parsed once)
# ---------------------------------------------------------------------------

@lru_cache(maxsize=None)
def _read_sheet(sheet_name: str) -> pd.DataFrame:
    """Read one sheet of PRC2_Screens_Counts_Scores_ZScores.xlsx in full.

    Cached per sheet name so repeated calls across notebook cells don't
    re-parse the workbook. Do not call this directly from notebook cells
    for exploratory purposes -- Calvin's data-access rule is head(50) only
    for ad hoc inspection of this file; use the sheet-specific loaders
    below, which return only the columns a given analysis needs.
    """
    return pd.read_excel(config.SCREEN_DATA_XLSX, sheet_name=sheet_name)


_ID_COLS = [
    "sgRNA ID", "sgRNA sequence", "Gene Symbol",
    "AtoG_pos", "AtoG_mutations", "AtoG_muttypes", "AtoG_muttype",
    "CtoT_pos", "CtoT_mutations", "CtoT_muttypes", "CtoT_muttype",
]


def _z(sheet_prefix: str, editor: str, comparison: str, screen: str) -> str:
    """Build the exact Z-score column name for (screen, editor, comparison),
    inserting the configured timepoint suffix only for the conditions that
    have more than one timepoint in the workbook (see config.TIMEPOINTS).
    """
    override = config.TIMEPOINT_OVERRIDES.get((screen, editor, comparison))
    if override is not None:
        return f"{sheet_prefix}-{editor}-{comparison}-{override}_Z"
    key = (screen, editor)
    if key in config.TIMEPOINTS:
        wk = config.TIMEPOINTS[key]
        return f"{sheet_prefix}-{editor}-{comparison}-{wk}_Z"
    return f"{sheet_prefix}-{editor}-{comparison}_Z"


# ---------------------------------------------------------------------------
# WT screen
# ---------------------------------------------------------------------------

def load_wt_screen() -> pd.DataFrame:
    """Long-format WT-screen table, equivalent to the old
    Table_S1_ScreenData_251214.csv for the WT-derived columns.

    Output columns (beyond the shared _ID_COLS):
      WT_LOF_dropout_ABE_zscore     <- WT-ABE-DMSO-D0-{Wk}_Z
      WT_LOF_dropout_CBE39_zscore   <- WT-CBE39-DMSO-D0-Wk3_Z   (single timepoint)
      WT_LOF_dropout_CBE6b_zscore   <- WT-CBE6b-DMSO-D0-{Wk}_Z
      WT_LOF_positiveselection_ABE_zscore   <- WT-ABE-PRC2i-D0-{Wk}_Z
      WT_LOF_positiveselection_CBE39_zscore <- WT-CBE39-PRC2i-D0-Wk8_Z (single)
      WT_LOF_positiveselection_CBE6b_zscore <- WT-CBE6b-PRC2i-D0-Wk6_Z
                                               (config.TIMEPOINT_OVERRIDES)

    See config.TIMEPOINTS for which timepoint {Wk} resolves to for ABE and
    CBE6b. Evidence so far (errors.md "Side result for E-01" and the
    2026-09-28 output-parity pass): the old outputs are
    reproduced by Wk3 for ABE DMSO-D0, ABE PRC2i-D0 and CBE6b DMSO-D0, and by
    Wk6 for CBE6b PRC2i-D0. The override in config.TIMEPOINT_OVERRIDES
    applies the Wk6 result.
    """
    df = _read_sheet("WT_Screens")
    out = df[_ID_COLS].copy()

    out["WT_LOF_dropout_ABE_zscore"] = df[_z("WT", "ABE", "DMSO-D0", "WT")]
    out["WT_LOF_dropout_CBE39_zscore"] = df["WT-CBE39-DMSO-D0-Wk3_Z"]
    out["WT_LOF_dropout_CBE6b_zscore"] = df[_z("WT", "CBE6b", "DMSO-D0", "WT")]

    out["WT_LOF_positiveselection_ABE_zscore"] = df[_z("WT", "ABE", "PRC2i-D0", "WT")]
    out["WT_LOF_positiveselection_CBE39_zscore"] = df["WT-CBE39-PRC2i-D0-Wk8_Z"]
    out["WT_LOF_positiveselection_CBE6b_zscore"] = df[_z("WT", "CBE6b", "PRC2i-D0", "WT")]

    return out


# ---------------------------------------------------------------------------
# Q575R screen
# ---------------------------------------------------------------------------

def load_q575r_screen() -> pd.DataFrame:
    """Long-format Q575R-screen table, equivalent to the old
    Table_S1_ScreenData_251214.csv for the Q575R-derived columns.

    CORRECTED 2026-09-28 (see errors.md, "Validation against
    previous_results", root-cause note and V-01/V-02/V-03/V-05): the
    original version of this loader mapped `positiveselection` to
    PRC2i-D0 and `dropout` to DMSO-D0, which does NOT match NB1's
    naming convention. Reconfirmed by exact reproduction of
    PRC2_Screens_KMeans_Clusters.xlsx, PRC2_Screens_PWES_Clusters.xlsx,
    and PRC2_Screens_AA_Mapping.xlsx: NB1's
    `Q575R_LOF_positiveselection_*` is the DMSO-PRC2i comparison, and
    `Q575R_LOF_dropout_*` is the PRC2i-D0 comparison.

    Output columns:
      Q575R_LOF_dropout_ABE_zscore          <- Q575R-ABE-PRC2i-D0_Z    (single tp)
      Q575R_LOF_dropout_CBE39_zscore        <- Q575R-CBE39-PRC2i-D0-{Wk}_Z
      Q575R_LOF_positiveselection_ABE_zscore   <- Q575R-ABE-DMSO-PRC2i_Z    (single tp)
      Q575R_LOF_positiveselection_CBE39_zscore <- Q575R-CBE39-DMSO-PRC2i-{Wk}_Z
      Q575R_LOF_dmso_d0_ABE_zscore           <- Q575R-ABE-DMSO-D0_Z (single tp)
      Q575R_LOF_dmso_d0_CBE39_zscore         <- Q575R-CBE39-DMSO-D0-{Wk}_Z

    (The DMSO-D0 comparison didn't have a name in NB1's convention --
    it wasn't one of the two primary "dropout"/"positiveselection"
    columns -- so it's given the new, self-describing name
    `Q575R_LOF_dmso_d0_*` here.)

    {Wk} for CBE39 is CONFIRMED as Wk3 for the DMSO-PRC2i comparison
    (reproduces PRC2_Screens_AA_Mapping.xlsx exactly -- see errors.md
    "Side result for E-01"). PRC2i-D0 and DMSO-D0 timepoints for CBE39
    are still unconfirmed placeholders.
    """
    df = _read_sheet("Q575R_Screens")
    out = df[_ID_COLS].copy()

    out["Q575R_LOF_dropout_ABE_zscore"] = df["Q575R-ABE-PRC2i-D0_Z"]
    out["Q575R_LOF_dropout_CBE39_zscore"] = df[_z("Q575R", "CBE39", "PRC2i-D0", "Q575R")]

    out["Q575R_LOF_positiveselection_ABE_zscore"] = df["Q575R-ABE-DMSO-PRC2i_Z"]
    out["Q575R_LOF_positiveselection_CBE39_zscore"] = df[_z("Q575R", "CBE39", "DMSO-PRC2i", "Q575R")]

    out["Q575R_LOF_dmso_d0_ABE_zscore"] = df["Q575R-ABE-DMSO-D0_Z"]
    out["Q575R_LOF_dmso_d0_CBE39_zscore"] = df[_z("Q575R", "CBE39", "DMSO-D0", "Q575R")]

    return out


# ---------------------------------------------------------------------------
# Stability screen (no timepoint ambiguity -- single lowGFP/highGFP split)
# ---------------------------------------------------------------------------

def load_stability_screen() -> pd.DataFrame:
    """Long-format Stability-screen table, equivalent to the old
    stability_{ABE,CBE6}_{EZH2,EED}_controlZscore.csv files, all in one
    frame.

    Output columns:
      Abundance_GFP-EZH2_ABE_zscore   <- GFP-EZH2-ABE-lowGFP-unsorted_Z
      Abundance_GFP-EZH2_CBE6b_zscore <- GFP-EZH2-CBE6-lowGFP-unsorted_Z
      Abundance_GFP-EED_ABE_zscore    <- GFP-EED-ABE-lowGFP-unsorted_Z
      Abundance_GFP-EED_CBE6b_zscore  <- GFP-EED-CBE6-lowGFP-unsorted_Z
    (high-GFP columns are also carried over with a _high suffix in case a
    later section wants the other tail of the sort.)
    """
    df = _read_sheet("Stability_Screens")
    out = df[_ID_COLS].copy()

    for gene in ("EZH2", "EED"):
        for editor, col_editor in (("ABE", "ABE"), ("CBE6b", "CBE6")):
            out[f"Abundance_GFP-{gene}_{editor}_zscore"] = df[
                f"GFP-{gene}-{col_editor}-lowGFP-unsorted_Z"
            ]
            out[f"Abundance_GFP-{gene}_{editor}_zscore_high"] = df[
                f"GFP-{gene}-{col_editor}-highGFP-unsorted_Z"
            ]

    return out


# ---------------------------------------------------------------------------
# MTF2 screen
# ---------------------------------------------------------------------------

def load_mtf2_screen() -> pd.DataFrame:
    """Long-format MTF2-screen table, equivalent to the old
    Table_S3_MTF2_ScreenData_251215.csv.

    CORRECTED 2026-09-28 (see errors.md, V-03): `aa_pos` is now taken
    from `AtoG_pos` (the amino-acid position used throughout the other
    three sheets), not `gene_pos` (a nucleotide-level coordinate that an
    earlier version of this loader used by mistake -- confirmed wrong
    because it does not reproduce PRC2_Screens_AA_Mapping.xlsx's
    MTF2_Mapping_Scores at all, while AtoG_pos reproduces all 4 columns
    there exactly).

    Output columns:
      WT_DMSO_D0_Wk3_zscore              <- WT-ABE-DMSO-D0-Wk3_Z
      Q575R_DMSO_D0_Wk3_zscore           <- Q575R-ABE-DMSO-D0-Wk3_Z
      Q575R_PRC2i_D0_Wk3_zscore          <- Q575R-ABE-PRC2i-D0-Wk3_Z
      Q575R_DMSO_PRC2i_Wk3_zscore        <- Q575R-ABE-DMSO-PRC2i-Wk3_Z
    Wk3 here is CONFIRMED (not a placeholder) from the timepoint encoded
    in previous_results/Figure3_Loess/*.svg filenames.
    """
    df = _read_sheet("MTF2_Screens")
    id_cols = [c for c in _ID_COLS if c in df.columns]  # MTF2 sheet has no "sgRNA ID" column
    out = df[id_cols].copy()
    out["Gene Symbol"] = df["Gene Symbol"]
    out["aa_pos"] = df["AtoG_pos"]

    wk = config.TIMEPOINTS[("MTF2", "ABE")]
    out["WT_DMSO_D0_Wk3_zscore"] = df[f"WT-ABE-DMSO-D0-{wk}_Z"]
    out["Q575R_DMSO_D0_Wk3_zscore"] = df[f"Q575R-ABE-DMSO-D0-{wk}_Z"]
    out["Q575R_PRC2i_D0_Wk3_zscore"] = df[f"Q575R-ABE-PRC2i-D0-{wk}_Z"]
    out["Q575R_DMSO_PRC2i_Wk3_zscore"] = df[f"Q575R-ABE-DMSO-PRC2i-{wk}_Z"]

    return out


# ---------------------------------------------------------------------------
# AA-level mapping / cluster tables (already residue-level; used directly)
# ---------------------------------------------------------------------------

def load_aa_mapping() -> pd.DataFrame:
    """PRC2_Screens_AA_Mapping.xlsx, sheet AA_Mapping_Scores -- this is
    already the per-residue LOF_Score / Stability_Score / Adjusted_Score
    table that NB1 built by hand in its PyMOL-export cells (filter_data +
    merge_full_df). Prefer this over re-deriving it from the sgRNA-level
    sheets when a residue-level score is all that's needed (PyMOL export,
    lollipop plots).
    """
    return pd.read_excel(config.AA_MAPPING_XLSX, sheet_name="AA_Mapping_Scores")


def load_mtf2_aa_mapping() -> pd.DataFrame:
    """PRC2_Screens_AA_Mapping.xlsx, sheet MTF2_Mapping_Scores."""
    return pd.read_excel(config.AA_MAPPING_XLSX, sheet_name="MTF2_Mapping_Scores")


def load_kmeans_clusters(editor: str) -> pd.DataFrame:
    """PRC2_Screens_KMeans_Clusters.xlsx, sheet {editor}_KMeans.

    editor: "ABE" or "CBE"
    """
    sheet = f"{editor}_KMeans"
    return pd.read_excel(config.KMEANS_CLUSTERS_XLSX, sheet_name=sheet)


def load_pwes_clusters(kind: str) -> pd.DataFrame:
    """PRC2_Screens_PWES_Clusters.xlsx.

    kind: "LOF", "Stability", or "LOF_Stab" (maps to LOF_Clusters,
    Stability_Clusters, LOF_Stab_Clusters).
    """
    sheet = {
        "LOF": "LOF_Clusters",
        "Stability": "Stability_Clusters",
        "LOF_Stab": "LOF_Stab_Clusters",
    }[kind]
    return pd.read_excel(config.PWES_CLUSTERS_XLSX, sheet_name=sheet)


# ---------------------------------------------------------------------------
# sgRNA-level structural-feature tables (for the structural-enrichment
# section ported from figures_ycm_xyh_v2.ipynb)
# ---------------------------------------------------------------------------

def load_sgrna_level(editor: str) -> pd.DataFrame:
    """sgRNA-level FoldX/ThermoMPNN/GFP-MAX table for one editor.

    editor: "ABE" or "CBE"

    Reads YoochanStuff/sgRNA_level_{editor}_251112.tsv, which has (nearly)
    the schema figures_ycm_xyh_v2.ipynb expects for sgRNA_level_ABE.tsv /
    sgRNA_level_CBE.tsv. See errors.md E-03 for the exact column
    comparison and any renames still needed.
    """
    path = config.SGRNA_LEVEL_ABE_TSV if editor == "ABE" else config.SGRNA_LEVEL_CBE_TSV
    return pd.read_csv(path, sep="\t", index_col=0)


def load_mutation_level_data() -> pd.DataFrame:
    """The residue-level structural-feature table (PSSM, ASA, secondary
    structure, contacting-residue counts, AlphaMissense, etc.) that the
    structural-enrichment section (odds ratios / forest plots / dist
    plots / substitution matrix / quantile boxplots) needs.

    NOT PRESENT in this folder as of the migration pass -- Calvin said
    he'll supply it. Raises a clear, actionable error rather than a bare
    FileNotFoundError so the notebook cell that calls this can catch it
    and skip the section cleanly (see master_analysis.ipynb, the
    Structural Enrichment section, and errors.md E-02).
    """
    if not Path(config.MUTATION_LEVEL_DATA_TSV).exists():
        raise FileNotFoundError(
            f"{config.MUTATION_LEVEL_DATA_TSV} not found. This file is required "
            "for the structural-enrichment section (odds ratios, forest plots, "
            "histograms, substitution matrix, quantile boxplots) ported from "
            "figures_ycm_xyh_v2.ipynb. Drop it into the RevisionAnalysis folder "
            "under this exact name/path, or edit config.MUTATION_LEVEL_DATA_TSV "
            "to point at it. See errors.md, item E-02, for details."
        )
    return pd.read_csv(config.MUTATION_LEVEL_DATA_TSV, sep="\t", index_col=0)
