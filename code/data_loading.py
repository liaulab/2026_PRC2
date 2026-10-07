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
    return pd.read_excel(config.SCREEN_DATA_XLSX, sheet_name=sheet_name)

_ID_COLS = [
    "sgRNA ID", "sgRNA sequence", "Gene Symbol",
    "AtoG_pos", "AtoG_mutations", "AtoG_muttypes", "AtoG_muttype",
    "CtoT_pos", "CtoT_mutations", "CtoT_muttypes", "CtoT_muttype",
]

def _z(sheet_prefix: str, editor: str, comparison: str, screen: str) -> str:
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
    """
    Output columns (beyond the shared _ID_COLS):
      WT_LOF_dropout_ABE_zscore     <- WT-ABE-DMSO-D0-{Wk}_Z
      WT_LOF_dropout_CBE39_zscore   <- WT-CBE39-DMSO-D0-Wk3_Z   (single timepoint)
      WT_LOF_dropout_CBE6b_zscore   <- WT-CBE6b-DMSO-D0-{Wk}_Z
      WT_LOF_positiveselection_ABE_zscore   <- WT-ABE-PRC2i-D0-{Wk}_Z
      WT_LOF_positiveselection_CBE39_zscore <- WT-CBE39-PRC2i-D0-Wk8_Z (single)
      WT_LOF_positiveselection_CBE6b_zscore <- WT-CBE6b-PRC2i-D0-Wk6_Z
                                               (config.TIMEPOINT_OVERRIDES)
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
    """
    Output columns:
      Q575R_LOF_dropout_ABE_zscore          <- Q575R-ABE-PRC2i-D0_Z    (single tp)
      Q575R_LOF_dropout_CBE39_zscore        <- Q575R-CBE39-PRC2i-D0-{Wk}_Z
      Q575R_LOF_positiveselection_ABE_zscore   <- Q575R-ABE-DMSO-PRC2i_Z    (single tp)
      Q575R_LOF_positiveselection_CBE39_zscore <- Q575R-CBE39-DMSO-PRC2i-{Wk}_Z
      Q575R_LOF_dmso_d0_ABE_zscore           <- Q575R-ABE-DMSO-D0_Z (single tp)
      Q575R_LOF_dmso_d0_CBE39_zscore         <- Q575R-CBE39-DMSO-D0-{Wk}_Z
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
# Stability screen
# ---------------------------------------------------------------------------

def load_stability_screen() -> pd.DataFrame:
    """
    Output columns:
      Abundance_GFP-EZH2_ABE_zscore   <- GFP-EZH2-ABE-lowGFP-unsorted_Z
      Abundance_GFP-EZH2_CBE6b_zscore <- GFP-EZH2-CBE6-lowGFP-unsorted_Z
      Abundance_GFP-EED_ABE_zscore    <- GFP-EED-ABE-lowGFP-unsorted_Z
      Abundance_GFP-EED_CBE6b_zscore  <- GFP-EED-CBE6-lowGFP-unsorted_Z
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
    """
    Output columns:
      WT_DMSO_D0_Wk3_zscore              <- WT-ABE-DMSO-D0-Wk3_Z
      Q575R_DMSO_D0_Wk3_zscore           <- Q575R-ABE-DMSO-D0-Wk3_Z
      Q575R_PRC2i_D0_Wk3_zscore          <- Q575R-ABE-PRC2i-D0-Wk3_Z
      Q575R_DMSO_PRC2i_Wk3_zscore        <- Q575R-ABE-DMSO-PRC2i-Wk3_Z
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
# AA-level mapping / cluster tables
# ---------------------------------------------------------------------------

def load_aa_mapping() -> pd.DataFrame:
    return pd.read_excel(config.AA_MAPPING_XLSX, sheet_name="AA_Mapping_Scores")


def load_mtf2_aa_mapping() -> pd.DataFrame:
    return pd.read_excel(config.AA_MAPPING_XLSX, sheet_name="MTF2_Mapping_Scores")


def load_kmeans_clusters(editor: str) -> pd.DataFrame:
    sheet = f"{editor}_KMeans"
    return pd.read_excel(config.KMEANS_CLUSTERS_XLSX, sheet_name=sheet)


def load_pwes_clusters(kind: str) -> pd.DataFrame:
    sheet = {
        "LOF": "LOF_Clusters",
        "Stability": "Stability_Clusters",
        "LOF_Stab": "LOF_Stab_Clusters",
    }[kind]
    return pd.read_excel(config.PWES_CLUSTERS_XLSX, sheet_name=sheet)


# ---------------------------------------------------------------------------
# sgRNA-level structural-feature tables
# ---------------------------------------------------------------------------

def load_sgrna_level(editor: str) -> pd.DataFrame:
    """
    """
    path = config.SGRNA_LEVEL_ABE_TSV if editor == "ABE" else config.SGRNA_LEVEL_CBE_TSV
    return pd.read_csv(path, sep="\t", index_col=0)


def load_mutation_level_data() -> pd.DataFrame:
    """
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
