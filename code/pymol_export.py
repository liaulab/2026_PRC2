from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

import numpy as np
import pandas as pd

from . import config


# ---------------------------------------------------------------------------
# Per-residue reshaping helpers
# ---------------------------------------------------------------------------

def filter_data(
    df: pd.DataFrame,
    value_col: str,
    pos_col: str,
    common_cols: Sequence[str] = ("Gene Symbol", "pos"),
) -> pd.DataFrame:
    """
    Keep positive-position rows and floor the position to an integer residue.

    Args:
        df: sgRNA-level table for one condition.
        value_col: the score column to carry through.
        pos_col: the source position column (e.g. ``AtoG_pos``/``CtoT_pos``/``aa_pos``).
        common_cols: identifier columns to retain (default ``Gene Symbol``, ``pos``).

    Returns:
        A copy with an integer ``pos`` column and only the requested columns.
    """
    common_cols = list(common_cols)
    df = df.copy()
    df = df[df[pos_col] > 0]
    df["pos"] = df[pos_col].apply(np.floor).astype(int)
    return df[common_cols + [value_col]]


def merge_full_df(
    df: pd.DataFrame,
    df_complete: pd.DataFrame,
    value_col: str,
    common_cols: Sequence[str] = ("Gene Symbol", "pos"),
    scale: Optional[float] = None,
) -> pd.DataFrame:
    """
    Collapse to one row per residue and fill the missing residues with 0.

    Args:
        df: filtered sgRNA-level table (see :func:`filter_data`).
        df_complete: the complete residue scaffold (see :func:`build_complete_scaffold`).
        value_col: the score column to aggregate.
        common_cols: the residue key columns to group/merge on.
        scale: if given, min-max scale the non-negative scores to ``[0, scale]``.

    Returns:
        The scaffold DataFrame with a filled, non-negative ``value_col``.
    """
    common_cols = list(common_cols)
    df = df.copy()
    result = df.groupby(common_cols, as_index=False)[value_col].max()
    result = df_complete.copy().merge(result, how="left", on=common_cols)
    result[value_col] = result[value_col].fillna(0)
    result[value_col] = result[value_col].apply(lambda x: x if x > 0 else 0)

    if scale:
        lo, hi = result[value_col].min(), result[value_col].max()
        if hi > lo:
            result[value_col] = (result[value_col] - lo) / (hi - lo) * scale
    return result


def build_complete_scaffold(
    genes: Sequence[str],
    protein_lengths: Optional[Dict[str, int]] = None,
) -> pd.DataFrame:
    """
    The complete ``Gene Symbol`` x 1..N residue scaffold for a set of genes.

    Args:
        genes: the genes to include.
        protein_lengths: gene -> length map (defaults to ``config.PROTEIN_LENGTHS``).

    Returns:
        A concatenated DataFrame with ``Gene Symbol`` and integer ``pos`` columns.
    """
    if protein_lengths is None:
        protein_lengths = config.PROTEIN_LENGTHS
    frames = []
    for gene in genes:
        n = protein_lengths[gene]
        frames.append(pd.DataFrame({
            "Gene Symbol": [gene] * n,
            "pos": [int(i + 1) for i in range(n)],
        }))
    return pd.concat(frames, ignore_index=True)


# ---------------------------------------------------------------------------
# Unified isoform / truncation remapping (errors.md E-11)
# ---------------------------------------------------------------------------

def remap_to_pdb_numbering(
    df: pd.DataFrame,
    gene_col: str = "Gene Symbol",
    pos_col: str = "pos",
    shift: bool = False,
) -> pd.DataFrame:
    """
    Remap screen residue numbering onto the 6WKR PDB numbering.

    Args:
        df: a table with a gene column and an integer position column.
        gene_col: name of the gene column.
        pos_col: name of the position column.
        shift: whether to shift EZH2 positions > 303 down by 5 (see above).

    Returns:
        A filtered (and optionally renumbered) copy of ``df``.
    """
    df = df.copy()

    # EZH2 isoform-2 insertion at 299-303: always drop it.
    ezh2_insert = (
        (df[gene_col] == "EZH2")
        & (df[pos_col] > 298)
        & (df[pos_col] <= 303)
    )
    df = df[~ezh2_insert]

    if shift:
        after = (df[gene_col] == "EZH2") & (df[pos_col] > 303)
        df.loc[after, pos_col] -= 5

    # AEBP2 truncation to <= 503.
    aebp2_drop = (df[gene_col] == "AEBP2") & (df[pos_col] > 503)
    df = df[~aebp2_drop]

    return df


# ---------------------------------------------------------------------------
# Export driver
# ---------------------------------------------------------------------------

def export_scores_to_txt(
    df: pd.DataFrame,
    value_col: str,
    genes: Sequence[str],
    out_dir: Union[str, Path],
    name_template: str,
    gene_col: str = "Gene Symbol",
    pos_col: str = "pos",
    apply_pdb_remap: bool = True,
) -> Dict[str, Path]:
    """
    Write one ``.txt`` of per-residue values per gene, in residue order.

    Args:
        df: a residue-level scaffold table (see :func:`merge_full_df`).
        value_col: the score column to write.
        genes: the genes to export.
        out_dir: directory to write into (created if needed).
        name_template: filename template taking a ``{gene}`` field, e.g.
            ``"stability-{gene}-ABE.txt"``.
        gene_col: gene column name.
        pos_col: position column name.
        apply_pdb_remap: whether to drop the EZH2 insertion / AEBP2 tail first.

    Returns:
        A gene -> written path map.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    written: Dict[str, Path] = {}
    for gene in genes:
        df_gene = df[df[gene_col] == gene]
        if apply_pdb_remap:
            df_gene = remap_to_pdb_numbering(
                df_gene, gene_col=gene_col, pos_col=pos_col, shift=False
            )
        df_gene = df_gene.sort_values(pos_col)
        path = out_dir / name_template.format(gene=gene)
        df_gene[value_col].to_csv(path, index=False, header=False)
        written[gene] = path
    return written
