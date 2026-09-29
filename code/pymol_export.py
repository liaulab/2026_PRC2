"""
Residue-score -> PyMOL b-factor ``.txt`` export pipeline.

Ports the PyMOL-mapping cells of ``251122_bescan_Plots_Github.ipynb`` (NB1
cells 60-76) into reusable functions, styled to match ``code/data_loading.py``.

Each exported ``.txt`` file holds one value per line, in ascending residue
order, with no header -- exactly what ``load_bfact.py`` expects to read back
into a PyMOL session and write into the b-factor column of ``6wkr.pdb`` (or a
per-gene structure) for surface coloring.

Two pieces of bookkeeping used to be copy-pasted several times in NB1 and are
consolidated here:

* ``build_complete_scaffold`` -- the "complete 1..N residue" DataFrame that NB1
  rebuilt inline in cells 63, 66, 69 and 72.
* ``remap_to_pdb_numbering`` -- the EZH2 isoform-2->isoform-1 offset and the
  AEBP2 truncation, which NB1 implemented three slightly different ways across
  the PyMOL-export cells (64, 67) and the PWES cell (81). See errors.md E-11.
"""

from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

import numpy as np
import pandas as pd

from . import config


# ---------------------------------------------------------------------------
# Per-residue reshaping helpers (NB1 cell 61)
# ---------------------------------------------------------------------------

def filter_data(
    df: pd.DataFrame,
    value_col: str,
    pos_col: str,
    common_cols: Sequence[str] = ("Gene Symbol", "pos"),
) -> pd.DataFrame:
    """Keep positive-position rows and floor the position to an integer residue.

    Equivalent to NB1's ``filter_data``: drops rows whose ``pos_col`` is <= 0,
    floors the (possibly fractional) position into an integer ``pos`` column,
    and returns only ``common_cols`` plus ``value_col``.

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
    """Collapse to one row per residue and fill the missing residues with 0.

    Equivalent to NB1's ``merge_full_df``: takes the per-residue maximum of
    ``value_col``, left-merges it onto the complete 1..N scaffold so every
    residue is present, replaces NaN and negative values with 0, and
    optionally min-max scales the result to ``[0, scale]``.

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
    """The complete ``Gene Symbol`` x 1..N residue scaffold for a set of genes.

    Replaces the four inline copies of this construction in NB1 (cells 63, 66,
    69, 72). One row per residue, positions 1..``protein_lengths[gene]``.

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
    """Remap screen residue numbering onto the 6WKR PDB numbering.

    Consolidates the three slightly different inline implementations NB1 used
    (PyMOL-export cells 64/67 and PWES cell 81) into one helper. Two genes need
    adjustment:

    * **EZH2** was screened as isoform 2, which carries a 5-residue insertion
      (positions 299-303) absent from the isoform-1 structure. Those five
      positions are always dropped. When ``shift`` is True, every position
      after the insertion (``pos`` > 303) is additionally shifted down by 5 so
      the numbers line up with the PDB chain -- this is required when residue
      numbers are used as a key (the PWES pipeline, cell 81). For the dense
      "complete scaffold" ``.txt`` export (cells 64/67), only the ordered value
      column is written, so ``shift`` is left False and dropping the five rows
      is sufficient.
    * **AEBP2** is truncated to residues <= 503 for PDB purposes.

    Other genes pass through unchanged.

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
    """Write one ``.txt`` of per-residue values per gene, in residue order.

    For each gene, optionally applies :func:`remap_to_pdb_numbering` (drop-only,
    matching NB1 cells 64/67/76), then writes the ``value_col`` column -- one
    value per line, ascending residue order, no header/index -- to
    ``out_dir/<name_template.format(gene=gene)>``.

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
