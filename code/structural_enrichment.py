"""
Structural-enrichment analyses ported from ``figures_ycm_xyh_v2.ipynb`` (NB2).

These are the analyses that ``be_scan`` has no equivalent for and that NB1 did
not cover: the Fisher/odds-ratio enrichment pipeline and its forest plots, the
hit-vs-non-hit histogram (dist) plots, the 20x20 substitution matrix, the
quantile-stratified boxplot/stripplot, and the split-violin plot of stability
predictors. Restyled to match ``be_scan.figure_plot`` conventions (type hints,
docstrings, option-object driven where practical) and parametrized so nothing
writes to a hard-coded path.

Data access
-----------
Everything that needs the residue-level structural-feature table calls
:func:`code.data_loading.load_mutation_level_data`, which raises a clear
``FileNotFoundError`` when the (not-yet-supplied) ``mutation_level_data.tsv`` is
absent (errors.md E-02). That error is deliberately **not** swallowed here --
the notebook layer wraps the affected subsection in ``try/except
FileNotFoundError`` and prints a skip message. The violin and quantile-boxplot
functions instead use ``load_sgrna_level`` (present in the repo) and run for
real.

Notes on the port (see errors.md)
----------------------------------
* E-13 (resolved): NB2 defined ``plot_distplot_on_wild_to_mutant_merged_v4``
  twice, in cells 21 and 23. Cell 22's call (the only call site in the
  notebook) runs *before* cell 23's redefinition, so it always executed the
  cell-21 version; the cell-23 version was dead code, never actually run.
  This port keeps one definition, matching cell 21's behavior (no legend,
  single ``.svg`` output, 0.7cm panel width, ``Others`` sorted last) by
  default, with ``add_legend``/``save_png`` kwargs for the cell-23 look if
  ever wanted. There is no old figure anywhere in ``previous_results/`` for
  this function (or the rest of this odds-ratio/forest-plot subsection), so
  neither version could be checked against a reference either way.
* E-15: the duplicate ``("Structural region", "Interface")`` key in
  ``main_order_map`` is removed.
* E-16: the dead ``mCSM-PPI`` scatter branch is not ported.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, norm, mannwhitneyu

from be_scan.figure_plot.figure_classes import AxisLabelOpts, AXIS


# ---------------------------------------------------------------------------
# Font (Arial.ttf if present, graceful fallback otherwise)
# ---------------------------------------------------------------------------

def _load_font_prop(size: int = 6, filename: str = "Arial.ttf") -> fm.FontProperties:
    for candidate in (Path.cwd() / filename, Path(__file__).parent.parent / filename):
        if candidate.is_file():
            return fm.FontProperties(fname=str(candidate), size=size)
    try:
        found = fm.findfont(fm.FontProperties(family="Arial"), fallback_to_default=False)
        return fm.FontProperties(fname=found, size=size)
    except Exception:
        return fm.FontProperties(family="sans-serif", size=size)


arial_font6 = _load_font_prop(size=6)

AA_ORDER = ['K', 'R', 'H', 'D', 'E', 'S', 'T', 'N', 'Q', 'C',
            'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']


# ---------------------------------------------------------------------------
# Odds-ratio / enrichment pipeline
# ---------------------------------------------------------------------------

def calculate_odds_ratio_with_ci(data, alpha: float = 0.05) -> Dict[str, Optional[float]]:
    """Odds ratio, 95% CI and Fisher's-exact p-value for a 2x2 table.

    Args:
        data: a 2x2 contingency table (list-of-lists or ``np.ndarray``).
        alpha: significance level for the confidence interval.

    Returns:
        A dict with ``odds_ratio``, ``ci_low``, ``ci_high`` and ``p_value``
        (values are ``None`` when the computation fails, e.g. a zero cell).
    """
    if not isinstance(data, (list, np.ndarray)):
        raise TypeError("Data must be a list or numpy array")
    if np.shape(data) != (2, 2):
        raise ValueError("Data must be a 2x2 contingency table")
    try:
        odds_ratio, p_value = fisher_exact(data, alternative="two-sided")
        a, b, c, d = np.ravel(data)
        log_or = np.log(odds_ratio)
        se_log_or = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
        z = norm.ppf(1 - alpha / 2)
        ci_low = np.exp(log_or - z * se_log_or)
        ci_high = np.exp(log_or + z * se_log_or)
        return {"odds_ratio": odds_ratio, "ci_low": ci_low,
                "ci_high": ci_high, "p_value": p_value}
    except Exception as e:  # noqa: BLE001 - preserve NB2 behavior
        print(f"An error occurred during calculation: {e}")
        return {"odds_ratio": None, "ci_low": None, "ci_high": None, "p_value": None}


def convert_to_binary(
    input_pd: pd.DataFrame,
    column_name: str,
    continuous=False,
    range_list=False,
) -> pd.DataFrame:
    """One-hot / range-binned binary indicator columns for ``column_name``.

    Mirrors NB2: with ``continuous`` a list of thresholds, produces ``< t``
    indicators; with ``range_list`` a list of ``"lo~hi"`` strings, produces
    half-open ``[lo, hi)`` indicators; otherwise one column per unique value.
    """
    content_list = input_pd[column_name].to_list()

    if continuous:
        binary_dict = {
            f"{column_name}_<{cat}": [
                1 if x is not None and float(x) < float(cat) else 0
                for x in content_list if x is not None
            ]
            for cat in continuous
        }
        return pd.DataFrame(binary_dict)

    if range_list:
        binary_dict = {
            f"{column_name}_{each_range}": [
                1 if x is not None and lower <= float(x) < upper else 0
                for x in content_list if x is not None
            ]
            for each_range in range_list
            if (parts := each_range.split("~")) and len(parts) == 2 and parts[0] and parts[1]
            if (lower := float(parts[0])) is not None and (upper := float(parts[1])) is not None
        }
        return pd.DataFrame(binary_dict)

    category_list = list(set(content_list))
    binary_dict = {f"{cat}": [1 if x == cat else 0 for x in content_list]
                   for cat in category_list}
    return pd.DataFrame(binary_dict)


def run_odds_ratio_test_on_all_conditions(
    result_pd: pd.DataFrame,
    column_name: str,
    zvalue: float,
    continuous=False,
    range_list=False,
    skip_binary: bool = False,
    label: str = "lowGFP-unsorted_Zscore_controls",
) -> list:
    """One-vs-rest odds-ratio tests for every category of ``column_name``.

    Splits rows into hits (``label`` > ``zvalue``) and non-hits, builds a 2x2
    table of each category against all others, and runs
    :func:`calculate_odds_ratio_with_ci` with a Bonferroni-adjusted p-value.
    Returns a list of ``[category, table, stats_dict]`` entries.
    """
    result_pd = result_pd[result_pd[column_name].notna()]

    if skip_binary:
        feats_pd = result_pd[[column_name]]
    else:
        feats_pd = convert_to_binary(result_pd, column_name, continuous, range_list)

    feats_pd.columns = feats_pd.columns.str.replace(" ", "_", regex=True)
    feats_pd.columns = feats_pd.columns.str.replace(r"\.", "_", regex=True)

    all_condition_list = feats_pd.columns.to_list()

    result_pd = result_pd.reset_index(drop=True)
    merged_pd = pd.concat([result_pd[[label]], feats_pd], axis=1)
    merged_pd = merged_pd.replace("-", None)
    merged_pd[label] = merged_pd[label].astype(float)

    hits_pd = merged_pd.query(f"`{label}` > @zvalue").copy()
    no_hits_pd = merged_pd.query(f"`{label}` <= @zvalue").copy()
    print(f"size_of_hits_pd:{column_name},{len(hits_pd)}")
    print(f"size_of_non-hits_pd:{column_name},{len(no_hits_pd)}")

    final_list = []
    for each_case in all_condition_list:
        non_case_list = all_condition_list.copy()
        non_case_list.remove(each_case)

        a_pd = hits_pd[hits_pd[each_case] == 1]
        b_pd = no_hits_pd[no_hits_pd[each_case] == 1]
        c_pd = pd.concat([hits_pd[hits_pd[nc] == 1] for nc in non_case_list]) \
            if non_case_list else pd.DataFrame()
        d_pd = pd.concat([no_hits_pd[no_hits_pd[nc] == 1] for nc in non_case_list]) \
            if non_case_list else pd.DataFrame()

        table = np.array([[len(a_pd), len(b_pd)], [len(c_pd), len(d_pd)]])
        odds_result_dict = calculate_odds_ratio_with_ci(table)
        if odds_result_dict["p_value"]:
            odds_result_dict["adjusted_p_value"] = np.minimum(
                odds_result_dict["p_value"] * len(all_condition_list), 1)
        else:
            odds_result_dict["adjusted_p_value"] = 1
        final_list.append([each_case, table, odds_result_dict])

    return final_list


def _process_results(result_list: list, result: list, feature_name: str) -> None:
    """Flatten one feature's OR test entries into ``result_list`` (NB2 helper)."""
    if not result:
        return
    for entry in result:
        stats = entry[2]
        result_list.append({
            "FeatureName": feature_name,
            "Category": entry[0],
            "odds_ratio": stats["odds_ratio"],
            "ci_low": stats["ci_low"],
            "ci_high": stats["ci_high"],
            "p_value": stats["p_value"],
            "adjusted_p_value": stats["adjusted_p_value"],
        })


def get_mutation_level_OR_table_tsv(
    mutation_level_pd: pd.DataFrame,
    out_tsv: Union[str, Path] = "odds_ratio_test_2.0_ALL.tsv",
    label: str = "max_lowGFP-unsorted_Zscore_controls",
) -> pd.DataFrame:
    """Run the full odds-ratio test battery and write the result TSV.

    Ports NB2's ``get_mutation_level_OR_table_tsv``: categorical features,
    ASA/AlphaMissense/energy range features, contacting-residue interaction
    features and the PSSM conservation feature, all one-vs-rest at a z=2.0 hit
    threshold. Writes ``out_tsv`` and returns the assembled DataFrame.
    """
    sequence_feature = ['wild_to_mutant', 'wild_to_mutant_property_change',
                        'wild_amino_acid', 'wild_amino_acid_properties',
                        'mutant_amino_acid', 'mutant_amino_acid_properties']
    energy_based_feature = ['mcsm_ppi_ddG_AF3', 'thermompnn_ddg_AF3',
                            'FoldX_AVG_total_energy_AF3']
    conservation_feature = ['PSSM_BLOSUM62']

    cat = '6WKR'
    interaction_columns = [
        'contacting_residues', f'Num_of_SUZ12_contacting_residues_{cat}',
        f'Num_of_EZH2_contacting_residues_{cat}', f'Num_of_EED_contacting_residues_{cat}',
        f'Num_of_RBBP4_contacting_residues_{cat}', f'Num_of_JARID2_contacting_residues_{cat}',
        f'Num_of_AEBP2_contacting_residues_{cat}',
    ]
    structural_feature = [
        '3_region1', '3_region2', 'core_surface_monomer_TF', 'core_surface_complex_TF',
        'PPI_residue_TF', 'PPI_residue', f'Secondary Structure {cat} (DSSP 3 states)',
        f'Normalized_ASA_{cat}_Monomer(%)', f'Normalized_ASA_PRC2_{cat}(%)',
        'contacting_residues',
    ] + interaction_columns

    feature_categorizer_dict = {
        'sequence_feature': sequence_feature,
        'conservation_feature': conservation_feature,
        'energy_based_feature': energy_based_feature,
        'structural_feature': structural_feature,
    }
    categorical_feature_list = [
        '3_region1', 'PPI_residue_TF', 'core_surface_monomer_TF', 'core_surface_complex_TF',
        'PPI_residue', 'wild_to_mutant', 'wild_to_mutant_property_change',
        'wild_amino_acid', 'wild_amino_acid_properties', 'mutant_amino_acid',
        'mutant_amino_acid_properties', 'Secondary Structure 6WKR (DSSP 3 states)',
    ]
    zval_threshold = 2.0
    result_list: list = []

    for each_column in categorical_feature_list:
        result = run_odds_ratio_test_on_all_conditions(
            mutation_level_pd, each_column, zval_threshold, label=label)
        _process_results(result_list, result, each_column)

    range_tests = {
        'Normalized_ASA_6WKR_Monomer(%)': ['0~5', '5~25', '25~50', '50~75', '75~100.1'],
        'Normalized_ASA_PRC2_6WKR(%)': ['0~5', '5~25', '25~50', '50~75', '75~100.1'],
        'AlphaMissense': ['0~0.34', '0.34~0.56', '0.56~1.01'],
        'mcsm_ppi_ddG_AF3': ['-100~-1', '-1~1.0001', '1.0001~100'],
        'thermompnn_ddg_AF3': ['-100~-1', '-1~1.0001', '1.0001~100'],
        'FoldX_AVG_total_energy_AF3': ['-100~-1', '-1~1.0001', '1.0001~100'],
    }
    for column, ranges in range_tests.items():
        result = run_odds_ratio_test_on_all_conditions(
            mutation_level_pd, column, zval_threshold, range_list=ranges, label=label)
        _process_results(result_list, result, column)

    for column in interaction_columns:
        result = run_odds_ratio_test_on_all_conditions(
            mutation_level_pd, column, zval_threshold,
            range_list=['0~1.0001', '1.0001~100'], label=label)
        _process_results(result_list, result, column)

    for column in conservation_feature:
        result = run_odds_ratio_test_on_all_conditions(
            mutation_level_pd, column, zval_threshold,
            range_list=['-100~-1', '-1~1.0001', '1.0001~100'], label=label)
        _process_results(result_list, result, column)

    df = pd.DataFrame(result_list)
    expected_columns = ['FeatureName', 'Category', 'odds_ratio', 'ci_low',
                        'ci_high', 'p_value', 'adjusted_p_value']
    df = df.reindex(columns=expected_columns)
    df = df.sort_values(by=['FeatureName', 'Category'], ascending=True)
    for ind, row in df.iterrows():
        for featureCategory, featureList in feature_categorizer_dict.items():
            if row['FeatureName'] in featureList:
                df.loc[ind, 'FeatureCategory'] = featureCategory

    out_tsv = Path(out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_tsv, sep='\t', index_label='index')
    return df


# ---------------------------------------------------------------------------
# Forest-plot remap dictionaries (NB2 cell 14; E-15 duplicate key removed)
# ---------------------------------------------------------------------------

FEATURE_NAME_REMAP = {
    '3_region1': 'Structural region',
    'FoldX_AVG_total_energy_AF3': 'FoldX',
    'Normalized_ASA_6WKR_Monomer(%)': 'Normalized ASA (6WKR, Monomer)',
    'Normalized_ASA_PRC2_6WKR(%)': 'Normalized ASA (6WKR, Complex)',
    'Num_of_AEBP2_contacting_residues_6WKR': 'Contacting residue',
    'Num_of_EED_contacting_residues_6WKR': 'Contacting residue',
    'Num_of_EZH2_contacting_residues_6WKR': 'Contacting residue',
    'Num_of_JARID2_contacting_residues_6WKR': 'Contacting residue',
    'Num_of_RBBP4_contacting_residues_6WKR': 'Contacting residue',
    'Num_of_SUZ12_contacting_residues_6WKR': 'Contacting residue',
    'PPI_residue': 'Contacting residue',
    'PPI_residue_TF': 'Contacting residue',
    'PSSM_BLOSUM62': 'PSSM conservation',
    'Secondary Structure 6WKR (DSSP 3 states)': 'Secondary structure',
    'mcsm_ppi_ddG_AF3': 'mCSM-PPI',
    'thermompnn_ddg_AF3': 'ThermoMPNN',
    'wild_amino_acid': 'Wild amino acid',
    'mutant_amino_acid': 'Mutant amino acid',
    'contacting_residues': 'Contacting residue',
}

CATEGORY_REMAP = {
    'core': 'Core', 'interface': 'Interface', 'surface': 'Surface',
    'unstructured/noASA': 'Unstructured/noASA',
    'FoldX_AVG_total_energy_AF3_-100~-1': 'Stabilizing (< -1)',
    'FoldX_AVG_total_energy_AF3_-1~1': 'Neutral (-1 to 1)',
    'FoldX_AVG_total_energy_AF3_1~100': 'Destabilizing (> 1)',
    'Normalized_ASA_6WKR_Monomer(%)_0~5': 'Core (<5)%',
    'Normalized_ASA_6WKR_Monomer(%)_5~25': 'Buried (5-25)%',
    'Normalized_ASA_6WKR_Monomer(%)_25~50': 'Medium-buried (25-50)%',
    'Normalized_ASA_6WKR_Monomer(%)_50~75': 'Medium-exposed (50-75)%',
    'Normalized_ASA_6WKR_Monomer(%)_75~100_1': 'Exposed (>75)%',
    'Normalized_ASA_PRC2_6WKR(%)_0~5': 'Core (<5)%',
    'Normalized_ASA_PRC2_6WKR(%)_5~25': 'Buried (5-25)%',
    'Normalized_ASA_PRC2_6WKR(%)_25~50': 'Medium-buried (25-50)%',
    'Normalized_ASA_PRC2_6WKR(%)_50~75': 'Medium-exposed (50-75)%',
    'Normalized_ASA_PRC2_6WKR(%)_75~100_1': 'Exposed (>75)%',
    'Num_of_AEBP2_contacting_residues_6WKR_1~100': 'Residues contacting with AEBP2',
    'Num_of_EED_contacting_residues_6WKR_1~100': 'Residues contacting with EED',
    'Num_of_EZH2_contacting_residues_6WKR_1~100': 'Residues contacting with EZH2',
    'Num_of_JARID2_contacting_residues_6WKR_1~100': 'Residues contacting with JARID2',
    'Num_of_RBBP4_contacting_residues_6WKR_1~100': 'Residues contacting with RBBP4',
    'Num_of_SUZ12_contacting_residues_6WKR_1~100': 'Residues contacting with SUZ12',
    'Num_of_AEBP2_contacting_residues_6WKR_0~1': 'Residues not contacting with AEBP2',
    'Num_of_EED_contacting_residues_6WKR_0~1': 'Residues not contacting with EED',
    'Num_of_EZH2_contacting_residues_6WKR_0~1': 'Residues not contacting with EZH2',
    'Num_of_JARID2_contacting_residues_6WKR_0~1': 'Residues not contacting with JARID2',
    'Num_of_RBBP4_contacting_residues_6WKR_0~1': 'Residues not contacting with RBBP4',
    'Num_of_SUZ12_contacting_residues_6WKR_0~1': 'Residues not contacting with SUZ12',
    'PSSM_BLOSUM62_-100~-1': 'Unfavorable',
    'PSSM_BLOSUM62_-1~1_0001': 'Neutral',
    'PSSM_BLOSUM62_1_0001~100': 'Favorable',
    'B': 'B (strand)', 'C': 'C (loop/coil)', 'H': 'H (helix)',
    'contacting_residues_0~1': 'Non-contacting residues',
    'contacting_residues_1~100': 'All contacting residues',
    'mcsm_ppi_ddG_AF3_-1~1_0001': 'Neutral',
    'mcsm_ppi_ddG_AF3_-100~-1': 'Stabilizing',
    'mcsm_ppi_ddG_AF3_1_0001~100': 'Destabilizing',
    'thermompnn_ddg_AF3_-1~1_0001': 'Neutral',
    'thermompnn_ddg_AF3_-100~-1': 'Stabilizing',
    'thermompnn_ddg_AF3_1_0001~100': 'Destabilizing',
    'Aspartic_Acid': 'Aspartic Acid', 'Glutamic_Acid': 'Glutamic Acid',
    'EZH2': 'Contacting residues in EZH2', 'SUZ12': 'Contacting residues in SUZ12',
    'EED': 'Contacting residues in EED',
}

# E-15: duplicate ("Structural region", "Interface") key removed.
MAIN_ORDER_MAP = {
    ("Wild amino acid", "Alanine"): 1,
    ("Wild amino acid", "Cysteine"): 2,
    ("Wild amino acid", "Phenylalanine"): 3,
    ("Wild amino acid", "Tryptophan"): 4,
    ("Mutant amino acid", "Arginine"): 5,
    ("Mutant amino acid", "Asparagine"): 6,
    ("Mutant amino acid", "Proline"): 7,
    ("PSSM conservation", "Unfavorable"): 8,
    ("PSSM conservation", "Neutral"): 9,
    ("PSSM conservation", "Favorable"): 10,
    ("Structural region", "Core"): 11,
    ("Structural region", "Interface"): 12,
    ("Structural region", "Surface"): 13,
    ("Structural region", "Unstructured/noASA"): 14,
}

SUPPL_WILD_ORDER_MAP = {("Wild amino acid", aa): i for i, aa in enumerate([
    "Alanine", "Arginine", "Asparagine", "Aspartic_Acid", "Cysteine",
    "Glutamic_Acid", "Glutamine", "Glycine", "Histidine", "Isoleucine",
    "Leucine", "Lysine", "Methionine", "Phenylalanine", "Proline",
    "Serine", "Threonine", "Tryptophan", "Tyrosine", "Valine"], start=1)}

SUPPL_MUTANT_ORDER_MAP = {("Mutant amino acid", aa): i for i, aa in enumerate([
    "Alanine", "Arginine", "Asparagine", "Aspartic_Acid", "Cysteine",
    "Glutamic_Acid", "Glutamine", "Glycine", "Histidine", "Isoleucine",
    "Leucine", "Lysine", "Methionine", "Phenylalanine", "Proline",
    "Serine", "Threonine", "Tryptophan", "Tyrosine", "Valine"], start=1)}

SUPPL_STRUCTURAL_FEATURES_ORDER_MAP = {
    ("Secondary structure", "B (strand)"): 1,
    ("Secondary structure", "C (loop/coil)"): 2,
    ("Secondary structure", "H (helix)"): 3,
    ("Structural region", "Core"): 4,
    ("Structural region", "Interface"): 5,
    ("Structural region", "Surface"): 6,
    ("Structural region", "Unstructured/noASA"): 7,
    ("Normalized ASA (6WKR, Monomer)", "Core (<5)%"): 8,
    ("Normalized ASA (6WKR, Monomer)", "Buried (5-25)%"): 9,
    ("Normalized ASA (6WKR, Monomer)", "Medium-buried (25-50)%"): 10,
    ("Normalized ASA (6WKR, Monomer)", "Medium-exposed (50-75)%"): 11,
    ("Normalized ASA (6WKR, Monomer)", "Exposed (>75)%"): 12,
    ("Normalized ASA (6WKR, Complex)", "Core (<5)%"): 13,
    ("Normalized ASA (6WKR, Complex)", "Buried (5-25)%"): 14,
    ("Normalized ASA (6WKR, Complex)", "Medium-buried (25-50)%"): 15,
    ("Normalized ASA (6WKR, Complex)", "Medium-exposed (50-75)%"): 16,
    ("Normalized ASA (6WKR, Complex)", "Exposed (>75)%"): 17,
}


def _prep_forest_df(df: pd.DataFrame, feature_name_remap, category_remap,
                    or_pthreshold: float) -> pd.DataFrame:
    """Shared forest-plot preprocessing: remap, log2-transform, color-code."""
    df = df.replace(feature_name_remap).replace(category_remap)
    df['name_for_plot'] = '(' + df['FeatureName'] + ') ' + df['Category']
    df['log_OR'] = np.log2(df['odds_ratio'])
    df['log_ci_low'] = np.log2(df['ci_low'])
    df['log_ci_high'] = np.log2(df['ci_high'])
    df['log_error_low'] = df['log_OR'] - df['log_ci_low']
    df['log_error_high'] = df['log_ci_high'] - df['log_OR']
    df.replace([np.inf, -np.inf], np.nan, inplace=True)

    def assign_color(row):
        if pd.notna(row['odds_ratio']) and pd.notna(row['adjusted_p_value']):
            if row['odds_ratio'] < 1 and row['adjusted_p_value'] < or_pthreshold:
                return "#005AB5"
            if row['odds_ratio'] > 1 and row['adjusted_p_value'] < or_pthreshold:
                return "#DC3220"
            return "#B5ADAD"
        return "white"

    df["marker_color"] = df.apply(assign_color, axis=1)
    return df


def plot_all_odds_ratio_manuscript_mainfigure(
    odds_ratio_test_result_tsv: Union[str, Path],
    out_svg: Union[str, Path] = "Forest_plot_Main.svg",
    feature_name_remap: Dict = FEATURE_NAME_REMAP,
    category_remap: Dict = CATEGORY_REMAP,
    main_order_map: Dict = MAIN_ORDER_MAP,
    or_pthreshold: float = 0.05,
) -> None:
    """Multi-panel main-figure forest plot of log2 odds ratios with CIs."""
    df = pd.read_csv(odds_ratio_test_result_tsv, sep='\t')
    df = _prep_forest_df(df, feature_name_remap, category_remap, or_pthreshold)

    selected = ['Wild amino acid', 'Mutant amino acid', 'PSSM conservation', 'Structural region']
    df = df[df["FeatureName"].isin(selected)].reset_index(drop=True)
    df['Order'] = df.apply(lambda r: main_order_map.get((r['FeatureName'], r['Category']), np.nan), axis=1)
    df = df[df['Order'].notna()].sort_values(by='Order').reset_index(drop=True)

    unique_features = df['FeatureName'].unique()
    n_features = len(unique_features)
    height_ratios = [len(df[df['FeatureName'] == f]) for f in unique_features]
    fig_height = 0.25 * len(df) + 1

    fig, axs = plt.subplots(nrows=n_features, ncols=1, figsize=(2.5, fig_height),
                            sharex=True, height_ratios=height_ratios, dpi=100)
    if n_features == 1:
        axs = [axs]
    plt.subplots_adjust(hspace=0.5)

    new_feature_name_dict = {
        'Wild amino acid': 'Initial residue',
        'Mutant amino acid': 'Resulting residue',
        'PSSM conservation': 'PSSM score of resulting residue',
        'Structural region': 'Structural region',
    }
    for j, feature in enumerate(unique_features):
        gcont = df[df['FeatureName'] == feature].reset_index(drop=True)
        ax = axs[j]
        import seaborn as sns
        sns.pointplot(data=gcont, x='log_OR', y='name_for_plot', dodge=True,
                      capsize=0.2, err_kws={'linewidth': 1.5}, linestyle='none',
                      markersize=4, ax=ax)
        for idx, row in gcont.iterrows():
            if pd.notna(row['log_OR']):
                ax.scatter(x=row['log_OR'], y=idx, color=row["marker_color"], s=30, zorder=3)
            if pd.notna(row['log_error_low']) and pd.notna(row['log_error_high']):
                ax.errorbar(x=row['log_OR'], y=idx,
                            xerr=[[row['log_error_low']], [row['log_error_high']]],
                            elinewidth=0.5, fmt='none', color='black', capsize=1.0, capthick=0.3)
        ax.axvline(x=0, linestyle="--", color="gray", linewidth=0.7)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_title(new_feature_name_dict[feature], fontsize=7, pad=2)
        ax.set_xticks([-4, -2, 0, 2, 4])
        ax.margins(y=0.15)

    out_svg = Path(out_svg)
    out_svg.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_svg, bbox_inches='tight')
    plt.close(fig)


def plot_group(df_group: pd.DataFrame, order_map: Dict, out_svg: Union[str, Path],
               figsize=(4, 3.5), xlim=(-6, 6)) -> None:
    """One forest-plot panel for a subset of features (NB2 ``plot_group``)."""
    import seaborn as sns
    df_group = df_group.copy()
    df_group['Order'] = df_group.apply(
        lambda r: order_map.get((r['FeatureName'], r['Category']), np.nan), axis=1)
    df_group = df_group[df_group['Order'].notna()].sort_values(by='Order').reset_index(drop=True)

    fig, ax = plt.subplots(figsize=figsize, dpi=1000)
    sns.pointplot(data=df_group, x='log_OR', y='name_for_plot', dodge=True,
                  capsize=0.2, err_kws={'linewidth': 1.5}, linestyle='none',
                  markersize=4, ax=ax)
    for i, row in df_group.iterrows():
        if pd.notna(row['log_OR']):
            ax.scatter(x=row['log_OR'], y=i, color=row["marker_color"], s=30, zorder=3)
        if pd.notna(row['log_error_low']) and pd.notna(row['log_error_high']):
            ax.errorbar(x=row['log_OR'], y=i,
                        xerr=[[row['log_error_low']], [row['log_error_high']]],
                        elinewidth=0.5, fmt='none', color='black', capsize=1.0, capthick=0.3)
    ax.axvline(x=0, linestyle="--", color="gray", linewidth=0.7)
    ax.set_xlim(*xlim)
    ax.set_xticks([-4, -2, 0, 2, 4])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("")
    ax.tick_params(axis='y', labelsize=6)
    plt.tight_layout()
    out_svg = Path(out_svg)
    out_svg.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_svg)
    plt.close(fig)


def plot_all_odds_ratio_manuscript_suppl(
    odds_ratio_test_result_tsv: Union[str, Path],
    out_dir: Union[str, Path] = ".",
    feature_name_remap: Dict = FEATURE_NAME_REMAP,
    category_remap: Dict = CATEGORY_REMAP,
    suppl_wild_order_map: Dict = SUPPL_WILD_ORDER_MAP,
    suppl_mutant_order_map: Dict = SUPPL_MUTANT_ORDER_MAP,
    suppl_structural_features_order_map: Dict = SUPPL_STRUCTURAL_FEATURES_ORDER_MAP,
    or_pthreshold: float = 0.05,
) -> None:
    """Supplemental forest plots: wild-AA, mutant-AA and structural-feature groups."""
    df = pd.read_csv(odds_ratio_test_result_tsv, sep='\t', quotechar='"')
    df = _prep_forest_df(df, feature_name_remap, category_remap, or_pthreshold)
    selected = ['Wild amino acid', 'Mutant amino acid', 'Structural region',
                'Normalized ASA (6WKR, Monomer)', 'Normalized ASA (6WKR, Complex)',
                'Contacting residue', 'Secondary structure', 'FoldX', 'ThermoMPNN',
                'mCSM-PPI', 'PSSM conservation']
    df = df[df["FeatureName"].isin(selected)].copy().reset_index(drop=True)

    out_dir = Path(out_dir)
    plot_group(df[df['FeatureName'] == 'Wild amino acid'], suppl_wild_order_map,
               out_dir / 'Forest_plot_Wild_AA.svg')
    plot_group(df[df['FeatureName'] == 'Mutant amino acid'], suppl_mutant_order_map,
               out_dir / 'Forest_plot_Mutant_AA.svg')
    struct_features = ["Secondary structure", "Structural region",
                       "Normalized ASA (6WKR, Monomer)", "Normalized ASA (6WKR, Complex)"]
    plot_group(df[df['FeatureName'].isin(struct_features)],
               suppl_structural_features_order_map,
               out_dir / 'Forest_plot_Structural_features.svg', figsize=(2, 3.5))


# ---------------------------------------------------------------------------
# Histogram / dist plots (E-13: the second, cell-23 definitions)
# ---------------------------------------------------------------------------

def _compute_bin_edges(df: pd.DataFrame, col: str, bins, axis: AxisLabelOpts) -> np.ndarray:
    if isinstance(bins, (np.ndarray, list)):
        return np.asarray(bins)
    if isinstance(bins, str):
        return np.histogram_bin_edges(df[col].dropna().to_numpy(), bins=bins)
    if axis.xlim is not None:
        return np.linspace(axis.xlim[0], axis.xlim[1], int(bins) + 1)
    lo, hi = np.nanmin(df[col]), np.nanmax(df[col])
    return np.linspace(lo, hi, int(bins) + 1)


def _style_axes_grid(g, axis: AxisLabelOpts, bottom_only_xlabel: bool = True,
                     vline_x: Optional[float] = 2.0) -> None:
    if axis.title:
        g.fig.suptitle(axis.title, fontproperties=arial_font6, y=1.02)
    if axis.ylabel:
        g.set_ylabels(axis.ylabel)
    axes = g.axes.flat
    last_ax = axes[-1]
    for ax in axes:
        if axis.xlim:
            ax.set_xlim(*axis.xlim)
        if axis.ylim:
            ax.set_ylim(*axis.ylim)
        if axis.yticks is not None:
            ax.set_yticks(axis.yticks)
        if axis.xticks is not None:
            ax.xaxis.set_major_locator(mticker.FixedLocator(axis.xticks))
        ax.tick_params(axis='x', which='both', bottom=True)
        for lbl in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
            lbl.set_fontproperties(arial_font6)
        ax.spines['left'].set_position(('outward', 3))
        ax.spines['left'].set_linewidth(axis.linewidth or 0.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if vline_x is not None:
            ax.axvline(x=vline_x, color='black', linestyle='--', linewidth=0.6)
    if bottom_only_xlabel:
        for ax in axes[:-1]:
            ax.tick_params(axis='x', labelbottom=False)
            ax.set_xlabel('')
        if axis.xlabel:
            last_ax.set_xlabel(axis.xlabel, fontproperties=arial_font6)
        last_ax.tick_params(axis='x', labelbottom=True)
    elif axis.xlabel:
        g.set_xlabels(axis.xlabel)


def plot_distplot_on_wild_to_mutant_merged_v4(
    df: pd.DataFrame,
    label: str = 'max_lowGFP-unsorted_Zscore_controls',
    axis: AxisLabelOpts = AXIS,
    out_svg: Optional[Union[str, Path]] = 'Dist_plot_Wild_to_Mutant_v2',
    bins=40, facet_by_series: bool = True, histtype: str = 'step',
    annotate_stats: bool = True, min_group_n: int = 5,
    add_legend: bool = False, save_png: bool = False,
    width_cm: float = 0.7,
    col_order_pref=('Cys to Arg', 'Leu to Pro', 'Trp to Arg', 'Others'),
) -> None:
    """Hit vs non-hit histograms faceted by wild->mutant substitution class.

    NB2 defined this function twice (errors.md E-13): a cell-21 version that
    was actually called (cell 22, the only call site in the notebook) and a
    cell-23 redefinition that came *after* that call and was therefore never
    executed for any real output. There is no old figure to check either
    version against (this whole odds-ratio/forest-plot subsection has no
    counterpart anywhere in ``previous_results/``). This port keeps the
    cell-21 behavior as the default (no legend, one output file, the narrower
    0.7cm panel width, ``Others`` sorted last) since that is what NB2 would
    actually have produced if run top to bottom; ``add_legend``/``save_png``
    let a caller opt into the cell-23 look if ever wanted.
    """
    import seaborn as sns
    df = df[[label, 'wild_to_mutant']].copy()
    remapper = {'CtoR': 'Cys to Arg', 'LtoP': 'Leu to Pro', 'WtoR': 'Trp to Arg'}
    df['wild_to_mutant'] = df['wild_to_mutant'].map(lambda x: remapper.get(x, 'Others'))
    df['hit_line'] = np.where(df[label] > 2, 'hit', 'non-hit')
    be = _compute_bin_edges(df, label, bins, axis)

    hue_order = ['non-hit', 'hit']
    palette = {'non-hit': '#005AB5', 'hit': '#DC3220'}
    linestyle_map = {'non-hit': '--', 'hit': '-'}

    def hist_step(data, x, color=None, **kwargs):
        ax = plt.gca()
        arr = data[x].dropna().to_numpy()
        if arr.size == 0:
            return
        key = data['hit_line'].iloc[0] if 'hit_line' in data and len(data) else 'hit'
        ax.hist(arr, bins=be, density=True, histtype=histtype, linewidth=0.7,
                linestyle=linestyle_map.get(key, '-'), color=color, alpha=1.0)

    col_order = [c for c in col_order_pref if c in df['wild_to_mutant'].unique()]
    for c in df['wild_to_mutant'].unique():
        if c not in col_order:
            col_order.append(c)
    width_in, height_in = width_cm / 2.54, 8.4 / 2.54

    if facet_by_series:
        g = sns.FacetGrid(df, col='wild_to_mutant', hue='hit_line', col_order=col_order,
                          hue_order=hue_order, height=height_in, aspect=0.1, sharex=True,
                          sharey=True, palette=palette, col_wrap=1, legend_out=False)
        g.map_dataframe(hist_step, x=label)
        if add_legend:
            g.add_legend(title='')
        if annotate_stats:
            for series, subdf in df.groupby('wild_to_mutant'):
                if series not in col_order:
                    continue
                ax_idx = col_order.index(series)
                if ax_idx >= len(g.axes.flat):
                    continue
                ax = g.axes.flat[ax_idx]
                hv = subdf.loc[subdf['hit_line'] == 'hit', label].dropna()
                nv = subdf.loc[subdf['hit_line'] == 'non-hit', label].dropna()
                if len(hv) >= min_group_n and len(nv) >= min_group_n:
                    _, p = mannwhitneyu(hv, nv, alternative='two-sided')
                    ax.text(0.98, 0.98,
                            f"hit:{hv.median():.2f}, non-hit:{nv.median():.2f}\n"
                            f"Δmedian={hv.median() - nv.median():.2f}\nM-W p={p:.1e}",
                            ha='right', va='top', transform=ax.transAxes,
                            fontsize=5, fontproperties=arial_font6)
    else:
        g = sns.FacetGrid(df, hue='hit_line', hue_order=hue_order, height=height_in,
                          aspect=0.2, sharex=True, sharey=True, palette=palette,
                          col_wrap=1, legend_out=False)
        g.map_dataframe(hist_step, x=label)
        if add_legend:
            g.add_legend(title='')

    _style_axes_grid(g, axis, bottom_only_xlabel=True, vline_x=2.0)
    g.fig.set_size_inches(width_in * max(1, len(col_order) if facet_by_series else 1),
                          height_in, forward=True)
    plt.tight_layout()
    if out_svg:
        out_svg = Path(out_svg)
        out_svg.parent.mkdir(parents=True, exist_ok=True)
        g.fig.savefig(f'{out_svg}.svg', format='svg', dpi=1200, bbox_inches='tight')
        if save_png:
            g.fig.savefig(f'{out_svg}.png', format='png', dpi=1200, bbox_inches='tight')
    plt.close(g.fig)


def plot_distplot_on_structural_regions_v4(
    df: pd.DataFrame,
    label: str = 'max_lowGFP-unsorted_Zscore_controls',
    axis: AxisLabelOpts = AXIS,
    out_svg: Optional[Union[str, Path]] = 'Histogram_structural_region',
    bins=40, histtype: str = 'step', annotate_stats: bool = True, min_group_n: int = 5,
) -> None:
    """Hit vs non-hit histograms faceted by structural region (core/interface/...)."""
    import seaborn as sns
    df = df[[label, '3_region1']].copy()
    region_remap = {'core': 'Core', 'interface': 'Interface', 'surface': 'Surface',
                    'unstructured/noASA': 'Unstructured/noASA'}
    df['3_region1'] = df['3_region1'].map(lambda x: region_remap.get(x, 'Unstructured/noASA'))
    df['hit_line'] = np.where(df[label] > 2, 'hit', 'non-hit')
    be = _compute_bin_edges(df, label, bins, axis)

    hue_order = ['non-hit', 'hit']
    palette = {'non-hit': '#005AB5', 'hit': '#DC3220'}
    linestyle_map = {'non-hit': '--', 'hit': '-'}

    def hist_step(data, x, color=None, **kwargs):
        ax = plt.gca()
        arr = data[x].dropna().to_numpy()
        if arr.size == 0:
            return
        key = data['hit_line'].iloc[0] if 'hit_line' in data and len(data) else 'hit'
        ax.hist(arr, bins=be, density=True, histtype=histtype, linewidth=0.7,
                linestyle=linestyle_map.get(key, '-'), color=color, alpha=1.0)

    col_order_pref = ['Core', 'Interface', 'Surface', 'Unstructured/noASA']
    col_order = [c for c in col_order_pref if c in df['3_region1'].unique()]
    for c in df['3_region1'].unique():
        if c not in col_order:
            col_order.append(c)
    width_in, height_in = .7 / 2.54, 8.4 / 2.54

    g = sns.FacetGrid(df, col='3_region1', hue='hit_line', col_order=col_order,
                      hue_order=hue_order, height=height_in, aspect=0.1, sharex=True,
                      sharey=True, palette=palette, col_wrap=1, legend_out=False)
    g.map_dataframe(hist_step, x=label)
    if annotate_stats:
        for region, subdf in df.groupby('3_region1'):
            if region not in col_order:
                continue
            ax_idx = col_order.index(region)
            if ax_idx >= len(g.axes.flat):
                continue
            ax = g.axes.flat[ax_idx]
            hv = subdf.loc[subdf['hit_line'] == 'hit', label].dropna()
            nv = subdf.loc[subdf['hit_line'] == 'non-hit', label].dropna()
            if len(hv) >= min_group_n and len(nv) >= min_group_n:
                _, p = mannwhitneyu(hv, nv, alternative='two-sided')
                ax.text(0.98, 0.98,
                        f"hit:{hv.median():.2f}, non-hit:{nv.median():.2f}\n"
                        f"Δmedian={hv.median() - nv.median():.2f}\nM-W p={p:.1e}",
                        ha='right', va='top', transform=ax.transAxes,
                        fontsize=5, fontproperties=arial_font6)

    _style_axes_grid(g, axis, bottom_only_xlabel=True, vline_x=2.0)
    g.fig.set_size_inches(width_in * max(1, len(col_order)), height_in, forward=True)
    plt.tight_layout()
    if out_svg:
        out_svg = Path(out_svg)
        out_svg.parent.mkdir(parents=True, exist_ok=True)
        g.fig.savefig(f'{out_svg}.svg', format='svg', dpi=1200, bbox_inches='tight')
    plt.close(g.fig)


# ---------------------------------------------------------------------------
# Substitution matrix
# ---------------------------------------------------------------------------

def plot_sub_matrix_on_all_merged(
    odds_ratio_test_tsv: Union[str, Path],
    out_svg: Union[str, Path] = 'Substitution_matrix.svg',
) -> None:
    """Side-by-side 20x20 substitution-matrix heatmaps of odds ratio and significance."""
    import seaborn as sns
    odds_ratio_test_pd = pd.read_csv(odds_ratio_test_tsv, sep='\t')
    df = odds_ratio_test_pd[odds_ratio_test_pd['FeatureName'] == 'wild_to_mutant'].copy()
    df = df[['Category', 'odds_ratio', 'adjusted_p_value']].rename(
        columns={'Category': 'wild_to_mutant', 'odds_ratio': 'OR', 'adjusted_p_value': 'pval'})
    df['pval'] = df['pval'].astype(float).apply(lambda x: 1 if x < 0.05 else 0)
    df['initial_residue'] = df['wild_to_mutant'].str[0]
    df['resulting_residue'] = df['wild_to_mutant'].str[-1]

    pivot_or = df.pivot(index='initial_residue', columns='resulting_residue',
                        values='OR').reindex(index=AA_ORDER, columns=AA_ORDER)
    pivot_pval = df.pivot(index='initial_residue', columns='resulting_residue',
                          values='pval').reindex(index=AA_ORDER, columns=AA_ORDER)

    fig, ax = plt.subplots(figsize=(8, 4), ncols=2, dpi=300)
    white_to_red = LinearSegmentedColormap.from_list("white_red", ["white", "red"])
    sns.heatmap(pivot_or, square=True, linewidths=0.3, linecolor='gray', cmap=white_to_red,
                vmin=0, xticklabels=True, yticklabels=True, cbar_kws={'shrink': 0.3}, ax=ax[0])
    sns.heatmap(pivot_pval, square=True, linewidths=0.3, linecolor='gray', cmap="Blues",
                vmin=0, vmax=1, xticklabels=True, yticklabels=True,
                cbar_kws={'shrink': 0.3}, ax=ax[1])
    for a in ax:
        a.set_xticklabels(a.get_xticklabels(), rotation=0, fontsize=7)
        a.set_yticklabels(a.get_yticklabels(), rotation=0, fontsize=7)
    ax[0].set_title("All merged, OR", fontsize=8)
    ax[1].set_title("All merged, $p$-value", fontsize=8)

    highlight_cells = [
        ('K', 'R'), ('K', 'E'), ('K', 'G'), ('R', 'G'), ('H', 'R'), ('D', 'G'), ('E', 'G'),
        ('S', 'G'), ('S', 'P'), ('T', 'A'), ('N', 'D'), ('N', 'S'), ('N', 'G'), ('Q', 'R'),
        ('C', 'R'), ('V', 'A'), ('I', 'T'), ('I', 'V'), ('I', 'M'), ('L', 'S'), ('L', 'P'),
        ('M', 'T'), ('M', 'V'), ('F', 'S'), ('F', 'P'), ('F', 'L'), ('Y', 'H'), ('Y', 'C'),
        ('W', 'R'), ('R', 'K'), ('R', 'H'), ('R', 'Q'), ('R', 'C'), ('R', 'W'), ('H', 'Y'),
        ('D', 'N'), ('E', 'K'), ('S', 'N'), ('S', 'L'), ('S', 'F'), ('T', 'I'), ('T', 'M'),
        ('C', 'Y'), ('G', 'K'), ('G', 'R'), ('G', 'D'), ('G', 'E'), ('G', 'S'), ('G', 'N'),
        ('P', 'S'), ('P', 'L'), ('P', 'F'), ('A', 'T'), ('A', 'V'), ('V', 'I'), ('V', 'M'),
        ('L', 'F'), ('M', 'I'),
    ]
    for subplot, pivot in zip(ax, [pivot_or, pivot_pval]):
        for wild, mut in highlight_cells:
            if wild in pivot.index and mut in pivot.columns:
                i = pivot.index.get_loc(wild)
                j = pivot.columns.get_loc(mut)
                subplot.add_patch(patches.Rectangle((j, i), 1, 1, linewidth=0.8,
                                                    edgecolor='black', facecolor='none', clip_on=False))
        for i in range(len(pivot)):
            subplot.add_patch(patches.Rectangle((i, i), 1, 1, linewidth=0, edgecolor=None,
                                                facecolor='lightgray', zorder=3, alpha=0.3, clip_on=False))
        for i in range(pivot.shape[0]):
            for j in range(pivot.shape[1]):
                if i != j and pd.isna(pivot.iloc[i, j]):
                    subplot.text(j + 0.5, i + 0.5, '·', ha='center', va='center',
                                 color='gray', fontsize=6)
        subplot.set_xlim(0, len(pivot.columns))
        subplot.set_ylim(len(pivot.index), 0)

    plt.tight_layout()
    out_svg = Path(out_svg)
    out_svg.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_svg)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Quantile-stratified boxplot/stripplot
# ---------------------------------------------------------------------------

def get_quantile_cutoffs(df: pd.DataFrame, column_name: str, threshold: float = 1) -> Tuple:
    df = df.replace('-', None)
    filtered = df[df[column_name] > threshold][column_name]
    if filtered.empty:
        s = df[column_name].dropna()
        return tuple(s.quantile(q) if len(s) else np.nan for q in [0, 0.25, 0.5, 0.75, 1.0])
    return tuple(filtered.quantile(q).round(2) for q in [0.00, 0.25, 0.50, 0.75, 1.00])


def categorize_quantile(value, q0, q25, q50, q75, q100) -> Optional[str]:
    if pd.isna(value):
        return None
    if value < -1:
        return 'Stabilizing'
    if -1 <= value <= q0:
        return 'Neutral'
    if q0 <= value < q25:
        return 'Q1'
    if q25 <= value < q50:
        return 'Q2'
    if q50 <= value < q75:
        return 'Q3'
    return 'Q4'


def pval_to_stars(pval: float) -> str:
    if pval < 1e-4:
        return '****'
    if pval < 1e-3:
        return '***'
    if pval < 1e-2:
        return '**'
    if pval < 0.05:
        return '*'
    return 'n.s.'


def plot_model_quantile_boxplot_v2(
    all_pd: pd.DataFrame,
    label: str,
    output_filename: Union[str, Path],
    axis: AxisLabelOpts = AXIS,
    width_cm: float = 6.6, height_cm: float = 3.3,
    name_dict: Optional[Dict[str, str]] = None,
    model_order: Optional[List[str]] = None,
    model_order_aa: Optional[List[str]] = None,
    quantile_order: Optional[List[str]] = None,
    stripplot_kws: Optional[Dict] = None,
    boxplot_kws: Optional[Dict] = None,
) -> Optional[pd.DataFrame]:
    """Quantile-stratified boxplots with overlaid per-gene stripplots per model.

    Returns the Mann-Whitney comparison DataFrame (with star annotations).
    """
    import seaborn as sns
    if name_dict is None:
        name_dict = {'FoldX_AVG_total_energy_AF3': 'FoldX', 'thermompnn_ddg_AF3': 'ThermoMPNN'}
    if model_order is None:
        model_order = ['FoldX', 'ThermoMPNN']
    if model_order_aa is None:
        model_order_aa = ['FoldX_aa', 'ThermoMPNN_aa']
    if quantile_order is None:
        quantile_order = ['Q1', 'Q2', 'Q3', 'Q4']
    if stripplot_kws is None:
        stripplot_kws = {'hue': 'Gene Symbol',
                         'palette': {'EED': '#7CCDF4', 'EZH2': '#B2B0B0', 'SUZ12': '#C19EE0'},
                         'alpha': 0.5, 'size': 3.0, 'jitter': 0.25, 'dodge': False, 'legend': False}
    if boxplot_kws is None:
        boxplot_kws = {'boxprops': dict(color='#3f3f3f'), 'medianprops': dict(color='#3f3f3f'),
                       'whiskerprops': dict(color='#3f3f3f'), 'capprops': dict(color='#3f3f3f'),
                       'fill': False, 'showcaps': True, 'showfliers': False, 'whis': 1.5, 'linewidth': 0.5}

    cols_needed = ['Gene Symbol', label, 'editor'] + list(name_dict.values()) + model_order_aa
    all_pd = all_pd[[c for c in cols_needed if c in all_pd.columns]].copy().reset_index(drop=True)

    combined = []
    for _, model_name in name_dict.items():
        if model_name not in all_pd.columns:
            continue
        q0, q25, q50, q75, q100 = get_quantile_cutoffs(all_pd, model_name)
        temp = all_pd[['Gene Symbol', label, f'{model_name}_aa']].copy()
        temp['Quantile'] = all_pd[model_name].apply(lambda x: categorize_quantile(x, q0, q25, q50, q75, q100))
        temp['Model'] = model_name
        temp['Model_Quantile'] = temp['Model'] + '-' + temp['Quantile'].astype(str)
        temp['editor'] = all_pd['editor']
        combined.append(temp)
    if not combined:
        raise ValueError("None of the configured model columns were found in the dataframe.")

    combined_df = pd.concat(combined, ignore_index=True)
    combined_df = combined_df[combined_df['Quantile'].isin(quantile_order + ['Neutral', 'Stabilizing'])]
    combined_df['Amino_acid_edit'] = combined_df[model_order_aa].bfill(axis=1).iloc[:, 0]
    combined_df = combined_df.drop(model_order_aa, axis=1)

    x_order = []
    for m in model_order:
        x_order.append(f'{m}-Stabilizing')
        x_order.append(f'{m}-Neutral')
        for q in quantile_order:
            x_order.append(f'{m}-{q}')
    x_order = [x for x in x_order if x in combined_df['Model_Quantile'].unique()]

    comparison_results = []
    for m in model_order:
        dfm = combined_df[combined_df['Model'] == m]
        q4 = dfm[dfm['Quantile'] == 'Q4'][label].dropna()
        for grp in ['Q1', 'Q2', 'Q3', 'Neutral', 'Stabilizing']:
            gvals = dfm[dfm['Quantile'] == grp][label].dropna()
            if len(q4) >= 3 and len(gvals) >= 3:
                _, pval = mannwhitneyu(q4, gvals, alternative='two-sided')
                comparison_results.append({'Model': m, 'Group_vs_Q4': grp, 'Q4_N': len(q4),
                                           f'{grp}_N': len(gvals), 'P-value': pval})
        dest = dfm[dfm['Quantile'].isin(quantile_order)][label].dropna()
        neu = dfm[dfm['Quantile'] == 'Neutral'][label].dropna()
        sta = dfm[dfm['Quantile'] == 'Stabilizing'][label].dropna()
        if len(dest) >= 3 and len(neu) >= 3:
            _, pval = mannwhitneyu(dest, neu, alternative='two-sided')
            comparison_results.append({'Model': m, 'Group_vs_Destabilizing': 'Neutral',
                                       'Destabilizing_N': len(dest), 'Neutral_N': len(neu), 'P-value': pval})
        if len(dest) >= 3 and len(sta) >= 3:
            _, pval = mannwhitneyu(dest, sta, alternative='two-sided')
            comparison_results.append({'Model': m, 'Group_vs_Destabilizing': 'Stabilizing',
                                       'Destabilizing_N': len(dest), 'Stabilizing_N': len(sta), 'P-value': pval})
    comparison_df = pd.DataFrame(comparison_results) if comparison_results else None

    width_in, height_in = width_cm / 2.54, height_cm / 2.54
    fig, ax = plt.subplots(figsize=(width_in, height_in))
    for m in model_order:
        for marker, ed in (("o", "ABE"), ("D", "CBE")):
            sdf = combined_df[(combined_df['Model'] == m) & (combined_df['editor'] == ed)]
            if not sdf.empty:
                sdf = sdf.sample(frac=1, random_state=0).reset_index(drop=True)
                ax = sns.stripplot(data=sdf, x='Model_Quantile', y=label, order=x_order,
                                   zorder=1, marker=marker, ax=ax, **stripplot_kws)
            for c in ax.collections:
                c.set_rasterized(True)

    sns.boxplot(data=combined_df, x='Model_Quantile', y=label, order=x_order, ax=ax,
                zorder=2, **boxplot_kws)

    if axis.xlim:
        ax.set_xlim(*axis.xlim)
    if axis.ylim:
        ax.set_ylim(*axis.ylim)
    if axis.yticks is not None:
        ax.set_yticks(axis.yticks)

    def _pretty_xtick(xlab: str) -> str:
        model, quant = xlab.split('-', 1)
        return quant if quant.startswith('Q') else model
    ax.set_xticklabels([_pretty_xtick(x) for x in x_order], rotation=0)
    ax.set_xlabel(axis.xlabel or "")
    ax.set_ylabel(axis.ylabel or "")
    for lbl in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
        lbl.set_fontproperties(arial_font6)
    ax.spines['left'].set_position(('outward', 3))
    ax.spines['left'].set_linewidth(axis.linewidth or 0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if axis.title:
        fig.suptitle(axis.title, fontproperties=arial_font6, y=1.02)

    plt.tight_layout()
    output_filename = Path(output_filename)
    output_filename.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f'{output_filename}.svg', format='svg', dpi=1200,
                bbox_inches='tight', transparent=True)
    plt.close(fig)

    if comparison_df is not None:
        comparison_df['P-value'] = comparison_df['P-value'].apply(pval_to_stars)
    return comparison_df


# ---------------------------------------------------------------------------
# Split-violin plot of stability predictors (no be_scan equivalent)
# ---------------------------------------------------------------------------

def with_alpha(color, alpha):
    rgba = list(mcolors.to_rgba(color))
    rgba[-1] = alpha
    return tuple(rgba)


def get_violin_plot_sgrna_level_v2(
    mapped_pd: pd.DataFrame,
    editor_type: Optional[str] = None,
    axis: AxisLabelOpts = AXIS,
    out_svg: Optional[Union[str, Path]] = 'Violin_plot.svg',
    width_cm: float = 4, height_cm: float = 3,
    annotate_p: bool = True,
    STABILIZING_COLOR: str = "#005AB5",
    DESTABILIZING_COLOR: str = "#DC3220",
    predictors: Optional[List[Tuple[str, str]]] = None,
) -> None:
    """Split-violin plots of dDDG predictors for hit vs non-hit variants per model.

    Colored by stability direction (blue stabilizing / red destabilizing) with
    per-predictor Mann-Whitney U tests. This function has no ``be_scan``
    equivalent and is preserved here.
    """
    import seaborn as sns
    if predictors is None:
        predictors = [('FoldX_AVG_total_energy_AF3', 'FoldX'),
                      ('thermompnn_ddg_AF3', 'ThermoMPNN')]

    df = mapped_pd.copy()
    df['hit'] = np.where(df['GFP-MAX'] > 2, 'hit', 'non-hit')
    if editor_type is not None:
        df = df[df['editor_type'] == editor_type]

    width_in, height_in = width_cm / 2.54, height_cm / 2.54
    fig, axes = plt.subplots(nrows=1, ncols=len(predictors), figsize=(width_in, height_in))
    if len(predictors) == 1:
        axes = [axes]

    for ax, (y_col, title) in zip(axes, predictors):
        if title not in df.columns:
            ax.text(0.5, 0.5, f'Missing:\n{title}', ha='center', va='center', fontsize=6)
            ax.set_axis_off()
            continue
        plot_df = df[['hit', title]].copy()
        plot_df['model'] = ''
        plot_df = plot_df.reset_index(drop=True)
        sns.violinplot(data=plot_df, x='model', y=title, hue='hit', split=True,
                       inner='quart', linewidth=0.7,
                       palette={'hit': with_alpha(DESTABILIZING_COLOR, 0.95),
                                'non-hit': with_alpha(STABILIZING_COLOR, 0.80)},
                       cut=0, hue_order=['non-hit', 'hit'], ax=ax)
        g_hit = plot_df[plot_df['hit'] == 'hit'][title].dropna()
        g_non = plot_df[plot_df['hit'] == 'non-hit'][title].dropna()
        med_destab = float(np.median(g_hit)) if len(g_hit) else np.nan
        med_stab = float(np.median(g_non)) if len(g_non) else np.nan
        print(f"{title}: median_destabilizing (hit) = {med_destab:.3f} [n={len(g_hit)}], "
              f"median_stabilizing (non-hit) = {med_stab:.3f} [n={len(g_non)}]")
        if len(g_hit) >= 3 and len(g_non) >= 3:
            _, pval = mannwhitneyu(g_hit, g_non, alternative='two-sided')
            stars = pval_to_stars(pval)
            print(f"    Mann-Whitney U: p = {pval:.3e} ({stars})")
            if annotate_p:
                ymax = np.nanpercentile(plot_df[title], 98)
                ytxt = ymax + 0.03 * (np.nanmax(plot_df[title]) - np.nanmin(plot_df[title]) + 1e-9)
                ax.text(0, ytxt, stars, ha='center', va='bottom', fontproperties=arial_font6)
        else:
            print("    Mann-Whitney U: insufficient data for testing.")

        ax.set_title(title, fontproperties=arial_font6)
        ax.set_xlabel("")
        if axis.ylabel and ax is axes[0]:
            ax.set_ylabel(axis.ylabel, fontproperties=arial_font6)
        else:
            ax.set_ylabel("")
        if title == "FoldX":
            ax.set_yticks([-10, 0, 10, 20, 30])
        elif title in ("ThermoMPNN", "mCSM-PPI"):
            ax.set_yticks([-3, 0, 3, 6])
        if axis.xlim:
            ax.set_xlim(*axis.xlim)
        if axis.ylim:
            ax.set_ylim(*axis.ylim)
        if axis.yticks is not None:
            ax.set_yticks(axis.yticks)
        ax.set_xticks([])
        for lbl in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
            lbl.set_fontproperties(arial_font6)
        ax.spines['left'].set_position(('outward', 3))
        ax.spines['left'].set_linewidth(axis.linewidth or 0.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        leg = ax.get_legend()
        if leg:
            leg.remove()

    if axis.title:
        fig.suptitle(axis.title, fontproperties=arial_font6, y=1.02)
    plt.tight_layout()
    if out_svg:
        out_svg = Path(out_svg)
        out_svg.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_svg, format='svg', dpi=1200, bbox_inches='tight', transparent=True)
    plt.close(fig)
