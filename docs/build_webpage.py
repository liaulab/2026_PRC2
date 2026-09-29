"""
Build docs/data + docs/manifest.json + the embedded MANIFEST_SEGMENTS
in docs/index.html from new_outputs_interactive/, with each sidebar
category split into "segments": a group of files that share one legend and
scheme, plus (for the group's files only) their own inline Plotly legend
patched off (layout.showlegend = False) so the shared header is the only
legend shown for that group. Segments with legend=None are left completely
untouched (byte-identical copy) -- this is how ClusteredHeatmaps, PWES, and
most Boxplots pass through unchanged.

Source of truth for "which files share a scheme" was determined by
parsing every file's embedded Plotly JSON and comparing legend entries
(name+color) -- see the accompanying investigation. Re-run this whenever
interactive_analysis.ipynb regenerates new_outputs_interactive/.
"""
import json
import os
import re
import shutil
import sys

sys.path.insert(0, os.path.dirname(__file__))
import plotly_json as pj

PROJECT = "/home/claude/project"
SRC = os.path.join(PROJECT, "new_outputs_interactive")
DOCS = os.path.join(PROJECT, "docs")
DATA = os.path.join(DOCS, "data")

# ---------------------------------------------------------------------------
# Legend schemes (label -> css color), determined empirically from the JSON.
# ---------------------------------------------------------------------------
MUT_TYPE_LEGEND = [
    {"label": "Missense", "color": "#fdb462"},
    {"label": "No Mutation/Silent", "color": "#80b1d3"},
    {"label": "Nonsense", "color": "#c03221"},
]
CORR_GENE_LEGEND = [
    {"label": "Control", "color": "#a7d5b9"},
    {"label": "EED", "color": "#7ccdf4"},
    {"label": "EZH2", "color": "#b2b0b0"},
    {"label": "Essential", "color": "#ffcc33"},
    {"label": "SUZ12", "color": "#c19ee0"},
]
BOXPLOT_GENE_LEGEND = [
    {"label": "AEBP2", "color": "#1fc1c1"},
    {"label": "Control", "color": "#f0c93d"},
    {"label": "EED", "color": "#7ccdf4"},
    {"label": "EZH2", "color": "#b2b0b0"},
    {"label": "Essential", "color": "#91caa9"},
    {"label": "JARID2", "color": "#ff9999"},
    {"label": "SUZ12", "color": "#c19ee0"},
]
LOF_ABE_LEGEND = [
    {"label": "Q575R", "color": "#ff5a5a"},
    {"label": "WT", "color": "#feb1c6"},
    {"label": "GFP-EED", "color": "#afffb9"},
    {"label": "GFP-EZH2", "color": "#5aa9e6"},
]
LOF_CBE_LEGEND = [
    {"label": "Q575R CBE39", "color": "#ff5a5a"},
    {"label": "WT CBE39", "color": "#ffd55a"},
    {"label": "WT CBE6", "color": "#feb1c6"},
    {"label": "GFP-EED", "color": "#afffb9"},
    {"label": "GFP-EZH2", "color": "#5aa9e6"},
]
Q575R_ABE_LOL_LEGEND = [{"label": "ABE", "color": "#feb1c6"}]
Q575R_CBE_LOL_LEGEND = [
    {"label": "CBE39", "color": "#feb1c6"},
    {"label": "CBE6", "color": "#5aa9e6"},
]

GENE_ORDER = ["EZH2", "EED", "SUZ12", "AEBP2", "JARID2"]


def lol(prefix, editor):
    return [f"{prefix}_{editor}_{g}_Lollipop_Split.html" for g in GENE_ORDER]


# ---------------------------------------------------------------------------
# Category -> list of segments. Each segment: {legend: [...] | None, files: [...]}
# folder = the new_outputs_interactive/<folder> subfolder the files live in.
# ---------------------------------------------------------------------------
CATEGORIES = {
    "WT Scatterplots": {
        "folder": "Scatterplots",
        "segments": [
            {"legend": MUT_TYPE_LEGEND, "files": [
                "EZH2-Karpas422_WT-ABE-CBE6b-DMSO-D0-Z-Scatterplot.html",
                "EED-Karpas422_WT-ABE-CBE6b-DMSO-D0-Z-Scatterplot.html",
                "SUZ12-Karpas422_WT-ABE-CBE6b-DMSO-D0-Z-Scatterplot.html",
                "AEBP2-Karpas422_WT-ABE-CBE6b-DMSO-D0-Z-Scatterplot.html",
                "JARID2-Karpas422_WT-ABE-CBE6b-DMSO-D0-Z-Scatterplot.html",
                "MTF2-Karpas422_WT-ABE-DMSO-D0-Z-Scatterplot.html",
            ]},
        ],
    },
    "Q575R Scatterplots": {
        "folder": "Scatterplots",
        "segments": [
            {"legend": MUT_TYPE_LEGEND, "files": [
                "EZH2-Karpas422_Q575R-ABE-CBE3-DMSO-PRC2i-Z-Scatterplot.html",
                "EZH2-Karpas422_Q575R-ABE-CBE3-PRC2i-D0-Z-Scatterplot5.html",
                "EED-Karpas422_Q575R-ABE-CBE3-DMSO-PRC2i-Z-Scatterplot.html",
                "EED-Karpas422_Q575R-ABE-CBE3-PRC2i-D0-Z-Scatterplot5.html",
                "SUZ12-Karpas422_Q575R-ABE-CBE3-DMSO-PRC2i-Z-Scatterplot.html",
                "SUZ12-Karpas422_Q575R-ABE-CBE3-PRC2i-D0-Z-Scatterplot5.html",
                "AEBP2-Karpas422_Q575R-ABE-CBE3-DMSO-PRC2i-Z-Scatterplot.html",
                "AEBP2-Karpas422_Q575R-ABE-CBE3-PRC2i-D0-Z-Scatterplot5.html",
                "JARID2-Karpas422_Q575R-ABE-CBE3-DMSO-PRC2i-Z-Scatterplot.html",
                "JARID2-Karpas422_Q575R-ABE-CBE3-PRC2i-D0-Z-Scatterplot5.html",
                "MTF2-Karpas422_Q575R-ABE-DMSO-D0-Z-Scatterplot.html",
                "MTF2-Karpas422_Q575R-ABE-DMSO-PRC2i-Z-Scatterplot.html",
                "MTF2-Karpas422_Q575R-ABE-PRC2i-D0-Z-Scatterplot.html",
            ]},
        ],
    },
    "EZH2 Stability Scatterplots": {
        "folder": "Scatterplots",
        "segments": [
            {"legend": MUT_TYPE_LEGEND, "files": [
                "EZH2-GFP_EZH2_EED_K562-ABE-lowGFP-unsorted-Z-Scatterplot.html",
                "EZH2-GFP_EZH2_EED_K562-CBE-lowGFP-unsorted-Z-Scatterplot.html",
                "EED-GFP_EZH2_EED_K562-ABE-lowGFP-unsorted-Z-Scatterplot.html",
                "EED-GFP_EZH2_EED_K562-CBE-lowGFP-unsorted-Z-Scatterplot.html",
                "SUZ12-GFP_EZH2_EED_K562-ABE-lowGFP-unsorted-Z-Scatterplot.html",
                "SUZ12-GFP_EZH2_EED_K562-CBE-lowGFP-unsorted-Z-Scatterplot.html",
                "AEBP2-GFP_EZH2_EED_K562-ABE-lowGFP-unsorted-Z-Scatterplot.html",
                "AEBP2-GFP_EZH2_EED_K562-CBE-lowGFP-unsorted-Z-Scatterplot.html",
                "JARID2-GFP_EZH2_EED_K562-ABE-lowGFP-unsorted-Z-Scatterplot.html",
                "JARID2-GFP_EZH2_EED_K562-CBE-lowGFP-unsorted-Z-Scatterplot.html",
            ]},
        ],
    },
    "EED Stability Scatterplots": {
        "folder": "Scatterplots",
        "segments": [
            {"legend": MUT_TYPE_LEGEND, "files": [
                "EZH2-GFP_EZH2_EED_K562-ABE-lowGFP-unsorted-Z-Scatterplot.html",
                "EZH2-GFP_EZH2_EED_K562-CBE-lowGFP-unsorted-Z-Scatterplot.html",
                "EED-GFP_EZH2_EED_K562-ABE-lowGFP-unsorted-Z-Scatterplot.html",
                "EED-GFP_EZH2_EED_K562-CBE-lowGFP-unsorted-Z-Scatterplot.html",
                "SUZ12-GFP_EZH2_EED_K562-ABE-lowGFP-unsorted-Z-Scatterplot.html",
                "SUZ12-GFP_EZH2_EED_K562-CBE-lowGFP-unsorted-Z-Scatterplot.html",
                "AEBP2-GFP_EZH2_EED_K562-ABE-lowGFP-unsorted-Z-Scatterplot.html",
                "AEBP2-GFP_EZH2_EED_K562-CBE-lowGFP-unsorted-Z-Scatterplot.html",
                "JARID2-GFP_EZH2_EED_K562-ABE-lowGFP-unsorted-Z-Scatterplot.html",
                "JARID2-GFP_EZH2_EED_K562-CBE-lowGFP-unsorted-Z-Scatterplot.html",
            ]},
        ],
    },
    "Correlation Scatterplots": {
        "folder": "Scatterplots",
        "segments": [
            {"legend": CORR_GENE_LEGEND, "files": [
                "DMSO-PRC2i-PRC2i-D0-ABE-correlation_scatterplot.html",
                "DMSO-PRC2i-PRC2i-D0-CBE-correlation_scatterplot.html",
                "GFP-EZH2-GFP-EED-ABE-correlation_scatterplot.html",
                "GFP-EZH2-GFP-EED-CBE-correlation_scatterplot.html",
            ]},
            {"legend": MUT_TYPE_LEGEND, "files": [
                "Parental-Q575R-ABE-correlation_scatterplot.html",
            ]},
        ],
    },
    "Boxplots": {
        "folder": "Boxplots",
        "segments": [
            {"legend": None, "files": [
                "WT_Boxplot_ByGene.html",
                "WT_Boxplot_ByMutation.html",
                "Q575R_Boxplot_ByGene.html",
                "Q575R_Boxplot_ByMutation.html",
                "Stability_Boxplot_ByGene.html",
                "Stability_Boxplot_ByMutation.html",
                "PRC2_Quartile_WT_Q575R_Density.html",
                "MTF2_WT_DMSOD0_Boxplot_ByMutation.html",
                "MTF2_Q575R_DMSOD0_Boxplot_ByMutation.html",
                "MTF2_Q575R_DMSOPRC2i_Boxplot_ByMutation.html",
            ]},
            {"legend": BOXPLOT_GENE_LEGEND, "files": [
                "PRC2_Quartile_WT_Q575R_Boxplot.html",
                "WT_Resistance_Boxplot.html",
            ]},
        ],
    },
    "Lollipops": {
        "folder": "Lollipops",
        "segments": [
            {"legend": LOF_ABE_LEGEND, "files": lol("LOF_Stability", "ABE")},
            {"legend": LOF_CBE_LEGEND, "files": lol("LOF_Stability", "CBE")},
            {"legend": Q575R_ABE_LOL_LEGEND, "files": lol("Q575R", "ABE")},
            {"legend": Q575R_CBE_LOL_LEGEND, "files": lol("Q575R", "CBE")},
        ],
    },
    "ClusteredHeatmaps": {
        "folder": "ClusteredHeatmaps",
        "segments": [
            {"legend": None, "files": [
                "ABE_ClusteredHeatmap_K6_Raw_Cutoff3.html",
                "CBE_ClusteredHeatmap_K6_Raw_Cutoff3.html",
            ]},
        ],
    },
    "PWES": {
        "folder": "PWES",
        "segments": [
            {"legend": None, "files": [
                "PWES-6WKR-LOF_Stability-PWES_heatmap.html",
                "PWES-6WKR-LOF_Stability-cluster_boxplots.html",
                "PWES-6WKR-LOF_Stability-cluster_clustermap.html",
                "PWES-6WKR-LOF_Stability-cluster_histogram.html",
                "PWES-6WKR-LOF_Stability-EZH2-cluster_scatterplots.html",
                "PWES-6WKR-LOF_Stability-EED-cluster_scatterplots.html",
                "PWES-6WKR-LOF_Stability-SUZ12-cluster_scatterplots.html",
            ]},
        ],
    },
}


_WRAPPER_DIV_RE = re.compile(r'<div style="height:[0-9.]+px; width:[0-9.]+px;">')


# The three PWES per-gene scatterplots (EED/EZH2/SUZ12-cluster_scatterplots)
# color points by cluster. EED and SUZ12 already export every cluster as the
# same neutral gray (#cccccc); only EZH2 has four clusters (31-34) exported
# with distinct colors, which reads as a meaningful highlight it isn't. Flatten
# every trace's marker color to that same neutral gray for all three files so
# none of them singles out clusters by color.
_PWES_SCATTER_RE = re.compile(r"PWES-6WKR-LOF_Stability-(EED|EZH2|SUZ12)-cluster_scatterplots\.html$")
_NEUTRAL_MARKER_COLOR = "#cccccc"


def _flatten_marker_colors(data):
    for trace in data:
        marker = trace.get("marker")
        if isinstance(marker, dict) and "color" in marker:
            marker["color"] = _NEUTRAL_MARKER_COLOR


# The Scatterplots and Lollipops exports draw the gray/colored protein-domain
# background bands as shapes with yref set to a bare data axis ("y", "y2", ...)
# and y0=0, y1=1 -- i.e. the rect only covers data-space y=0 to y=1, not the
# axis's actual range (e.g. [-10, 7.5] or [0, 12]). Depending on the axis, that
# leaves the band as a thin sliver instead of a full-height background fill.
# (The PWES per-gene scatterplots already export these correctly with
# yref="paper", so they're untouched by this.) Fix: for exactly that pattern,
# switch yref to the axis's "domain" variant ("y domain", "y2 domain", ...),
# which Plotly resolves to 0-1 of that axis's own plotting area regardless of
# its data range -- so the band always fills the full height of its subplot.
_BARE_Y_AXIS_RE = re.compile(r"^y(\d*)$")


def _fix_domain_rect_shapes(layout):
    for shape in layout.get("shapes", []):
        if shape.get("type") != "rect":
            continue
        yref = shape.get("yref", "")
        if shape.get("y0") == 0 and shape.get("y1") == 1 and _BARE_Y_AXIS_RE.match(yref):
            shape["yref"] = yref + " domain"


def process_figure(src_path, dst_path, hide_legend, flatten_colors=False):
    """Write src_path to dst_path with fixes applied to every figure
    (not just legend-suppressed ones):

    1. Every exported figure hardcodes layout.width/height to its native
       export pixel size (e.g. 325x212) and wraps the plot in an outer div
       with that same fixed pixel size baked into its inline style. The
       page's card/iframe is much bigger than that, so the plot only ever
       drew at its tiny native size in the corner of a mostly-empty card --
       config.responsive:true does NOT override an explicitly-set
       layout.width/height on its own. Fix: drop width/height from layout,
       set autosize=True, and rewrite the wrapper div's inline style to
       100%/100% so Plotly measures the actual (much bigger) card and
       draws to fill it from the first paint -- no runtime resize-dispatch
       hack needed.
    2. hide_legend=True also forces layout.showlegend=False, for files
       that belong to a shared-legend group (see CATEGORIES above).
    3. flatten_colors=True forces every trace's marker.color to a single
       neutral gray -- see _flatten_marker_colors above.
    4. Domain-background rect shapes with yref pinned to bare data-axis
       y=0..1 are switched to the "<axis> domain" yref so they fill the
       whole subplot height -- see _fix_domain_rect_shapes above.
    """
    html, div_id, data, layout, config, (start, end) = pj.load(src_path)
    layout.pop("width", None)
    layout.pop("height", None)
    layout["autosize"] = True
    if hide_legend:
        layout["showlegend"] = False
    if flatten_colors:
        _flatten_marker_colors(data)
    _fix_domain_rect_shapes(layout)
    new_call = "Plotly.newPlot(\"%s\", %s, %s, %s)" % (
        div_id, json.dumps(data), json.dumps(layout), json.dumps(config),
    )
    new_html = html[:start] + new_call + html[end:]
    new_html = _WRAPPER_DIV_RE.sub('<div style="height:100%; width:100%;">', new_html, count=1)
    os.makedirs(os.path.dirname(dst_path), exist_ok=True)
    with open(dst_path, "w", encoding="utf-8") as f:
        f.write(new_html)


def main():
    # folders actually touched (so we can clear+rebuild data/ cleanly for them)
    touched_folders = {cfg["folder"] for cfg in CATEGORIES.values()}
    for folder in touched_folders:
        dst_dir = os.path.join(DATA, folder)
        if os.path.isdir(dst_dir):
            shutil.rmtree(dst_dir)
        os.makedirs(dst_dir, exist_ok=True)

    manifest = {}
    n_suppressed = 0
    n_copied = 0
    for category, cfg in CATEGORIES.items():
        folder = cfg["folder"]
        seg_out = []
        seen_files_in_cat = set()
        for seg in cfg["segments"]:
            for fname in seg["files"]:
                src_path = os.path.join(SRC, folder, fname)
                dst_path = os.path.join(DATA, folder, fname)
                if not os.path.isfile(src_path):
                    raise FileNotFoundError(src_path)
                if fname not in seen_files_in_cat:
                    process_figure(
                        src_path, dst_path,
                        hide_legend=seg["legend"] is not None,
                        flatten_colors=bool(_PWES_SCATTER_RE.search(fname)),
                    )
                    if seg["legend"] is not None:
                        n_suppressed += 1
                    else:
                        n_copied += 1
                    seen_files_in_cat.add(fname)
            seg_out.append({"legend": seg["legend"], "files": seg["files"]})
        manifest[category] = seg_out

    with open(os.path.join(DOCS, "manifest.json"), "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)

    print(f"legend-suppressed + fill-fixed {n_suppressed} files, fill-fixed only {n_copied} files")
    total = sum(len(seg["files"]) for cfg in CATEGORIES.values() for seg in cfg["segments"])
    print(f"total file-slots across all segments: {total}")

    # sanity: every file that exists in new_outputs_interactive for these
    # folders should appear at least once in the manifest, and vice versa.
    for folder in touched_folders:
        on_disk = set(os.listdir(os.path.join(SRC, folder)))
        in_manifest = set()
        for category, cfg in CATEGORIES.items():
            if cfg["folder"] != folder:
                continue
            for seg in cfg["segments"]:
                in_manifest.update(seg["files"])
        missing = on_disk - in_manifest
        extra = in_manifest - on_disk
        if missing:
            print(f"WARNING: {folder} has files on disk not in manifest: {missing}")
        if extra:
            print(f"WARNING: {folder} manifest references files not on disk: {extra}")

    return manifest


if __name__ == "__main__":
    main()
