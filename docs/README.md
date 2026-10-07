# Interactive figures — browsable data page

Static site: `index.html` (single page, light/dark toggle, a left sidebar
listing each figure type; clicking one lays out every figure in that
category side by side in a responsive grid — no per-figure scrolling
required, just more rows as the window narrows) plus `data/<Category>/`
(one HTML file per interactive figure, copied straight from
`interactive_analysis.ipynb`'s `new_outputs_interactive/`).

The sidebar has 8 entries. `Boxplots`, `Lollipops`, `ClusteredHeatmaps` and
`PWES` map 1:1 to their `new_outputs_interactive/` subfolder. `Scatterplots`
is split into 4 sidebar categories — `WT Scatterplots`, `Q575R Scatterplots`,
`Stability Scatterplots` and `Correlation Scatterplots` — that all point at
the single `data/Scatterplots/` folder (see `FOLDER_OVERRIDE` in
`index.html`); the 10 combined per-gene stability scatterplots (each already
shows both the GFP-EZH2 and GFP-EED constructs together) live under the one
`Stability Scatterplots` entry (was split into separate `EZH2 Stability
Scatterplots` / `EED Stability Scatterplots` entries pointing at the same
files -- merged per Calvin's 2026-10-07 request, since the data was always
identical either way).

Per a 2026-09-28 request: `Lollipops` now only contains the screen-derived
lollipops (`Q575R_*`, `LOF_Stability_*`); the FoldX/ThermoMPNN
predicted-ddG stability-prediction lollipops were removed entirely (the
84th-percentile FoldX cutoff was passing 17-21% of guides per gene — too
dense on a small plot — and rather than re-tune it, the predictor-based
lollipops are dropped, keeping only screen data). Lollipops are ordered
EZH2, EED, SUZ12, AEBP2, JARID2 (was alphabetical). `ClusteredHeatmaps` now
only keeps the base `ABE_ClusteredHeatmap_K6` / `CBE_ClusteredHeatmap_K6`
(the cutoff-3 and per-cluster subset heatmaps were dropped). These same
changes were made in `master_analysis.ipynb` (the static figures).

Per a 2026-09-29 follow-up request: `ClusteredHeatmaps` now keeps the
cutoff-3, ARI-validated guide set instead — `ABE_ClusteredHeatmap_K6_Raw_Cutoff3`
/ `CBE_ClusteredHeatmap_K6_Raw_Cutoff3` (still 1 file per editor).
`Boxplots` dropped the MTF2-screen "PRC2i-D0" mutation-type plot and
renamed the 3 remaining MTF2-derived boxplots with an `MTF2_` prefix.
`Q575R Scatterplots` was trimmed per gene to just the large "DMSO-PRC2i"
and small "PRC2i-D0" plots (MTF2's 3 scatterplots, which have no
large/small size variants, are unchanged). `Correlation Scatterplots` now
keeps only `Parental-Q575R-ABE-correlation_scatterplot`, the 2
GFP-EZH2/GFP-EED plots, and the 2 DMSO-PRC2i/PRC2i-D0 plots (the
DMSO-PRC2i/DMSO-D0 pair and the MTF2-screen
`Q575R-DMSO-PRC2i-DMSO-D0-ABE` plot were dropped). See `errors.md`'s
2026-09-29 entry for the full breakdown. These same changes were made in
`master_analysis.ipynb` (the static figures).

Per a same-day (2026-09-29) follow-up: `Correlation Scatterplots` also
dropped the `Scatter_GFP-MAX_{ABE,CBE,merged}_{FoldX,ThermoMPNN}` predictor
plots (5 files left). `Lollipops` was trimmed to the split-axis version only
for both the Q575R and LOF-vs-Stability sets (20 files, was 40) — while
checking this, a `be_scan` library bug was found and fixed: the interactive
lollipops' negative/bottom-panel hit filter (`abs(value) >= threshold` with
a negative threshold) never actually filtered anything, so the bottom panel
of every interactive lollipop was showing every data point instead of just
the real hits; it's now monkey-patched in `interactive_analysis.ipynb`'s
imports cell to match the static notebook's filtering. `Boxplots`,
`ClusteredHeatmaps`, `PWES` and `WT/Q575R/Stability Scatterplots` are
unchanged from the entry above. See `errors.md` for the full breakdown.

### Shared legends (2026-09-29 follow-up)

Each sidebar page now shows one static (non-clickable) legend bar at the
top, and the plots underneath it have their own inline Plotly legend
hidden, instead of every plot carrying its own legend. This is implemented
in `build_webpage.py` (new — see "Regenerating" below), not in
`interactive_analysis.ipynb`: `be_scan.figure_plot.figure_interactive` has
no clean parameter to redirect a plot's legend the way the static
pipeline's `legend=replace(LEGEND, ..., path=...)` writes one to a
separate file, so the build script instead parses each candidate file's
embedded `Plotly.newPlot(...)` JSON, groups files that share an identical
color scheme, and writes a copy to `data/` with `layout.showlegend`
forced to `false` for files in a shared group (files with a unique or
nonexistent inline legend are copied through unchanged).

Groups found: all 29 WT/Q575R/Stability-gene scatterplots plus the 1
Parental correlation-scatterplot share one 3-color mutation-type palette;
4 of the 5 Correlation Scatterplots share a 5-color gene palette (the 5th,
`Parental-Q575R-ABE`, gets its own legend bar since it's on the
mutation-type scheme instead); exactly 2 Boxplots
(`PRC2_Quartile_WT_Q575R_Boxplot`, `WT_Resistance_Boxplot`) share a
7-color gene palette — the other 10 Boxplots have no inline legend at all
(mutation/gene categories are already printed on the x-axis) so nothing
changed for them; Lollipops split into 4 groups of 5 genes each
(`LOF_Stability_ABE/CBE`, `Q575R_ABE/CBE`), now reordered into 4
contiguous gene-ordered blocks (EZH2/EED/SUZ12/AEBP2/JARID2) instead of
interleaved, so each legend sits directly above the 5 plots it covers.
`ClusteredHeatmaps` and `PWES` are unchanged.

`manifest.json`'s shape changed from `{category: [files]}` to `{category:
[{legend: [...]|null, files: [...]}]}` — a "segment"; `index.html`'s
embedded `MANIFEST_SEGMENTS` mirrors it, and `renderCategory()` renders a
legend bar above each segment's own grid instead of one flat grid per
category. See `errors.md`'s 2026-09-29 follow-up entry for the full
breakdown of how the groups were determined.

### Figures filling their card, the real fix (2026-09-29 follow-up)

The 2026-09-28 fix below (stretch the wrapper div, dispatch a `resize`
event) turned out not to be reliable enough in practice -- figures were
still drawing at their tiny native export size inside much bigger cards,
because each figure's `layout` hardcodes its native pixel `width`/`height`
and `config.responsive: true` doesn't override that on its own.
`build_webpage.py` now fixes this properly for every figure in every
category (not just the legend-touched ones): it drops `layout.width`/
`height`, sets `layout.autosize = true`, and rewrites the wrapper div's
inline style to `100%`/`100%`, so the figure measures its real container
and fills it from the first paint. `index.html`'s `attachResponsiveFill`
is now just a harmless margin/overflow cleanup on iframe load -- the
resize-dispatch trick is gone.

While tracking this down, a real bug was also found and fixed in
`plotly_json.py`: its bracket-scanner located the end of the
`Plotly.newPlot(id, data, layout, config` arguments but not the call's own
closing `)`, so every legend-suppressed file got a stray extra `)` left
over from the original text -- a silent JS syntax error that meant
`Plotly.newPlot` never ran at all in any of those 66 files (they rendered
as pure blank space, not just badly-sized). See `errors.md`'s same-day
follow-up entry for the full story, including how this was actually
verified by rendering in headless Chromium with the exact Plotly build
the figures were exported with, rather than trusting a JSON re-parse.

### Figures filling their card (2026-09-28 follow-up fix, superseded above)

Each exported Plotly figure has a small **fixed-pixel** outer wrapper div
(e.g. ~325x212px for a scatterplot, ~423x605px for a ClusteredHeatmap)
around the actual `plotly-graph-div`. Just enlarging the iframe/card around
it did nothing, because that outer div never grows past its native
export size — this was the "mostly wasted space" bug. The fix
(`attachResponsiveFill` in `index.html`) forces that outer wrapper's
`width`/`height` to `100%` once the iframe loads, then fires a `resize`
event on the iframe's `contentWindow` (immediately, then again at 60ms and
250ms) so Plotly's own `{responsive: true}` config re-measures the now
full-size container and relayouts the plot to fill it.

Per-category sizing (`LAYOUT` in `index.html`, replacing one fixed card
size for everything):

- **ClusteredHeatmaps**: tall cards (`calc(100vh - 200px)`, min 560px) in
  an `auto-fit minmax(380px, 1fr)` grid, so both heatmaps fill the page
  height with no scrolling.
- **PWES**: a strict 2-column grid. The heatmap, boxplots, histogram and
  clustermap files each span the full row at `calc(100vh - 200px)` height
  (one long window per plot, per request); the 3 gene `cluster_scatterplots`
  files (EZH2, EED, SUZ12) stay at normal height, 2 across, in the same
  grid.
- Everything else (Scatterplots, Boxplots, Lollipops): `auto-fit
  minmax(460px, 1fr)` at 480px card height.

Also fixed this round: the gene-order sort for the split scatterplot
categories (`WT Scatterplots`, `Q575R Scatterplots`, `Stability
Scatterplots`) had silently fallen back to alphabetical order, because the
stability-scatterplot filenames all contain the substring `GFP_EZH2_EED_K562`
(the construct name), so a naive "does the filename contain EZH2"
substring check matched every file instead of the actual guide gene. The
sort now keys off the filename's **leading prefix** (the gene name before
the first `-`) instead of a substring search, and all scatterplot
categories are correctly ordered EZH2, EED, SUZ12, AEBP2, JARID2 (MTF2
last, since it's outside the requested gene list).

`index.html` has the figure list (which files exist per category) baked
into it directly, as an inline `<script>` block, rather than fetched from a
separate `manifest.json` at load time — a `fetch()` of a sibling file is
blocked by the browser when the page is opened as a plain `file://` page
(e.g. double-clicked from Finder or a Drive folder), which is how this page
will be opened most of the time. `manifest.json` still exists alongside it
as the source the build step reads from (see below), but the page itself
doesn't need it at runtime.

## Viewing it

Because the figure list is embedded, **you can just double-click
`index.html`** — no server, no GitHub Pages, no build step required. That
also means it can be opened straight from this Drive folder.

## Publishing via GitHub Pages (optional)

If you want a shareable URL instead of (or in addition to) opening the file
locally: this folder is already named `docs/` for exactly this purpose —
in the repo's Settings → Pages, set the source to the `main` branch,
`/docs` folder. Push, and everything here is client-side, same as viewing
it locally; no further build step runs on GitHub's end.

## Regenerating after re-running `interactive_analysis.ipynb`

Re-run the notebook, then from `docs/` run `python3 build_webpage.py`.
It reads `new_outputs_interactive/`, rebuilds `data/` (suppressing the
inline legend on files that belong to a shared-legend group, copying
everything else unchanged), and rewrites `manifest.json`. It prints a
warning if any on-disk figure isn't referenced by any category's file
list, or vice versa, so a change in `interactive_analysis.ipynb`'s output
(a renamed/added/removed file) doesn't silently fall out of the page.

`build_webpage.py` does **not** touch `index.html`'s embedded
`MANIFEST_SEGMENTS` — that's a separate step, since the segment groupings
(which files share a legend) were determined once by inspecting the
figures' actual Plotly JSON and are unlikely to change just because the
notebook re-ran with the same figures. If the notebook's output changes in
a way that could change a color scheme or grouping (a new gene added, a
palette recolored), re-derive the groups (`plotly_json.py` has the
JSON-parsing helper used to do this) before re-embedding
`manifest.json`'s contents into `index.html`.

## Note on file sizes

Most figures are a few KB to a few hundred KB (plotly.js itself is loaded
from a CDN, not embedded — see `interactive_analysis.ipynb`'s imports
cell). Two PWES figures are much larger because they're full residue x
residue matrices: `PWES-6WKR-LOF_Stability-PWES_heatmap.html` (~45 MB) and
`PWES-6WKR-LOF_Stability-cluster_clustermap.html` (~13 MB). Both are under
GitHub's 100 MB hard file-size limit, but GitHub will warn about the larger
one when it's pushed — that's expected, not an error.
