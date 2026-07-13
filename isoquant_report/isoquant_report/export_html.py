"""Static HTML export for headless pipeline runs.

Builds a comprehensive, navigable single-file report grouping figures and
summary tables into interpreted sections.
"""

from __future__ import annotations

import datetime as _dt
import json
from pathlib import Path
from typing import Any

import pandas as pd
from jinja2 import Template

from isoquant_report.report import (
    ReportData,
    build_all_figures,
    build_summary_tables,
    headline_metrics,
)

# --- Section layout ------------------------------------------------------
# Each section groups related figures (by key from build_all_figures) and
# tables (by title from build_summary_tables) with an interpretive blurb.
SECTION_LAYOUT: list[dict[str, Any]] = [
    {
        "id": "assignment",
        "title": "Read assignment QC",
        "desc": (
            "How IsoQuant assigned long reads to reference transcripts. A high "
            "proportion of <em>unique</em> assignments and a low proportion of "
            "<em>ambiguous</em>/<em>inconsistent</em> reads indicates confident "
            "isoform-level quantification. Structural categories follow the "
            "SQANTI-style scheme (FSM, ISM, NIC, NNIC, etc.)."
        ),
        "tables": ["Assignment types", "Structural categories"],
        "figs": [
            "assignment_type",
            "structural_category",
            "novel_known",
            "assignment_per_cell",
            "assignment_rate",
        ],
    },
    {
        "id": "cellqc",
        "title": "Cell calling & QC",
        "desc": (
            "Per-cell quality metrics derived from the gene count matrix. The "
            "barcode-rank (knee) plot separates real cells from empty droplets. "
            "Violin plots show the distribution of total UMIs, genes detected, "
            "and mitochondrial content; cells outside the configured thresholds "
            "are flagged in the summary table."
        ),
        "tables": ["Cell QC summary"],
        "figs": [
            "knee_plot",
            "total_counts_violin",
            "genes_detected_violin",
            "mito_violin",
            "depth_vs_genes",
        ],
    },
    {
        "id": "splice",
        "title": "Splice / unspliced QC",
        "desc": (
            "Molecule-level splicing status projected to each cell. The unspliced "
            "fraction reflects intron retention / nascent transcripts; unusually "
            "high values can indicate genomic-DNA contamination or immature "
            "nuclei. Spliced vs unspliced UMIs should correlate across cells."
        ),
        "tables": ["Splice summary"],
        "figs": [
            "splice_composition",
            "unspliced_fraction_hist",
            "spliced_vs_unspliced",
        ],
    },
    {
        "id": "genevsiso",
        "title": "Gene vs isoform detection",
        "desc": (
            "Comparison of gene-level and isoform-level detection per cell, on the "
            "intersection of barcodes present in both matrices. Cells detecting "
            "many genes but few isoforms may have insufficient read depth for "
            "confident isoform assignment."
        ),
        "tables": [],
        "figs": ["gene_vs_isoform_detection", "top_genes", "top_isoforms"],
    },
    {
        "id": "structure",
        "title": "Isoform structure",
        "desc": (
            "Structural properties of the transcript models IsoQuant reported, "
            "including how many are novel (not present in the reference "
            "annotation), the number of isoforms per gene, and the distribution "
            "of transcript lengths and exon counts."
        ),
        "tables": ["Isoform structure", "Transcript models (novel vs known)"],
        "figs": [
            "model_novelty",
            "isoforms_per_gene",
            "tx_length",
            "exon_count",
        ],
    },
    {
        "id": "saturation",
        "title": "Sequencing saturation",
        "desc": (
            "Rarefaction curves obtained by downsampling reads/UMIs. A curve that "
            "is still climbing steeply indicates that deeper sequencing would "
            "recover appreciably more isoforms or genes."
        ),
        "tables": [],
        "figs": ["isoform_saturation", "umi_saturation"],
    },
    {
        "id": "reconcile",
        "title": "Matrix reconciliation",
        "desc": (
            "Agreement between IsoQuant's native grouped count matrix and the "
            "derived spliced/unspliced AnnData matrices built by the pipeline. "
            "High cell/feature Jaccard and total-count correlation confirm the "
            "derived matrices faithfully represent the native output."
        ),
        "tables": ["Matrix reconciliation"],
        "figs": [],
    },
]

CARD_META = {
    "n_cells": ("Cells", "Called barcodes in the gene matrix"),
    "n_genes": ("Genes", "Features in the gene matrix"),
    "n_isoforms": ("Isoforms", "Features in the transcript matrix"),
    "median_umis": ("Median UMIs", "Median total UMIs per cell"),
    "median_unspliced_fraction": ("Median unspliced frac", "Per-cell median"),
    "pct_unique_assignments": ("% unique", "Reads uniquely assigned"),
    "cell_jaccard": ("Cell Jaccard", "Native vs derived cell overlap"),
    "count_correlation": ("Count corr.", "Native vs derived totals"),
}

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>IsoQuant QC Report — {{ sample }}</title>
  <script src="https://cdn.plot.ly/plotly-3.0.1.min.js" charset="utf-8"></script>
  <style>
    :root { --fg:#1f2933; --muted:#647079; --line:#e3e8ee; --accent:#2b6cb0;
            --bg:#ffffff; --panel:#f8fafc; --warn:#b7791f; --ok:#2f855a; }
    * { box-sizing: border-box; }
    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", system-ui, sans-serif;
           margin: 0; color: var(--fg); background: var(--bg); line-height: 1.5; }
    a { color: var(--accent); text-decoration: none; }
    .layout { display: flex; align-items: flex-start; }
    nav { position: sticky; top: 0; align-self: flex-start; width: 240px; height: 100vh;
          overflow-y: auto; border-right: 1px solid var(--line); padding: 1.5rem 1rem;
          background: var(--panel); flex-shrink: 0; }
    nav h2 { font-size: 0.8rem; text-transform: uppercase; letter-spacing: .05em;
             color: var(--muted); margin: 0 0 .5rem; }
    nav ol { list-style: none; padding: 0; margin: 0; }
    nav li { margin: .15rem 0; }
    nav a { display: block; padding: .3rem .5rem; border-radius: 6px; font-size: .9rem; }
    nav a:hover { background: #eef2f7; }
    main { flex: 1; padding: 2rem 2.5rem; max-width: 1100px; margin: 0 auto; }
    header h1 { margin: 0 0 .25rem; font-size: 1.6rem; }
    header .sub { color: var(--muted); font-size: .9rem; margin-bottom: 1.5rem; }
    .cards { display: grid; grid-template-columns: repeat(auto-fill, minmax(150px, 1fr));
             gap: .75rem; margin: 1rem 0 2rem; }
    .card { border: 1px solid var(--line); border-radius: 10px; padding: .85rem 1rem;
            background: var(--panel); }
    .card h3 { margin: 0; font-size: 1.4rem; }
    .card p { margin: .2rem 0 0; color: var(--muted); font-size: .78rem; }
    section { margin: 2.5rem 0; scroll-margin-top: 1rem; }
    section > h2 { font-size: 1.25rem; border-bottom: 2px solid var(--line);
                   padding-bottom: .4rem; }
    .desc { color: var(--muted); font-size: .92rem; margin: .5rem 0 1.25rem; }
    .grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(440px, 1fr));
            gap: 1rem; }
    .plot { border: 1px solid var(--line); border-radius: 10px; padding: .5rem; min-height: 320px; }
    table { border-collapse: collapse; width: 100%; margin: .5rem 0 1.25rem; font-size: .9rem; }
    caption { text-align: left; font-weight: 600; padding: .35rem 0; }
    th, td { border-bottom: 1px solid var(--line); padding: .45rem .6rem; text-align: left; }
    th { background: var(--panel); font-weight: 600; }
    tbody tr:hover { background: #fafcff; }
    .warn { color: var(--warn); }
    .pill { display:inline-block; padding:.1rem .5rem; border-radius: 999px; font-size:.75rem; }
    .pill.ok { background:#e6fffa; color: var(--ok); }
    .pill.no { background:#fff5f5; color:#c53030; }
    .empty { color: var(--muted); font-style: italic; }
    footer { color: var(--muted); font-size: .8rem; margin: 3rem 0 1rem; }
  </style>
</head>
<body>
<div class="layout">
  <nav>
    <h2>Contents</h2>
    <ol>
      <li><a href="#overview">Overview</a></li>
      <li><a href="#completeness">Data completeness</a></li>
      {% for s in sections %}<li><a href="#{{ s.id }}">{{ s.title }}</a></li>{% endfor %}
    </ol>
  </nav>
  <main>
    <header>
      <h1>IsoQuant single-cell QC report</h1>
      <div class="sub">Sample <strong>{{ sample }}</strong> &middot; generated {{ generated }}</div>
    </header>

    <section id="overview">
      <h2>Overview</h2>
      <div class="desc">Headline metrics summarising the run at a glance.</div>
      <div class="cards">
        {% for c in cards %}
        <div class="card"><h3>{{ c.value }}</h3><p title="{{ c.help }}">{{ c.label }}</p></div>
        {% endfor %}
      </div>
    </section>

    <section id="completeness">
      <h2>Data completeness</h2>
      <div class="desc">Which IsoQuant / pipeline artefacts were located and loaded.</div>
      <table>
        <thead><tr><th>Artefact</th><th>Status</th></tr></thead>
        <tbody>
        {% for name, ok in completeness.items() %}
          <tr><td>{{ name }}</td><td>
            {% if ok %}<span class="pill ok">loaded</span>{% else %}<span class="pill no">missing</span>{% endif %}
          </td></tr>
        {% endfor %}
        </tbody>
      </table>
      {% if errors %}<p class="warn">Warnings: {{ errors | join("; ") }}</p>{% endif %}
    </section>

    {% for s in sections %}
    <section id="{{ s.id }}">
      <h2>{{ s.title }}</h2>
      <div class="desc">{{ s.desc | safe }}</div>
      {% if not s.tables and not s.figs %}
        <p class="empty">No data available for this section.</p>
      {% endif %}
      {% for t in s.tables %}
      <table>
        <caption>{{ t.title }}</caption>
        <thead><tr>{% for col in t.columns %}<th>{{ col }}</th>{% endfor %}</tr></thead>
        <tbody>
          {% for row in t.rows %}<tr>{% for cell in row %}<td>{{ cell }}</td>{% endfor %}</tr>{% endfor %}
        </tbody>
      </table>
      {% endfor %}
      {% if s.figs %}
      <div class="grid">
        {% for f in s.figs %}<div class="plot" id="{{ f }}"></div>{% endfor %}
      </div>
      {% endif %}
    </section>
    {% endfor %}

    <footer>Generated by isoquant_report. Interactive figures rendered with Plotly.</footer>
  </main>
</div>
<script>
  var figures = {{ figures_json | safe }};
  for (var id in figures) {
    try {
      Plotly.newPlot(id, figures[id].data, figures[id].layout, {responsive: true});
    } catch (e) {
      var el = document.getElementById(id);
      if (el) { el.innerHTML = '<p class="warn">Failed to render figure: ' + e + '</p>'; }
      console.error("Failed to render", id, e);
    }
  }
</script>
</body>
</html>
"""


def _table_payload(title: str, df: pd.DataFrame) -> dict[str, Any]:
    """Serialise a DataFrame into template-friendly rows/columns."""
    return {
        "title": title,
        "columns": [str(c) for c in df.columns],
        "rows": df.astype(object).where(pd.notnull(df), "").values.tolist(),
    }


def export_html(data: ReportData, outpath: str | Path) -> None:
    """Write a self-contained, sectioned HTML report."""
    figs = build_all_figures(data)
    tables = build_summary_tables(data)
    metrics = headline_metrics(data)

    cards = [
        {
            "label": CARD_META.get(k, (k, ""))[0],
            "help": CARD_META.get(k, (k, ""))[1],
            "value": v,
        }
        for k, v in metrics.items()
        if k != "sample"
    ]

    figures_json: dict[str, Any] = {}
    for key, fig in figs.items():
        figures_json[key] = json.loads(fig.to_json())

    # Build the sections that actually have content.
    rendered_sections = []
    for spec in SECTION_LAYOUT:
        sec_tables = [
            _table_payload(t, tables[t])
            for t in spec["tables"]
            if t in tables and not tables[t].empty
        ]
        sec_figs = [f for f in spec["figs"] if f in figures_json]
        rendered_sections.append(
            {
                "id": spec["id"],
                "title": spec["title"],
                "desc": spec["desc"],
                "tables": sec_tables,
                "figs": sec_figs,
            }
        )

    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    html = Template(HTML_TEMPLATE).render(
        sample=data.config.sample,
        generated=_dt.datetime.now().strftime("%Y-%m-%d %H:%M"),
        cards=cards,
        completeness=data.paths.completeness(),
        errors=data.errors,
        sections=rendered_sections,
        figures_json=json.dumps(figures_json),
    )
    outpath.write_text(html)
