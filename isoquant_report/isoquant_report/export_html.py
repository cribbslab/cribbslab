"""Static HTML export for headless pipeline runs."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from jinja2 import Template

from isoquant_report.report import ReportData, build_all_figures, headline_metrics

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>IsoQuant QC Report — {{ sample }}</title>
  <script src="https://cdn.plot.ly/plotly-3.0.1.min.js" charset="utf-8"></script>
  <style>
    body { font-family: system-ui, sans-serif; margin: 2rem; max-width: 1200px; }
    .cards { display: flex; flex-wrap: wrap; gap: 1rem; margin: 1rem 0; }
    .card { border: 1px solid #ddd; border-radius: 8px; padding: 1rem; min-width: 140px; }
    .card h3 { margin: 0; font-size: 1.5rem; }
    .card p { margin: 0.25rem 0 0; color: #555; font-size: 0.85rem; }
    .section { margin: 2rem 0; }
    .plot { margin: 1rem 0; }
    table { border-collapse: collapse; width: 100%; }
    th, td { border: 1px solid #ddd; padding: 0.5rem; text-align: left; }
    th { background: #f5f5f5; }
    .warn { color: #a00; }
  </style>
</head>
<body>
  <h1>IsoQuant QC Report</h1>
  <p>Sample: <strong>{{ sample }}</strong></p>

  <div class="cards">
    {% for label, value in metrics.items() %}
    <div class="card"><h3>{{ value }}</h3><p>{{ label }}</p></div>
    {% endfor %}
  </div>

  <div class="section">
    <h2>Data completeness</h2>
    <table>
      <tr><th>Artifact</th><th>Loaded</th></tr>
      {% for name, ok in completeness.items() %}
      <tr><td>{{ name }}</td><td>{{ "yes" if ok else "no" }}</td></tr>
      {% endfor %}
    </table>
    {% if errors %}
    <p class="warn">Warnings: {{ errors | join("; ") }}</p>
    {% endif %}
  </div>

  {% for section, fig_id in sections %}
  <div class="section">
    <h2>{{ section }}</h2>
    <div class="plot" id="{{ fig_id }}"></div>
  </div>
  {% endfor %}

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


def export_html(data: ReportData, outpath: str | Path) -> None:
    """Write a self-contained HTML report with embedded Plotly figures."""
    figs = build_all_figures(data)
    metrics = headline_metrics(data)
    # Friendly labels for cards
    card_labels = {
        "n_cells": "Cells",
        "n_genes": "Genes",
        "n_isoforms": "Isoforms",
        "median_umis": "Median UMIs",
        "median_unspliced_fraction": "Median unspliced frac",
        "pct_unique_assignments": "% unique assignments",
        "cell_jaccard": "Cell Jaccard",
        "count_correlation": "Count correlation",
    }
    cards = {
        card_labels.get(k, k): v
        for k, v in metrics.items()
        if k != "sample"
    }

    import json

    figures_json: dict[str, Any] = {}
    sections = []
    for i, (name, fig) in enumerate(figs.items()):
        fid = f"fig_{i}"
        figures_json[fid] = json.loads(fig.to_json())
        sections.append((name.replace("_", " ").title(), fid))

    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    tmpl = Template(HTML_TEMPLATE)
    html = tmpl.render(
        sample=data.config.sample,
        metrics=cards,
        completeness=data.paths.completeness(),
        errors=data.errors,
        sections=sections,
        figures_json=json.dumps(figures_json),
    )
    outpath.write_text(html)
