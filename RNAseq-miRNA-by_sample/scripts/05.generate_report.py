#!/usr/bin/env python3
"""
05 - Functional enrichment HTML report generator (GO/KEGG + miRNA)

Description:
    Iterates over the summarized "Top" tables produced by the R script
    (04.continueR.R) for each comparison (results/RNAm&miRNA_top/), and
    builds one HTML report per comparison with two sections (GO and KEGG):
    a summary dotplot and a table with the most relevant terms/pathways,
    including the target genes and miRNAs involved, with their score and
    log2FoldChange values shown as tooltips. It also generates an
    index.html that links to every generated report.
    
Usage:
    python3 scripts/05.generate_report.py

Requirements:
    - Must be run after 04.continueR.sh, since it consumes the CSV/PNG
      files produced by that step (results/RNAm&miRNA_top,
      results/Activated_Suppressed_Graphics, etc.).
    - pandas is required; Pillow (PIL) is optional and only used to
      detect a blank pathview.png more accurately.
"""

import os
import re
import glob
import pandas as pd

try:
    from PIL import Image
    HAS_PIL = True
except ImportError:
    HAS_PIL = False

# Fixed GSEA parameter (adjust it if your pipeline uses a different maxSize)
MAX_SIZE = 100

# If a local pathview.png is smaller than this, it is suspected to be blank/broken
BLANK_SIZE_THRESHOLD_BYTES = 15000

script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(script_dir)
base = os.path.join(project_root, "results")
data_dir = os.path.join(project_root, "data")

top_csv_dir = os.path.join(base, "RNAm&miRNA_top_padj0.05")
dotplot_dir = os.path.join(base, "Activated_Suppressed_Graphics_padj0.05")
graphics_root = os.path.join(data_dir, "mrna_input", "GO&KEGG_Graphics")

reports = os.path.join(base, "reports")
os.makedirs(reports, exist_ok=True)

# Stylesheet used across all generated HTML reports
css = """
body { font-family: "Times New Roman", Times, serif; margin: 40px; color: black; line-height: 1.4; }
h1 { text-align: center; background-color: #fff3b0; padding: 15px; border-radius: 10px; }
.subtitle { text-align: center; color: #b00000; font-weight: bold; font-size: 18px; margin-bottom: 25px; }
#index { background-color: #f2f2f2; padding: 15px; border-radius: 8px; margin-bottom: 30px; }
#index a { color: black; text-decoration: none; }
#index a:hover { text-decoration: underline; }
.section-title { background-color: #ffe0b2; text-align: center; padding: 10px; border-radius: 8px; margin-top: 40px; }
.main-image { width: 650px; max-width: 90%; display: block; margin: auto; border: 1px solid #cccccc; }
.mini-plot { width: 70px; border: 1px solid #cccccc; margin: 0 auto; display: block; }
.plot-cell { text-align: center; }
.click-text { font-size: 12px; text-align: center; color: #555; }
.table-title { text-align: center; font-size: 18px; font-weight: bold; margin-top: 25px; margin-bottom: 10px; }
table { width: 100%; border-collapse: collapse; margin-bottom: 30px; }
th { background-color: #1e3a8a; color: white; padding: 8px; font-size: 18px; }
td { border: 1px solid #dddddd; padding: 8px; vertical-align: middle; font-size: 16px; }
tr:nth-child(even) { background-color: #e6f0ff; }
tr:hover { background-color: #d4e4ff; }
.note-text { font-size: 15px; font-style: italic; color: #333; margin-top: 5px; margin-bottom: 15px; }
.mirna-tag { cursor: help; border-bottom: 1px dotted blue; }
.gene-tag { cursor: help; border-bottom: 1px dotted #008000; }
.id-cell { font-weight: bold; }
"""

GENE_DETAIL_RE = re.compile(r'^(.+?) \(score=([^)]+)\)$')
MIRNA_DETAIL_RE = re.compile(r'^(.+?) \(L2FC=([^)]+)\)$')


def add_image(abs_path, html_rel_path, out, cls):
    """Writes a clickable image into the HTML output if the file exists."""
    if abs_path and os.path.exists(abs_path):
        out.write(f"""
<a href="{html_rel_path}" target="_blank">
<img src="{html_rel_path}" class="{cls}">
</a>
<p class="click-text">(Click on image to enlarge)</p>
""")
        return True
    return False


def display_name(comparison):
    """Strips the leading 'sample_' prefix from a comparison name, if present."""
    return comparison[len("sample_"):] if comparison.startswith("sample_") else comparison


def is_blank_image(path):
    """Checks whether a PNG image appears to be blank (or broken)."""
    if HAS_PIL:
        try:
            img = Image.open(path).convert("L")
            pixels = img.getdata()
            mean = sum(pixels) / len(pixels)
            return mean > 250
        except Exception:
            pass
    try:
        return os.path.getsize(path) < BLANK_SIZE_THRESHOLD_BYTES
    except OSError:
        return False


def discover_comparisons():
    """Automatically detects the available comparisons from the GO CSV files."""
    pattern = os.path.join(top_csv_dir, "GO_*_RNAm_miRNA_Top.csv")
    comparisons = []
    for path in sorted(glob.glob(pattern)):
        fname = os.path.basename(path)
        comp = fname[len("GO_"):-len("_RNAm_miRNA_Top.csv")]
        comparisons.append(comp)
    return comparisons


def parse_gene_scores(details_str):
    """gene -> score (a single value per gene)."""
    scores = {}
    if not isinstance(details_str, str) or not details_str.strip():
        return scores
    for entry in details_str.split("<br>"):
        entry = entry.strip()
        m = GENE_DETAIL_RE.match(entry)
        if not m:
            continue
        gene, score = m.groups()
        try:
            score = float(score)
        except ValueError:
            score = None
        scores[gene.strip()] = score
    return scores


def parse_mirna_l2fc(details_str):
    """miRNA -> L2FC (a single value per miRNA)."""
    l2fcs = {}
    if not isinstance(details_str, str) or not details_str.strip():
        return l2fcs
    for entry in details_str.split("<br>"):
        entry = entry.strip()
        m = MIRNA_DETAIL_RE.match(entry)
        if not m:
            continue
        mirna, l2fc = m.groups()
        try:
            l2fc = float(l2fc)
        except ValueError:
            l2fc = None
        l2fcs[mirna.strip()] = l2fc
    return l2fcs


def mirna_cell_html(mirnas_str, mirna_details_str):
    """Builds the miRNA table cell, with a per-miRNA L2FC tooltip."""
    if not isinstance(mirnas_str, str) or not mirnas_str.strip():
        return "-"

    l2fcs = parse_mirna_l2fc(mirna_details_str)
    mirna_names = [m.strip() for m in mirnas_str.split(",") if m.strip()]

    spans = []
    for name in mirna_names:
        l2fc = l2fcs.get(name)
        tooltip = f"L2FC: {l2fc:.2f}" if l2fc is not None else "No L2FC data available"
        spans.append(f'<span class="mirna-tag" title="{tooltip}">{name}</span>')

    return "<br>".join(spans)


def target_genes_cell_html(genes_str, gene_details_str):
    """Builds the target-genes table cell, with a per-gene score tooltip."""
    if not isinstance(genes_str, str) or not genes_str.strip():
        return "-"

    scores = parse_gene_scores(gene_details_str)
    gene_names = [g.strip() for g in genes_str.split(",") if g.strip()]

    spans = []
    for name in gene_names:
        score = scores.get(name)
        tooltip = f"Score: {score:.1f}" if score is not None else "No score data available"
        spans.append(f'<span class="gene-tag" title="{tooltip}">{name}</span>')

    return "<br>".join(spans)


def enrichment_score_html(value):
    """Colours the NES value: green (positive), red (negative), grey (zero)."""
    try:
        v = float(value)
    except (TypeError, ValueError):
        return value
    color = "#008000" if v > 0 else "#b00000" if v < 0 else "#555555"
    return f'<span style="color:{color}; font-weight:bold;">{v:.3f}</span>'


def render_db_section(f, db, comparison, graphics_dir):
    """Writes the HTML section (GO or KEGG) of a comparison: dotplot + table."""
    name = "GO terms" if db == "GO" else "KEGG pathways"
    plot_col_name = "Enrichtment plot" if db == "GO" else "Pathway map"

    f.write(f"""
<h2 id="{db.lower()}" class="section-title">{db} enrichment</h2>
<h3 id="{db.lower()}-graph">Summary plot</h3>
""")
    dotplot_path = os.path.join(dotplot_dir, f"{db}_{comparison}_dotplot.png")
    dotplot_rel = os.path.relpath(dotplot_path, reports)
    add_image(dotplot_path, dotplot_rel, f, "main-image")

    csv_path = os.path.join(top_csv_dir, f"{db}_{comparison}_RNAm_miRNA_Top.csv")
    table_title = "Top 20 GO terms (10 up / 10 down by NES)" if db == "GO" else "All significant KEGG pathways (padj < 0.05)"
    f.write(f'<div id="{db.lower()}-table" class="table-title">{table_title}</div><table>')

    if os.path.exists(csv_path):
        df = pd.read_csv(csv_path, sep="\t", quotechar='"')

        headers = [plot_col_name, "ID", "Description", "NES", "p.adjust", "Target Genes", "miRNAs"]
        f.write("<tr>" + "".join(f"<th>{h}</th>" for h in headers) + "</tr>")

        for _, row in df.iterrows():
            row = row.to_dict()
            pathway_id = row.get("ID", "")

            if db == "GO":
                go_id_us = str(pathway_id).replace(":", "_")
                matches = glob.glob(os.path.join(graphics_dir, f"cProf.GOgseaplot.*.{go_id_us}.png"))
                plot_abs = matches[0] if matches else None
                plot_src = os.path.relpath(plot_abs, reports) if plot_abs else None
                plot_title = None
            else:
                local_path = os.path.join(graphics_dir, f"{pathway_id}.pathview.png")
                if os.path.exists(local_path) and not is_blank_image(local_path):
                    plot_src = os.path.relpath(local_path, reports)
                    plot_title = None
                else:
                    org_match = re.match(r'^([A-Za-z]+)', str(pathway_id))
                    org = org_match.group(1) if org_match else ""
                    plot_src = f"https://www.kegg.jp/kegg/pathway/{org}/{pathway_id}.png"
                    plot_title = "Official KEGG image (the locally generated one did not render correctly)"

            if plot_src:
                title_attr = f' title="{plot_title}"' if plot_title else ""
                plot_html = (
                    f'<a href="{plot_src}" target="_blank">'
                    f'<img src="{plot_src}" class="mini-plot"{title_attr}></a>'
                )
            else:
                plot_html = "-"

            padj = row.get("p.adjust")
            padj_str = f"{padj:.4e}" if isinstance(padj, (int, float)) else padj

            values = [
                (plot_html, "plot-cell"),
                (pathway_id, "id-cell"),
                (row.get("Description", ""), ""),
                (enrichment_score_html(row.get("NES")), ""),
                (padj_str, ""),
                (target_genes_cell_html(row.get("Target Genes"), row.get("Target_Details")), ""),
                (mirna_cell_html(row.get("miRNAs"), row.get("miRNA_Details")), ""),
            ]
            f.write(
                "<tr>"
                + "".join(f'<td class="{cls}">{v}</td>' if cls else f"<td>{v}</td>" for v, cls in values)
                + "</tr>"
            )

    f.write("</table>")

    f.write(f'<p class="note-text"><strong>NOTE:</strong> Click on the "{plot_col_name}" image to enlarge it.</p>')
    f.write('<p class="note-text"><strong>NOTE:</strong> Hover over a Target Gene to see its Score, which estimates the strength of the interaction between a miRNA and its target gene.</p>')
    f.write('<p class="note-text"><strong>NOTE:</strong> Hover over a miRNA to see its L2FC (Log2FoldChange), the value obtained from DESeq2.</p>')
    f.write('<p class="note-text"><strong>NOTE:</strong> NES: Normalized Enrichment Score.</p>')


def render_report(comparison):
    """Generates the full HTML file (GO + KEGG) for a given comparison."""
    graphics_dir = os.path.join(graphics_root, comparison)
    file_name = f"{comparison}.html"
    file_path = os.path.join(reports, file_name)

    with open(file_path, "w", encoding="utf-8") as f:
        f.write(f"""
<html>
<head>
<meta charset="UTF-8">
<title>Functional Enrichment Report</title>
<style>{css}</style>
</head>
<body>
<h1>Functional Enrichment Analysis (GSEA)</h1>
<div class="subtitle">{display_name(comparison)} | maxSize = {MAX_SIZE}</div>
<p><a href="index.html">&larr; Back to index</a></p>
<div id="index">
<h3>Index</h3>
<ul>
<li><a href="#go">GO enrichment</a>
    <ul><li><a href="#go-graph">Summary plot</a></li><li><a href="#go-table">Top GOs</a></li></ul>
</li>
<li><a href="#kegg">KEGG enrichment</a>
    <ul><li><a href="#kegg-graph">Summary plot</a></li><li><a href="#kegg-table">Top KEGG pathways</a></li></ul>
</li>
</ul>
</div>
""")
        render_db_section(f, "GO", comparison, graphics_dir)
        render_db_section(f, "KEGG", comparison, graphics_dir)
        f.write("</body></html>")

    return file_name


def render_index(comparisons_and_files):
    """Generates the index.html linking to every comparison report."""
    index_path = os.path.join(reports, "index.html")
    with open(index_path, "w", encoding="utf-8") as f:
        f.write(f"""
<html>
<head>
<meta charset="UTF-8">
<title>Functional Enrichment Reports</title>
<style>{css}</style>
</head>
<body>
<h1>Functional Enrichment Analysis (GSEA)</h1>
<div id="index">
<h3>Comparisons</h3>
<ul>
""")
        for comparison, file_name in comparisons_and_files:
            f.write(f'<li><a href="{file_name}">{display_name(comparison)}</a></li>\n')
        f.write("</ul></div></body></html>")


def main():
    comparisons = discover_comparisons()
    if not comparisons:
        print(f"No GO tables found in {top_csv_dir}")
        return

    comparisons_and_files = []
    for comparison in comparisons:
        file_name = render_report(comparison)
        comparisons_and_files.append((comparison, file_name))
        print(f"  - {comparison} -> reports/{file_name}")

    render_index(comparisons_and_files)
    print(f"\nHTML reports generated in '{reports}' ({len(comparisons)} comparisons + index.html)")


if __name__ == "__main__":
    main()
