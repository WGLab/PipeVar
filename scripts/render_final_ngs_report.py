#!/usr/bin/env python3
import argparse
import csv
import gzip
import html
import json
from datetime import datetime
from pathlib import Path


MAX_CELL_CHARS = 600


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace", newline="")


def shorten(value):
    value = "" if value is None else str(value)
    if len(value) <= MAX_CELL_CHARS:
        return value
    return value[:MAX_CELL_CHARS] + "..."


def esc(value):
    return html.escape(shorten(value), quote=True)


def parse_info(info_text):
    parsed = {}
    if not info_text or info_text == ".":
        return parsed
    for item in info_text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            parsed[key] = value
        else:
            parsed[item] = "true"
    return parsed


def read_vcf(path):
    rows = []
    sample_names = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.lstrip("#").split("\t")
                sample_names = header[9:] if len(header) > 9 else []
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            info = parse_info(parts[7])
            row = {
                "Rank": info.get("PRIO_RANK", "."),
                "Gene": info.get("Gene", info.get("Gene.refGene", ".")),
                "Chrom": parts[0],
                "Pos": parts[1],
                "Ref": parts[3],
                "Alt": parts[4],
                "Filter": parts[6],
                "Priority": info.get("PRIO_CAT", "."),
                "Score": info.get("PRIO_MAX_SCORE", "."),
                "Group": info.get("PRIO_GROUP_SIZE", "."),
                "Source": info.get("SOURCE", "."),
                "Model": info.get("MODEL", "."),
                "Genotype": ".",
                "Details": line,
            }
            if len(parts) >= 10:
                fmt = parts[8].split(":")
                sample = parts[9].split(":")
                if "GT" in fmt:
                    gt_idx = fmt.index("GT")
                    if gt_idx < len(sample):
                        row["Genotype"] = sample[gt_idx]
                elif sample:
                    row["Genotype"] = sample[0]
            rows.append(row)
    meta = {"samples": ", ".join(sample_names) if sample_names else "."}
    return rows, meta


def read_tsv(path, comment_prefix=None):
    rows = []
    with open_text(path) as handle:
        lines = [line.rstrip("\n") for line in handle if line.rstrip("\n")]
    if not lines:
        return []

    header_idx = 0
    if comment_prefix:
        for idx, line in enumerate(lines):
            if line.startswith(comment_prefix):
                header_idx = idx
                break

    header = lines[header_idx].split("\t")
    header = [h.strip() for h in header]
    for line in lines[header_idx + 1:]:
        if not line or line.startswith("##"):
            continue
        parts = line.split("\t")
        row = {}
        for idx, col in enumerate(header):
            row[col] = parts[idx] if idx < len(parts) else ""
        if len(parts) > len(header):
            row["Extra"] = "\t".join(parts[len(header):])
        rows.append(row)
    return rows


def read_prio_gene_report(path):
    rows = []
    header = None
    with open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("#RANK"):
                header = [col.strip().lstrip("#") for col in line.split("\t")]
                continue
            if header is None or line.startswith("##"):
                continue
            parts = line.split("\t")
            row = {}
            for idx, col in enumerate(header):
                if col == "VCF_LINE":
                    row[col] = "\t".join(parts[idx:])
                    break
                row[col] = parts[idx] if idx < len(parts) else ""
            rows.append(row)
    return rows


def read_repeat_report(path):
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    nonempty_lines = [line for line in text.splitlines() if line.strip()]
    if nonempty_lines and "\t" in nonempty_lines[0] and not nonempty_lines[0].lstrip().startswith("["):
        return read_tsv(path)
    try:
        start = text.index("[")
        end = text.rindex("]") + 1
        data = json.loads(text[start:end])
        if isinstance(data, list):
            return [
                {
                    "LocusId": item.get("LocusId", "."),
                    "Genotype": item.get("Genotype", "."),
                    "GenotypeConfidenceInterval": item.get("GenotypeConfidenceInterval", "."),
                    "Pathogenic": "yes",
                }
                for item in data
            ]
    except (ValueError, json.JSONDecodeError):
        pass
    return [{"Status": line} for line in text.splitlines() if line.strip()]


def table_html(rows, empty_message):
    if not rows:
        return f'<div class="empty">{esc(empty_message)}</div>'
    columns = list(rows[0].keys())
    for row in rows[1:]:
        for key in row:
            if key not in columns:
                columns.append(key)
    head = "".join(f"<th>{esc(col)}</th>" for col in columns)
    body_rows = []
    for row in rows:
        cells = "".join(f"<td>{esc(row.get(col, ''))}</td>" for col in columns)
        body_rows.append(f"<tr>{cells}</tr>")
    return (
        '<div class="table-wrap">'
        '<input class="filter" type="search" placeholder="Filter rows" aria-label="Filter rows">'
        '<table><thead><tr>'
        + head
        + "</tr></thead><tbody>"
        + "".join(body_rows)
        + "</tbody></table></div>"
    )


def tab_button(tab_id, label, count, active=False):
    selected = "true" if active else "false"
    klass = "tab active" if active else "tab"
    return (
        f'<button class="{klass}" role="tab" aria-selected="{selected}" '
        f'aria-controls="{tab_id}" data-tab="{tab_id}">{esc(label)} '
        f'<span>{count}</span></button>'
    )


def tab_panel(tab_id, title, body, active=False):
    hidden = "" if active else " hidden"
    return (
        f'<section id="{tab_id}" class="panel"{hidden} role="tabpanel">'
        f"<h2>{esc(title)}</h2>{body}</section>"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--prio-vcf", required=True)
    parser.add_argument("--prio-gene-vcf", required=True)
    parser.add_argument("--repeat-tsv")
    parser.add_argument("--mito-tsv")
    parser.add_argument("--output-html", required=True)
    parser.add_argument("--title", default="PipeVar final variant report")
    parser.add_argument("--mito-min-af", default="0.5")
    parser.add_argument("--mito-min-apogee2", default="0.5")
    parser.add_argument("--mito-min-mitotip", default="12.66")
    args = parser.parse_args()

    prio_vcf = Path(args.prio_vcf)
    prio_gene = Path(args.prio_gene_vcf)
    if not prio_vcf.exists():
        raise SystemExit(f"Missing required prioritized variant VCF: {prio_vcf}")
    if not prio_gene.exists():
        raise SystemExit(f"Missing required prioritized gene report: {prio_gene}")

    variant_rows, vcf_meta = read_vcf(prio_vcf)
    gene_rows = read_prio_gene_report(prio_gene)
    mito_path = Path(args.mito_tsv) if args.mito_tsv else None
    repeat_path = Path(args.repeat_tsv) if args.repeat_tsv else None
    mito_rows = read_tsv(mito_path) if mito_path and mito_path.exists() else []
    repeat_rows = read_repeat_report(repeat_path) if repeat_path and repeat_path.exists() else []

    files = [
        {"Input": "Prioritized variants", "File": prio_vcf.name, "Status": "loaded"},
        {"Input": "Prio gene report", "File": prio_gene.name, "Status": "loaded"},
        {
            "Input": "Mitochondrial",
            "File": mito_path.name if mito_path else ".",
            "Status": "loaded" if mito_path and mito_path.exists() else "not run or no file provided",
        },
        {
            "Input": "Repeat expansion",
            "File": repeat_path.name if repeat_path else ".",
            "Status": "loaded" if repeat_path and repeat_path.exists() else "not run or no file provided",
        },
    ]
    provenance = [
        {"Field": "Sample", "Value": args.sample},
        {"Field": "Generated", "Value": datetime.now().isoformat(timespec="seconds")},
        {"Field": "VCF samples", "Value": vcf_meta.get("samples", ".")},
        {"Field": "Mito AF threshold", "Value": f">{args.mito_min_af}"},
        {"Field": "Mito APOGEE2 threshold", "Value": f">{args.mito_min_apogee2}"},
        {"Field": "Mito MitoTip threshold", "Value": f">{args.mito_min_mitotip}"},
        {
            "Field": "Ranking note",
            "Value": "Sections are displayed side by side; mitochondrial and repeat expansion findings are not cross-prioritized against nuclear variants.",
        },
    ] + files

    tabs = [
        ("prio-gene", "Gene/Scenario Prioritization", gene_rows, "No prioritized gene/scenario rows were found."),
        ("variants", "Prioritized Nuclear Variants", variant_rows, "No prioritized variant rows were found."),
        ("mitochondrial", "Mitochondrial Variants", mito_rows, "Mitochondrial analysis was not run or no mitochondrial rows passed filters."),
        ("repeat-expansion", "STR/Repeat Expansions", repeat_rows, "Repeat expansion analysis was not run or no pathogenic loci were reported."),
        ("provenance", "Provenance", provenance, "No provenance available."),
    ]

    nav = "".join(tab_button(tab_id, label, len(rows), active=(idx == 0)) for idx, (tab_id, label, rows, _) in enumerate(tabs))
    panels = "".join(
        tab_panel(tab_id, label, table_html(rows, empty), active=(idx == 0))
        for idx, (tab_id, label, rows, empty) in enumerate(tabs)
    )

    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(args.title)} - {esc(args.sample)}</title>
<style>
:root {{
  color-scheme: light;
  --ink: #18212f;
  --muted: #5c6675;
  --line: #d7dde5;
  --panel: #ffffff;
  --soft: #f5f7fa;
  --accent: #0f766e;
  --accent-2: #7c2d12;
}}
* {{ box-sizing: border-box; }}
body {{
  margin: 0;
  font-family: Arial, Helvetica, sans-serif;
  color: var(--ink);
  background: #eef2f6;
}}
header {{
  padding: 22px 28px 16px;
  background: #ffffff;
  border-bottom: 1px solid var(--line);
}}
h1 {{ margin: 0 0 8px; font-size: 26px; font-weight: 700; }}
.meta {{ color: var(--muted); font-size: 14px; }}
.counts {{ display: flex; flex-wrap: wrap; gap: 10px; margin-top: 14px; }}
.count {{ border: 1px solid var(--line); border-radius: 6px; padding: 7px 10px; background: var(--soft); font-size: 13px; }}
main {{ padding: 18px 24px 28px; }}
nav {{ display: flex; flex-wrap: wrap; gap: 8px; margin-bottom: 14px; }}
.tab {{
  border: 1px solid var(--line);
  border-radius: 6px;
  background: #fff;
  color: var(--ink);
  padding: 9px 12px;
  cursor: pointer;
  font-size: 14px;
}}
.tab span {{ color: var(--muted); margin-left: 5px; }}
.tab.active {{ border-color: var(--accent); box-shadow: inset 0 -3px 0 var(--accent); }}
.panel {{
  background: var(--panel);
  border: 1px solid var(--line);
  border-radius: 8px;
  padding: 16px;
}}
h2 {{ margin: 0 0 12px; font-size: 19px; }}
.empty {{ padding: 18px; background: var(--soft); border: 1px dashed var(--line); border-radius: 6px; color: var(--muted); }}
.table-wrap {{ overflow: auto; }}
.filter {{
  width: min(420px, 100%);
  margin-bottom: 10px;
  padding: 8px 10px;
  border: 1px solid var(--line);
  border-radius: 6px;
  font-size: 14px;
}}
table {{ border-collapse: collapse; width: 100%; min-width: 760px; font-size: 13px; }}
th, td {{ border: 1px solid var(--line); padding: 7px 9px; text-align: left; vertical-align: top; }}
th {{ position: sticky; top: 0; background: #f8fafc; cursor: pointer; }}
td {{ max-width: 520px; overflow-wrap: anywhere; }}
</style>
</head>
<body>
<header>
  <h1>{esc(args.title)}</h1>
  <div class="meta">Sample: <strong>{esc(args.sample)}</strong> · Generated: {esc(datetime.now().isoformat(timespec="seconds"))}</div>
  <div class="counts">
    <div class="count">Gene/Scenario: {len(gene_rows)}</div>
    <div class="count">Nuclear Variants: {len(variant_rows)}</div>
    <div class="count">Mitochondrial Variants: {len(mito_rows)}</div>
    <div class="count">STR/Repeats: {len(repeat_rows)}</div>
  </div>
</header>
<main>
<nav role="tablist">{nav}</nav>
{panels}
</main>
<script>
document.querySelectorAll('.tab').forEach(button => {{
  button.addEventListener('click', () => {{
    document.querySelectorAll('.tab').forEach(b => {{
      b.classList.remove('active');
      b.setAttribute('aria-selected', 'false');
    }});
    document.querySelectorAll('.panel').forEach(p => p.hidden = true);
    button.classList.add('active');
    button.setAttribute('aria-selected', 'true');
    document.getElementById(button.dataset.tab).hidden = false;
  }});
}});
document.querySelectorAll('.filter').forEach(input => {{
  input.addEventListener('input', () => {{
    const table = input.closest('.table-wrap').querySelector('tbody');
    const needle = input.value.toLowerCase();
    table.querySelectorAll('tr').forEach(row => {{
      row.hidden = !row.textContent.toLowerCase().includes(needle);
    }});
  }});
}});
document.querySelectorAll('th').forEach((th, idx) => {{
  th.addEventListener('click', () => {{
    const table = th.closest('table');
    const tbody = table.querySelector('tbody');
    const rows = Array.from(tbody.querySelectorAll('tr'));
    const asc = th.dataset.asc !== 'true';
    rows.sort((a, b) => {{
      const av = a.children[idx].textContent;
      const bv = b.children[idx].textContent;
      const an = Number(av);
      const bn = Number(bv);
      if (!Number.isNaN(an) && !Number.isNaN(bn)) return asc ? an - bn : bn - an;
      return asc ? av.localeCompare(bv) : bv.localeCompare(av);
    }});
    th.dataset.asc = asc ? 'true' : 'false';
    rows.forEach(row => tbody.appendChild(row));
  }});
}});
</script>
</body>
</html>
"""
    Path(args.output_html).write_text(document, encoding="utf-8")


if __name__ == "__main__":
    main()
