import csv
import json
import re
import sys
from pathlib import Path

if len(sys.argv) != 2:
    print("Usage: python add_assembly_stats_to_icarus.py <quast_output_dir>", file=sys.stderr)
    sys.exit(1)

outdir = Path(sys.argv[1]).resolve()
report_tsv = outdir / "report.tsv"

if not report_tsv.exists():
    print(f"ERROR: report.tsv not found: {report_tsv}", file=sys.stderr)
    sys.exit(1)

html_paths = list(outdir.rglob("all_chromosomes.html"))

if not html_paths:
    print(f"ERROR: all_chromosomes.html not found under: {outdir}", file=sys.stderr)
    sys.exit(1)

raw = {}

with report_tsv.open(newline="") as f:
    reader = csv.reader(f, delimiter="\t")
    for row in reader:
        if len(row) >= 2:
            key = row[0].strip()
            value = row[1].strip()
            raw[key] = value

stats = {
    "assembly": raw.get("Assembly", ""),
    "total_length": raw.get("Total length", ""),
    "num_contigs": raw.get("# contigs", ""),
    "n50": raw.get("N50", ""),
    "l50": raw.get("L50", ""),
    "largest_contig": raw.get("Largest contig", ""),
    "genome_fraction": (
        raw.get("Genome fraction (%)", "") + "%"
        if raw.get("Genome fraction (%)", "")
        else ""
    ),
}

script = (
    '<script id="assembly-stats-data">\n'
    'window.assemblyStats = '
    + json.dumps(stats, ensure_ascii=False, indent=2)
    + ';\n'
    '</script>\n'
)

for html_path in html_paths:
    html = html_path.read_text(errors="replace")

    if 'id="assembly-stats-data"' in html:
        html = re.sub(
            r'<script id="assembly-stats-data">.*?</script>\s*',
            script,
            html,
            flags=re.S,
        )
    elif "</body>" in html:
        html = html.replace("</body>", script + "</body>")
    else:
        html += "\n" + script

    html_path.write_text(html)
    print(f"Injected assembly stats into: {html_path}")
    print(json.dumps(stats, ensure_ascii=False, indent=2))
