"""Render cargo-modules DOT output into readable SVGs.

Graphviz is not required: the module graphs contain only nodes and two kinds
of edges, so a small layout is enough. Run after regenerating the .dot files:

    python scripts/render_module_graphs.py
"""

import re
from pathlib import Path
from xml.sax.saxutils import escape

ROOT = Path(__file__).resolve().parent.parent
GRAPHS = ROOT / "docs" / "graphs"

NODE = re.compile(r'^\s*"([^"]+)"\s+\[label="([^"]+)",\s*fillcolor="([^"]+)"')
EDGE = re.compile(r'^\s*"([^"]+)"\s+->\s+"([^"]+)"\s+\[label="([^"]+)"')

BOX_W = 210
BOX_H = 42
COL_GAP = 270
ROW_GAP = 66
MARGIN = 30
HEADER = 58


def parse(dot_text):
    nodes = {}
    edges = []
    for line in dot_text.splitlines():
        match = NODE.match(line)
        if match:
            node_id, label, color = match.groups()
            nodes[node_id] = (label.split("|")[-1], color)
            continue
        match = EDGE.match(line)
        if match:
            source, target, kind = match.groups()
            edges.append((source, target, kind))
    return nodes, edges


def depth(node_id):
    return node_id.count("::")


def render(dot_path, svg_path):
    nodes, edges = parse(dot_path.read_text(encoding="utf-8"))
    levels = {}
    for node_id in nodes:
        levels.setdefault(depth(node_id), []).append(node_id)
    for level in levels.values():
        level.sort()

    positions = {}
    for level, ids in levels.items():
        for index, node_id in enumerate(ids):
            positions[node_id] = (
                MARGIN + level * COL_GAP,
                HEADER + index * ROW_GAP,
            )

    max_level = max(levels)
    rows = max(len(ids) for ids in levels.values())
    width = MARGIN * 2 + max_level * COL_GAP + BOX_W
    height = HEADER + rows * ROW_GAP + 70
    title = next(iter(nodes)).split("::")[0] if nodes else "modules"

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{MARGIN}" y="34" font-family="Segoe UI, sans-serif" font-size="22" font-weight="600" fill="#1f2937">{escape(title)} module graph</text>',
        '<defs>',
        '<marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="3" orient="auto">',
        '<path d="M0,0 L7,3 L0,6 Z" fill="#6b7280"/>',
        '</marker>',
        '</defs>',
    ]

    for offset, (source, target, kind) in enumerate(edges):
        if source not in positions or target not in positions:
            continue
        x1, y1 = positions[source]
        x2, y2 = positions[target]
        shift = 9 if kind == "uses" else 0
        if x2 >= x1 + BOX_W:
            start = (x1 + BOX_W, y1 + BOX_H / 2 + shift)
            end = (x2, y2 + BOX_H / 2 + shift)
            middle = (start[0] + end[0]) / 2
            path = (
                f"M {start[0]},{start[1]} C {middle},{start[1]} "
                f"{middle},{end[1]} {end[0]},{end[1]}"
            )
        else:
            start = (x1 + BOX_W / 2, y1 + BOX_H)
            end = (x2 + BOX_W / 2, y2 + BOX_H + shift)
            dip = max(start[1], end[1]) + 26 + (offset % 3) * 12
            middle = (start[0] + end[0]) / 2
            path = f"M {start[0]},{start[1]} Q {middle},{dip} {end[0]},{end[1]}"
        dash = ' stroke-dasharray="5 4"' if kind == "uses" else ""
        parts.append(
            f'<path d="{path}" fill="none" stroke="#9ca3af" stroke-width="1.4"'
            f'{dash} marker-end="url(#arrow)" data-kind="{escape(kind)}">'
            f'<title>{escape(source)} --{escape(kind)}--&gt; {escape(target)}</title></path>'
        )

    for node_id, (label, color) in nodes.items():
        x, y = positions[node_id]
        parts.append(
            f'<rect x="{x}" y="{y}" width="{BOX_W}" height="{BOX_H}" rx="6"'
            f' fill="{escape(color)}" stroke="#ffffff" stroke-width="1.5"/>'
            f'<text x="{x + BOX_W / 2}" y="{y + BOX_H / 2 + 4}" text-anchor="middle"'
            f' font-family="Consolas, monospace" font-size="12" fill="#111827">'
            f'{escape(label)}</text>'
        )

    parts.append(
        '<text x="{}" y="{}" font-family="Segoe UI, sans-serif" font-size="12" fill="#4b5563">'
        'solid: owns  |  dashed: uses</text>'.format(MARGIN, height - 10)
    )
    parts.append("</svg>")
    svg_path.write_text("\n".join(parts) + "\n", encoding="utf-8")
    return len(nodes), len(edges)


def render_index(names):
    panels = []
    for name in names:
        svg = (GRAPHS / f"{name}-modules.svg").read_text(encoding="utf-8")
        svg = svg.replace("<svg ", '<svg class="graph" ', 1)
        panels.append(
            f'<section><h2>{escape(name)}</h2>{svg}</section>'
        )
    html = (
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        '<title>Nodris module graphs</title><style>'
        'body{margin:0;padding:32px;background:#f9fafb;color:#1f2937;'
        'font-family:Segoe UI,Arial,sans-serif}'
        'h1{font-size:26px;margin:0 0 8px}'
        'p{color:#4b5563;margin:0 0 28px}'
        'section{background:#fff;border:1px solid #e5e7eb;border-radius:10px;'
        'padding:20px;margin-bottom:24px;overflow:auto}'
        'h2{font-size:17px;margin:0 0 12px}'
        '.graph{display:block;max-width:100%;height:auto}'
        '</style></head><body>'
        '<h1>Nodris module graphs</h1>'
        '<p>Generated with cargo-modules 0.28.0. Solid edges are ownership, '
        'dashed edges are uses.</p>'
        + "".join(panels)
        + '</body></html>\n'
    )
    (GRAPHS / "index.html").write_text(html, encoding="utf-8")


if __name__ == "__main__":
    names = ("nodris-core", "nodris-io")
    for name in names:
        dot = GRAPHS / f"{name}-modules.dot"
        svg = GRAPHS / f"{name}-modules.svg"
        node_count, edge_count = render(dot, svg)
        print(f"{svg.name}: {node_count} nodes, {edge_count} edges")
    render_index(names)
    print("index.html written")
