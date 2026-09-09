import base64
import colorsys
import html
from pathlib import Path
import numpy as np
import pandas as pd


def distinct_colors(labels, palette, overrides=None, fixed=None):
    """Return deterministic, non-repeating colors for one categorical plot."""
    labels = [str(label) for label in labels]
    overrides = {str(k): str(v) for k, v in (overrides or {}).items()}
    fixed = {str(k): str(v) for k, v in (fixed or {}).items()}
    colors = {}
    used = set()
    for label in labels:
        if label in fixed or label in overrides:
            color = fixed.get(label, overrides.get(label))
            normalized = color.lower()
            if normalized in used:
                raise ValueError(f'Duplicate configured color {color} in one categorical plot')
            colors[label] = color
            used.add(normalized)
    candidates = list(dict.fromkeys(color for color in palette if color.lower() not in used))
    # Golden-ratio hue spacing with alternating saturation/value extends short
    # journal palettes without recycling a color for a different cell type.
    for i in range(max(64, len(labels) * 4)):
        hue = (0.11 + i * 0.618033988749895) % 1
        saturation = (0.62, 0.78, 0.48)[i % 3]
        value = (0.82, 0.68)[(i // 3) % 2]
        rgb = colorsys.hsv_to_rgb(hue, saturation, value)
        color = '#{:02X}{:02X}{:02X}'.format(*(round(255 * channel) for channel in rgb))
        if color.lower() not in used:
            candidates.append(color)
            used.add(color.lower())
    used = {color.lower() for color in colors.values()}
    for label in labels:
        if label in colors:
            continue
        while candidates and candidates[0].lower() in used:
            candidates.pop(0)
        if not candidates:
            raise ValueError('Could not generate enough distinct categorical colors')
        colors[label] = candidates.pop(0)
        used.add(colors[label].lower())
    return colors


def embedded_file(path, label, mime='application/octet-stream'):
    path = Path(path)
    payload = base64.b64encode(path.read_bytes()).decode()
    return f'<a download="{html.escape(path.name, quote=True)}" href="data:{mime};base64,{payload}">{html.escape(label)}</a>'


def table(frame):
    return '<div class="table-wrap">' + frame.to_html(index=False, escape=True, border=0, float_format=lambda v: f'{v:.4g}') + '</div>'


def design_notes(obs, batch_key):
    # Do not count sequencing libraries or individual cells as biological replicates.
    notes = []
    donors = obs[['donor', 'condition']].drop_duplicates().groupby('condition', observed=True)['donor'].nunique()
    for condition, n in donors.items():
        if n < 2:
            notes.append(f'{condition}: donor n={n}. Cell counts do not represent biological replication.')
    cross = pd.crosstab(obs[batch_key].astype(str), obs['condition'].astype(str))
    if cross.shape[0] > 1 and cross.shape[1] > 1 and ((cross > 0).sum(axis=1) == 1).all():
        notes.append(f'{batch_key} and condition are completely confounded. Batch correction may remove biological differences; review the experimental design and the embeddings before and after correction.')
    if cross.shape[0] == 1:
        notes.append('Only one batch level is present. Mixing between batches cannot be evaluated.')
    return notes
