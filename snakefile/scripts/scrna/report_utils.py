import base64
import html
from pathlib import Path
import numpy as np
import pandas as pd


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
