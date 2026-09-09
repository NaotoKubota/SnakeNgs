import numpy as np
import pandas as pd


def configure_training(cfg, threads):
    import torch
    import scvi
    import scanpy as sc
    scvi.settings.seed = cfg['seed']
    torch.set_num_threads(threads)
    sc.settings.n_jobs = threads
    if cfg['integration']['accelerator'] == 'gpu' and not torch.cuda.is_available():
        raise ValueError('GPU requested but unavailable; check --singularity-args --nv and scheduler allocation')


def save_history(model, path):
    records = []
    for metric, frame in model.history.items():
        for step, row in frame.iterrows():
            for column, value in row.items():
                records.append({'step': step, 'metric': metric, 'column': column, 'value': float(value)})
    pd.DataFrame(records, columns=['step', 'metric', 'column', 'value']).to_csv(path, sep='\t', index=False)
