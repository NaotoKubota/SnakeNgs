import sys
sys.path.insert(0, snakemake.scriptdir)
from pathlib import Path
import shutil
import celltypist
from runtime import logged, sha256, versions, write_json

with logged(snakemake):
    p = snakemake.params.cfg['celltypist']
    if snakemake.input:
        source = snakemake.input[0]
    else:
        celltypist.models.download_models(force_update=False, model=p['model'])
        source = celltypist.models.get_model_path(p['model'])
    digest = sha256(source)
    if p['model_sha256'] and p['model_sha256'].lower() != digest:
        raise ValueError('CellTypist model SHA256 differs from config')
    shutil.copyfile(source, snakemake.output.model)
    model = celltypist.models.Model.load(snakemake.output.model)
    write_json(snakemake.output.metadata, {'name': p['model'], 'source': str(source),
        'sha256': digest, 'n_genes': len(model.features), 'n_classes': len(model.cell_types), 'versions': versions()})
