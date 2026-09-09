import sys
sys.path.insert(0, snakemake.scriptdir)
from pathlib import Path
from runtime import logged

with logged(snakemake):
    out = Path(snakemake.output[0])
    out.mkdir(parents=True, exist_ok=True)
    # Only this dedicated rule-owned directory is regenerated on lane changes.
    for old in out.glob('*.fastq.gz'):
        if not old.is_symlink():
            raise ValueError(f'Refusing to replace non-symlink staged input: {old}')
        old.unlink()
    for read, files in [('R1', snakemake.input.r1), ('R2', snakemake.input.r2)]:
        for lane, path in enumerate(files, 1):
            (out / f'{snakemake.wildcards.sample}_S1_L{lane:03d}_{read}_001.fastq.gz').symlink_to(Path(path).resolve())
