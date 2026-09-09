import sys
sys.path.insert(0, snakemake.scriptdir)
from datetime import datetime, timezone
from pathlib import Path
import subprocess
from runtime import logged, sha256, versions, write_json

with logged(snakemake):
    source = Path(snakemake.input.snakefile).parent.parent
    try:
        commit = subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip()
        dirty = bool(subprocess.check_output(['git', '-C', str(source), 'status', '--porcelain'], text=True).strip())
    except (FileNotFoundError, subprocess.CalledProcessError):
        commit, dirty = 'unavailable', None
    write_json(snakemake.output[0], {
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'git_commit': commit, 'git_dirty': dirty,
        'config': snakemake.params.cfg, 'samples': snakemake.params.samples,
        'source_sha256': {str(p): sha256(p) for p in snakemake.input},
        'versions': versions(), 'count_source': 'kb only',
    })
