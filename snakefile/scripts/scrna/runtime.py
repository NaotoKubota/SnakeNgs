"""Small stdlib helpers also usable in the kb container."""
import contextlib
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys
import traceback


def write_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w') as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False, allow_nan=False, default=str)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def versions():
    result = {'python': sys.version.split()[0]}
    for name in ['kb-python', 'anndata', 'scanpy', 'scvi-tools', 'celltypist', 'scrublet',
                 'torch', 'numpy', 'scipy', 'pandas', 'matplotlib', 'scikit-learn']:
        try:
            result[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            pass
    return result


def command(argv):
    print('Executing:', repr([str(a) for a in argv]), flush=True)
    subprocess.run([str(a) for a in argv], check=True, stdout=sys.stdout, stderr=sys.stderr)


@contextlib.contextmanager
def logged(sm):
    Path(sm.log[0]).parent.mkdir(parents=True, exist_ok=True)
    for out in sm.output:
        Path(out).parent.mkdir(parents=True, exist_ok=True)
    with open(sm.log[0], 'w', buffering=1) as handle:
        with contextlib.redirect_stdout(handle), contextlib.redirect_stderr(handle):
            try:
                print('Versions:', versions())
                yield
            except Exception:
                traceback.print_exc()
                raise
