import os
from pathlib import Path
import shutil
import subprocess
import sys
import yaml
import pytest
from make_fixture import make_fixture

ROOT = Path(__file__).resolve().parents[2]


def test_dag_counts_and_fastq(tmp_path):
    executable = shutil.which('snakemake') or str(Path(sys.executable).parent / 'snakemake')
    if not Path(executable).exists():
        pytest.skip('snakemake not installed')
    path = make_fixture(tmp_path)
    cfg = yaml.safe_load(path.read_text())
    for mode, build, kb_workflow in [('counts', False, 'nac'), ('fastq', True, 'nac'),
                                    ('fastq', False, 'nac'), ('fastq', True, 'standard'),
                                    ('fastq', False, 'standard')]:
        cfg['input_mode'] = mode
        cfg['kb']['workflow'] = kb_workflow
        cfg['reference']['build'] = build
        if mode == 'fastq':
            for key in ['fasta', 'gtf', 'index', 't2g', 'cdna_t2c', 'intron_t2c']:
                ref = tmp_path / key
                ref.touch()
                cfg['reference'][key] = str(ref)
            refdir = tmp_path / 'cr_reference'
            refdir.mkdir(exist_ok=True)
            cfg['cellranger']['transcriptome'] = str(refdir)
            lines = ['sample\tR1\tR2\tbatch\tdonor\tcondition']
            for sample in ['a', 'b']:
                for read in ['R1', 'R2']: (tmp_path / f'{sample}_{read}.fastq.gz').touch()
                lines.append(f'{sample}\t{tmp_path}/{sample}_R1.fastq.gz\t{tmp_path}/{sample}_R2.fastq.gz\t{sample}\t{sample}\tcontrol')
            (tmp_path / 'samples.tsv').write_text('\n'.join(lines) + '\n')
        path.write_text(yaml.safe_dump(cfg))
        result = subprocess.run([executable, '-s', str(ROOT / 'snakefile/scrna_seq.smk'),
                                  '--configfile', str(path), '--cores', '2', '--dry-run'],
                                  text=True, capture_output=True, cwd=ROOT)
        assert result.returncode == 0, result.stdout + result.stderr
        for stage in ['sample_qc', 'soupx', 'scvi', 'celltypist', 'scanvi', 'report']:
            assert stage in result.stdout
        if mode == 'fastq':
            assert 'kb_count' in result.stdout and 'cellranger_metrics' in result.stdout


def test_all_rules_declare_container_and_no_host_run():
    source = (ROOT / 'snakefile/scrna_seq.smk').read_text()
    import re
    rules = re.split(r'^\s*rule \w+:', source, flags=re.MULTILINE)[1:]
    assert len(rules) >= 10
    assert all(re.search(r'^\s*container:', body, flags=re.MULTILINE) for body in rules)
    assert not re.search(r'^\s*run:', source, flags=re.MULTILINE)
