import copy
from pathlib import Path
import pytest
import yaml
from workflow_config import prepare_config
from report_utils import design_notes, table
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture
def config(tmp_path):
    defaults = yaml.safe_load((ROOT / 'examples/scrna_seq/config.yaml').read_text())
    sheet = tmp_path / 'samples.tsv'
    sheet.write_text('sample\tR1\tR2\tbatch\tdonor\tcondition\na\tr1.fastq.gz\tr2.fastq.gz\tb1\td1\tc\n')
    cfg = copy.deepcopy(defaults)
    cfg['sample_table'] = str(sheet)
    return cfg, defaults, tmp_path


def test_paths_resolved_before_workdir(config):
    cfg, defaults, base = config
    result, samples = prepare_config(cfg, defaults, base)
    assert samples['a']['R1'] == [str(base / 'r1.fastq.gz')]
    assert result['workdir'] == str(base / 'results/scrna_seq')


def test_duplicate_library_and_fastq_rejected(config):
    cfg, defaults, base = config
    path = Path(cfg['sample_table'])
    text = path.read_text()
    path.write_text(text + text.splitlines()[-1] + '\n')
    with pytest.raises(ValueError, match='Duplicate library'): prepare_config(cfg, defaults, base)
    path.write_text(text + 'b\tr1.fastq.gz\tr3.fastq.gz\tb1\td2\tc\n')
    with pytest.raises(ValueError, match='FASTQ reused'): prepare_config(cfg, defaults, base)


@pytest.mark.parametrize('section,key,value', [('depth', 'model_layer', 'lognorm'),
    ('soupx', 'contamination_fraction', 1.2), ('scanvi', 'posterior_samples', 1),
    ('scrublet', 'enabled', 'false'), ('integration', 'devices', 2)])
def test_invalid_parameters(config, section, key, value):
    cfg, defaults, base = config
    cfg[section][key] = value
    with pytest.raises(ValueError): prepare_config(cfg, defaults, base)


def test_unknown_setting_rejected(config):
    cfg, defaults, base = config
    cfg['scrublet']['threshhold'] = .3
    with pytest.raises(ValueError, match='Unknown setting'): prepare_config(cfg, defaults, base)


def test_design_and_escaping():
    obs = pd.DataFrame({'sample': ['s1', 's2'], 'donor': ['d1', 'd2'],
                        'condition': ['Ctrl', 'KO'], 'batch': ['a', 'b']})
    notes = design_notes(obs, 'batch')
    assert any('confounded' in n for n in notes)
    assert len([n for n in notes if 'donor n=1' in n]) == 2
    assert '<script>' not in table(pd.DataFrame({'x': ['<script>']}))
