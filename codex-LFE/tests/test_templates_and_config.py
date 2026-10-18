import numpy as np
import pytest
import yaml
from pathlib import Path

from codex_lfe.pilot import validate_config
from codex_lfe.templates import build_station_template


def test_train_calibration_test_overlap_rejected():
    cfg = yaml.safe_load((Path(__file__).parents[1] / 'configs/cpu_pilot.yaml').read_text())
    validate_config(cfg)
    cfg['test_dates'] = cfg['train_dates'][:1]
    with pytest.raises(ValueError, match='overlap'):
        validate_config(cfg)


def test_independent_half_stack_recovers_injected_template():
    rng = np.random.default_rng(7)
    pulse = rng.normal(size=(3, 40))
    windows = rng.normal(size=(100, 3, 1400))
    windows[:, :, 600:640] += pulse * 2
    template, qc = build_station_template(windows, 20)
    assert qc['usable']
    assert qc['split_stack_cc'] > 0.9
    assert 19 < qc['peak_offset_from_published_s'] < 23
    assert template.shape == (3, 160)


def test_incoherent_or_too_few_training_examples_rejected():
    rng = np.random.default_rng(8)
    windows = rng.normal(size=(10, 3, 1400))
    template, qc = build_station_template(windows, 20)
    assert template is None and not qc['usable']


def test_sinc_template_cut_preserves_high_frequency_and_rejects_neighbor_gap():
    from types import SimpleNamespace
    from codex_lfe.templates import cut
    fs=20; t=np.arange(600)/fs
    x=np.tile(np.sin(2*np.pi*8*t),(3,1))
    wave=SimpleNamespace(x=x,valid=np.ones_like(x,dtype=bool),fs=fs,start_s=0.)
    first=100.5; count=200
    expected=np.tile(np.sin(2*np.pi*8*(first+np.arange(count))/fs),(3,1))
    actual=cut(wave,first/fs,count,alignment='sinc')
    assert np.max(np.abs(actual-expected)) < 3e-4
    linear=cut(wave,first/fs,count)
    assert np.sqrt(np.mean(linear**2)) < .4*np.sqrt(np.mean(expected**2))
    wave.valid[:,90]=False
    assert cut(wave,first/fs,count,alignment='sinc') is None
    assert cut(wave,first/fs,count) is not None
