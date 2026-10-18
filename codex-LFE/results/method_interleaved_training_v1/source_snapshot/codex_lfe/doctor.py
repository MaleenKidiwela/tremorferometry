"""Read-only input/hardware checks and an optional CUDA numerical parity check."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import yaml

from .matched_filter import normalized_cc
from .pilot import ROOT, validate_config
from .provenance import environment


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', default='configs/cpu_pilot.yaml')
    ap.add_argument('--cuda-parity', action='store_true')
    a = ap.parse_args()
    cfg = yaml.safe_load(Path(a.config).read_text())
    validate_config(cfg)
    try:
        import torch
    except ImportError:
        torch = None
    missing, total = [], 0
    for date in cfg['train_dates'] + cfg['calibration_dates'] + cfg['test_dates']:
        d = pd.Timestamp(date)
        for station in cfg['stations']:
            p = ROOT / cfg['raw_root'] / station / str(d.year) / f'{d.dayofyear:03d}.mseed'
            if not p.exists():
                missing.append(str(p))
            else:
                total += p.stat().st_size
    report = {'environment': environment(), 'cuda_available': bool(torch and torch.cuda.is_available()),
              'requested_raw_bytes': total, 'missing_raw_paths': missing,
              'reference_exists': (ROOT / cfg['reference_csv']).exists(),
              'site_count': len({m['site_id'] for m in cfg['stations'].values()})}
    if a.cuda_parity:
        if torch is None or not torch.cuda.is_available():
            raise RuntimeError('CUDA parity cannot run on a CPU-only pod')
        rng = np.random.default_rng(12)
        t = rng.normal(size=(3, 160))
        x = rng.normal(size=(3, 20000))
        x[:, 12000:12160] = t * 5
        valid = np.ones_like(x, bool)
        valid[:, 3000:4000] = False
        cpu, cm = normalized_cc(x, t, valid)
        gpu, gm = normalized_cc(x, t, valid, backend='torch', device='cuda')
        np.testing.assert_array_equal(cm, gm)
        np.testing.assert_allclose(cpu[cm], gpu[cm], atol=1e-8, rtol=1e-8)
        assert np.nanargmax(gpu) == 12000
        report['cuda_parity'] = 'passed'
    print(json.dumps(report, indent=2))
    if missing or not report['reference_exists']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
