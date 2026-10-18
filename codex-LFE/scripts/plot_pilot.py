"""Standalone diagnostic figure from a completed, immutable pilot report."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument('run')
a = ap.parse_args()
run = Path(a.run)
report = json.loads((run / 'report.json').read_text())
fig, axs = plt.subplots(2, 2, figsize=(11, 7.4), constrained_layout=True)
colors = {'calibration': '#2166ac', 'test': '#b35806'}
for split, scores in report['splits'].items():
    thresholds = sorted(float(t) for t in scores)
    rows = [scores[f'{t:g}'] for t in thresholds]
    obs = [r['reference']['n_detections'] for r in rows]
    null = [np.mean(r['null']['control_counts']) if r['null'] else np.nan for r in rows]
    recall = [r['reference']['same_family_recall'] for r in rows]
    agreement = [r['reference']['reference_agreement_fraction'] for r in rows]
    axs[0, 0].plot(thresholds, obs, 'o-', color=colors[split], label=f'{split}: observed')
    axs[0, 0].plot(thresholds, null, 'o--', color=colors[split], alpha=.6, label=f'{split}: mean control')
    axs[0, 1].plot(thresholds, recall, 'o-', color=colors[split], label=split)
    axs[1, 0].plot(thresholds, agreement, 'o-', color=colors[split], label=split)
axs[0, 0].set(yscale='symlog', ylabel='Final candidates (symlog)', title='Observed versus site-shift controls')
axs[0, 1].set(ylabel='Correct-family reference recall', ylim=(0, 1), title='Recovery after final family assignment')
axs[1, 0].set(ylabel='Fraction matching reference', ylim=(0, 1), title='Reference agreement — not measured precision')
for ax in [axs[0, 0], axs[0, 1], axs[1, 0]]:
    ax.set_xlabel('Network threshold (MAD units)')
    ax.grid(alpha=.2)
    ax.legend(fontsize=8)
    ax.axvline(report['selection']['diagnostic_threshold'], color='gray', linewidth=1, alpha=.5)
qc = pd.read_csv(run / 'template_qc.csv')
for k, (family, g) in enumerate(qc.groupby('family_id')):
    axs[1, 1].plot(g.station, g.split_stack_cc, 'o-', label=f'Source group {k + 1}')
axs[1, 1].axhline(.25, linestyle='--', color='gray', label='Training QC threshold')
axs[1, 1].set(ylabel='Alternating training-half stack correlation', ylim=(0, 1), title='Templates built from original raw data')
axs[1, 1].tick_params(axis='x', rotation=20)
axs[1, 1].legend(fontsize=8)
axs[1, 1].grid(alpha=.2)
fig.suptitle('codex-LFE: independent northern CPU pilot\nExperimental candidates; calibration gate ' +
             ('passed' if report['selection']['calibration_gate_passed'] else 'not passed'), fontsize=14)
fig.savefig(run / 'pilot_diagnostics.png', dpi=170)
fig.savefig(run / 'pilot_diagnostics.pdf')
print(run / 'pilot_diagnostics.png')
