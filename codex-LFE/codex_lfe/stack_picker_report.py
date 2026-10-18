"""Summarize saved stack diagnostics with family-cluster uncertainty."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .provenance import write_json,sha256


def main():
    root=Path('results/picker_labeled_stacks_tail_v1')
    f=pd.read_parquet(root/'arrival_results.parquet')
    f=f[f.variant.eq('native_Z12')]
    rng=np.random.default_rng(20260914)
    rows=[]
    fig,axes=plt.subplots(2,2,figsize=(11,7))
    for j,phase in enumerate(['P','S']):
        phase_rows=f[f.phase.eq(phase)]
        for tol in [.25,.5,1.0]:
            points=[]
            for threshold,g in phase_rows.groupby('threshold'):
                g=g.copy();g['matched']=g.covered&(abs(g.residual_s)<=tol)
                family=g.groupby('source_family_id').agg(n=('covered','sum'),k=('matched','sum'))
                n=family.n.to_numpy();k=family.k.to_numpy()
                samples=rng.integers(0,len(family),size=(2000,len(family)))
                den=n[samples].sum(axis=1);num=k[samples].sum(axis=1)
                ci=np.quantile(num[den>0]/den[den>0],[.025,.975])
                recall=k.sum()/n.sum()
                rows.append(dict(phase=phase,threshold=threshold,tolerance_s=tol,n_labels=len(g),n_covered=int(n.sum()),n_matched=int(k.sum()),reference_recall=float(recall),family_bootstrap_ci95_low=float(ci[0]),family_bootstrap_ci95_high=float(ci[1])))
                points.append((threshold,recall))
            x,y=zip(*points);axes[0,j].plot(x,y,marker='o',label=f'±{tol:g} s')
        axes[0,j].set(title=f'{phase}: published-stack label recovery',xlabel='Model score threshold',ylabel='Fraction of covered labels recovered',ylim=(0,1))
        axes[0,j].legend();axes[0,j].grid(alpha=.2)
        sub=phase_rows[phase_rows.threshold.eq(.3)&phase_rows.covered]
        residual=sub.residual_s.dropna()
        within=residual[abs(residual)<=1]
        axes[1,j].hist(within,bins=np.arange(-1.025,1.051,.05),color='#337c92')
        axes[1,j].axvline(0,color='black',lw=1)
        axes[1,j].set(title=f'{phase}: residuals within ±1 s at score 0.3\n{len(within)}/{len(sub)} covered labels; misses excluded',xlabel='Predicted minus published label (s)',ylabel='Stack arrivals')
    fig.suptitle('Picker diagnostic on published stacks — not continuous-event precision')
    fig.tight_layout();fig.savefig(root/'timing_diagnostic.png',dpi=160);plt.close(fig)
    table=pd.DataFrame(rows);table.to_csv(root/'family_cluster_intervals.csv',index=False)
    write_json(root/'report_verification.json',dict(status='saved_results_summarized',source_sha256=sha256(root/'arrival_results.parquet'),report_source_sha256=sha256(__file__),bootstrap='2000 resamples of source families; descriptive uncertainty on existing stack population, not natural-event completeness',production_picker_validated=False))
    chosen=table[(table.threshold==.3)&table.tolerance_s.isin([.5,1.0])]
    lines=['# Published-stack picker timing diagnostic','', 'This missing diagnostic uses local jms weights (documented training excludes Cascadia), existing published nominal P/S labels and no waveform downloads. Earlier continuous and injection experiments remain available and were not rerun.','', 'The initial fixed-grid scan excluded most S labels because short records did not supply the next full 60-second window. The tail check evaluates only previously uncovered labels in a final unpadded window placed by record length; already covered results are reused unchanged. Production inference code has not been changed.','', '| Phase | Tolerance | Covered labels | Matches | Recovery | Family-bootstrap 95% interval |','|---|---:|---:|---:|---:|---:|']
    for r in chosen.itertuples():lines.append(f'| {r.phase} | ±{r.tolerance_s:g} s | {r.n_covered} | {r.n_matched} | {100*r.reference_recall:.1f}% | {100*r.family_bootstrap_ci95_low:.1f}–{100*r.family_bootstrap_ci95_high:.1f}% |')
    lines+=['', 'All table rows use score 0.3 as a diagnostic, not a selected production threshold. Other predeclared thresholds and tolerances are in the CSV.','', 'Limits: stacks average many events; all P labels are aligned at sample 800; source families overlap the original row-level splits; label uncertainties are absent; unmatched peaks are not proven false LFEs. The original PO station identifiers are preserved rather than silently relabeled as current roster stations. CPU/GPU verification on one stack found identical validity masks and maximum probability difference 0.000858; this is not exhaustive numerical validation.','', 'The timing figure and intervals describe stack behavior only. Natural recordings with independent phase/event review, continuous usable exposure, instrument transfer and independent-site association remain necessary before full-history production.']
    (root/'README.md').write_text('\n'.join(lines)+'\n')
    print(chosen.to_string(index=False))

if __name__=='__main__':main()
