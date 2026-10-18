"""Summarize corrected development and fixed next-day panel without truth inflation."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .provenance import sha256,write_json

PANELS={'development_2010_11_22':'ducellier_multifamily_november_edges_v2','next_day_withheld_B039':'ducellier_multifamily_november23_holdout_edges_v2'}

def main():
    out=Path('results/ducellier_validation_summary_v1');out.mkdir(exist_ok=True)
    rows=[];held_rows=[];sources={}
    for label,name in PANELS.items():
        root=Path('results')/name
        data=pd.read_csv(root/'summary.csv');data=data[data.threshold_mad==8]
        for family,table in data.groupby('source_template_id'):
            observed=table[table.control==-1].iloc[0];controls=table[table.control>=0]
            rows.append(dict(panel=label,source_template_id=family,published_family_id=observed.published_family_id,candidates=int(observed.candidates),reference_matches=int(observed.matched_reference),eligible_reference=int(observed.eligible_reference),exposure_hours=observed.exposure_s/3600,mean_control_candidates=float(controls.candidates.mean()),max_control_candidates=int(controls.candidates.max()),n_controls=len(controls)))
        grouped=pd.read_csv(root/'within_family_review_candidates.csv')
        if 'withheld_cc' in grouped:
            for (family,control),table in grouped.groupby(['published_family_id','control']):
                held_rows.append(dict(panel=label,published_family_id=family,control=int(control),candidate_rows=len(table),median_withheld_cc=float(table.withheld_cc.median()),median_withheld_score_mad=float(table.withheld_score_mad.median()),positive_withheld_fraction=float(table.withheld_cc.gt(0).mean())))
        for f in ['summary.csv','within_family_review_candidates.csv','design.json','report.json']:sources[str(root/f)]=sha256(root/f)
    frame=pd.DataFrame(rows);frame.to_csv(out/'per_template_comparison.csv',index=False)
    pd.DataFrame(held_rows).to_csv(out/'withheld_station_summary.csv',index=False)
    fig,axes=plt.subplots(2,1,figsize=(11,7),layout='constrained')
    for ax,(label,table) in zip(axes,frame.groupby('panel',sort=False)):
        x=np.arange(len(table));ax.bar(x-.22,table.candidates,.22,label='Observed candidate rows')
        ax.bar(x,table.reference_matches,.22,label='Published reference matches')
        ax.bar(x+.22,table.mean_control_candidates,.22,label='Mean shifted-control rows')
        ax.set_xticks(x,table.source_template_id,rotation=20,ha='right');ax.set_ylabel('Rows per template')
        ax.set_title(label.replace('_',' ')+' — 8 MAD');ax.legend(fontsize=8)
    fig.suptitle('Southern development and next-day checks: not unique-event or precision estimates')
    fig.savefig(out/'comparison.png',dpi=160);plt.close(fig)
    write_json(out/'provenance.json',dict(status='complete_bounded_validation_summary',source_sha256=sha256(__file__),inputs_sha256=sources,full_margin_method_validated=False,phase_accuracy_measured=False,cross_family_duplicates_resolved=False))
    (out/'README.md').write_text('''# Southern bounded validation summary

These tables compare one development day with a next-day panel whose 8-MAD threshold and withheld station were specified before examining its event records. No configuration is approved for full-margin production.

Seven template groups represent five published Plourde families. The comparison figure counts rows per template; do not sum them as unique LFEs. A training-template-only clock diagnostic supports provisional grouping within family 026. Cross-family duplicates remain unresolved, and template clocks are not physical origins. Published locations are borrowed family assignments, not new hypocentral solutions.

All observed/control comparisons use common valid exposure, including complete-template exclusion of response taper and trace boundaries. The extra withheld-station mask uses data availability only. PB.B039 is excluded from next-day detection scores; its correlation is evaluated afterward. Control ratios and withheld correlations are diagnostics, not independently reviewed precision. An unmatched candidate may be a new event, an ordinary earthquake, or an artifact.

Raw recordings, response metadata, published templates, source hashes, per-station correlations and control shifts are retained in the linked input directories. No full-history production scan has started. Next steps are independent candidate checks, broader regional validation and monthly runner integration.
''')
    print(frame.to_string(index=False));print(pd.DataFrame(held_rows).to_string(index=False))

if __name__=='__main__':main()
