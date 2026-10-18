"""Compare final catalogs on shared eligible reference IDs, without rescanning."""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import pandas as pd
from codex_lfe.evaluate import match_events

ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--baseline',required=True)
ap.add_argument('--variant',required=True)
ap.add_argument('--output',required=True)
ap.add_argument('--separate-training-anchors',action='store_true',help='Score each run against its own training-derived anchor; verify published source times agree')
args=ap.parse_args()
a,b=Path(args.baseline),Path(args.variant)
ra=pd.read_parquet(a/'calibration_eligible_reference.parquet')
rb=pd.read_parquet(b/'calibration_eligible_reference.parquet')
common=set(ra.reference_id)&set(rb.reference_id)
r=ra[ra.reference_id.isin(common)].sort_values('reference_id').reset_index(drop=True)
other=rb[rb.reference_id.isin(common)].sort_values('reference_id').reset_index(drop=True)
if args.separate_training_anchors:
 pd.testing.assert_frame_equal(r[['reference_id','family_id']],other[['reference_id','family_id']])
 anchors_a=json.loads((a/'training.json').read_text())['anchors_from_published_s']
 anchors_b=json.loads((b/'training.json').read_text())['anchors_from_published_s']
 import numpy as np
 np.testing.assert_allclose(r.time_s-r.family_id.map(anchors_a),other.time_s-other.family_id.map(anchors_b),rtol=0,atol=1e-6)
else:
 pd.testing.assert_frame_equal(r[['reference_id','time_s','time_basis','family_id']],other[['reference_id','time_s','time_basis','family_id']])
cfg=json.loads((a/'manifest.json').read_text())['config']
rows=[]
for label,path in [('baseline',a),('variant',b)]:
 target=other if label=='variant' and args.separate_training_anchors else r
 frames=[pd.read_parquet(p) for p in sorted(path.glob('calibration_????-??-??_candidates.parquet'))]
 c=pd.concat(frames,ignore_index=True)
 for threshold in cfg['thresholds']:
  group=c[(c.realization==0)&(c.threshold==threshold)]
  matched=match_events(group,target,tolerance_s=cfg['reference_tolerance_s'],same_family=True)
  rows.append({'variant':label,'threshold':threshold,'n_shared_eligible':len(r),'n_matched':len(matched),'recall':len(matched)/len(r) if len(r) else None})
pd.DataFrame(rows).to_csv(args.output,index=False)
