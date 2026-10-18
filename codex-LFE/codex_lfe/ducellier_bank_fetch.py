"""Fetch a bounded published comparison bank on the existing seven-station panel."""
from pathlib import Path
import json,urllib.request
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
from .provenance import sha256,write_json

def main():
    root=Path('results/southern_benchmark_inventory_v1');tree=json.loads((root/'github_tree.json').read_text())
    overlap=pd.read_csv(root/'family_roster_overlap.csv').fillna({'allowed_triplets':''})
    allowed=set(pd.read_csv(root/'targeted_reference_packet_plan.csv').station)
    overlap=overlap[(overlap.n_allowed_triplets>=5)&overlap.allowed_triplets.map(lambda s:set(s.split(',')).issubset(allowed))]
    paths=[]
    for row in overlap.itertuples():
        stations={s.split('.')[1] for s in row.allowed_triplets.split(',')}
        for item in tree['tree']:
            p=item['path']
            if p.startswith('src/templates_v1/'+row.family_id+'/') and p.endswith('.pkl') and Path(p).stem.split('_')[0] in stations:paths.append(p)
        paths.append('data/Ducellier/catalogs/'+row.family_id+'/catalog_2004_2011.pkl')
    def fetch(p):
        path=root/'source'/p
        if not path.exists():
            url='https://raw.githubusercontent.com/ArianeDucellier/catalog/'+tree['sha']+'/'+p
            with urllib.request.urlopen(url,timeout=30) as response:data=response.read(2_000_001)
            if len(data)>2_000_000:raise ValueError('Source exceeds 2 MB cap')
            path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
        return dict(path=p,sha256=sha256(path),bytes=path.stat().st_size)
    write_json(root/'comparison_bank_design.json',dict(commit=tree['sha'],families=overlap.family_id.tolist(),files=len(paths),maximum_bytes=len(paths)*2_000_000,selection='At least five roster triplets, subset of the existing seven-station packet panel; no new waveform stations'))
    with ThreadPoolExecutor(max_workers=3) as pool:records=list(pool.map(fetch,paths))
    write_json(root/'comparison_bank_manifest.json',dict(commit=tree['sha'],records=records,total_bytes=sum(r['bytes'] for r in records)))
    print('Comparison bank complete:',len(overlap),'families,',len(records),'files,',sum(r['bytes'] for r in records),'bytes',flush=True)

if __name__=='__main__':main()
