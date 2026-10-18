from pathlib import Path
import json
from tempfile import TemporaryDirectory
import pandas as pd
import pytest
from codex_lfe.catalog_products import verify_unit
from codex_lfe.provenance import write_json

@pytest.fixture
def products():
    base=Path(__file__).resolve().parents[1]/'results'/'product_verification_tests';base.mkdir(exist_ok=True)
    with TemporaryDirectory(dir=base) as name:
        root=Path(name);paths={k:root/(k+'.parquet') for k in ['events','families','coverage','controls']}
        pd.DataFrame([dict(event_id='event1',time_s=5.,time_basis='anchor',family_id='f1',quality_status='experimental_unverified',supporting_stations='XX.A',support_windows=json.dumps([dict(station='XX.A',start_s=4.,end_s=8.)]))]).to_parquet(paths['events'])
        pd.DataFrame([dict(family_id='f1')]).to_parquet(paths['families'])
        pd.DataFrame([dict(station='XX.A',start_s=0.,end_s=10.,status='processed',evidence_id='input1'),dict(station='XX.B',start_s=0.,end_s=10.,status='provider_no_data',evidence_id='input2')]).to_parquet(paths['coverage'])
        pd.DataFrame([dict(realization=1,count=0,exposure_s=10.)]).to_parquet(paths['controls'])
        paths['provenance']=root/'provenance.json'
        write_json(paths['provenance'],dict(time_basis='anchor',coverage_evidence={k:dict(description='test receipt',source_sha256='a'*64) for k in ['input1','input2']}))
        yield paths

def test_integrity_does_not_certify_science(products):
    report=verify_unit(products,0.,10.,{'XX.A','XX.B'})
    assert report['passed'] and not report['scientific_validation_performed'] and not report['production_validated']

@pytest.mark.parametrize('failure',['missing_station','gap','overlap','failed_request','missing_evidence','unsupported_event','end_boundary','unknown_family','support_window_gap'])
def test_invalid_products_block_cleanup_verification(products,failure):
    c=pd.read_parquet(products['coverage']);e=pd.read_parquet(products['events'])
    if failure=='missing_station':c=c.iloc[:1]
    elif failure=='gap':c.loc[0,'start_s']=1.
    elif failure=='overlap':c=pd.concat([c,c.iloc[:1]])
    elif failure=='failed_request':c.loc[1,'status']='request_failed'
    elif failure=='missing_evidence':c.loc[1,'evidence_id']='absent'
    elif failure=='unsupported_event':e.loc[0,'supporting_stations']='XX.B'
    elif failure=='end_boundary':e.loc[0,'time_s']=10.
    elif failure=='unknown_family':e.loc[0,'family_id']='unknown'
    elif failure=='support_window_gap':
        c=pd.concat([c,pd.DataFrame([dict(station='XX.A',start_s=6.,end_s=10.,status='waveform_gap',evidence_id='input1')])]);c.loc[0,'end_s']=6.
    c.to_parquet(products['coverage']);e.to_parquet(products['events'])
    with pytest.raises(ValueError):verify_unit(products,0.,10.,{'XX.A','XX.B'})

def test_verified_lifecycle_retains_products_and_deletes_only_owned_copy(products):
    import shutil
    from codex_lfe.monthly_storage import MonthlyStorage
    root=products['events'].parent
    s=MonthlyStorage(root/'monthly',100);s.begin_month('2026-09')
    original=root/'original_waveform';original.write_bytes(b'waveform')
    staged=s.stage('unit1','raw',original.read_bytes())
    retained=s.root/'retained';retained.mkdir()
    paths={k:retained/p.name for k,p in products.items()}
    for k,p in products.items():shutil.copyfile(p,paths[k])
    paths['verification']=retained/'verification.json'
    report=s.commit_verified_unit('unit1',paths,0.,10.,{'XX.A','XX.B'})
    assert report['passed'] and not report['production_validated']
    s=MonthlyStorage(s.root,100);s.cleanup('unit1');s.finish_month(['unit1'])
    assert original.read_bytes()==b'waveform' and not staged.exists()
    assert all(p.exists() for p in paths.values())


def add_padding(products):
    products['support_coverage']=products['coverage'].with_name('support_coverage.parquet')
    pd.DataFrame([dict(station='XX.A',start_s=-2.,end_s=0.,status='processed',evidence_id='input1'),
                  dict(station='XX.A',start_s=10.,end_s=12.,status='processed',evidence_id='input1')]).to_parquet(products['support_coverage'])
    p=json.loads(products['provenance'].read_text());p['support_coverage_bounds']=dict(start_s=-2.,end_s=12.)
    write_json(products['provenance'],p)
    e=pd.read_parquet(products['events']);e.loc[0,'support_windows']=json.dumps([dict(station='XX.A',start_s=-2.,end_s=12.)])
    e.to_parquet(products['events'])


def test_padding_supports_boundary_event_without_changing_event_core(products):
    add_padding(products)
    report=verify_unit(products,0.,10.,{'XX.A','XX.B'})
    assert report['n_padding_coverage_intervals']==2
    assert 'support_coverage' in report['products']
    assert not report['production_validated']


@pytest.mark.parametrize('failure',['gap','core_override','unknown_station','overlap','missing_bounds','missing_evidence','outside_bounds'])
def test_padding_cannot_hide_missing_or_conflicting_support(products,failure):
    add_padding(products);c=pd.read_parquet(products['support_coverage'])
    if failure=='gap':c.loc[0,'end_s']=-1.
    elif failure=='core_override':c.loc[0,'end_s']=1.
    elif failure=='unknown_station':c.loc[0,'station']='XX.C'
    elif failure=='overlap':c=pd.concat([c,c.iloc[:1]])
    elif failure=='missing_bounds':
        p=json.loads(products['provenance'].read_text());p.pop('support_coverage_bounds');write_json(products['provenance'],p)
    elif failure=='missing_evidence':c.loc[0,'evidence_id']='absent'
    elif failure=='outside_bounds':c.loc[0,'start_s']=-3.
    c.to_parquet(products['support_coverage'])
    with pytest.raises(ValueError):verify_unit(products,0.,10.,{'XX.A','XX.B'})
