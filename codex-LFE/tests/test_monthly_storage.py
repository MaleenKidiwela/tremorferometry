from pathlib import Path
from tempfile import TemporaryDirectory
import pytest
from codex_lfe.monthly_storage import MonthlyStorage,REQUIRED
from codex_lfe.provenance import write_json

@pytest.fixture
def storage():
    base=Path(__file__).resolve().parents[1]/'results'/'storage_lifecycle_tests'
    base.mkdir(exist_ok=True)
    with TemporaryDirectory(dir=base) as directory:
        yield MonthlyStorage(directory,10)

def artifacts(s):
    folder=s.root/'retained';folder.mkdir()
    result={}
    for kind in REQUIRED:
        p=folder/(kind+'.json');write_json(p,dict(passed=True));result[kind]=p
    return result

def test_budget_and_active_month(storage):
    s=storage;s.begin_month('2026-09');p=s.stage('day14','raw',b'123456')
    with pytest.raises(ValueError,match='budget'):s.stage('day14','extra',b'12345')
    with pytest.raises(ValueError,match='unfinished'):s.begin_month('2026-08')
    with pytest.raises(ValueError,match='not committed'):s.cleanup('day14')
    assert p.read_bytes()==b'123456'

def test_changed_retained_data_preserves_raw(storage):
    s=storage;s.begin_month('2026-09');p=s.stage('day14','raw',b'123')
    a=artifacts(s);s.commit('day14',a);a['events'].write_text('changed')
    with pytest.raises(ValueError,match='output changed'):s.cleanup('day14')
    assert p.exists()

def test_changed_raw_preserved(storage):
    s=storage;s.begin_month('2026-09');p=s.stage('day14','raw',b'123')
    s.commit('day14',artifacts(s));p.write_bytes(b'other')
    with pytest.raises(ValueError,match='Staged data changed'):s.cleanup('day14')
    assert p.read_bytes()==b'other'

def test_verified_cleanup_restart_and_reverse_month(storage):
    s=storage;s.begin_month('2026-09');p=s.stage('day14','raw',b'123')
    a=artifacts(s);s.commit('day14',a)
    s=MonthlyStorage(s.root,10);s.cleanup('day14')
    assert not p.exists() and all(p.exists() for p in a.values())
    with pytest.raises(ValueError,match='coverage incomplete'):s.finish_month(['day14','day13'])
    s.finish_month(['day14'])
    with pytest.raises(ValueError,match='backward'):s.begin_month('2026-10')
    s.begin_month('2026-08')

def test_external_symlink_is_never_deleted(storage):
    s=storage;s.begin_month('2026-09');p=s.stage('day14','raw',b'123')
    s.commit('day14',artifacts(s))
    external=s.root/'preserved';external.write_bytes(b'123');p.unlink();p.symlink_to(external)
    with pytest.raises(ValueError,match='owned'):s.cleanup('day14')
    assert external.exists()

def test_empty_unit_requires_verified_products_and_survives_restart(storage):
    s=storage;s.begin_month('2010-01');s.begin_unit('day31')
    with pytest.raises(ValueError,match='not committed'):s.cleanup('day31')
    s=MonthlyStorage(s.root,10)
    with pytest.raises(ValueError,match='Missing retained'):s.commit('day31',{})
    a=artifacts(s);s.commit('day31',a);s.cleanup('day31');s.finish_month(['day31'])
    assert s.state['units']['2010-01/day31']['files']=={}
    assert all(p.exists() for p in a.values())
