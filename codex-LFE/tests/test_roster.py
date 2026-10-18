import pytest
from codex_lfe.roster import read_roster, validate_membership


def test_network_identity_and_outside_station_rejection(tmp_path):
    p=tmp_path/'roster.csv';p.write_text('network,station\nCN,PGC\nUW,LRIV\n')
    assert read_roster(p)=={'CN.PGC','UW.LRIV'}
    assert validate_membership(['CN.PGC'],p)['n_allowed']==2
    with pytest.raises(ValueError,match='outside target roster'):
        validate_membership(['UW.PGC'],p)


@pytest.mark.parametrize('text',[
    'network,station\nCN,PGC\nCN,PGC\n',
    'station\nPGC\n',
    'network,station\nCN,*\n',
    'network,station\n',
])
def test_bad_rosters_rejected(tmp_path,text):
    p=tmp_path/'roster.csv';p.write_text(text)
    with pytest.raises(ValueError): read_roster(p)
