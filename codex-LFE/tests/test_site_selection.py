import pandas as pd
from codex_lfe.site_selection import representative_sites

def test_colocated_instruments_do_not_multiply_site_support():
    frame=pd.DataFrame([dict(id='PB.B032',latitude=43.668,longitude=-123.3923,score=.99),dict(id='PB.B031',latitude=43.6643,longitude=-123.3967,score=.01),dict(id='UO.IRYS',latitude=43.8856,longitude=-123.257,score=.5)])
    selected,groups=representative_sites(frame)
    assert selected.id.tolist()==['PB.B031','UO.IRYS']
    assert groups['PB.B032']==groups['PB.B031']
    assert groups['UO.IRYS']!=groups['PB.B031']
