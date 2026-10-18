from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np
from scipy.io import savemat
from codex_lfe.chestler_reference import matlab_nominal_seconds,parse_reference

def test_serial_dates_preserve_fractional_seconds_and_epoch():
    np.testing.assert_allclose(matlab_nominal_seconds([719529,719529.5,719529+1.25/86400]),[0,43200,1.25],atol=1e-5)

def test_structured_mat_retains_duplicates_and_uncertain_clock():
    base=Path(__file__).resolve().parents[1]/'results'/'chestler_import_tests';base.mkdir(exist_ok=True)
    with TemporaryDirectory(dir=base) as directory:
        data=np.empty((2,1),dtype=[('name','O'),('peakTimesMo','O'),('lat','O'),('lon','O')])
        data[0,0]=('A',np.array([[719529.5],[719529.5]]),48.,-123.)
        data[1,0]=('B',np.array([[719530.]]),47.8,-123.2)
        path=Path(directory)/'sample.mat';savemat(path,{'LFEs':data})
        families,events=parse_reference(path)
        assert families.source_family_id.tolist()==['A','B']
        assert len(events)==3 and events.reference_id.is_unique
        assert events.nominal_time_s_unverified.tolist()==[43200.,43200.,86400.]
        assert not families.final_published_selection_verified.any()
        assert events.time_basis.str.endswith('_unverified').all()
