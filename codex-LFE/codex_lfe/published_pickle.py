"""Read this published numerical dataset using explicit data-class allowlists."""
import importlib
import pickle

ALLOWED = {
 ('obspy.core.trace','Trace'),('obspy.core.trace','Stats'),
 ('obspy.core.utcdatetime','UTCDateTime'),('obspy.core.util.attribdict','AttribDict'),
 ('obspy.core.inventory.response','Response'),('obspy.core.inventory.response','InstrumentSensitivity'),
 ('obspy.core.inventory.response','PolesZerosResponseStage'),('obspy.core.inventory.response','ResponseStage'),
 ('obspy.core.inventory.response','CoefficientsTypeResponseStage'),('obspy.core.inventory.response','CoefficientWithUncertainties'),
 ('obspy.core.inventory.util','Frequency'),('obspy.core.inventory.util','Azimuth'),('obspy.core.inventory.util','Dip'),
 ('obspy.core.util.obspy_types','ComplexWithUncertainties'),('obspy.core.util.obspy_types','_ComplexUncertainty'),
 ('obspy.core.util.obspy_types','FloatWithUncertaintiesAndUnit'),
 ('numpy.core.multiarray','_reconstruct'),('numpy.core.multiarray','scalar'),('numpy','ndarray'),('numpy','dtype'),
 ('pandas.core.frame','DataFrame'),('pandas.core.internals.managers','BlockManager'),
 ('pandas.core.indexes.base','_new_Index'),('pandas.core.indexes.base','Index'),('pandas.core.indexes.range','RangeIndex'),
 ('obspy.core.inventory.response','FIRResponseStage'),('obspy.core.inventory.response','FilterCoefficient'),
 ('builtins','slice')}

class DataUnpickler(pickle.Unpickler):
    def find_class(self,module,name):
        if (module,name) not in ALLOWED:
            raise pickle.UnpicklingError(f'Unapproved global in numerical data: {module}.{name}')
        if module=='numpy.core.multiarray':module='numpy._core.multiarray'
        return getattr(importlib.import_module(module),name)


def read_published_pickle(path):
    with open(path,'rb') as stream:return DataUnpickler(stream).load()
