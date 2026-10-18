import numpy as np
import pytest
from codex_lfe.training_artifact import load_training_artifact
from codex_lfe.provenance import sha256, write_json


def test_training_reuse_rejects_changed_inputs(tmp_path):
    cfg={'stations':{'CN.PGC':{}}}
    write_json(tmp_path/'training.json',{'families_usable':['f'],'anchors_from_published_s':{'f':20.}})
    np.savez(tmp_path/'templates.npz',f__CN_PGC=np.array([]),
             **{'f__CN.PGC__template':np.ones((3,160)),'f__CN.PGC__offset_s':np.array(-2.)})
    (tmp_path/'template_qc.csv').write_text('usable\nTrue\n')
    write_json(tmp_path/'input_manifest.json',{})
    names=['training.json','templates.npz','template_qc.csv','input_manifest.json']
    write_json(tmp_path/'training_artifact.json',{'config':cfg,'reference_sha256':'ref',
                                              'files':{n:sha256(tmp_path/n) for n in names}})
    templates,anchors,_=load_training_artifact(tmp_path,cfg,'ref')
    assert templates['f']['CN.PGC'][0].shape==(3,160)
    assert anchors['f']==20.
    with pytest.raises(ValueError,match='mismatch'):
        load_training_artifact(tmp_path,{'stations':{}},'ref')
    (tmp_path/'template_qc.csv').write_text('changed')
    with pytest.raises(ValueError,match='changed'):
        load_training_artifact(tmp_path,cfg,'ref')


def test_sealed_training_manifest_survives_runtime_updates(tmp_path):
    cfg={'stations':{'CN.PGC':{}}}
    write_json(tmp_path/'training.json',{'families_usable':['f'],'anchors_from_published_s':{'f':20.}})
    np.savez(tmp_path/'templates.npz',**{'f__CN.PGC__template':np.ones((3,160)),'f__CN.PGC__offset_s':np.array(-2.)})
    (tmp_path/'template_qc.csv').write_text('usable\nTrue\n')
    write_json(tmp_path/'training_input_manifest.json',{'training':'sealed'})
    names=['training.json','templates.npz','template_qc.csv','training_input_manifest.json']
    write_json(tmp_path/'training_artifact.json',{'config':cfg,'reference_sha256':'ref','files':{n:sha256(tmp_path/n) for n in names}})
    write_json(tmp_path/'input_manifest.json',{'training':'sealed','calibration':'added later'})
    _,_,inputs=load_training_artifact(tmp_path,cfg,'ref')
    assert inputs=={'training':'sealed'}
    write_json(tmp_path/'training_input_manifest.json',{'training':'changed'})
    with pytest.raises(ValueError,match='changed'):
        load_training_artifact(tmp_path,cfg,'ref')
