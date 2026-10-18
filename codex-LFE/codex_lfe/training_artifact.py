"""Load reusable training-only templates with explicit configuration/content checks."""
import json
from pathlib import Path
import numpy as np
from .provenance import sha256


def load_training_artifact(path, config, reference_sha256):
    root = Path(path)
    record = json.loads((root/'training_artifact.json').read_text())
    if record['config'] != config or record['reference_sha256'] != reference_sha256:
        raise ValueError('training artifact configuration/reference mismatch')
    input_name = 'training_input_manifest.json' if 'training_input_manifest.json' in record['files'] else 'input_manifest.json'
    expected = {'training.json','templates.npz','template_qc.csv',input_name}
    if set(record['files']) != expected:
        raise ValueError('incomplete training artifact')
    for name in expected:
        if sha256(root/name) != record['files'][name]:
            raise ValueError(f'training artifact changed: {name}')
    training = json.loads((root/'training.json').read_text())
    templates = {f:{} for f in training['families_usable']}
    with np.load(root/'templates.npz',allow_pickle=False) as archive:
        for family in templates:
            for station in config['stations']:
                key = f'{family}__{station}'
                if key+'__template' in archive:
                    templates[family][station] = (archive[key+'__template'].copy(),float(archive[key+'__offset_s']))
    return templates, training['anchors_from_published_s'], json.loads((root/input_name).read_text())
