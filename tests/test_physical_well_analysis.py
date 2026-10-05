"""Physical plate-position analysis, distinct from the library-order test."""
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / 'analysis/well-adjacency-20261005'
META = ROOT / 'analysis/batch-adjacency-2026-09-07/batch_meta.tsv'


def test_plate_inputs_and_current_coordinates_are_bound_to_the_replay():
    provenance = json.loads((ANALYSIS/'outputs/input_provenance.json').read_text())
    for key, path in [('positions_sha256', ANALYSIS/'well_positions.csv'),
                      ('metadata_sha256', META)]:
        assert hashlib.sha256(path.read_bytes()).hexdigest() == provenance[key]
    assert provenance['features_sha256'] == '129f47d8f0db8d9afd6f8c67b8d80b5bead90155d8c91d2f43ea6a4139b0cb12'
    wells = pd.read_csv(ANALYSIS/'outputs/wells.tsv', sep='\t')
    metadata = pd.read_csv(META, sep='\t').set_index('profile')
    assert len(wells) == wells.sid.nunique() == 288
    assert not wells.duplicated(['Plate', 'Position']).any()
    for column in ['lat', 'lon']:
        np.testing.assert_allclose(wells[column], metadata.loc[wells.sid, column], atol=1e-12, rtol=0)
    assert provenance['source_position_rows'] == 344
    assert provenance['unique_source_profiles'] == 298
    assert provenance['repeated_profile_rows_removed'] == 46
    assert provenance['unmapped_profiles'] == 10


def test_full_position_permutation_replay_and_manuscript_rounding(tmp_path):
    for name in ['wells.tsv', 'bc_matrix.tsv']:
        shutil.copy2(ANALYSIS/'outputs'/name, tmp_path/name)
    subprocess.run([sys.executable, str(ANALYSIS/'02_neighbor_permtest.py'),
                    '--metadata', str(META), '--output', str(tmp_path)],
                   check=True, capture_output=True, text=True)
    expected = json.loads((ANALYSIS/'outputs/results.json').read_text())
    actual = json.loads((tmp_path/'results.json').read_text())
    assert actual == expected
    for name in ['null_distribution.tsv', 'strata_summary.tsv']:
        assert (tmp_path/name).read_bytes() == (ANALYSIS/'outputs'/name).read_bytes()
    primary = actual['results'][0]
    assert primary['trip'] == 1 and primary['neighbours'] == 'ortho+diag'
    assert primary['samples'] == 259 and primary['adjacent_pairs'] == 441
    assert primary['strata'] == 20 and primary['permutations'] == 999
    assert f"{primary['delta']:+.3f}" == '+0.005'
    assert f"{primary['p_one_sided']:.2f}" == '0.87'
    supplement = (ROOT/'empty-quarter-amplicon/supplement.tex').read_text()
    assert 'same or different soil compartments' in supplement
    assert 'first listed position' in supplement
    assert 'random seed 42' in supplement
