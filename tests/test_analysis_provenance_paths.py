import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VERDICTS = (
    "analysis/v3/compartment_composition/claim_verdict.json",
    "analysis/v3/depth_extraction/claim_verdict.json",
    "analysis/v3/distance_decay_turnover/claim_verdict.json",
    "analysis/v3/evenness_decomposition/claim_verdict.json",
    "analysis/v3/geographic_prediction/claim_verdict.json",
    "analysis/v3/pma_endpoint_results/pma_summary.json",
    "analysis/v3/xrf_community_clr/claim_verdict.json",
)


def test_canonical_verdicts_use_portable_input_paths():
    for relative in VERDICTS:
        path = ROOT / relative
        payload = json.loads(path.read_text(encoding="utf-8"))
        serialized = json.dumps(payload)
        assert "/home/" not in serialized, relative
        assert "/ibex/" not in serialized, relative
        assert "/mnt/" not in serialized, relative
        assert "/datawaha/" not in serialized, relative
