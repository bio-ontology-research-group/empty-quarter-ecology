"""Physical sample labels must not filter a group-linked pH exposure."""
import importlib.util
from pathlib import Path

import pandas as pd

spec = importlib.util.spec_from_file_location("ph_group_linkage", Path(__file__).parents[1] / "analysis/v3/ph_ecology_analysis.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_group_join_keeps_other_field_replicates_and_counts_profiles_once(tmp_path):
    cache = tmp_path / "analysis/v2/review/cache"
    cache.mkdir(parents=True)
    pd.DataFrame({"Trip": [4, 4], "Site": [4, 4], "Type": ["Surface", "Surface"],
                  "shannon": [1., 3.], "depth": [10, 20]}, index=["sequenced_A", "sequenced_B"]).to_csv(cache / "alpha.tsv", sep="\t")
    pd.DataFrame({"sequenced_A": [5, 5], "sequenced_B": [8, 12]}, index=["GenusA", "GenusB"]).to_csv(cache / "genus_counts.tsv", sep="\t")
    accepted = tmp_path / "accepted.tsv"
    pd.DataFrame({"sample_id": ["unused_field_A", "unused_field_B"], "trip": [4, 4],
                  "site": [4, 4], "compartment": ["Surface", "Surface"],
                  "ph_value": [8., 9.], "disposition": [module.ADMITTED] * 2}).to_csv(accepted, sep="\t", index=False)
    joined, groups = module.load_sample_level_join(tmp_path, accepted)
    assert len(joined) == 2
    assert groups.iloc[0]["ph"] == 8.5
    assert groups.iloc[0]["shannon"] == 2.
    assert groups.iloc[0]["sequencing_depth"] == 30
    assert groups.iloc[0]["n_profiles"] == 2
    assert set(joined["physical_specimen_identity"]) == {"unverified"}
    assert set(joined["archived_material_relation"]) == {"other_field_replicate_reported"}
    counts, _, info = module.load_grouped_counts(tmp_path, accepted, groups, 1)
    assert counts.to_numpy().sum() == 30
    assert info["selected_profile_columns"] == 2
