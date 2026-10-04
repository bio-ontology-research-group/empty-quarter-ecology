from pathlib import Path
import sys
import unittest
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis" / "v3"))
from ph_context_reanalysis import site_ph


class PHContextTests(unittest.TestCase):
    def test_groups_have_equal_weight_despite_different_assay_counts(self):
        groups = pd.DataFrame([(1, 1, "Surface", 7., 3), (2, 1, "Surface", 9., 1)],
                              columns=["trip", "site", "compartment", "ph", "n_ph_specimens"])
        self.assertEqual(site_ph(groups).loc[1], 8.)

    def test_duplicate_group_keys_rejected(self):
        groups = pd.DataFrame([(1, 1, "Surface", 7.), (1, 1, "Surface", 9.)],
                              columns=["trip", "site", "compartment", "ph"])
        with self.assertRaisesRegex(ValueError, "Repeated assay-group"):
            site_ph(groups)


if __name__ == "__main__":
    unittest.main()
