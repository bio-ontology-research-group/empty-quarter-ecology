import importlib.util
from pathlib import Path
import sys
import unittest

import pandas as pd

SCRIPTS = Path(__file__).resolve().parents[1] / "analysis" / "v3"
sys.path.insert(0, str(SCRIPTS))
import paired_alpha_sensitivity as paired


class PairedAlphaTests(unittest.TestCase):
    def test_unmatched_campaign_does_not_enter_matched_difference(self):
        frame = pd.DataFrame([
            (1, 1, "Surface", 10), (1, 1, "Deep", 20),
            (2, 1, "Surface", 100),
            (1, 2, "Surface", 20), (1, 2, "Deep", 25),
        ], columns=["Trip", "Site", "Type", "value"])
        matched, count = paired.site_differences(frame, "value", "Deep", "Surface", True)
        marginal, _ = paired.site_differences(frame, "value", "Deep", "Surface", False)
        self.assertEqual(count, 2)
        self.assertEqual(matched.to_dict(), {1: 10, 2: 5})
        self.assertEqual(marginal.to_dict(), {1: -35, 2: 5})

    def test_balanced_campaigns_agree_and_sites_have_equal_weight(self):
        frame = pd.DataFrame([
            (1, 1, "Surface", 10), (1, 1, "Deep", 20),
            (2, 1, "Surface", 100), (2, 1, "Deep", 120),
            (1, 2, "Surface", 20), (1, 2, "Deep", 25),
        ], columns=["Trip", "Site", "Type", "value"])
        matched, count = paired.site_differences(frame, "value", "Deep", "Surface", True)
        marginal, _ = paired.site_differences(frame, "value", "Deep", "Surface", False)
        self.assertEqual(count, 3)
        pd.testing.assert_series_equal(matched, marginal)
        self.assertEqual(matched.mean(), 10)


if __name__ == "__main__":
    unittest.main()
