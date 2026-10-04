"""Reproduce the reported T4 correlations from the compact frozen inputs."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
root = Path(__file__).resolve().parent
before = pd.read_csv(root / 'before.tsv', sep='\t', index_col=0)
removed = pd.read_csv(root / 'removed_asv_counts.tsv', sep='\t', index_col=0).loc[before.index]
assert before.shape == (95, 3) and removed.shape == (95, 7)
n = before.depth.to_numpy()
h = before.shannon.to_numpy()
c = removed.to_numpy()
remaining = n - c.sum(axis=1)
removed_clogc = (c * np.log(np.maximum(c, 1))).sum(axis=1)
after_h = np.log(remaining) - (n * (np.log(n) - h) - removed_clogc) / remaining
after_richness = before.richness_raw.to_numpy() - (c > 0).sum(axis=1)
shannon_rho = float(spearmanr(h, after_h).statistic)
richness_rho = float(spearmanr(before.richness_raw, after_richness).statistic)
expected = json.loads((root / 'summary.json').read_text())
assert abs(shannon_rho - expected['shannon_spearman']) < 1e-12
assert abs(richness_rho - expected['raw_richness_spearman']) < 1e-12
print(json.dumps({'n': 95, 'shannon_spearman': shannon_rho, 'raw_richness_spearman': richness_rho}, indent=2))
