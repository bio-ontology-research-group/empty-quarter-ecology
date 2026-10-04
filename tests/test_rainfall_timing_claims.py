"""Keep the product-specific timing limits visible in the reported rain analysis."""
from pathlib import Path
import pandas as pd
import manuscript_text as mt

ROOT = Path(__file__).resolve().parents[1] / 'analysis/v3'


def test_fixed_interval_peaks_and_positive_open_meteo_placebo_are_reported():
    peaks = []
    for directory in ('rain_pulse_response', 'rain_pulse_response_open_meteo'):
        bins = pd.read_csv(ROOT / directory / 'disjoint_lag_bin_sensitivity.tsv', sep='\t')
        richness = bins[bins.endpoint == 'richness_hurlbert_25000']
        peak = richness.loc[richness.estimate_per_mm.idxmax()]
        peaks.append((int(peak.lag_start_complete_days), int(peak.lag_end_complete_days)))
    assert peaks == [(3, 4), (1, 1)]
    future = pd.read_csv(ROOT / 'rain_pulse_response_open_meteo/future_rain_placebo_scan.tsv', sep='\t')
    antecedent = pd.read_csv(ROOT / 'rain_pulse_response_open_meteo/pulse_peak_scan.tsv', sep='\t')
    assert future.classical_t.max() > antecedent.classical_t.max()
    assert round(future.partial_r2.max(), 3) == .043
    results = mt.section('main', 'Three additional checks', r'\begintable')
    assert '3–4 days with NASA POWER and at 1 day with Open-Meteo' in results
    assert 'comparable richness association with Open-Meteo' in results
    supplement = mt.section('supplement', r'\subsubsectionRobustness checks', r'\subsubsectionWet spell attribution')
    assert '0.043' in supplement and 'uncertain timing' in supplement
    for unsupported in ('at 3–4 days in both products', 'negative control produced no comparable signal'):
        assert unsupported not in supplement
