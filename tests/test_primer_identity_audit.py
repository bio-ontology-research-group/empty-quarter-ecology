import gzip

from analysis.v3.primer_identity_audit import audit_file


def test_primer_audit_distinguishes_reverse_primers(tmp_path):
    path = tmp_path / "reads_R2_001.fastq.gz"
    sequences = [
        "GACTACAGGGGTATCTAATCCAAA",
        "GACTACCVGGGTATCTAATCCAAA".replace("V", "G"),
        "GGACTACAGGGGTATCTAATAAA",
        "CCTACGGGAGGCTGCAGAAA",
    ]
    with gzip.open(path, "wt", encoding="ascii") as handle:
        for index, sequence in enumerate(sequences):
            handle.write(f"@read{index}\n{sequence}\n+\n")
            handle.write("I" * len(sequence) + "\n")

    result = audit_file(path, limit=10)
    assert result["n_reads"] == 4
    assert result["matches_bakt_785r"] == 2
    assert result["matches_apprill_806rb"] == 1
    assert result["matches_bakt_341f"] == 1
