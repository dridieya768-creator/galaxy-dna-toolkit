"""Test d'intégration du script CLI (nécessite Biopython)."""
from pathlib import Path

import pytest

pytest.importorskip("Bio")

import dna_analyzer  # noqa: E402

FASTA = Path(__file__).resolve().parents[1] / "tools/dna_analyzer/test-data/example.fasta"


def run(tmp_path, **extra):
    out = {k: tmp_path / f"{k}.tsv" for k in ("gc", "kmer", "orf")}
    argv = ["--input", str(FASTA), "--out-gc", str(out["gc"]),
            "--out-kmer", str(out["kmer"]), "--out-orf", str(out["orf"])]
    for k, v in extra.items():
        argv += [f"--{k}", str(v)]
    return dna_analyzer.main(argv), out


def test_cli_outputs(tmp_path):
    code, out = run(tmp_path)
    assert code == 0
    gc_lines = out["gc"].read_text().splitlines()
    assert gc_lines[0] == "seq_id\tstart\tend\tgc_percent"
    assert len(gc_lines) == 15  # en-tête + 7 fenêtres x 2 séquences
    orf_text = out["orf"].read_text()
    assert "seq1_test\t+1\t151\t276\t126" in orf_text
    assert "MAAAAAAA" in orf_text  # protéine traduite (M + 20 A + 20 K)


def test_cli_empty_fasta(tmp_path):
    empty = tmp_path / "empty.fasta"
    empty.write_text("")
    argv = ["--input", str(empty), "--out-gc", str(tmp_path / "a"),
            "--out-kmer", str(tmp_path / "b"), "--out-orf", str(tmp_path / "c")]
    assert dna_analyzer.main(argv) == 1
