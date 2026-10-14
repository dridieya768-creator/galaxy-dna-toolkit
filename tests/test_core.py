import pytest

import dna_core as core


def test_reverse_complement():
    assert core.reverse_complement("ATGC") == "GCAT"
    assert core.reverse_complement("AAANCC") == "GGNTTT"


def test_gc_content():
    assert core.gc_content("GGCC") == 100.0
    assert core.gc_content("ATAT") == 0.0
    assert core.gc_content("ATGC") == 50.0
    assert core.gc_content("NNNN") is None


def test_gc_sliding_window_basic():
    rows = core.gc_sliding_window("GGGGAAAA", window=4, step=4)
    assert [(r["start"], r["end"]) for r in rows] == [(1, 4), (5, 8)]
    assert rows[0]["gc_percent"] == 100.0
    assert rows[1]["gc_percent"] == 0.0


def test_gc_sliding_window_overlap_and_short_sequence():
    rows = core.gc_sliding_window("ATGCATGC", window=4, step=2)
    assert [(r["start"], r["end"]) for r in rows] == [(1, 4), (3, 6), (5, 8)]
    # séquence plus courte que la fenêtre : une seule fenêtre tronquée
    rows = core.gc_sliding_window("GC", window=10, step=5)
    assert len(rows) == 1 and rows[0]["gc_percent"] == 100.0


def test_gc_ignores_n_and_all_n_window():
    rows = core.gc_sliding_window("GCNN" "NNNN", window=4, step=4)
    assert rows[0]["gc_percent"] == 100.0
    assert rows[1]["gc_percent"] is None


def test_gc_invalid_parameters():
    with pytest.raises(ValueError):
        core.gc_sliding_window("ACGT", window=0, step=1)
    with pytest.raises(ValueError):
        core.gc_sliding_window("ACGT", window=4, step=0)


def test_kmer_counts():
    counts = core.kmer_counts("ATATAT", 2)
    assert counts == {"AT": 3, "TA": 2}
    assert core.kmer_counts("ACNGT", 2) == {"AC": 1, "GT": 1}  # N ignoré
    assert core.kmer_counts("AC", 3) == {}
    with pytest.raises(ValueError):
        core.kmer_counts("ACGT", 0)


def test_orf_forward_strand():
    orfs = core.find_orfs("ATGAAATAG", min_len_nt=9)
    assert len(orfs) == 1
    o = orfs[0]
    assert (o["frame"], o["start"], o["end"]) == ("+1", 1, 9)
    assert o["length_nt"] == 9 and o["length_aa"] == 2


def test_orf_offset_frames():
    assert core.find_orfs("CATGAAATAG", 9)[0]["frame"] == "+2"
    assert core.find_orfs("CCATGAAATAG", 9)[0]["frame"] == "+3"


def test_orf_min_length_filter():
    assert core.find_orfs("ATGAAATAG", min_len_nt=12) == []


def test_orf_reverse_strand_coordinates():
    fwd = "CCATGAAATAGGG"
    rc = core.reverse_complement(fwd)
    orfs = core.find_orfs(rc, 9)
    # même ORF, vu sur le brin opposé, mêmes coordonnées que sur fwd
    reverse = [o for o in orfs if o["frame"].startswith("-")]
    assert len(reverse) == 1
    o = reverse[0]
    assert o["sequence"] == "ATGAAATAG"
    assert (o["start"], o["end"]) == (3, 11)
    assert rc[o["start"] - 1:o["end"]] == core.reverse_complement("ATGAAATAG")


def test_orf_no_stop_not_reported():
    assert core.find_orfs("ATGAAAGCT", 3) == []
