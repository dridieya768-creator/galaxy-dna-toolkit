"""Logique pure (sans dépendance externe) pour l'analyse de séquences ADN.

Fonctions : GC par fenêtre glissante, comptage de k-mers, ORFs sur 6 cadres.
"""
from collections import Counter
from typing import Dict, Iterator, List, Optional

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
START_CODON = "ATG"
STOP_CODONS = {"TAA", "TAG", "TGA"}


def clean_sequence(seq: str) -> str:
    """Met en majuscules et supprime espaces / retours à la ligne."""
    return "".join(seq.split()).upper()


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def gc_sliding_window(seq: str, window: int, step: int) -> List[Dict]:
    """GC (%) sur fenêtre glissante.

    Le dénominateur ne compte que les bases A/C/G/T de la fenêtre ;
    une fenêtre sans base valide donne gc_percent=None.
    Positions 1-based, bornes incluses.
    """
    if window <= 0 or step <= 0:
        raise ValueError("window et step doivent être > 0")
    seq = clean_sequence(seq)
    n = len(seq)
    # sommes préfixes : O(n) au lieu de O(n * window)
    gc_pref = [0] * (n + 1)
    acgt_pref = [0] * (n + 1)
    for i, base in enumerate(seq):
        gc_pref[i + 1] = gc_pref[i] + (base in "GC")
        acgt_pref[i + 1] = acgt_pref[i] + (base in "ACGT")
    rows = []
    last_start = max(n - window, 0)
    for start in range(0, last_start + 1, step):
        end = min(start + window, n)
        valid = acgt_pref[end] - acgt_pref[start]
        gc = None if valid == 0 else 100.0 * (gc_pref[end] - gc_pref[start]) / valid
        rows.append({"start": start + 1, "end": end, "gc_percent": gc})
    return rows


def gc_content(seq: str) -> Optional[float]:
    seq = clean_sequence(seq)
    valid = sum(seq.count(b) for b in "ACGT")
    if valid == 0:
        return None
    return 100.0 * (seq.count("G") + seq.count("C")) / valid


def kmer_counts(seq: str, k: int) -> Counter:
    """Compte les k-mers ; ignore ceux contenant autre chose que A/C/G/T."""
    if k <= 0:
        raise ValueError("k doit être > 0")
    seq = clean_sequence(seq)
    counts: Counter = Counter()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if set(kmer) <= set("ACGT"):
            counts[kmer] += 1
    return counts


def _scan_frame(seq: str, frame: int, min_len_nt: int) -> Iterator[tuple]:
    """Génère (début0, fin0_exclu) des ORFs ATG...STOP dans un cadre (0,1,2)."""
    start = None
    for i in range(frame, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        if start is None:
            if codon == START_CODON:
                start = i
        elif codon in STOP_CODONS:
            end = i + 3
            if end - start >= min_len_nt:
                yield start, end
            start = None


def find_orfs(seq: str, min_len_nt: int = 100) -> List[Dict]:
    """ORFs (ATG -> codon stop inclus) sur les 6 cadres de lecture.

    Cadres : +1,+2,+3 (brin direct) et -1,-2,-3 (brin complémentaire inverse).
    Coordonnées 1-based sur la séquence d'origine (start <= end).
    """
    seq = clean_sequence(seq)
    length = len(seq)
    orfs: List[Dict] = []
    for strand, s in (("+", seq), ("-", reverse_complement(seq))):
        for frame in range(3):
            for a, b in _scan_frame(s, frame, min_len_nt):
                if strand == "+":
                    start, end = a + 1, b
                else:
                    start, end = length - b + 1, length - a
                orfs.append({
                    "frame": f"{strand}{frame + 1}",
                    "start": start,
                    "end": end,
                    "length_nt": b - a,
                    "length_aa": (b - a) // 3 - 1,
                    "sequence": s[a:b],
                })
    orfs.sort(key=lambda o: (o["start"], o["frame"]))
    return orfs
