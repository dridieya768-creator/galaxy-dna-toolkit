#!/usr/bin/env python3
"""Analyse de séquences ADN : GC par fenêtre glissante, k-mers, ORFs (6 cadres).

Produit trois tableaux TSV (consommés ensuite par plots.R).
"""
import argparse
import csv
import sys

from Bio import SeqIO
from Bio.Seq import Seq

import dna_core as core


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="Fichier FASTA")
    p.add_argument("--window", type=int, default=100, help="Taille de fenêtre GC")
    p.add_argument("--step", type=int, default=50, help="Pas de la fenêtre GC")
    p.add_argument("--kmer", type=int, default=3, help="Taille des k-mers")
    p.add_argument("--min-orf", type=int, default=100,
                   help="Longueur minimale d'un ORF (nucléotides, codon stop inclus)")
    p.add_argument("--out-gc", required=True)
    p.add_argument("--out-kmer", required=True)
    p.add_argument("--out-orf", required=True)
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    records = list(SeqIO.parse(args.input, "fasta"))
    if not records:
        print("Erreur : aucune séquence FASTA valide trouvée.", file=sys.stderr)
        return 1

    with open(args.out_gc, "w", newline="") as f_gc, \
         open(args.out_kmer, "w", newline="") as f_k, \
         open(args.out_orf, "w", newline="") as f_o:
        w_gc = csv.writer(f_gc, delimiter="\t")
        w_k = csv.writer(f_k, delimiter="\t")
        w_o = csv.writer(f_o, delimiter="\t")
        w_gc.writerow(["seq_id", "start", "end", "gc_percent"])
        w_k.writerow(["seq_id", "kmer", "count", "frequency"])
        w_o.writerow(["seq_id", "frame", "start", "end", "length_nt",
                      "length_aa", "protein"])

        for rec in records:
            seq = str(rec.seq)
            for r in core.gc_sliding_window(seq, args.window, args.step):
                gc = "NA" if r["gc_percent"] is None else f"{r['gc_percent']:.2f}"
                w_gc.writerow([rec.id, r["start"], r["end"], gc])

            counts = core.kmer_counts(seq, args.kmer)
            total = sum(counts.values())
            for kmer, c in sorted(counts.items(), key=lambda x: (-x[1], x[0])):
                w_k.writerow([rec.id, kmer, c, f"{c / total:.6f}"])

            for o in core.find_orfs(seq, args.min_orf):
                protein = str(Seq(o["sequence"]).translate(to_stop=True))
                w_o.writerow([rec.id, o["frame"], o["start"], o["end"],
                              o["length_nt"], o["length_aa"], protein])
    return 0


if __name__ == "__main__":
    sys.exit(main())
