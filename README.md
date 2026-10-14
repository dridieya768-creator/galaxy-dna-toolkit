# Outil Galaxy d'analyse de séquences ADN

Outil [Galaxy](https://galaxyproject.org/) d'analyse de séquences ADN, développé en Python (Biopython) et R (ggplot2), intégré via un wrapper XML et testé avec Planemo.

## Fonctionnalités

| Analyse | Détail |
|---|---|
| **GC par fenêtre glissante** | taille de fenêtre et pas configurables, bases ambiguës (N) ignorées |
| **k-mers** | fréquences pour k = 1 à 8 |
| **ORFs sur 6 cadres** | ATG → codon stop, brins + et −, protéine traduite (Biopython) |
| **Tracés R / ggplot2** | profil GC, top 20 k-mers, carte des ORFs |

## Structure

```
tools/dna_analyzer/
├── dna_analyzer.xml   # wrapper Galaxy
├── dna_analyzer.py    # CLI (Biopython)
├── dna_core.py        # logique pure, testable sans dépendance
├── plots.R            # tracés ggplot2
└── test-data/example.fasta
tests/                 # pytest
```

## Utilisation en ligne de commande

```bash
pip install -r requirements.txt
python tools/dna_analyzer/dna_analyzer.py \
  --input tools/dna_analyzer/test-data/example.fasta \
  --window 100 --step 50 --kmer 3 --min-orf 100 \
  --out-gc gc.tsv --out-kmer kmer.tsv --out-orf orf.tsv

Rscript tools/dna_analyzer/plots.R gc.tsv kmer.tsv orf.tsv gc.png kmer.png orf.png
```

(R et le paquet `ggplot2` sont requis pour les tracés.)

## Tests

```bash
pip install -r requirements-dev.txt
pytest --cov=tools/dna_analyzer          # tests unitaires et CLI
planemo lint tools/dna_analyzer          # validation du wrapper
planemo test --conda_auto_init --conda_auto_install tools/dna_analyzer
```

## Tester sur une instance Galaxy locale

```bash
planemo serve tools/dna_analyzer
```

Ouvrir ensuite http://127.0.0.1:9090, charger un fichier FASTA et lancer l'outil « Analyse de séquences ADN ».

## Format des sorties

* `gc_table` : `seq_id, start, end, gc_percent`
* `kmer_table` : `seq_id, kmer, count, frequency`
* `orf_table` : `seq_id, frame, start, end, length_nt, length_aa, protein` (positions en base 1)

## Licence

MIT
