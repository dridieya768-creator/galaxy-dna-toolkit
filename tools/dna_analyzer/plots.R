#!/usr/bin/env Rscript
# Tracés ggplot2 : profil GC, k-mers les plus fréquents, ORFs par cadre.
# Usage : Rscript plots.R gc.tsv kmer.tsv orf.tsv out_gc.png out_kmer.png out_orf.png
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) stop("Usage: plots.R gc kmer orf out_gc out_kmer out_orf")

gc   <- read.delim(args[1], stringsAsFactors = FALSE, na.strings = "NA")
kmer <- read.delim(args[2], stringsAsFactors = FALSE)
orf  <- read.delim(args[3], stringsAsFactors = FALSE)

# 1. Profil GC (fenêtre glissante)
gc$mid <- (gc$start + gc$end) / 2
p1 <- ggplot(gc, aes(x = mid, y = gc_percent)) +
  geom_line(colour = "#2c7fb8") +
  facet_wrap(~ seq_id, scales = "free_x", ncol = 1) +
  labs(x = "Position (pb)", y = "GC (%)", title = "GC par fenêtre glissante") +
  theme_minimal()
ggsave(args[4], p1, width = 8, height = max(3, 2.5 * length(unique(gc$seq_id))),
       dpi = 120, limitsize = FALSE)

# 2. Top 20 k-mers (par séquence)
top <- do.call(rbind, lapply(split(kmer, kmer$seq_id), function(d) head(d[order(-d$count), ], 20)))
p2 <- ggplot(top, aes(x = reorder(kmer, count), y = count)) +
  geom_col(fill = "#41ab5d") +
  coord_flip() +
  facet_wrap(~ seq_id, scales = "free") +
  labs(x = "k-mer", y = "Occurrences", title = "k-mers les plus fréquents") +
  theme_minimal()
ggsave(args[5], p2, width = 8, height = 6, dpi = 120)

# 3. ORFs : segments par cadre (message si aucun ORF)
if (nrow(orf) > 0) {
  p3 <- ggplot(orf, aes(x = start, xend = end, y = frame, yend = frame)) +
    geom_segment(linewidth = 4, colour = "#d95f0e") +
    facet_wrap(~ seq_id, scales = "free_x", ncol = 1) +
    labs(x = "Position (pb)", y = "Cadre de lecture", title = "ORFs détectés") +
    theme_minimal()
  h <- max(3, 2.5 * length(unique(orf$seq_id)))
} else {
  p3 <- ggplot() + annotate("text", x = 0, y = 0, label = "Aucun ORF détecté") +
    theme_void()
  h <- 3
}
ggsave(args[6], p3, width = 8, height = h, dpi = 120, limitsize = FALSE)
