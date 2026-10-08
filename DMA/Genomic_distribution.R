#########Genomic Distribution Annotation#########
#definir anotación del genoma 
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#anotar mis DMls usando el objeto GRanges (gr_dml)
peakAnno_all <- annotatePeak(
  gr_dml,
  tssRegion = c(-2000,500),
  TxDb = txdb,
  level = "gene",
  assignGenomicAnnotation = TRUE,
  genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                "Downstream", "Intergenic"),
  annoDb = "org.Hs.eg.db", 
  addFlankGeneInfo = FALSE,
  flankDistance = 5000,
  sameStrand = FALSE,
  ignoreOverlap = FALSE,
  ignoreUpstream = FALSE,
  ignoreDownstream = FALSE,
  overlap = "all", #Pprobar con TSS 
  verbose = TRUE,
  columns = c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME")
)

anno_df_all <- as.data.frame(peakAnno_all)
anno_df_tss <- as.data.frame(peakAnno_tss) #  


# anotar las DMRs

peakAnno_dmr_all <- annotatePeak(
  gr_dmr,
  tssRegion = c(-3000,3000),
  TxDb = txdb,
  level = "gene",
  assignGenomicAnnotation = TRUE,
  genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                "Downstream", "Intergenic"),
  annoDb = "org.Hs.eg.db", 
  addFlankGeneInfo = FALSE,
  flankDistance = 5000,
  sameStrand = FALSE,
  ignoreOverlap = FALSE,
  ignoreUpstream = FALSE,
  ignoreDownstream = FALSE,
  overlap = "all", #Pprobar con TSS y all
  verbose = TRUE,
  columns = c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME")
)
anno_df_dmr_all <- as.data.frame(peakAnno_dmr_all)

#############CpG Context Annotation###########


# 2. Construir anotaciones CpG para hg38
annots <- build_annotations(
  genome = "hg38",
  annotations = c(
    "hg38_cpg_islands",
    "hg38_cpg_shores",
    "hg38_cpg_shelves",
    "hg38_cpg_inter"
  )
)

# 3. Anotar DMLs
dml_annot <- annotate_regions(
  regions = gr_dml,
  annotations = annots,
  ignore.strand = TRUE,
  quiet = FALSE
)

dml_annot_df <- as.data.frame(dml_annot)

# 4. Resumir categorías
dml_annot_df <- dml_annot_df %>%
  mutate(
    meth_group = case_when(
      diff > 0 ~ "Hyper",
      diff < 0 ~ "Hypo",
      TRUE ~ NA_character_
    ),
    cpg_context = case_when(
      str_detect(annot.type, "islands") ~ "Island",
      str_detect(annot.type, "shores") ~ "Shore",
      str_detect(annot.type, "shelves") ~ "Shelf",
      str_detect(annot.type, "inter") ~ "Open Sea",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(meth_group), !is.na(cpg_context))

# 5. Si una DML sale repetida, quedarse con una sola categoría
#    Prioridad: Island > Shore > Shelf > Open Sea
# ========================
dml_annot_df <- dml_annot_df %>%
  mutate(
    cpg_context = factor(
      cpg_context,
      levels = c("Island", "Shore", "Shelf", "Open Sea")
    )
  ) %>%
  arrange(seqnames, start, cpg_context) %>%
  distinct(seqnames, start, end, diff, .keep_all = TRUE)

plot_sum <- dml_annot_df %>%
  count(meth_group, cpg_context) %>%
  group_by(meth_group) %>%
  mutate(perc = 100 * n / sum(n)) %>%
  ungroup()

order_regions <- c("Island", "Shore", "Shelf", "Open Sea")

plot_sum <- expand.grid(
  meth_group = c("Hypo", "Hyper"),
  cpg_context = order_regions,
  stringsAsFactors = FALSE
) %>%
  left_join(plot_sum, by = c("meth_group", "cpg_context")) %>%
  mutate(
    n = ifelse(is.na(n), 0, n),
    perc = ifelse(is.na(perc), 0, perc),
    cpg_context = factor(cpg_context, levels = rev(order_regions)),
    perc_plot = ifelse(meth_group == "Hypo", -perc, perc),
    label = sprintf("%.2f%%", perc)
  )


ggplot(plot_sum, aes(x = perc_plot, y = cpg_context, fill = meth_group)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_text(
    aes(
      label = label,
      hjust = ifelse(meth_group == "Hypo", 1.08, -0.08)
    ),
    size = 3.8
  ) +
  scale_x_continuous(
    labels = function(x) paste0(abs(x), "%"),
    expand = expansion(mult = c(0.18, 0.18))
  ) +
  scale_fill_manual(
    values = c("Hyper" = "pink3", "Hypo" = "thistle"),
    name = NULL
  ) +
  labs(
    x = "Percentage of DMLs",
    y = "CpG region",
    title = "CpG context distribution of hypo- and hypermethylated DMLs"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom"
  )

ggsave("results/genomic_region_annotatr_DMLs.pdf", width = 7, height = 3.8, dpi = 300, units = "in" ) 

###################Enrichment DMLs #####################

library(GenomicFeatures)
library(GenomicRanges)
library(ChIPseeker)
library(dplyr)
library(stringr)
library(ggplot2)

# 1. SIMPLIFICAR ANOTACIONES OBSERVADAS
simplify_annot <- function(x) {
  case_when(
    str_detect(x, "^Promoter") ~ "Promoter",
    str_detect(x, "5' UTR") ~ "5'UTR",
    str_detect(x, "3' UTR") ~ "3'UTR",
    str_detect(x, "^Exon") ~ "Exon",
    str_detect(x, "^Intron") ~ "Intron",
    str_detect(x, "^Downstream")       ~ "Distal Intergenic",
    str_detect(x, "Distal Intergenic") ~ "Distal Intergenic"
  )
}

# observed: DMLs ya anotadas
obs_df <- dml_annot_df %>%
  mutate(region = simplify_annot(annotation)) %>%
  count(region, name = "observed_n") %>%
  mutate(observed_pct = observed_n / sum(observed_n))

# 2. CONSTRUIR REGIONES ESPERADAS DESDE TxDb

# Promotores: mismo criterio que ChIPseeker
prom_gr <- promoters(genes(txdb), upstream = 2000, downstream = 500)

# 5'UTR y 3'UTR
utr5_gr <- unlist(fiveUTRsByTranscript(txdb), use.names = FALSE)
utr3_gr <- unlist(threeUTRsByTranscript(txdb), use.names = FALSE)

# Exones e intrones
exon_gr <- exons(txdb)
intron_gr <- unlist(intronsByTranscript(txdb), use.names = FALSE)

# Downstream: aproximación simétrica respecto al gen
gene_gr <- genes(txdb)

# Reducir para evitar contar solapamientos internos repetidos
prom_gr <- GenomicRanges::reduce(prom_gr)
utr5_gr <- GenomicRanges::reduce(utr5_gr)
utr3_gr <- GenomicRanges::reduce(utr3_gr)
exon_gr <- GenomicRanges::reduce(exon_gr)
intron_gr <-GenomicRanges::reduce(intron_gr)

# 3. DEFINIR ESPACIO GENÓMICO TOTA
# sacar longitudes cromosómicas desde seqinfo del txdb
si <- seqinfo(txdb)
valid_seqs <- names(seqlengths(si))[!is.na(seqlengths(si))]
valid_seqs <- paste0("chr", 1:22)
genome_gr <- GRanges(
  seqnames = valid_seqs,
  ranges = IRanges(start = 1, end = seqlengths(si)[valid_seqs])
)

# 4. HACER PARTICIONES JERÁRQUICAS para no contar dos veces la misma base
prom_final <- prom_gr
utr5_final <- setdiff(utr5_gr, prom_final)
utr3_final <- setdiff(utr3_gr, reduce(c(prom_final, utr5_final)))
exon_final <- setdiff(exon_gr, reduce(c(prom_final, utr5_final, utr3_final)))
intron_final <- setdiff(intron_gr, reduce(c(prom_final, utr5_final, utr3_final, exon_final)))

used_gr <- reduce(c(prom_final, utr5_final, utr3_final, exon_final, intron_final))
intergenic_final <- setdiff(genome_gr, used_gr)
intergenic_final <- reduce(intergenic_final)

# 5. CALCULAR LONGITUDES ESPERADAS
expected_df <- data.frame(
  region = c("Promoter", "5'UTR", "3'UTR", "Exon", "Intron","Distal Intergenic"),
  feature_length = c(
    sum(width(prom_final)),
    sum(width(utr5_final)),
    sum(width(utr3_final)),
    sum(width(exon_final)),
    sum(width(intron_final)),
    sum(width(intergenic_final))
  )
) %>%
  mutate(expected_pct = feature_length / sum(feature_length))

# 6. UNIR OBSERVED + EXPECTED
order_regions <- c("Promoter", "5'UTR", "3'UTR", "Exon",
                   "Intron", "Distal Intergenic")

fold_df <- full_join(obs_df, expected_df, by = "region") %>%
  mutate(
    observed_n = ifelse(is.na(observed_n), 0, observed_n),
    observed_pct = ifelse(is.na(observed_pct), 0, observed_pct),
    feature_length = ifelse(is.na(feature_length), 0, feature_length),
    expected_pct = ifelse(is.na(expected_pct), 0, expected_pct)
  ) %>%
  filter(region != "Other", expected_pct > 0) %>%
  mutate(
    fold_enrichment = observed_pct / expected_pct,
    log2_FE = log2(fold_enrichment),
    region = factor(region, levels = order_regions)
  ) %>%
  arrange(log2_FE)

# 7. CLEVELAND PLOT
palette_regions <- c(
  "Promoter" = "#4E79A7",
  "5'UTR" = "#A0CBE8",
  "3'UTR" = "#F28E2B",
  "Exon" = "#E15759",
  "Intron" = "#76B7B2",
  "Distal Intergenic" = "#EDC948"
)

ggplot(fold_df, aes(x = log2_FE, y = region, color = region)) +
  geom_point(size = 5, shape = 19) +
  scale_y_discrete(limits = rev(levels(fold_df$region)), position = "right") +
  scale_colour_manual(values = palette_regions) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.2) +
  labs(
    x = "log2(Fold enrichment observed / expected)",
    y = "",
    title = "Genomic region enrichment of DMLs"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 9), 
    axis.text.y = element_text(size = 9),
    axis.title.x.bottom = element_text(size= 9)
  )

ggsave("genomic_region_enrichment_DML.pdf", width = 7, height = 3.8, dpi = 300,

####### Enrichemnt DMRs #############
# observed: DMRs ya anotadas
obs_df_dmr <- anno_df_dmr_all %>%
  mutate(region = simplify_annot(annotation)) %>%
  count(region, name = "observed_n") %>%
  mutate(observed_pct = observed_n / sum(observed_n))


# 2. CONSTRUIR REGIONES ESPERADAS DESDE TxDb
# Promotores: mismo criterio que ChIPseeker
prom_gr <- promoters(genes(txdb), upstream = 2000, downstream = 500)

# 5'UTR y 3'UTR
utr5_gr <- unlist(fiveUTRsByTranscript(txdb), use.names = FALSE)
utr3_gr <- unlist(threeUTRsByTranscript(txdb), use.names = FALSE)

# Exones e intrones
exon_gr <- exons(txdb)
intron_gr <- unlist(intronsByTranscript(txdb), use.names = FALSE)

# Downstream: aproximación simétrica respecto al gen
gene_gr <- genes(txdb)

# Reducir para evitar contar solapamientos internos repetidos
prom_gr <- reduce(prom_gr)
utr5_gr <- reduce(utr5_gr)
utr3_gr <- reduce(utr3_gr)
exon_gr <- reduce(exon_gr)
intron_gr <- reduce(intron_gr)

# 3. DEFINIR ESPACIO GENÓMICO TOTAL
# sacar longitudes cromosómicas desde seqinfo del txdb
si <- seqinfo(txdb)
valid_seqs <- names(seqlengths(si))[!is.na(seqlengths(si))]
genome_gr <- GRanges(
  seqnames = valid_seqs,
  ranges = IRanges(start = 1, end = seqlengths(si)[valid_seqs])
)

# 4. HACER PARTICIONES JERÁRQUICAS  para no contar dos veces la misma base

prom_final <- prom_gr

utr5_final <- setdiff(utr5_gr, prom_final)
utr3_final <- setdiff(utr3_gr, reduce(c(prom_final, utr5_final)))
exon_final <- setdiff(exon_gr, reduce(c(prom_final, utr5_final, utr3_final)))
intron_final <- setdiff(intron_gr, reduce(c(prom_final, utr5_final, utr3_final, exon_final)))

used_gr <- reduce(c(prom_final, utr5_final, utr3_final, exon_final, intron_final))
intergenic_final <- setdiff(genome_gr, used_gr)
intergenic_final <- reduce(intergenic_final)

# 5. CALCULAR LONGITUDES ESPERADAS
expected_df_dmr <- data.frame(
  region = c("Promoter", "5'UTR", "3'UTR", "Exon", "Intron","Distal Intergenic"),
  feature_length = c(
    sum(width(prom_final)),
    sum(width(utr5_final)),
    sum(width(utr3_final)),
    sum(width(exon_final)),
    sum(width(intron_final)),
    sum(width(intergenic_final))
  )
) %>%
  mutate(expected_pct = feature_length / sum(feature_length))

# 6. UNIR OBSERVED + EXPECTED
order_regions <- c("Promoter", "5'UTR", "3'UTR", "Exon",
                   "Intron", "Distal Intergenic")

fold_df_dmr <- full_join(obs_df_dmr, expected_df_dmr, by = "region") %>%
  mutate(
    observed_n = ifelse(is.na(observed_n), 0, observed_n),
    observed_pct = ifelse(is.na(observed_pct), 0, observed_pct),
    feature_length = ifelse(is.na(feature_length), 0, feature_length),
    expected_pct = ifelse(is.na(expected_pct), 0, expected_pct)
  ) %>%
  filter(region != "Other", expected_pct > 0) %>%
  mutate(
    fold_enrichment = observed_pct / expected_pct,
    log2_FE = log2(fold_enrichment),
    region = factor(region, levels = order_regions)
  ) %>%
  arrange(log2_FE)

# 7. CLEVELAND PLOT
palette_regions <- c(
  "Promoter" = "",
  "5'UTR" = "#6D597A",
  "3'UTR" = "#F28E2B",
  "Exon" = "#E15759",
  "Intron" = "#76B7B2",
  "Distal Intergenic" = "#EDC948"
)

# 8 Plot
ggplot(fold_df_dmr, aes(x = log2_FE, y = region, color = region)) +
  geom_point(size = 5, shape = 19) +
  scale_y_discrete(limits = rev(levels(fold_df$region)), position = "right") +
  scale_colour_manual(values = cols) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.2) +
  labs(
    x = "log2(Fold enrichment observed / expected)",
    y = "",
    title = "Genomic region enrichment of DMRs"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 9), 
    axis.text.y = element_text(size = 9),
    axis.title.x.bottom = element_text(size= 9)
  )

ggsave("genomic_region_enrichment_DMR.pdf", width = 7, height = 3.8, dpi = 300, units = "in" ) 
