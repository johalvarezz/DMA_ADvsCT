
BiocManager::install(c("rGREAT", "GenomicRanges", "IRanges", "org.Hs.eg.db", "AnnotationDbi"))

library(rGREAT)
library(GenomicRanges)
library(IRanges)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Anotación para todas las DMLs (3384) DMRs (96): obtenemos GO terms (BP,cc,,MF), KEGG (C2:CP:KEGG)

#Preparar SOLO DMLs significativas como GRanges. Tener en cuenta que en DSS, DML suele tener chr y pos. Una DML es un punto (1 bp).

# DML_sig_smooth es la tabla de DMLs significativas filtradas pot FDR

gr_dml <- GRanges(
  seqnames = DML_sig_smooth$chr,
  ranges   = IRanges(start = DML_sig_smooth$pos, end = DML_sig_smooth$pos),
  strand = "*",
  diff = DML_sig_smooth$diff,
  pval = DML_sig_smooth$pval,
  fdr = DML_sig_smooth$fdr
)
# opcional, pero útil para trazar
names(gr_dml) <- paste0(as.character(seqnames(gr_dml)), ":", start(gr_dml))
#convertimos nuestro DMLtest_BSobj_noSNP_smoth a GRange y será nuestro background

bg <- GRanges(
  seqnames = DMLtest_BSobj_noSNP_smooth$chr,
  ranges   = IRanges(start = DMLtest_BSobj_noSNP_smooth$pos,
                     end   = DMLtest_BSobj_noSNP_smooth$pos)
)
# 4) Correr rGREAT local (hg38) para GO y KEGG 
obj_bp <- great(gr_dml, gene_sets = "GO:BP", tss_source = "hg38", background = bg)
obj_cc <- great(gr_dml, gene_sets = "GO:CC", tss_source = "hg38", background = bg)
obj_mf <- great(gr_dml, gene_sets = "GO:MF", tss_source = "hg38", background = bg)
obj_kegg<-great(gr_dml, gene_sets = "msigdb:C2:CP:KEGG",  tss_source = "hg38", background = bg)

# Sacar las tablas de enriquecimiento con la función getEnrichment() 

tb_bp   <- getEnrichmentTable(obj_bp,   min_region_hits = 5 )
tb_cc   <- getEnrichmentTable(obj_cc,   min_region_hits = 5)
tb_mf   <- getEnrichmentTable(obj_mf,   min_region_hits = 5)
tb_kegg <- getEnrichmentTable(obj_kegg, min_region_hits = 5)

#Filtrar las significativas

sig_bp   <- subset(tb_bp,   p_adjust < 0.05)
sig_cc   <- subset(tb_cc,   p_adjust < 0.05)
sig_mf   <- subset(tb_mf,   p_adjust < 0.05)
sig_kegg <- subset(tb_kegg, p_adjust < 0.05)


write.csv(sig_bp, "go_bp.csv", row.names = FALSE)
write.csv(sig_cc, "go_cc.csv", row.names = FALSE)
write.csv(sig_mf, "go_mf.csv", row.names = FALSE)
write.csv(sig_kegg, "sig_kegg_dml.csv", row.names =FALSE)

# Genes cercanos / dominios 
assoc_bp <- getRegionGeneAssociations(obj_mf, use_symbols = TRUE)
mcols(assoc_bp)
# para adjuntar a mi DMLs originales de DSS, primero convierto a data.frame 
assoc_df <- as.data.frame(assoc_bp)
## como puede que en alguna DML caiga en más de un gen se convierte a string la columna de genes para que uego pueda hacerse el merge sin problemas
assoc_df$annotated_genes <- sapply(assoc_df$annotated_genes, function(x) paste(x, collapse=";"))
assoc_df$dist_to_TSS <- sapply(assoc_df$dist_to_TSS, function(x) paste(x, collapse=";"))
head(assoc_df)
a <- DML_sig_smooth 
a$id   <- paste0(a$chr, ":", a$pos)
assoc_df$id <- paste0(assoc_df$seqnames, ":", assoc_df$start)
final_df <- merge(
  a,
  assoc_df[, c("id","annotated_genes", "dist_to_TSS")],
  by = "id",
  all.x = TRUE
)
write.csv2(final_df,"results/DMLs_sig_smooth_rGREAT.csv")
HACER un DOTPLOT DE LOS 20 GO:BPs SELECCIONADOS: 
library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(forcats)


go_bp <- fread("go_bp.csv")
go_cc <- fread("go_cc.csv")
go_mf <- fread("go_mf.csv")

terms_bp <- c(
  "regulation of fatty acid metabolic process",
  "regulation of lipid storage",
  "production of molecular mediator involved in inflammatory response",
  "oxidative phosphorylation",
  "receptor localization to synapse",
  "endocytosis",
  "glial cell migration",
  "forebrain neuron development",
  "innate immune response",
  "antigen processing and presentation",
  "protein ubiquitination",
  "autophagy",
  "cell adhesion molecule production",
  "cell-cell adhesion mediated by integrin",
  "leukocyte cell-cell adhesion",
  "cellular response to reactive oxygen species",
  "vesicle-mediated transport in synapse",
  "positive regulation of neuron projection development",
  "lipid homeostasis",
  "synaptic vesicle endocytosis",
  "presynaptic endocytosis"
)


terms_cc <- c(
  "recycling endosome membrane",
  "early endosome membrane",
  "lipid droplet",
  "mitochondrial membrane",
  "immunological synapse",
  "synaptic membrane",
  "neuron to neuron synapse",
  "immunological synapse",
  "postsynaptic membrane",
  "protein complex involved in cell-matrix adhesion",
  "dendritic spine"
)

terms_mf <- c(
  "sterol ester esterase activity",
  "ubiquitin conjugating enzyme binding",
  "ligand-modulated transcription factor activity",
  "nuclear receptor activity",
  "ubiquitin ligase activator activity",
  "neurotransmitter transmembrane transporter activity",
  "cell-matrix adhesion mediator activity",
  "lipid binding",
  "cell adhesion mediator activity",
  "protein binding involved in heterotypic cell-cell adhesion",
  "cholesterol transfer activity",
  "sterol transfer activity"
) 

  prepare_go <- function(df, selected_terms, ontology_name) {
  
  df %>%
    filter(description %in% selected_terms) %>%
    mutate(
      minusLog10FDR = -log10(p_adjust),
      DML_count = observed_region_hits,
      description = str_wrap(description, width = 40),
      GeneRatio = observed_gene_hits / gene_set_size,
      ontology = ontology_name
    ) %>%
    arrange(p_adjust) %>%
    mutate(description = fct_reorder(description, minusLog10FDR))
}

bp_plot <- prepare_go(go_bp, terms_bp, "BP")
write.csv2(bp_plot,"GO_BP_figura.csv")
cc_plot <- prepare_go(go_cc, terms_cc, "CC")
write.csv2(cc_plot,"GO_CC_figura.csv")
mf_plot <- prepare_go(go_mf, terms_mf, "MF")
write.csv2(mf_plot,"GO_MF_Dfigura.csv")

bp_plot<- read.csv2("GO_BP_figura.csv")

 p_bp<- ggplot(bp_plot, aes(x = GeneRatio, y = description,
                            size = DML_count, color = minusLog10FDR)) +
  geom_point(alpha = 0.9) +
   scale_color_gradient(low = "skyblue", high = "navy") +
  labs(
    title = "GO Biological Process of DMLs",
    x = "GeneRatio",
    y = NULL,
    size = "DMLs count",
    color = expression("- log10(adjusted p−value)")
  ) +
   guides(
  size = guide_legend(order = 1),
  color = guide_colorbar(order = 2)
)+
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 9),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
p_bp

ggsave("results/GO_BP_DML_gratio.pdf", width = 7, height = 9, dpi = 300 )

p_cc <- ggplot(cc_plot, aes(x = minusLog10FDR, y = description,
                            size = DML_count, color = GeneRatio)) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(low = "orange", high = "orange4") +
  labs(
    title = "GO Cellular Component of DMLs",
    x = expression(-log[10]("adjusted p-value")),
    y = NULL,
    size = "DMLs count",
    color = "Gene Ratio"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 9),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
p_cc

ggsave("results/GO_CC_DML_gratio.pdf", width = 7, height = 5, dpi = 300 )
p_mf <- ggplot(mf_plot, aes(x = minusLog10FDR, y = description,
                            size = DML_count, color = GeneRatio )) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(low = "coral", high = "chocolate4") +
  labs(
    title = "GO Molecular function of DMLs",
    x = expression(-log[10]("adjusted p-value")),
    y = NULL,
    size = "DMLs count",
    color = "Gene Ratio"
  ) + guides(
  size = guide_legend(order = 1),
  color = guide_colorbar(order = 2)
)  +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 9),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
p_mf
ggsave("results/GO_MF_DML_gratio.pdf", width = 7, height = 5, dpi = 300 )

                               
# Filtrar significativos si no lo hiciste
sig_kegg <- tb_kegg %>% filter(p_adjust < 0.05)

terms_kegg <- c(
  "Adipocytokine Signaling Pathway",
  "Ppar Signaling Pathway",
  "Alzheimers Disease",
  "Purine Metabolism",
  "Endocytosis",
  "Peroxisome",
  "Oxidative Phosphorylation"
)


kegg <- kegg %>%
  mutate(
    description = id %>%
      str_remove("^KEGG_") %>%      # quita "KEGG_" del inicio
      str_replace_all("_", " ") %>%
      str_to_lower() %>%            # todo a minúsculas
      str_to_title()                # primera letra de cada palabra en mayúscula
  )

kegg$description

kegg_plot <- kegg %>%
  filter(description %in% terms_kegg) %>%
  mutate(
    minusLog10FDR = -log10(p_adjust),
    DML_count = observed_region_hits,
    GeneRatio = observed_gene_hits / gene_set_size,
    description = str_wrap(description, width = 35)
  ) %>%
  arrange(p_adjust) %>%
  mutate(description = fct_reorder(description, minusLog10FDR))

p_kegg<- ggplot(kegg_plot,
                 aes(x = minusLog10FDR, y = description , size = DML_count, color= minusLog10FDR)) +
  geom_point(alpha = 0.9)+
  scale_color_gradient(low = "mediumpurple1", high="mediumpurple4")  +
  labs(
    title = "KEGG pathway enrichment analysis of DMLs",
    x = expression(-log[10]("adjusted p-value")),
    y = NULL,
    size = "DML count",
    color= "Gene Ratio"
  ) +
  guides(
  size = guide_legend(order = 1),
  color = guide_colorbar(order = 2)
)+
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )
p_kegg

ggsave("results/kegg_DMLs_gratio.pdf", width = 7, height = 3.8, dpi = 300, units = "in")
                               
#########GO terms de las DMRs ############## 

#anotación de las DMRs(96) con rGREAT

BiocManager::install("GO.db")

DMR_smooth <-DMRs_BSobj_noSNP_smoth

gr_dmr <- GRanges(
  seqnames = DMR_smooth$chr,
  ranges   = IRanges(start = DMR_smooth$start, end = DMR_smooth$end),
  strand = "*",
  diff.Methy = DMR_smooth$diff.Methy,
  areaStat= DMR_smooth$areaStat
  
)

# opcional, pero útil para trazar
names(gr_dmr) <- paste0(as.character(seqnames(gr_dmr)), ":", start(gr_dmr))
#convertimos nuestro DMLtest_BSobj_noSNP_smoth a GRange y será nuestro background

bg <- GRanges(
  seqnames = DMLtest_BSobj_noSNP_smoth$chr,
  ranges   = IRanges(start = DMLtest_BSobj_noSNP_smoth$pos,
                     end   = DMLtest_BSobj_noSNP_smoth$pos)
)

# 4) Correr rGREAT local (hg38) para GO y KEGG 
obj_bp_dmr <- great(gr_dmr, gene_sets = "GO:BP", tss_source = "hg38")
obj_cc_dmr <- great(gr_dmr, gene_sets = "GO:CC", tss_source = "hg38")
obj_mf_dmr <- great(gr_dmr, gene_sets = "GO:MF", tss_source = "hg38")
obj_kegg_dmr<-great(gr_dmr, gene_sets = "msigdb:C2:CP:KEGG",  tss_source = "hg38")
obj_react_dmr<-great(gr_dmr, gene_sets = "msigdb:C2:CP:REACTOME", tss_source = "hg38")
obj_inmuno_dmr <- great(gr_dmr, gene_sets = "msigdb:C7:IMMUNESIGDB", tss_source= "hg38") 

# Sacar las tablas de enriquecimiento con la función getEnrichment() 

tb_bp_dmr   <- getEnrichmentTable(obj_bp_dmr,   min_region_hits = 3 )
tb_cc_dmr   <- getEnrichmentTable(obj_cc_dmr,   min_region_hits = 3)
tb_mf_dmr   <- getEnrichmentTable(obj_mf_dmr,   min_region_hits = 3)
tb_kegg_dmr <- getEnrichmentTable(obj_kegg_dmr, min_region_hits = 5)
tb_reactome_dmr <- getEnrichmentTable(obj_react_dmr, min_region_hits = 1)
tb_inmuno_dmr <-  getEnrichmentTable(obj_inmuno_dmr, min_region_hits = 5)

#Filtrar las significativas

sig_bp_dmr   <- subset(tb_bp_dmr,   p_adjust < 0.05)
sig_cc_dmr   <- subset(tb_cc_dmr,   p_adjust < 0.05)
sig_mf_dmr   <- subset(tb_mf_dmr,   p_adjust < 0.05)@no sale nada 
sig_kegg_dmr <- subset(tb_kegg_dmr, p_adjust < 0.05)
sig_react_dmr <- subset(tb_reactome_dmr, p_adjust < 0.05) #no sale nada 
sig_inmuno_dmr <- subset(tb_inmuno_dmr, p_adjust < 0.05) # no sale nada 

write.csv(sig_bp_dmr, "sig_GO_bp_dmr.csv", row.names = FALSE)
bp_dmr <- fread("sig_GO_bp_dmr.csv")
write.csv(sig_cc_dmr, "sig_GO_cc_dmr.csv", row.names = FALSE)
cc_dmr <- fread("sig_GO_cc_dmr.csv")
write.csv2(tb_kegg_dmr,"results/tb_kegg_dmr.csv")

# Asociaciones
assoc_bp_dmr <- getRegionGeneAssociations(obj_cc_dmr, use_symbols = TRUE)

assoc_df_dmr <- as.data.frame(assoc_bp_dmr)

# Convertir listas a texto
assoc_df_dmr$annotated_genes <- sapply(assoc_df_dmr$annotated_genes, paste, collapse=";")
assoc_df_dmr$dist_to_TSS     <- sapply(assoc_df_dmr$dist_to_TSS, paste, collapse=";")

# Crear IDs iguales en ambos dataframes (chr:start-end)

b <- DMR_smooth
b$id <- paste0(b$chr, ":", b$start, "-", b$end)

assoc_df_dmr$id <- paste0(
  as.character(assoc_df_dmr$seqnames), ":",
  assoc_df_dmr$start, "-",
  assoc_df_dmr$end
)

# Merge
final_df_dmr <- merge(
  b,
  assoc_df_dmr[, c("id","annotated_genes","dist_to_TSS")],
  by = "id",
  all.x = TRUE
)

# Comprobación
nrow(final_df_dmr)
nrow(b)

write.csv2(final_df_dmr, "results/DMRs_sig_smooth_anotadas_rGREAT.csv", row.names = FALSE)                               
                               
#KEGG plot
kegg_dmr <- sig_kegg_dmr %>%
  mutate(
    description = id %>%
      str_remove("^KEGG_") %>%      # Quita "KEGG_" al inicio
      str_replace_all("_", " ") %>% # Cambia guiones bajos por espacios
      str_to_lower() %>%            # Todo a minúsculas
      str_to_title()                # Primera letra de cada palabra en mayúscula
  )


kegg_dmr_plot <- kegg_dmr %>%
  mutate(
    DML_count = observed_region_hits,   # <-- cambia si tu columna se llama distinto
    neglog10FDR = -log10(p_adjust),
    GeneRatio = observed_gene_hits / gene_set_size,
    term = str_wrap(description, 45)
  ) %>%
  arrange(p_adjust) %>%
  mutate(term = factor(term, levels = rev(term)))

p_kegg_dmr <-ggplot(kegg_dmr_plot, aes(x = GeneRatio, y =reorder(description,GeneRatio),
                      size = DML_count, color = neglog10FDR)) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(low = "darkorchid1", high = "mediumpurple4") +
  labs(
    title = "Diferential Methylation Region KEGG pathway enrichment",
    x = "Gene Ratio",
    y = NULL,
    size = "DML count",
    color = expression(-log[10]("FDR (p_adjust)"))
  ) +
  guides(
  size = guide_legend(order = 1),
  color = guide_colorbar(order = 2)
)+
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 10)
  )

p_kegg_dmr
ggsave("results/kegg_DMRs.pdf", width = 9, height = 5, dpi = 300, units = "in") 
                               
