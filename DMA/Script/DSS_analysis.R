#required_packages



BACW_43<-fread("C:/Users/alvarez.166446/Desktop/DMA_ADvsCT/coverage.dss.gz/BACW_43_1_val_1_bismark_bt2_pe.deduplicated.bismark.dss.gz")
BACW_44<-fread("C:/Users/alvarez.166446/Desktop/DMA_ADvsCT/coverage.dss.gz/BACW_44_1_val_1_bismark_bt2_pe.deduplicated.bismark.dss.gz")
BACW_46<-fread("C:/Users/alvarez.166446/Desktop/DMA_ADvsCT/coverage.dss.gz/BACW_46_1_val_1_bismark_bt2_pe.deduplicated.bismark.dss.gz")
BACW_48<-fread("C:/Users/alvarez.166446/Desktop/DMA_ADvsCT/coverage.dss.gz/BACW_48_1_val_1_bismark_bt2_pe.deduplicated.bismark.dss.gz")
BACW_49<-fread("C:/Users/alvarez.166446/Desktop/DMA_ADvsCT/coverage.dss.gz/BACW_49_1_val_1_bismark_bt2_pe.deduplicated.bismark.dss.gz")
BACW_51<-fread("C:/Users/alvarez.166446/Desktop/DMA_ADvsCT/coverage.dss.gz/BACW_51_1_val_1_bismark_bt2_pe.deduplicated.bismark.dss.gz")
BACW_52<-fread("C:/Users/alvarez.166446/Desktop/DMA_ADvsCT/coverage.dss.gz/BACW_52_1_val_1_bismark_bt2_pe.deduplicated.bismark.dss.gz")
BACW_54<-fread("C:/Users/alvarez.166446/Desktop/DMA_ADvsCT/coverage.dss.gz/BACW_54_1_val_1_bismark_bt2_pe.deduplicated.bismark.dss.gz")
BACW_55<-fread("C:/Users/alvarez.166446/Desktop/DMA_ADvsCT/coverage.dss.gz/BACW_55_1_val_1_bismark_bt2_pe.deduplicated.bismark.dss.gz")
BACW_57<-fread("C:/Users/alvarez.166446/Desktop/DMA_ADvsCT/coverage.dss.gz/BACW_57_1_val_1_bismark_bt2_pe.deduplicated.bismark.dss.gz")
BACW_59<-fread("C:/Users/alvarez.166446/Desktop/DMA_ADvsCT/coverage.dss.gz/BACW_59_1_val_1_bismark_bt2_pe.deduplicated.bismark.dss.gz")
BACW_60<-fread("C:/Users/alvarez.166446/Desktop/DMA_ADvsCT/coverage.dss.gz/BACW_60_1_val_1_bismark_bt2_pe.deduplicated.bismark.dss.gz")
BACW_61<-fread("C:/Users/alvarez.166446/Desktop/DMA_ADvsCT/coverage.dss.gz/BACW_61_1_val_1_bismark_bt2_pe.deduplicated.bismark.dss.gz")
BACW_62<-fread("C:/Users/alvarez.166446/Desktop/DMA_ADvsCT/coverage.dss.gz/BACW_62_1_val_1_bismark_bt2_pe.deduplicated.bismark.dss.gz")
BACW_63<-fread("C:/Users/alvarez.166446/Desktop/DMA_ADvsCT/coverage.dss.gz/BACW_63_1_val_1_bismark_bt2_pe.deduplicated.bismark.dss.gz")
BACW_64<-fread("C:/Users/alvarez.166446/Desktop/DMA_ADvsCT/coverage.dss.gz/BACW_64_1_val_1_bismark_bt2_pe.deduplicated.bismark.dss.gz")

colnames(BACW_43)<-c("chr",     "pos",     "N",       "X")
colnames(BACW_44)<-c("chr",     "pos",     "N",       "X")
colnames(BACW_46)<-c("chr",     "pos",     "N",       "X")
colnames(BACW_48)<-c("chr",     "pos",     "N",       "X")

colnames(BACW_49)<-c("chr",     "pos",     "N",       "X")
colnames(BACW_51)<-c("chr",     "pos",     "N",       "X")
colnames(BACW_52)<-c("chr",     "pos",     "N",       "X")
colnames(BACW_54)<-c("chr",     "pos",     "N",       "X")

colnames(BACW_55)<-c("chr",     "pos",     "N",       "X")
colnames(BACW_57)<-c("chr",     "pos",     "N",       "X")
colnames(BACW_59)<-c("chr",     "pos",     "N",       "X")
colnames(BACW_60)<-c("chr",     "pos",     "N",       "X")

colnames(BACW_61)<-c("chr",     "pos",     "N",       "X")
colnames(BACW_62)<-c("chr",     "pos",     "N",       "X")
colnames(BACW_63)<-c("chr",     "pos",     "N",       "X")
colnames(BACW_64)<-c("chr",     "pos",     "N",       "X")

# Lista de cromosomas válidos (1 a 22
valid_chr <- paste0("chr", c(1:22))

# Filtrar los datos de cada objeto cargado
BACW_43 <- BACW_43[BACW_43$chr %in% valid_chr, ]
BACW_44 <- BACW_44[BACW_44$chr %in% valid_chr, ]
BACW_46 <- BACW_46[BACW_46$chr %in% valid_chr, ]
BACW_48 <- BACW_48[BACW_48$chr %in% valid_chr, ]
BACW_49 <- BACW_49[BACW_49$chr %in% valid_chr, ]
BACW_51 <- BACW_51[BACW_51$chr %in% valid_chr, ]
BACW_52 <- BACW_52[BACW_52$chr %in% valid_chr, ]
BACW_54 <- BACW_54[BACW_54$chr %in% valid_chr, ]
BACW_55 <- BACW_55[BACW_55$chr %in% valid_chr, ]
BACW_57 <- BACW_57[BACW_57$chr %in% valid_chr, ]
BACW_59 <- BACW_59[BACW_59$chr %in% valid_chr, ]
BACW_60 <- BACW_60[BACW_60$chr %in% valid_chr, ]
BACW_61 <- BACW_61[BACW_61$chr %in% valid_chr, ]
BACW_62 <- BACW_62[BACW_62$chr %in% valid_chr, ]
BACW_63 <- BACW_63[BACW_63$chr %in% valid_chr, ]
BACW_64 <- BACW_64[BACW_64$chr %in% valid_chr, ]

# Ahora los datos filtrados pueden ser usados para crear objetos de DSS o continuar con el análisis

BSobj = makeBSseqData( list(BACW_43, 
                            BACW_44, 
                            BACW_46, 
                            BACW_48, 
                            BACW_49, 
                            BACW_51, 
                            BACW_52, 
                            BACW_54, 
                            BACW_55, 
                            BACW_57, 
                            BACW_59, 
                            BACW_60, 
                            BACW_61, 
                            BACW_62, 
                            BACW_63, 
                            BACW_64), 
                       c("BACW_43_AD",
                         "BACW_44_CT",
                         "BACW_46_AD",
                         "BACW_48_CT",
                         "BACW_49_AD",
                         "BACW_51_AD",
                         "BACW_52_CT",
                         "BACW_54_AD",
                         "BACW_55_CT",
                         "BACW_57_CT",
                         "BACW_59_CT",
                         "BACW_60_AD",
                         "BACW_61_CT",
                         "BACW_62_AD",
                         "BACW_63_AD",
                         "BACW_64_CT"))

head(BSobj)

saveRDS(BSobj, "./BSobj.rds")
BSobj <- readRDS("./BSobj.rds")

samps <- sampleNames(BSobj) 

group <- sub(".*_", "", samps)   # extrae AD o CT

table(group)
cov <- getCoverage(BSobj, type = "Cov")

ad_samples <- which(group == "AD")
ct_samples <- which(group == "CT")

# Keep CpGs with coverage >=10 in at least 4 AD and 4 CT samples
keep <- rowSums(cov[, ad_samples, drop = FALSE] >= 10) >= 4 &
        rowSums(cov[, ct_samples, drop = FALSE] >= 10) >= 4

# Filter BSobj
BSobj_filt <- BSobj[keep, ]

table(is.na(BSobj_filt))

saveRDS( BSobj_filt, "./BSobj_filtered.rds")

# 2 Quality Control 
library(bsseq)
library(ggplot2)
library(reshape2)

 ##############Global Methylation################  

beta_mean <- colMeans(beta_values, na.rm = TRUE)

df_mean <- data.frame(
  Sample = samps,
  Group = group,
  GlobalMeth = beta_mean
)

by(df_mean$GlobalMeth, df_mean$Group, summary)
test <- wilcox.test(GlobalMeth ~ Group, data = df_mean)
pval <- signif(test$p.value)

ggplot(df_mean, aes(x = Group, y = GlobalMeth, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_jitter(width = 0.05, size = 2) +
  scale_fill_manual(values = c("AD" = "palevioletred4",
                               "CT" = "bisque2")) +
  theme_minimal() +
  stat_compare_means(method = "wilcox.test", 
  comparisons = list(c("AD", "CT")),
  label.y = 0.74) +
  labs(title = "Global Methylation",
      y = expression("Average Methylation (" * beta * ")"),
       x = "") +
theme(legend.position = "none")

ggsave("results/Global_methylation_violin.pdf", width = 6, height = 5, units = "in" )

##############PCA#################

# beta_sample: matriz2B CpG x muestra (submuestreada)
BSobj_filt<- readRDS("BSobj_filtered.rds")
beta_values_pca <- bsseq::getMeth(BSobj_filt, type = "raw")

beta_values_pca<-na.omit(beta_values_pca) 

summary(beta_values_pca)

# Transformación a M-values
eps <- 1e-3  # más robusto que 1e-6 para evitar Inf numéricos
M <- log2((beta_values_pca + eps) / (1 - beta_values_pca + eps))

#  PCA
pca <- prcomp(t(M), center = TRUE, scale. = FALSE)
# Variance
var_prcomp <- pca$sdev^2
# Plot PCs variance
pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), pc=c(1:length(var_prcomp)))
head(round(pcvar, 3))
ggplot(pcvar[1:15,], aes(x = pc)) + geom_line(aes(y=var)) + 
  labs(x="Principal Component", y="Explained Variance") + 
  geom_point(aes(y=var)) + 
  ggtitle("Principal Components variance")

ggsave("Principal Components variance.pdf", width = 6, height = 5, units = "in" )

############Filtrar SNPS##############

) Cargar SNPs filtrados
snps <- fread("snp151_filtered_MAF10.txt", header = FALSE)
setnames(snps, c("chr","start","end","rsID","refUCSC","alleles","alleleFreqs","MAF"))

# 2) SNPs como puntos en start (1-based)
snp_gr <- GRanges(
  seqnames = snps$chr,
  ranges   = IRanges(start = snps$start, end = snps$start)
)

# 3) CpGs del BSobj (locus CpG)
cpg_gr <- rowRanges(BSobj_filt)

# 4) Overlap (solo misma coordenada) y filtrado
hits <- findOverlaps(cpg_gr, snp_gr, ignore.strand = TRUE)
idx_remove <- unique(queryHits(hits))
BSobj_noSNP <- BSobj_filt[-idx_remove, ]

# 5) Resumen y guardado
cat("CpGs antes:", nrow(BSobj_filt), "\n")
cat("CpGs eliminados (SNP MAF>=0.10 en el locus):", length(idx_remove), "\n")
cat("CpGs después:", nrow(BSobj_noSNP), "\n")

saveRDS(BSobj_noSNP, "BSobj_noSNP_MAF10.rds")

BSobj_noSNP_MAF10 <-readRDS("BSobj_noSNP_MAF10.rds")

############Analysisi diferencial de metilación###############
sampleNames(BSobj_noSNP_MAF10)
group1<-c("BACW_43_AD", "BACW_46_AD","BACW_49_AD", "BACW_51_AD","BACW_54_AD", "BACW_60_AD","BACW_62_AD", "BACW_63_AD")
group2<-c("BACW_44_CT", "BACW_48_CT","BACW_52_CT", "BACW_55_CT","BACW_57_CT", "BACW_59_CT","BACW_61_CT", "BACW_64_CT")

#OBTENGO LOS DML/DMR CON SMOOTH FUNCTION 

DMLtest_BSobj_noSNP_smooth <- DMLtest(BSobj_noSNP_MAF10, group1= group1, group2= group2,
        smoothing=TRUE, smoothing.span=200)

DMLtest_BSobj_noSNP_smooth$fdr <- p.adjust(DMLtest_BSobj_noSNP_smooth$pval, method="BH")

DMLs_BSobj_noSNP_smooth <- callDML(DMLtest_BSobj_noSNP_smooth, delta = 0.1, p.threshold = 0.05)

DML_sig_smooth <- subset(DMLs_BSobj_noSNP_smooth, fdr <= 0.05)

DMRs_BSobj_noSNP_smoth <- callDMR(DMLtest_BSobj_noSNP_smooth,  delta = 0.1,
    p.threshold = 0.05,
    minlen = 50,
    minCG = 3,
    dis.merge = 100,
    pct.sig = 0.5
  )

write.csv2(DML_sig_smooth, "DML_sig_smooth.csv", row.names= FALSE)
write.csv2(DMRs_BSobj_noSNP_smoth, "DMRs_BSobj_noSNP_smoth", row.names = FALSE)
DMLtest_BSobj_noSNP_smooth<- readRDS("DMLTest_BSobj_noSNP_smoth.rds")

########Anotar las DMLs ###########

gr_dml_smooth <- GRanges(
  seqnames = DML_sig_smooth$chr,
  ranges   = IRanges(start = DML_sig_smooth$pos,
                     end   = DML_sig_smooth$pos)
)

gr_dml_smooth
length(gr_dml_smooth)

# si no lo tienes:
BiocManager::install("rGREAT")
library(rGREAT)

job <- submitGreatJob(
  gr = gr_dml_smooth,
  species = "hg38",
  rule = "basalPlusExt",
  adv_upstream = 5,     # kb (basal upstream)
  adv_downstream = 1,   # kb (basal downstream)
  adv_span = 1000       # kb (máxima extensión)
)

assoc <- getRegionGeneAssociations(job)
assoc

top_gene <- vapply(seq_along(assoc$annotated_genes), function(i) {
  g <- as.character(assoc$annotated_genes[[i]])
  d <- as.numeric(assoc$dist_to_TSS[[i]])
  if (length(g) == 0) return(NA_character_)
  g[which.min(abs(d))]
}, character(1))

annot <- data.frame(
  chr = as.character(seqnames(assoc)),
  pos = start(assoc),
  top_gene = top_gene
)
head(annot)

DMLs_sig_smooth_annotadas <- merge(
  DML_sig_smooth,
  annot,
  by = c("chr", "pos"),
  all.x = TRUE
)
write.csv2(DMLs_sig_smooth_annotadas,"DMLs_sig_smooth_annotadas.csv", row.names = FALSE)
