##############################################################################
## Transcriptomic Analysis of MYC Target Genes as a Proxy for MYC
## Transcriptional Activity
##
## Supplementary R script for the book chapter:
##   Arnal Segura, M.  "Transcriptomic Analysis of MYC Target Genes as a Proxy
##   for MYC Transcriptional Activity". In: The MYC Gene: Methods and
##   Protocols (3rd ed.).
##
## Worked example: bulk RNA-seq from GEO accession GSE309250 (MDA-MB-231 cells
## treated with 30 uM Omomyc, indicated as DP30 in the sample names, or
## vehicle for 120 h; two biological replicates per condition).
##
## Author:  M. Arnal Segura (Peptomyc S.L., Barcelona, Spain)
## Tested:  R 4.4.2
##############################################################################


## ---- 1. Software requirements ----------------------------------------------

library(data.table)            # Fast import of the table of counts
library(dplyr)                 # Pipe operator (%>%) and tidy data manipulation
library(edgeR)                 # TMM normalization of RNA-seq counts
library(limma)                 # DE analysis with voom + lmFit + eBayes
library(msigdbr)               # Access to the MSigDB gene set collections
library(clusterProfiler)       # Gene Set Enrichment Analysis (GSEA)
library(enrichplot)            # dotplot() for gseaResult objects
library(GSEABase)              # GeneSet / GeneSetCollection containers
library(GSVA)                  # Gene Set Variation Analysis (GSVA)
library(SummarizedExperiment)  # assay() accessor for the GSVA result
library(ComplexHeatmap)        # Heatmap visualization
library(grid)                  # gpar() for heatmap graphical parameters
library(fastcluster)           # Fast hierarchical clustering
library(decoupleR)             # CollecTRI regulon retrieval + run_ulm()

# Note: this script was prepared with msigdbr v7.5.1. The msigdbr API was
# changed in v9.0.0 (the `category` argument was renamed to `collection`);
# adapt the msigdbr() calls below if a newer version of the package is used.


## ---- 2. Import of the table of counts (GEO GSE309250) ----------------------

counts_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE309nnn/GSE309250/",
  "suppl/GSE309250%5Fcounts%2Etxt%2Egz"
)

dt  <- fread(counts_url)
mat <- as.matrix(dt, rownames = 1)

# Four samples, two replicates per condition:
#   MDAMB231_120h_DP30_B3, MDAMB231_120h_DP30_B4  -> Omomyc (30 uM)
#   MDAMB231_120h_veh_B3,  MDAMB231_120h_veh_B4   -> vehicle
print(colnames(mat))


## ---- 3. Pre-processing -----------------------------------------------------

dim(mat)  # 26364 genes x 4 samples

# Retain genes with more than 10 counts in all 4 samples
keep     <- rowSums(mat > 10) >= 4
counts_f <- mat[keep, ]
dim(counts_f)  # 13367 genes x 4 samples


## ---- 4. Normalization: TMM and log2-CPM ------------------------------------

d <- DGEList(counts = counts_f)
d <- calcNormFactors(d, method = "TMM")

# log2-CPM matrix used as input for GSVA; the prior count stabilizes the
# variance near zero
cpm_matrix <- cpm(d, log = TRUE, prior.count = 1)


## ---- 5. Differential expression analysis (limma + voom) --------------------

# Define the condition from the column names
cond <- ifelse(grepl("veh", colnames(counts_f)), "vehicle", "omomyc")
table(cond)

# A no-intercept design simplifies the specification of pairwise contrasts
design <- model.matrix(~ 0 + cond)
rownames(design) <- colnames(counts_f)
colnames(design) <- gsub("cond", "", colnames(design))

# Contrast: Omomyc-treated samples versus vehicle-treated samples
cont_matrix <- makeContrasts(OMOvsVEH = omomyc - vehicle, levels = design)

# voom() is called on the matrix of raw counts. Quantile normalization
# is an optional between-array correction and should be evaluated empirically.
voom_res <- voom(counts_f, design, plot = TRUE, normalize = "quantile")

fit      <- lmFit(voom_res, design)
fit_main <- contrasts.fit(fit, cont_matrix)
fit_main <- eBayes(fit_main)

# Full DE table (FDR-adjusted)
top_diff <- topTable(fit_main, coef = "OMOvsVEH", number = Inf, adjust.method = "fdr")
head(top_diff)


## ---- 6. MYC-related gene sets from MSigDB ----------------------------------

# Hallmark collection (H)
H <- msigdbr(species = "Homo sapiens", category = "H")
H_symbol <- as.data.frame(H %>% dplyr::select(gs_name, gene_symbol))
H_symbol$collection <- "Hallmark"
dim(H_symbol)   # 8209 rows x 3 columns

# C2 curated collection
C2 <- msigdbr(species = "Homo sapiens", category = "C2")
C2_symbol <- as.data.frame(C2 %>% dplyr::select(gs_name, gene_symbol))
C2_symbol$collection <- "C2"
dim(C2_symbol)  # ~594903 rows x 3 columns (size depends on the MSigDB version)

# Combine both collections and retain MYC-related gene sets
all_gs  <- rbind(H_symbol, C2_symbol)
myc_gsc <- all_gs[grepl("MYC_|MYCN_|MYCL_", all_gs$gs_name), ]
length(unique(myc_gsc$gs_name))  # 51 MYC-related gene sets


## ---- 7. Gene Set Enrichment Analysis (GSEA) --------------------------------

# Ranking metric: -log10(p) * sign(log2FC)
signL  <- sign(top_diff$logFC)
logP   <- -log10(top_diff$P.Value)
metric <- logP * signL
names(metric) <- rownames(top_diff)
metric <- sort(metric, decreasing = TRUE)

# GSEA relies on a permutation step; fix the random seed for reproducibility
set.seed(2)
gsea_res <- GSEA(
  geneList     = metric,
  TERM2GENE    = myc_gsc[, c("gs_name", "gene_symbol")],
  verbose      = FALSE,
  eps          = 0,
  pvalueCutoff = 0.05
)
gsea_tbl <- gsea_res@result

# HALLMARK_MYC_TARGETS_V1 and DANG_MYC_TARGETS_UP are negatively enriched,
# consistent with MYC inhibition in the Omomyc-vs-vehicle contrast.
print(gsea_tbl[, c("NES", "pvalue", "p.adjust")])

# Dotplot of the enriched gene sets, sorted by NES
print(dotplot(gsea_res, x = "NES"))


## ---- 8. Gene Set Variation Analysis (GSVA) ---------------------------------

# Split the MYC-related gene sets by gs_name
gene_lists <- split(myc_gsc$gene_symbol, myc_gsc$gs_name)

# Build a GeneSetCollection (the format required by GSVA)
gsc <- GeneSetCollection(
  lapply(names(gene_lists), function(name) {
    GeneSet(
      unique(as.character(gene_lists[[name]])),
      setName    = name,
      geneIdType = SymbolIdentifier()  # change if Entrez or Ensembl IDs are used
    )
  })
)

# Run GSVA on the log2-CPM matrix
par_obj  <- gsvaParam(cpm_matrix, gsc, minSize = 5, maxSize = 500)
gsva_res <- gsva(par_obj)

# Extract the score matrix
gsva_m=gsva_res[1:nrow(gsva_res),1:ncol(gsva_res)]
dim(gsva_m)#47  4

# Heatmap of GSVA enrichment scores
fh <- function(x) fastcluster::hclust(dist(x), method = "ward.D2")

draw(
  Heatmap(
    gsva_m,
    name            = " ",
    cluster_rows    = fh,
    cluster_columns = fh,
    row_names_gp    = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 10)
  )
)
# GSVA scores cleanly separate Omomyc-treated from vehicle-treated samples;
# HALLMARK_MYC_TARGETS_V1 is among the most discriminant gene sets.


## ---- 9. Transcription-factor activity with CollecTRI + ULM -----------------

# Retrieve the CollecTRI regulon (signed TF-target interactions)
net <- get_collectri(organism = "human", split_complexes = FALSE)

# Reshape the gene-level metric as a one-column matrix so that the contrast
# name is preserved in the `condition` column of the run_ulm() output.
metric_mat <- matrix(
  metric,
  ncol     = 1,
  dimnames = list(names(metric), "OMOvsVEH")
)

contrast_acts <- run_ulm(
  mat     = metric_mat,
  net     = net,
  .source = "source",
  .target = "target",
  .mor    = "mor",
  minsize = 5
)
contrast_acts <- as.data.frame(contrast_acts)

# MYC activity in the Omomyc-vs-vehicle contrast: negative ULM score with a
# significant p-value, consistent with MYC inhibition.
print(contrast_acts[contrast_acts$source == "MYC", ])


## ---- 10. Session information -----------------------------------------------

sessionInfo()
