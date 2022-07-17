#########################################################
#
# PBMC data analyzed without TCR
# This code is broken up into 5 main parts
# 1. Processing
# 2. Pseudobulk analysis
# 3. GSEA analysis
# 4. Counting cells in cell types
# 5. Figures from paper
#
#########################################################



######################################
# 1. Processing
######################################


library(Seurat)
library(SeuratDisk)
library(SeuratData)

# load PBMC reference
reference <- LoadH5Seurat("E:/pbmc_multimodal.h5seurat")

tcr<-c('TRA','TRAV1-1','TRAV1-2','TRAV2','TRAV3','TRAV4','TRAV5','TRAV6','TRAV7','TRAV8-1','TRAV8-2','TRAV8-3','TRAV8-4','TRAV8-5','TRAV8-6','TRAV8-7','TRAV9-1','TRAV9-2','TRAV10','TRAV11','TRAV12-1','TRAV12-2','TRAV12-3','TRAV13-1','TRAV13-2','TRAV14DV4','TRAV15','TRAV16','TRAV17','TRAV18','TRAV19','TRAV20','TRAV21','TRAV22','TRAV23DV6','TRAV24','TRAV25','TRAV26-1','TRAV26-2','TRAV27','TRAV28','TRAV29DV5','TRAV30','TRAV31','TRAV32','TRAV33','TRAV34','TRAV35','TRAV36DV7','TRAV37','TRAV38-1','TRAV38-2DV8','TRAV39','TRAV40','TRAV41','TRAJ1','TRAJ2','TRAJ3','TRAJ4','TRAJ5','TRAJ6','TRAJ7','TRAJ8','TRAJ9','TRAJ10','TRAJ11','TRAJ12','TRAJ13','TRAJ14','TRAJ15','TRAJ16','TRAJ17','TRAJ18','TRAJ19','TRAJ20','TRAJ21','TRAJ22','TRAJ23','TRAJ24','TRAJ25','TRAJ26','TRAJ27','TRAJ28','TRAJ29','TRAJ30','TRAJ31','TRAJ32','TRAJ33','TRAJ34','TRAJ35','TRAJ36','TRAJ37','TRAJ38','TRAJ39','TRAJ40','TRAJ41','TRAJ42','TRAJ43','TRAJ44','TRAJ45','TRAJ46','TRAJ47','TRAJ48','TRAJ49','TRAJ50','TRAJ51','TRAJ52','TRAJ53','TRAJ54','TRAJ55','TRAJ56','TRAJ57','TRAJ58','TRAJ59','TRAJ60','TRAJ61','TRAC','TRB','TRBV1','TRBV2','TRBV3-1','TRBV3-2','TRBV4-1','TRBV4-2','TRBV4-3','TRBV5-1','TRBV5-2','TRBV5-3','TRBV5-4','TRBV5-5','TRBV5-6','TRBV5-7','TRBV5-8','TRBV6-1','TRBV6-2','TRBV6-3','TRBV6-4','TRBV6-5','TRBV6-6','TRBV6-7','TRBV6-8','TRBV6-9','TRBV7-1','TRBV7-2','TRBV7-3','TRBV7-4','TRBV7-5','TRBV7-6','TRBV7-7','TRBV7-8','TRBV7-9','TRBV8-1','TRBV8-2','TRBV9','TRBV10-1','TRBV10-2','TRBV10-3','TRBV11-1','TRBV11-2','TRBV11-3','TRBV12-1','TRBV12-2','TRBV12-3','TRBV12-4','TRBV12-5','TRBV13','TRBV14','TRBV15','TRBV16','TRBV17','TRBV18','TRBV19','TRBV20-1','TRBV21-1','TRBV22-1','TRBV23-1','TRBV24-1','TRBV25-1','TRBV26','TRBV27','TRBV28','TRBV29-1','TRBV30','TRBVA','TRBVB','TRBD1','TRBVC','TRBD2','TRBJ1-1','TRBJ1-2','TRBJ1-3','TRBJ1-4','TRBJ1-5','TRBJ1-6','TRBJ2-1','TRBJ2-2','TRBJ2-2P','TRBJ2-3','TRBJ2-4','TRBJ2-5','TRBJ2-6','TRBJ2-7','TRBC1','TRBC2','TRD','TRDV1','TRDV2','TRDV3','TRDD1','TRDD2','TRDD3','TRDJ1','TRDJ2','TRDJ3','TRDJ4','TRDC','TRG','TRGV1','TRGV2','TRGV3','TRGV4','TRGV5','TRGV5P','TRGV6','TRGV7','TRGV8','TRGV9','TRGV10','TRGV11','TRGVA','TRGVB','TRGJ1','TRGJ2','TRGJP','TRGJP1','TRGJP2','TRGC1','TRGC2')

for (i in c('02','04','05','07','11')){
  # Reading in each sample as raw data
  assign(paste0("hNB13",i,"_data"),Read10X_h5(paste0("E:/Research/Gennaro_sc/Cellranger_output_files/hNB13",i,"_outs/filtered_feature_bc_matrix.h5")))
  # Removing the TCR gene counts as they aren't informative and make interpretation 
  # complicated
  assign(paste0("hNB13",i,"_data_no_tcr"),get(paste0("hNB13",i,"_data"))[!(row.names(get(paste0("hNB13",i,"_data"))) %in% tcr),])
  
  # Creating a seurat object for each sample
  assign(paste0("hNB13",i,"_no_tcr"),CreateSeuratObject(get(paste0("hNB13",i,"_data_no_tcr"))))
  # Adding metadata so when they're combined they can be differentiated
  assign(paste0("hNB13",i,"_no_tcr"),AddMetaData(get(paste0("hNB13",i,"_no_tcr")), i, col.name = "Sample"))
  
}

# This will no longer work as rscrublet has changed. The new version requires
# a parameter if it can't be determined and the input is a transpose, but can take
# sparse matrices now.
library(rscrublet)
for (i in c('hNB1302_no_tcr', 'hNB1304_no_tcr', 'hNB1305_no_tcr', 'hNB1307_no_tcr', 'hNB1311_no_tcr')) {
  
  doublets <- scrubDoublets(as.matrix(get(i)@assays[["RNA"]]@counts))
  cellnames <- colnames(as.matrix(get(i)@assays[["RNA"]]@counts))
  bool.scrub <- as.matrix(doublets[["scrubDoublets"]])
  rownames(bool.scrub)<-cellnames
  assign(i,AddMetaData(get(i), bool.scrub, col.name = "scrublets"))
  
}

# Here's a new rscrublet bit of code that should work, 
# but the above code was originally used.
# for (i in c('hNB1302_no_tcr', 'hNB1304_no_tcr', 'hNB1305_no_tcr', 'hNB1307_no_tcr', 'hNB1311_no_tcr')) {
#   
# doublets <- scrub_doublets(t(as.matrix(get(i)@assays[["RNA"]]@counts)))
# doublets <- call_doublets(doublets)
# cellnames <- colnames(get(i)@assays[["RNA"]]@counts)
# bool.scrub <- as.matrix(doublets[["predicted_doublets"]])
# rownames(bool.scrub)<-cellnames
# assign(i,AddMetaData(get(i),bool.scrub, col.name = 'scrublets'))
# }

hNB13_no_tcr <- merge(hNB1302_no_tcr, c(hNB1304_no_tcr, hNB1305_no_tcr, hNB1307_no_tcr, hNB1311_no_tcr))


rm(hNB1302_data, hNB1304_data, hNB1305_data, hNB1307_data, hNB1311_data)
rm(hNB1302_data_no_tcr, hNB1304_data_no_tcr, hNB1305_data_no_tcr, hNB1307_data_no_tcr, hNB1311_data_no_tcr)
rm(hNB1302_no_tcr, hNB1304_no_tcr, hNB1305_no_tcr, hNB1307_no_tcr, hNB1311_no_tcr)

# Filtering
hNB13_no_tcr[["percent.mt"]] <-
  PercentageFeatureSet(hNB13_no_tcr, pattern = "^MT-")
VlnPlot(
  hNB13_no_tcr,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "Sample",
  ncol = 3
)

# Didn't filter out cells with low MT percent or too few or many genes
# Not much convincing data to do that for modern techniques, 
# but here's how to do it
#
# hNB13 <- subset(hNB13, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)


# Normalize via pearson residuals
# Don't forget to cite the paper for SCTransform
hNB13_no_tcr <- SCTransform(hNB13_no_tcr)

# Run PCA dimensionality reduction
hNB13_no_tcr <- RunPCA(hNB13_no_tcr)

# Run UMAP on PCA
hNB13_no_tcr <- RunUMAP(hNB13_no_tcr, dims = 1:50)

# Unsupervised clustering techniques
hNB13_no_tcr <- FindNeighbors(hNB13_no_tcr, dims = 1:50)
hNB13_no_tcr <- FindClusters(hNB13_no_tcr, resolution = 0.1) 
# adjust resolution to find however many clusters you want, 0.001 for 2-3 0.8 for ~12-15

# Note, that due to the nature of UMAP, the graphs will probably look
# noticeably different if you run the analysis yourself, but the local
# structure should be preserved, and the analyses done should be identical


# Plot the data by sample and by cluster
p1 <- DimPlot(hNB13_no_tcr, 
              reduction = "umap", 
              group.by = "Sample")+ RotatedAxis()
p2 <- DimPlot(hNB13_no_tcr,
              reduction = "umap",
              label = TRUE,
              repel = TRUE)
p1 + p2

anchors <- FindTransferAnchors(
  reference = reference,
  query = hNB13_no_tcr,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50,
  # This is here to fix a bug they have, I'm not sure of it's impact
  recompute.residuals = FALSE
)

# Map cell types and predicted Ab from Reference to hNB13_no_tcr
hNB13_no_tcr <- MapQuery(
  anchorset = anchors,
  query = hNB13_no_tcr,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    celltype.l3 = "celltype.l3",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca",
  reduction.model = "wnn.umap"
)

p2 <- DimPlot(hNB13_no_tcr,
              reduction = "umap",
              label = TRUE,
              repel = TRUE,
              group.by = "predicted.celltype.l2")+ RotatedAxis()
p1 + p2

save(hNB13_no_tcr, file = "E:/Research/Gennaro_sc/hNB13_no_tcr_analyzed.Rdata")


######################################
# 2. Pseudobulk analysis of genesets
######################################

monocytes<-subset(hNB13_no_tcr, subset = predicted.celltype.l1 == "Mono")

pseudo_bulk_mono<-AggregateExpression(monocytes, group.by = "Sample")

write.csv(pseudo_bulk_mono[["SCT"]], "E:/Research/Gennaro_sc/pseudo_bulk_results/pseudo_bulk_monocyte_sct.csv")
write.csv(pseudo_bulk_mono[["RNA"]], "E:/Research/Gennaro_sc/pseudo_bulk_results/pseudo_bulk_monocyte_rna.csv")

# Showing only genes in gene sets

# Reads in gene sets and formats them as vectors removing the top two lines
alpha_set<-as.vector(t(read.table("E:/Research/Gennaro_sc/ifna_response_geneset.txt", sep = '\n')))[-2][-1]
gamma_set<-as.vector(t(read.table("E:/Research/Gennaro_sc/ifng_response_geneset.txt", sep = '\n')))[-2][-1]
tnfa_by_nfkb_set<-as.vector(t(read.table("E:/Research/Gennaro_sc/tnfa_signaling_via_nfkb_geneset.txt", sep = '\n')))[-2][-1]
# These files can be provided upon request, but are simply a text file with gene names
# from the MSigDB files


# Keep only the genes in the gene lists 
alpha_set_pseudo_bulk_mono_sct<-pseudo_bulk_mono[["SCT"]][rownames(pseudo_bulk_mono[["SCT"]]) %in% alpha_set, ]
gamma_set_pseudo_bulk_mono_sct<-pseudo_bulk_mono[["SCT"]][rownames(pseudo_bulk_mono[["SCT"]]) %in% gamma_set, ]
tnfa_by_nfkb_set_pseudo_bulk_mono_sct<-pseudo_bulk_mono[["SCT"]][rownames(pseudo_bulk_mono[["SCT"]]) %in% tnfa_by_nfkb_set, ]

alpha_set_pseudo_bulk_mono_rna<-pseudo_bulk_mono[["RNA"]][rownames(pseudo_bulk_mono[["RNA"]]) %in% alpha_set, ]
gamma_set_pseudo_bulk_mono_rna<-pseudo_bulk_mono[["RNA"]][rownames(pseudo_bulk_mono[["RNA"]]) %in% gamma_set, ]
tnfa_by_nfkb_set_pseudo_bulk_mono_rna<-pseudo_bulk_mono[["RNA"]][rownames(pseudo_bulk_mono[["RNA"]]) %in% tnfa_by_nfkb_set, ]

write.csv(alpha_set_pseudo_bulk_mono_sct, "E:/Research/Gennaro_sc/pseudo_bulk_results/pseudo_bulk_monocyte_alpha_set_sct.csv")
write.csv(gamma_set_pseudo_bulk_mono_sct, "E:/Research/Gennaro_sc/pseudo_bulk_results/pseudo_bulk_monocyte_gamma_set_sct.csv")
write.csv(tnfa_by_nfkb_set_pseudo_bulk_mono_sct, "E:/Research/Gennaro_sc/pseudo_bulk_results/pseudo_bulk_monocyte_tnfa_by_nfkb_set_sct.csv")

write.csv(alpha_set_pseudo_bulk_mono_rna, "E:/Research/Gennaro_sc/pseudo_bulk_results/pseudo_bulk_monocyte_alpha_set_rna.csv")
write.csv(gamma_set_pseudo_bulk_mono_rna, "E:/Research/Gennaro_sc/pseudo_bulk_results/pseudo_bulk_monocyte_gamma_set_rna.csv")
write.csv(tnfa_by_nfkb_set_pseudo_bulk_mono_rna, "E:/Research/Gennaro_sc/pseudo_bulk_results/pseudo_bulk_monocyte_tnfa_by_nfkb_set_rna.csv")

######################################
# 3. GSEA analysis
######################################

library(fgsea)
library(data.table)
library(ggplot2)

#library(clusterProfiler)
#library(GOSemSim)

# Since looking through all genes is quite time consuming future is used
# to speed up the process through parallelization 
library(future)
plan(strategy = "multicore", workers = 16)  # change to however many threads/cores you have

# Load the object if not already present 
load("E:/Research/Gennaro_sc/hNB13_analyzed.Rdata")

# Change idents to do DE
Idents(hNB13) <- "days.post.infection"

# Only look at differences across monocytes
mono_hNB13 <- subset(hNB13, subset = predicted.celltype.l1 == "Mono")

# Find markers, which will take a while. The future library may help
day2_mono_markers <-
  FindMarkers(
    mono_hNB13,
    ident.1 = 2,
    logfc.threshold = 0,
    min.pct = 0,
    test.use = "MAST"
  )
# Since the threshold is 0 it will take a long time!

# Export it so that we don't need to do this again due to 
# how long it takes
write.csv(day2_mono_markers, file = "E:/Research/Gennaro_sc/day2_mono_markers.csv")

# This is how to read it back in, the row.names argument is used 
# since when exported the rownames are added to a column labeled by "X"
day2_mono_markers<-read.csv("E:/Research/Gennaro_sc/day2_mono_markers.csv", row.names = "X")

# Reformat into a ranking of log2FC
ranks <- day2_mono_markers$avg_log2FC
names(ranks) <- rownames(day2_mono_markers)

# Use this to reset the plots
dev.off(dev.list()["RStudioGD"])

# Onboard the Pathways or Gene Sets
kegg_sigdb_pathways <- gmtPathways("E:/Research/Gennaro_sc/GSEA/c2.cp.kegg.v7.4.symbols.gmt")
# Run fGSEA (fast Gene Set Enrichment Analysis)
# minSize and maxSize exclude pathways or gene sets with less than minSize
# genes or more than maxSize genes. Not sure why this is done.
fgseaRes_kegg <- fgsea(pathways = kegg_sigdb_pathways, 
                       stats = ranks,
                       minSize=10,
                       maxSize=500)

# Pick most significantly upregulated
topPathwaysUp <- fgseaRes_kegg[ES > 0][head(order(pval), n=10), pathway]
# Pick most significantly downregulated
topPathwaysDown <- fgseaRes_kegg[ES < 0][head(order(pval), n=10), pathway]
# Combine the two
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

# Plot them all via ranking metric
plotGseaTable(kegg_sigdb_pathways[topPathways], ranks, fgseaRes_kegg, 
              gseaParam=0.5)

# Can look at main branches to remove redundancy by collapsing pathways
# the plotting just the top 20 of these
collapsedPathways <- collapsePathways(fgseaRes_kegg[order(pval)][padj < 0.01], 
                                      kegg_sigdb_pathways, ranks)
mainPathways <- fgseaRes_kegg[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]
plotGseaTable(kegg_sigdb_pathways[mainPathways], ranks, fgseaRes_kegg, 
              gseaParam = 0.5)

# Can also plot specific running Enrichment plots like this
plotEnrichment(kegg_sigdb_pathways[["KEGG_TIGHT_JUNCTION"]],
               ranks) + labs(title="Genes related to Tight Junctions")

# Use this to reset the plots
dev.off(dev.list()["RStudioGD"])


# Repeat the same procedure for other pathways or gene sets as shown below.
# KEGG and GO seem the most useful since there are a lot of extraneous
# sets in the other curated .cmt files
# May try WikiPathways as well



# Use this to reset the plots
dev.off(dev.list()["RStudioGD"])

GOBP_sigdb_pathways <- gmtPathways("E:/Research/Gennaro_sc/GSEA/c5.go.bp.v7.4.symbols.gmt")
fgseaRes_GOBP <- fgsea(pathways = GOBP_sigdb_pathways, 
                       stats = ranks,
                       minSize=10,
                       maxSize=500)
topPathwaysUp <- fgseaRes_GOBP[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes_GOBP[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(GOBP_sigdb_pathways[topPathways], ranks, fgseaRes_GOBP, 
              gseaParam=0.5)

# Use this to reset the plots
dev.off(dev.list()["RStudioGD"])
GOCC_sigdb_pathways <- gmtPathways("E:/Research/Gennaro_sc/GSEA/c5.go.cc.v7.4.symbols.gmt")
fgseaRes_GOCC <- fgsea(pathways = GOCC_sigdb_pathways, 
                       stats = ranks,
                       minSize=10,
                       maxSize=500)
topPathwaysUp <- fgseaRes_GOCC[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes_GOCC[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(GOCC_sigdb_pathways[topPathways], ranks, fgseaRes_GOCC, 
              gseaParam=0.5)

# Use this to reset the plots
dev.off(dev.list()["RStudioGD"])
GOMF_sigdb_pathways <- gmtPathways("E:/Research/Gennaro_sc/GSEA/c5.go.mf.v7.4.symbols.gmt")
fgseaRes_GOMF <- fgsea(pathways = GOMF_sigdb_pathways, 
                       stats = ranks,
                       minSize=10,
                       maxSize=500)
topPathwaysUp <- fgseaRes_GOMF[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes_GOMF[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(GOMF_sigdb_pathways[topPathways], ranks, fgseaRes_GOMF, 
              gseaParam=0.5)

dev.off(dev.list()["RStudioGD"])
c8type_sigdb_pathways <- gmtPathways("E:/Research/Gennaro_sc/GSEA/c8.all.v7.4.symbols.gmt")
fgseaRes_typec8 <- fgsea(pathways = c8type_sigdb_pathways, 
                         stats = ranks,
                         minSize=10,
                         maxSize=500)
topPathwaysUp <- fgseaRes_typec8[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes_typec8[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(c8type_sigdb_pathways[topPathways], ranks, fgseaRes_typec8, 
              gseaParam=0.5)

# Use this to reset the plots
dev.off(dev.list()["RStudioGD"])
h_sigdb_pathways <- gmtPathways("E:/Research/Gennaro_sc/GSEA/h.all.v7.3.symbols.gmt")
fgseaRes_hall <- fgsea(pathways = h_sigdb_pathways, 
                       stats = ranks,
                       minSize=10,
                       maxSize=500)
topPathwaysUp <- fgseaRes_hall[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes_hall[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(h_sigdb_pathways[topPathways], ranks, fgseaRes_hall, 
              gseaParam=0.5)

# Use this to reset the plots
dev.off(dev.list()["RStudioGD"])
c7_all_pathways <- gmtPathways("E:/Research/Gennaro_sc/GSEA/c7.all.v7.3.symbols.gmt")
fgseaRes_c7all <- fgsea(pathways = c7_all_pathways, 
                        stats = ranks,
                        minSize=10,
                        maxSize=500)
c7_sigdb_pathways <- gmtPathways("E:/Research/Gennaro_sc/GSEA/c7.immunesigdb.v7.3.symbols.gmt")
fgseaRes_c7sigdb <- fgsea(pathways = c7_sigdb_pathways, 
                          stats = ranks,
                          minSize=10,
                          maxSize=500)
plotEnrichment(c7_sigdb_pathways[["GSE42724_NAIVE_BCELL_VS_PLASMABLAST_UP"]],
               ranks) + labs(title="B v Plasmablast")
topPathwaysUp <- fgseaRes_c7sigdb[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes_c7sigdb[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(c7_sigdb_pathways[topPathways], ranks, fgseaRes_c7sigdb, 
              gseaParam=0.5)
collapsedPathways <- collapsePathways(fgseaRes_c7sigdb[order(pval)][padj < 0.01], 
                                      c7_sigdb_pathways, ranks)
mainPathways <- fgseaRes_c7sigdb[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]
plotGseaTable(head(c7_sigdb_pathways[mainPathways], n=20), ranks, fgseaRes_c7sigdb, 
              gseaParam = 0.5)

# Use this to reset the plots
dev.off(dev.list()["RStudioGD"])


##############################################
# Comparing Sample 2 values of each cell type cluster to all the rest

# Loading pathways to test via fGSEA (from MSigDB)
kegg_sigdb_pathways <- gmtPathways("E:/Research/Gennaro_sc/GSEA/c2.cp.kegg.v7.4.symbols.gmt")
GOMF_sigdb_pathways <- gmtPathways("E:/Research/Gennaro_sc/GSEA/c5.go.mf.v7.4.symbols.gmt")
GOBP_sigdb_pathways <- gmtPathways("E:/Research/Gennaro_sc/GSEA/c5.go.bp.v7.4.symbols.gmt")
GOCC_sigdb_pathways <- gmtPathways("E:/Research/Gennaro_sc/GSEA/c5.go.cc.v7.4.symbols.gmt")


# For l1

# An automated way to do this for each cell type

Idents(hNB13) <- "predicted.celltype.l1"
l1_idents <- levels(Idents(hNB13))
Idents(hNB13) <- "days.post.infection"

for (i in 1:7) {
  
  cat(paste0("Starting analysis on ",l1_idents[i]," (",i,"/7)\n"))
  
  # Using MAST and Seurat to rank all genes in terms of differential expression log2fc
  assign(
    paste0("l1_",l1_idents[i], "_ranks"),
    FindMarkers(
      subset(hNB13, subset = predicted.celltype.l1 == l1_idents[i]),
      ident.1 = 2,         #for sample 2
      logfc.threshold = 0,
      min.pct = 0,
      test.use = "MAST"
    )
  )
  
  cat(paste0("Running fGSEA on ",l1_idents[i]," (",i,"/7)\n"))
  
  # Moving the currently calculated genes into a ranked list
  ranks <- get(paste0("l1_",l1_idents[i], "_ranks"))$avg_log2FC
  names(ranks) <- rownames(get(paste0("l1_",l1_idents[i], "_ranks")))
  
  # Running multiple fGSEAs and saving the outputs
  fgseaRes_kegg <- fgsea(pathways = kegg_sigdb_pathways, 
                         stats = ranks,
                         minSize=10,
                         maxSize=500)
  fgseaRes_kegg$leadingEdge<-as.character(fgseaRes_kegg$leadingEdge)
  write.csv(fgseaRes_kegg, file = paste0("E:/Research/Gennaro_sc/GSEA/l1 enrichments/",l1_idents[i],"_l1_KEGG_pathway_enrichment.csv"))
  
  fgseaRes_GOMF <- fgsea(pathways = GOMF_sigdb_pathways, 
                         stats = ranks,
                         minSize=10,
                         maxSize=500)
  fgseaRes_GOMF$leadingEdge<-as.character(fgseaRes_GOMF$leadingEdge)
  write.csv(fgseaRes_GOMF, file = paste0("E:/Research/Gennaro_sc/GSEA/l1 enrichments/",l1_idents[i],"_l1_GOMF_pathway_enrichment.csv"))
  
  fgseaRes_GOBP <- fgsea(pathways = GOBP_sigdb_pathways, 
                         stats = ranks,
                         minSize=10,
                         maxSize=500)
  fgseaRes_GOBP$leadingEdge<-as.character(fgseaRes_GOBP$leadingEdge)
  write.csv(fgseaRes_GOBP, file = paste0("E:/Research/Gennaro_sc/GSEA/l1 enrichments/",l1_idents[i],"_l1_GOBP_pathway_enrichment.csv"))
  
  fgseaRes_GOCC <- fgsea(pathways = GOCC_sigdb_pathways, 
                         stats = ranks,
                         minSize=10,
                         maxSize=500)
  fgseaRes_GOCC$leadingEdge<-as.character(fgseaRes_GOCC$leadingEdge)
  write.csv(fgseaRes_GOCC, file = paste0("E:/Research/Gennaro_sc/GSEA/l1 enrichments/",l1_idents[i],"_l1_GOCC_pathway_enrichment.csv"))
  
  cat(paste0("Finished analysis on ",l1_idents[i]," (",i,"/7)\n"))
  
}

# For l2

# An automated way to do this for each cell type
Idents(hNB13) <- "predicted.celltype.l2"
l2_idents <- levels(Idents(hNB13))
Idents(hNB13) <- "days.post.infection"

for (i in 1:28) {
  
  cat(paste0("Starting analysis on ",l2_idents[i]," (",i,"/28)\n"))
  
  assign(
    paste0("l2_",l2_idents[i], "_ranks"),
    FindMarkers(
      subset(hNB13, subset = predicted.celltype.l2 == l2_idents[i]),
      ident.1 = 2,         #for sample 2
      logfc.threshold = 0,
      min.pct = 0,
      test.use = "MAST"
    )
  )
  
  cat(paste0("Running fGSEA on ",l2_idents[i]," (",i,"/28)\n"))
  
  ranks <- get(paste0("l2_",l2_idents[i], "_ranks"))$avg_log2FC
  names(ranks) <- rownames(get(paste0("l2_",l2_idents[i], "_ranks")))
  fgseaRes_kegg <- fgsea(pathways = kegg_sigdb_pathways, 
                         stats = ranks,
                         minSize=10,
                         maxSize=500)
  fgseaRes_kegg$leadingEdge<-as.character(fgseaRes_kegg$leadingEdge)
  write.csv(fgseaRes_kegg, file = paste0("E:/Research/Gennaro_sc/GSEA/l2 enrichments/",l2_idents[i],"_l2_KEGG_pathway_enrichment.csv"))
  
  fgseaRes_GOMF <- fgsea(pathways = GOMF_sigdb_pathways, 
                         stats = ranks,
                         minSize=10,
                         maxSize=500)
  fgseaRes_GOMF$leadingEdge<-as.character(fgseaRes_GOMF$leadingEdge)
  write.csv(fgseaRes_GOMF, file = paste0("E:/Research/Gennaro_sc/GSEA/l2 enrichments/",l2_idents[i],"_l2_GOMF_pathway_enrichment.csv"))
  
  fgseaRes_GOBP <- fgsea(pathways = GOBP_sigdb_pathways, 
                         stats = ranks,
                         minSize=10,
                         maxSize=500)
  fgseaRes_GOBP$leadingEdge<-as.character(fgseaRes_GOBP$leadingEdge)
  write.csv(fgseaRes_GOBP, file = paste0("E:/Research/Gennaro_sc/GSEA/l2 enrichments/",l2_idents[i],"_l2_GOBP_pathway_enrichment.csv"))
  
  fgseaRes_GOCC <- fgsea(pathways = GOCC_sigdb_pathways, 
                         stats = ranks,
                         minSize=10,
                         maxSize=500)
  fgseaRes_GOCC$leadingEdge<-as.character(fgseaRes_GOCC$leadingEdge)
  write.csv(fgseaRes_GOCC, file = paste0("E:/Research/Gennaro_sc/GSEA/l2 enrichments/",l2_idents[i],"_l2_GOCC_pathway_enrichment.csv"))
  
  cat(paste0("Finished analysis on ",l2_idents[i]," (",i,"/28)\n"))
}
# 
# # l1 and l2 without the fGSEA so that it can run without errors hopefully 
# for (i in 1:7) {
#   assign(
#     paste0(l1_idents[i], "_ranks"),
#     FindMarkers(
#       subset(hNB13, subset = predicted.celltype.l1 == l1_idents[i]),
#       ident.1 = 2,         #for sample 2
#       logfc.threshold = 0,
#       min.pct = 0,
#       test.use = "MAST"
#     )
#   )
#   ranks <- get(paste0(l1_idents[i], "_ranks"))$avg_log2FC
#   names(ranks) <- rownames(get(paste0(l1_idents[i], "_ranks")))
# }
# 
# for (i in 1:28) {
#   assign(
#     paste0(l2_idents[i], "_ranks"),
#     FindMarkers(
#       subset(hNB13, subset = predicted.celltype.l2 == l2_idents[i]),
#       ident.1 = 2,         #for sample 2
#       logfc.threshold = 0,
#       min.pct = 0,
#       test.use = "MAST"
#     )
#   )
#   ranks <- get(paste0(l2_idents[i], "_ranks"))$avg_log2FC
#   names(ranks) <- rownames(get(paste0(l2_idents[i], "_ranks")))
# 
# }
# 
# # Now that I've figured it out we can run the fGSEA and save the results for each l1 and l2
# for (i in 1:7) {
#   ranks <- get(paste0(l1_idents[i], "_ranks"))$avg_log2FC
#   names(ranks) <- rownames(get(paste0(l1_idents[i], "_ranks")))
#   
#   fgseaRes_kegg <- fgsea(pathways = kegg_sigdb_pathways, 
#                          stats = ranks,
#                          minSize=10,
#                          maxSize=500)
#   fgseaRes_kegg$leadingEdge<-as.character(fgseaRes_kegg$leadingEdge)
#   write.csv(fgseaRes_kegg, file = paste0("E:/Research/Gennaro_sc/GSEA/",l1_idents[i],"_l1_KEGG_pathway_enrichment.csv"))
#   
#   fgseaRes_GOMF <- fgsea(pathways = GOMF_sigdb_pathways, 
#                          stats = ranks,
#                          minSize=10,
#                          maxSize=500)
#   fgseaRes_GOMF$leadingEdge<-as.character(fgseaRes_GOMF$leadingEdge)
#   write.csv(fgseaRes_GOMF, file = paste0("E:/Research/Gennaro_sc/GSEA/",l1_idents[i],"_l1_GOMF_pathway_enrichment.csv"))
#   
#   fgseaRes_GOBP <- fgsea(pathways = GOBP_sigdb_pathways, 
#                          stats = ranks,
#                          minSize=10,
#                          maxSize=500)
#   fgseaRes_GOBP$leadingEdge<-as.character(fgseaRes_GOBP$leadingEdge)
#   write.csv(fgseaRes_GOBP, file = paste0("E:/Research/Gennaro_sc/GSEA/",l1_idents[i],"_l1_GOBP_pathway_enrichment.csv"))
#   
#   fgseaRes_GOCC <- fgsea(pathways = GOCC_sigdb_pathways, 
#                          stats = ranks,
#                          minSize=10,
#                          maxSize=500)
#   fgseaRes_GOCC$leadingEdge<-as.character(fgseaRes_GOCC$leadingEdge)
#   write.csv(fgseaRes_GOCC, file = paste0("E:/Research/Gennaro_sc/GSEA/",l1_idents[i],"_l1_GOCC_pathway_enrichment.csv"))
#   
#   
# }
# 
# for (i in 1:28) {
#   ranks <- get(paste0(l2_idents[i], "_ranks"))$avg_log2FC
#   names(ranks) <- rownames(get(paste0(l2_idents[i], "_ranks")))
#   
#   fgseaRes_kegg <- fgsea(pathways = kegg_sigdb_pathways, 
#                          stats = ranks,
#                          minSize=10,
#                          maxSize=500)
#   fgseaRes_kegg$leadingEdge<-as.character(fgseaRes_kegg$leadingEdge)
#   write.csv(fgseaRes_kegg, file = paste0("E:/Research/Gennaro_sc/GSEA/l2 enrichments/",l2_idents[i],"_l2_KEGG_pathway_enrichment.csv"))
#   
#   fgseaRes_GOMF <- fgsea(pathways = GOMF_sigdb_pathways, 
#                          stats = ranks,
#                          minSize=10,
#                          maxSize=500)
#   fgseaRes_GOMF$leadingEdge<-as.character(fgseaRes_GOMF$leadingEdge)
#   write.csv(fgseaRes_GOMF, file = paste0("E:/Research/Gennaro_sc/GSEA/l2 enrichments/",l2_idents[i],"_l2_GOMF_pathway_enrichment.csv"))
#   
#   fgseaRes_GOBP <- fgsea(pathways = GOBP_sigdb_pathways, 
#                          stats = ranks,
#                          minSize=10,
#                          maxSize=500)
#   fgseaRes_GOBP$leadingEdge<-as.character(fgseaRes_GOBP$leadingEdge)
#   write.csv(fgseaRes_GOBP, file = paste0("E:/Research/Gennaro_sc/GSEA/l2 enrichments/",l2_idents[i],"_l2_GOBP_pathway_enrichment.csv"))
#   
#   fgseaRes_GOCC <- fgsea(pathways = GOCC_sigdb_pathways, 
#                          stats = ranks,
#                          minSize=10,
#                          maxSize=500)
#   fgseaRes_GOCC$leadingEdge<-as.character(fgseaRes_GOCC$leadingEdge)
#   write.csv(fgseaRes_GOCC, file = paste0("E:/Research/Gennaro_sc/GSEA/l2 enrichments/",l2_idents[i],"_l2_GOCC_pathway_enrichment.csv"))
#   
#   
# }

# Run again for marker genes and save the output

for (i in 1:7) {
  
  cat(paste0("Starting analysis on ",l1_idents[i]," (",i,"/7)\n"))
  
  assign(
    paste0("l1_",l1_idents[i], "_markers"),
    FindMarkers(
      subset(hNB13, subset = predicted.celltype.l1 == l1_idents[i]),
      ident.1 = 2,         #for sample 2
      logfc.threshold = 0.25,
      test.use = "MAST"
    )
  )
  write.csv(get(paste0("l1_",l1_idents[i], "_markers")), file = paste0("E:/Research/Gennaro_sc/GSEA/marker genes/",l1_idents[i],"_l1_sample_2_marker_genes.csv"))
}


for (i in 1:28) {
  
  cat(paste0("Starting analysis on ",l2_idents[i]," (",i,"/28)\n"))
  
  assign(
    paste0("l2_",l2_idents[i], "_markers"),
    FindMarkers(
      subset(hNB13, subset = predicted.celltype.l2 == l2_idents[i]),
      ident.1 = 2,         #for sample 2
      logfc.threshold = 0.25,
      test.use = "MAST"
    )
  )
  write.csv(get(paste0("l2_",l2_idents[i], "_markers")), file = paste0("E:/Research/Gennaro_sc/GSEA/marker genes/",l2_idents[i],"_l2_sample_2_marker_genes.csv"))
  
}

# Run again for marker genes between sample 2 and 11 and save the output

for (i in 1:7) {
  
  cat(paste0("Starting analysis on ",l1_idents[i]," (",i,"/7)\n"))
  
  assign(
    paste0("l1_",l1_idents[i], "_markers_2to11"),
    FindMarkers(
      subset(hNB13, subset = predicted.celltype.l1 == l1_idents[i]),
      ident.1 = 2,         #for sample 2
      ident.2 = 11,        #compared to sample 11
      logfc.threshold = 0.25,
      test.use = "MAST"
    )
  )
  write.csv(get(paste0("l1_",l1_idents[i], "_markers_2to11")), file = paste0("E:/Research/Gennaro_sc/GSEA/marker genes/",l1_idents[i],"_l1_sample_2to11_marker_genes.csv"))
}


for (i in 1:28) {
  
  cat(paste0("Starting analysis on ",l2_idents[i]," (",i,"/28)\n"))
  
  assign(
    paste0("l2_",l2_idents[i], "_markers_2to11"),
    FindMarkers(
      subset(hNB13, subset = predicted.celltype.l2 == l2_idents[i]),
      ident.1 = 2,         #for sample 2
      ident.2 = 11,        #compared to sample 11
      logfc.threshold = 0.25,
      test.use = "MAST"
    )
  )
  write.csv(get(paste0("l2_",l2_idents[i], "_markers_2to11")), file = paste0("E:/Research/Gennaro_sc/GSEA/marker genes/",l2_idents[i],"_l2_sample_2to11_marker_genes.csv"))
  
}


# An automated way to do GSEA for each cell type comparing sample 2 to 11
# For l1

Idents(hNB13) <- "predicted.celltype.l1"
l1_idents <- levels(Idents(hNB13))
Idents(hNB13) <- "days.post.infection"

for (i in 1:7) {
  
  cat(paste0("Starting analysis on ",l1_idents[i]," (",i,"/7)\n"))
  
  #Adding tryCatch
  tryCatch( 
    # Using MAST and Seurat to rank all genes in terms of differential expression log2fc
    expr = {
      assign(
        paste0("l1_",l1_idents[i], "_ranks_2to11"),
        FindMarkers(
          subset(hNB13, subset = predicted.celltype.l1 == l1_idents[i]),
          ident.1 = 2,         #for sample 2
          ident.2 = 11,        #compared to sample 11
          logfc.threshold = 0,
          min.pct = 0,
          test.use = "MAST"
        )
      ) 
      
      
      cat(paste0("Running fGSEA on ",l1_idents[i]," (",i,"/7)\n"))
      
      # Moving the currently calculated genes into a ranked list
      ranks <- get(paste0("l1_",l1_idents[i], "_ranks_2to11"))$avg_log2FC
      names(ranks) <- rownames(get(paste0("l1_",l1_idents[i], "_ranks_2to11")))
      
      # Running multiple fGSEAs and saving the outputs
      fgseaRes_kegg <- fgsea(pathways = kegg_sigdb_pathways, 
                             stats = ranks,
                             minSize=10,
                             maxSize=500)
      fgseaRes_kegg$leadingEdge<-as.character(fgseaRes_kegg$leadingEdge)
      write.csv(fgseaRes_kegg, file = paste0("E:/Research/Gennaro_sc/GSEA/2to11/l1 enrichments/",l1_idents[i],"_l1_KEGG_pathway_enrichment.csv"))
      
      fgseaRes_GOMF <- fgsea(pathways = GOMF_sigdb_pathways, 
                             stats = ranks,
                             minSize=10,
                             maxSize=500)
      fgseaRes_GOMF$leadingEdge<-as.character(fgseaRes_GOMF$leadingEdge)
      write.csv(fgseaRes_GOMF, file = paste0("E:/Research/Gennaro_sc/GSEA/2to11/l1 enrichments/",l1_idents[i],"_l1_GOMF_pathway_enrichment.csv"))
      
      fgseaRes_GOBP <- fgsea(pathways = GOBP_sigdb_pathways, 
                             stats = ranks,
                             minSize=10,
                             maxSize=500)
      fgseaRes_GOBP$leadingEdge<-as.character(fgseaRes_GOBP$leadingEdge)
      write.csv(fgseaRes_GOBP, file = paste0("E:/Research/Gennaro_sc/GSEA/2to11/l1 enrichments/",l1_idents[i],"_l1_GOBP_pathway_enrichment.csv"))
      
      fgseaRes_GOCC <- fgsea(pathways = GOCC_sigdb_pathways, 
                             stats = ranks,
                             minSize=10,
                             maxSize=500)
      fgseaRes_GOCC$leadingEdge<-as.character(fgseaRes_GOCC$leadingEdge)
      write.csv(fgseaRes_GOCC, file = paste0("E:/Research/Gennaro_sc/GSEA/2to11/l1 enrichments/",l1_idents[i],"_l1_GOCC_pathway_enrichment.csv"))
      message("Completed without error")},
    
    error = function(e){
      message("There was an error")},
    
    warning = function(w){
      message("There was a warning")}
    
  )
  cat(paste0("Finished analysis on ",l1_idents[i]," (",i,"/7)\n"))
  
}

# For l2

for (i in 1:28) {
  
  cat(paste0("Starting analysis on ",l2_idents[i]," (",i,"/28)\n"))
  #Adding tryCatch
  tryCatch( 
    # Using MAST and Seurat to rank all genes in terms of differential expression log2fc
    expr = {
      assign(
        paste0("l2_",l2_idents[i], "_ranks_2to11"),
        FindMarkers(
          subset(hNB13, subset = predicted.celltype.l2 == l2_idents[i]),
          ident.1 = 2,         #for sample 2
          ident.2 = 11,        #compared to sample 11
          logfc.threshold = 0,
          min.pct = 0,
          test.use = "MAST"
        )
      )
      cat(paste0("Running fGSEA on ",l2_idents[i]," (",i,"/28)\n"))
      
      ranks <- get(paste0("l2_",l2_idents[i], "_ranks_2to11"))$avg_log2FC
      names(ranks) <- rownames(get(paste0("l2_",l2_idents[i], "_ranks_2to11")))
      
      
      fgseaRes_kegg <- fgsea(pathways = kegg_sigdb_pathways, 
                             stats = ranks,
                             minSize=10,
                             maxSize=500)
      fgseaRes_kegg$leadingEdge<-as.character(fgseaRes_kegg$leadingEdge)
      write.csv(fgseaRes_kegg, file = paste0("E:/Research/Gennaro_sc/GSEA/2to11/l2 enrichments/",l2_idents[i],"_l2_KEGG_pathway_enrichment.csv"))
      
      
      fgseaRes_GOMF <- fgsea(pathways = GOMF_sigdb_pathways, 
                             stats = ranks,
                             minSize=10,
                             maxSize=500)
      fgseaRes_GOMF$leadingEdge<-as.character(fgseaRes_GOMF$leadingEdge)
      write.csv(fgseaRes_GOMF, file = paste0("E:/Research/Gennaro_sc/GSEA/2to11/l2 enrichments/",l2_idents[i],"_l2_GOMF_pathway_enrichment.csv"))
      
      
      
      fgseaRes_GOBP <- fgsea(pathways = GOBP_sigdb_pathways, 
                             stats = ranks,
                             minSize=10,
                             maxSize=500)
      fgseaRes_GOBP$leadingEdge<-as.character(fgseaRes_GOBP$leadingEdge)
      write.csv(fgseaRes_GOBP, file = paste0("E:/Research/Gennaro_sc/GSEA/2to11/l2 enrichments/",l2_idents[i],"_l2_GOBP_pathway_enrichment.csv"))
      
      fgseaRes_GOCC <- fgsea(pathways = GOCC_sigdb_pathways, 
                             stats = ranks,
                             minSize=10,
                             maxSize=500)
      fgseaRes_GOCC$leadingEdge<-as.character(fgseaRes_GOCC$leadingEdge)
      write.csv(fgseaRes_GOCC, file = paste0("E:/Research/Gennaro_sc/GSEA/2to11/l2 enrichments/",l2_idents[i],"_l2_GOCC_pathway_enrichment.csv"))
      message("Completed without error")},
    
    error = function(e){
      message("There was an error")},
    
    warning = function(w){
      message("There was a warning")}
    
  )
  
  cat(paste0("Finished analysis on ",l2_idents[i]," (",i,"/28)\n"))
}

######################################
# 4. Cell type counts and markers
######################################


library(Seurat)
library(SeuratDisk)
library(openxlsx)

# Finding cell counts for MapQuery cell typing

# Showing healthy PBMC cell amounts and proportions

Idents(reference) <- "celltype.l1"
control.l1 <- table(Idents(reference))
prop.control.l1 <- prop.table(table(Idents(reference)))

Idents(reference) <- "celltype.l2"
control.l2 <- table(Idents(reference))
prop.control.l2 <- prop.table(table(Idents(reference)))

Idents(reference) <- "celltype.l3"
control.l3 <- table(Idents(reference))
prop.control.l3 <- prop.table(table(Idents(reference)))

# Showing amounts of each cell type by granularity level
Idents(hNB13) <- "predicted.celltype.l1"
sample.count.l1 <- table(Idents(hNB13), hNB13$days.post.infection)
sample.prop.count.l1 <-
  prop.table(table(Idents(hNB13), hNB13$days.post.infection), margin = 2)

Idents(hNB13) <- "predicted.celltype.l2"
sample.count.l2 <- table(Idents(hNB13), hNB13$days.post.infection)
sample.prop.count.l2 <-
  prop.table(table(Idents(hNB13), hNB13$days.post.infection), margin = 2)

Idents(hNB13) <- "predicted.celltype.l3"
sample.count.l3 <- table(Idents(hNB13), hNB13$days.post.infection)
sample.prop.count.l3 <-
  prop.table(table(Idents(hNB13), hNB13$days.post.infection), margin = 2)

# Make a workbook object via openxlsx
cell_counts <- createWorkbook("cell count")

# Add all the worksheets by name
addWorksheet(cell_counts, "control l1")
addWorksheet(cell_counts, "proportions control l1")
addWorksheet(cell_counts, "control l2")
addWorksheet(cell_counts, "proportions control l2")
addWorksheet(cell_counts, "control l3")
addWorksheet(cell_counts, "proportions control l3")

addWorksheet(cell_counts, "sample l1")
addWorksheet(cell_counts, "proportions sample l1")
addWorksheet(cell_counts, "sample l2")
addWorksheet(cell_counts, "proportions sample l2")
addWorksheet(cell_counts, "sample l3")
addWorksheet(cell_counts, "proportions sample l3")

# Write data to worksheets
writeData(cell_counts, "control l1", control.l1)
writeData(cell_counts, "proportions control l1", prop.control.l1)
writeData(cell_counts, "control l2", control.l2)
writeData(cell_counts, "proportions control l2", prop.control.l2)
writeData(cell_counts, "control l3", control.l3)
writeData(cell_counts, "proportions control l3", prop.control.l3)

writeData(cell_counts, "sample l1", sample.count.l1)
writeData(cell_counts, "proportions sample l1", sample.prop.count.l1)
writeData(cell_counts, "sample l2", sample.count.l2)
writeData(cell_counts, "proportions sample l2", sample.prop.count.l2)
writeData(cell_counts, "sample l3", sample.count.l3)
writeData(cell_counts, "proportions sample l3", sample.prop.count.l3)

# Save workbook object to xlsx file
saveWorkbook(cell_counts, "E:/Gennaro_sc/cell_amounts.xlsx")


# Find markers across days between cell types l1
Idents(hNB13) <- "days.post.infection"

mono_markers <- FindAllMarkers(
  subset(hNB13, subset = predicted.celltype.l1 == "Mono"),
  min.pct = 0.25,
  test.use = "MAST"
)

CD8_markers <- FindAllMarkers(
  subset(hNB13, subset = predicted.celltype.l1 == "CD8 T"),
  min.pct = 0.25,
  test.use = "MAST"
)

CD4_markers <- FindAllMarkers(
  subset(hNB13, subset = predicted.celltype.l1 == "CD4 T"),
  min.pct = 0.25,
  test.use = "MAST"
)

NK_markers <- FindAllMarkers(
  subset(hNB13, subset = predicted.celltype.l1 == "NK"),
  min.pct = 0.25,
  test.use = "MAST"
)

DC_markers <- FindAllMarkers(
  subset(hNB13, subset = predicted.celltype.l1 == "DC"),
  min.pct = 0.25,
  test.use = "MAST"
)

otherT_markers <- FindAllMarkers(
  subset(hNB13, subset = predicted.celltype.l1 == "other T"),
  min.pct = 0.25,
  test.use = "MAST"
)

# Output markers across days between cell types l1 into excel

# JUST REALIZED B CELLS ARE MISSING, TO LOOK INTO DATA MORE

l1_markers_across_days <-
  createWorkbook("Markers Across days between cell types l1")

# Add all the worksheets by name
addWorksheet(l1_markers_across_days, "Monocyte markers")
addWorksheet(l1_markers_across_days, "CD8 markers")
addWorksheet(l1_markers_across_days, "NK markers")
addWorksheet(l1_markers_across_days, "CD4 markers")
addWorksheet(l1_markers_across_days, "other T cell markers")
addWorksheet(l1_markers_across_days, "DC markers")

# Add data to worksheets
writeData(l1_markers_across_days, "Monocyte markers", mono_markers)
writeData(l1_markers_across_days, "CD8 markers", CD8_markers)
writeData(l1_markers_across_days, "NK markers", NK_markers)
writeData(l1_markers_across_days, "CD4 markers", CD4_markers)
writeData(l1_markers_across_days, "other T cell markers", otherT_markers)
writeData(l1_markers_across_days, "DC markers", DC_markers)

# Save workbook object to xlsx file
saveWorkbook(
  l1_markers_across_days,
  "E:/Gennaro_sc/markers_across_days_by_l1_cell_type.xlsx"
)


# An automated way to do this for l2
Idents(hNB13) <- "predicted.celltype.l2"
l2_idents <- levels(Idents(hNB13))
Idents(hNB13) <- "days.post.infection"

for (i in 1:28) {
  assign(
    paste0(l2_idents[i], "_markers"),
    FindAllMarkers(
      subset(hNB13, subset = predicted.celltype.l2 == l2_idents[i]),
      min.pct = 0.25,
      test.use = "MAST"
    )
  )
}

l2_markers_across_days <- createWorkbook("Markers Across days between cell types l2")

for (i in 1:28) {
  addWorksheet(l2_markers_across_days, paste0(l2_idents[i], " markers"))
  writeData(l2_markers_across_days,
            paste0(l2_idents[i], " markers"),
            get(paste0(l2_idents[i], "_markers")))
}

saveWorkbook(
  l2_markers_across_days,
  "E:/Gennaro_sc/markers_across_days_by_l2_cell_type.xlsx"
)

# Same marker identification done for l3
Idents(hNB13) <- "predicted.celltype.l3"
l3_idents <- levels(Idents(hNB13))
Idents(hNB13) <- "days.post.infection"

for (i in 1:28) {
  assign(
    paste0(l3_idents[i], "_markers"),
    FindAllMarkers(
      subset(hNB13, subset = predicted.celltype.l3 == l3_idents[i]),
      min.pct = 0.25,
      test.use = "MAST"
    )
  )
}

l3_markers_across_days <- createWorkbook("Markers Across days between cell types l3")

for (i in 1:28) {
  addWorksheet(l3_markers_across_days, paste0(l3_idents[i], " markers"))
}

for (i in 1:28) {
  writeData(l3_markers_across_days,
            paste0(l3_idents[i], " markers"),
            get(paste0(l3_idents[i], "_markers")))
}

saveWorkbook(
  l3_markers_across_days,
  "E:/Gennaro_sc/markers_across_days_by_l3_cell_type.xlsx"
)

# Make a workbook object via openxlsx
blish_cell_counts <- createWorkbook("Blish data cell count")

# Huge for loop which does everything for finding cell types and
# outputing the data on them into an xlsx
for (i in 1:3) {
  Idents(vent_blish) <- paste0("predicted.celltype.l", i)
  assign(paste0("blish.count.l", i),
         table(Idents(vent_blish), vent_blish$Donor.full))
  assign(paste0("blish.prop.count.l", i),
         prop.table(table(
           Idents(vent_blish), vent_blish$Donor.full
         ), margin = 2))
  addWorksheet(blish_cell_counts, paste0("blish counts l", i))
  addWorksheet(blish_cell_counts, paste0("blish cell proportions l", i))
  writeData(blish_cell_counts, paste0("blish counts l", i), 
            get(paste0("blish.count.l", i)))
  writeData(blish_cell_counts, paste0("blish cell proportions l", i),
            get(paste0("blish.prop.count.l", i)))
}

# Save workbook object to xlsx file
saveWorkbook(blish_cell_counts, "E:/Gennaro_sc/Blish_cell_counts.xlsx")


# Finding cell counts for original cell typing
# Make a workbook object via openxlsx
blish_cell_counts_orig <-
  createWorkbook("Blish data cell count original typing")

# Add all the worksheets by name

Idents(vent_blish) <- "cell.type.coarse"
blish.count.coarse <- table(Idents(vent_blish), vent_blish$Donor.full)
blish.prop.count.coarse <-
  prop.table(table(Idents(vent_blish), vent_blish$Donor.full), margin = 2)
Idents(vent_blish) <- "cell.type.fine"
blish.count.fine <- table(Idents(vent_blish), vent_blish$Donor.full)
blish.prop.count.fine <-
  prop.table(table(Idents(vent_blish), vent_blish$Donor.full), margin = 2)

addWorksheet(blish_cell_counts_orig, "blish orig count, coarse")
addWorksheet(blish_cell_counts_orig, "blish orig proportions, coarse")
addWorksheet(blish_cell_counts_orig, "blish orig count, fine")
addWorksheet(blish_cell_counts_orig, "blish orig proportions, fine")

writeData(blish_cell_counts_orig,
          "blish orig count, coarse",
          blish.count.coarse)
writeData(blish_cell_counts_orig,
          "blish orig proportions, coarse",
          blish.prop.count.coarse)
writeData(blish_cell_counts_orig,
          "blish orig count, fine",
          blish.count.fine)
writeData(blish_cell_counts_orig,
          "blish orig proportions, fine",
          blish.prop.count.fine)

# Save workbook object to xlsx file
saveWorkbook(blish_cell_counts_orig,
             "E:/Gennaro_sc/Blish_cell_counts_orig.xlsx")


# Markers for other cell types for sample 2
Idents(hNB13)<-"days.post.infection"

CD8_2markers <- FindMarkers(
  subset(hNB13, subset = predicted.celltype.l1 == "CD8 T"),
  ident.1 = 2,
  min.pct = 0.25,
  test.use = "MAST"
)

CD4_2markers <- FindMarkers(
  subset(hNB13, subset = predicted.celltype.l1 == "CD4 T"),
  ident.1 = 2,
  min.pct = 0.25,
  test.use = "MAST"
)

NK_2markers <- FindMarkers(
  subset(hNB13, subset = predicted.celltype.l1 == "NK"),
  ident.1 = 2,
  min.pct = 0.25,
  test.use = "MAST"
)


otherT_2markers <- FindMarkers(
  subset(hNB13, subset = predicted.celltype.l1 == "other T"),
  ident.1 = 2,
  min.pct = 0.25,
  test.use = "MAST"
)

write.csv(CD8_2markers, file = "E:/Gennaro_sc/CD8_sample2_markers.csv")
write.csv(CD4_2markers, file = "E:/Gennaro_sc/CD4_sample2_markers.csv")
write.csv(NK_2markers, file = "E:/Gennaro_sc/NK_sample2_markers.csv")
write.csv(otherT_2markers, file = "E:/Gennaro_sc/otherT_sample2_markers.csv")


######################################
# 5. Figures from paper:
######################################


library(Seurat)
library(SeuratDisk)
library(ggplot2)

if (!require(devtools)) {
  install.packages("devtools")
}

devtools::install_github("elliefewings/DoMultiBarHeatmap")
library(DoMultiBarHeatmap)
# Couldn't save images from this


devtools::install_github("xmc811/Scillus", ref = "development")
library(Scillus)

# You'll need to change this to the file path on your machine
load("E:/Research/Gennaro_sc/hNB13_no_tcr_analyzed.Rdata")

# Saved svg with default settings height at 825 pixels for large and 525 for small
# Did this using RStudio, but this can be done in R using 
DimPlot(hNB13_no_tcr, group.by = 'Sample')
DimPlot(hNB13_no_tcr, group.by = 'predicted.celltype.l1')
DimPlot(hNB13_no_tcr, group.by = 'predicted.celltype.l2')
DimPlot(hNB13_no_tcr, group.by = 'Sample') + NoLegend()
DimPlot(hNB13_no_tcr, group.by = 'predicted.celltype.l1') + NoLegend()
DimPlot(hNB13_no_tcr, group.by = 'predicted.celltype.l2') + NoLegend()


# Heatmaps split by cell 
DoHeatmap(subset(hNB13_no_tcr, subset = predicted.celltype.l1 == "Mono"), 
          features = c('IFI27','IFITM1','IFITM3','IFI6','ISG15','MT2A','IFI44L','LY6E','MX1'),
          group.by = 'Sample', combine = T) + scale_fill_gradientn(colors = c("green", "black", "red"))
DoHeatmap(subset(hNB13_no_tcr, subset = predicted.celltype.l1 == "CD4 T"), 
          features = c('IFI27','IFITM1','IFITM3','IFI6','ISG15','MT2A','IFI44L','LY6E','MX1'),
          group.by = 'Sample', combine = T) + scale_fill_gradientn(colors = c("green", "black", "red"))
DoHeatmap(subset(hNB13_no_tcr, subset = predicted.celltype.l1 == "CD8 T"), 
          features = c('IFI27','IFITM1','IFITM3','IFI6','ISG15','MT2A','IFI44L','LY6E','MX1'),
          group.by = 'Sample', combine = T) + scale_fill_gradientn(colors = c("green", "black", "red"))
DoHeatmap(subset(hNB13_no_tcr, subset = predicted.celltype.l1 == "NK"), 
          features = c('IFI27','IFITM1','IFITM3','IFI6','ISG15','MT2A','IFI44L','LY6E','MX1'),
          group.by = 'Sample', combine = T) + scale_fill_gradientn(colors = c("green", "black", "red"))

# Heatmaps combined per group
DoHeatmap(subset(hNB13_no_tcr, subset = predicted.celltype.l1 == "Mono"), 
          features = c('IFI27','IFITM1','IFITM3','IFI6','ISG15','MT2A','IFI44L','LY6E','MX1'),
          group.by = 'Sample', combine = T) + scale_fill_gradientn(colors = c("green", "black", "red"))
DoHeatmap(subset(hNB13_no_tcr, subset = predicted.celltype.l1 == "CD4 T"), 
          features = c('IFI27','IFITM1','IFITM3','IFI6','ISG15','MT2A','IFI44L','LY6E','MX1'),
          group.by = 'Sample', combine = T) + scale_fill_gradientn(colors = c("green", "black", "red"))
DoHeatmap(subset(hNB13_no_tcr, subset = predicted.celltype.l1 == "CD8 T"), 
          features = c('IFI27','IFITM1','IFITM3','IFI6','ISG15','MT2A','IFI44L','LY6E','MX1'),
          group.by = 'Sample', combine = T) + scale_fill_gradientn(colors = c("green", "black", "red"))
DoHeatmap(subset(hNB13_no_tcr, subset = predicted.celltype.l1 == "NK"), 
          features = c('IFI27','IFITM1','IFITM3','IFI6','ISG15','MT2A','IFI44L','LY6E','MX1'),
          group.by = 'Sample', combine = T) + scale_fill_gradientn(colors = c("green", "black", "red"))

ISG_markers = c('IFI27','IFITM1','IFITM3','IFI6','ISG15','MT2A','IFI44L','LY6E','MX1')
library(scales)

plot_heatmap(dataset = subset(hNB13_no_tcr, subset = predicted.celltype.l1 == c("CD8 T","Mono","CD4 T","NK")), 
             markers = ISG_markers,
             sort_var = c("predicted.celltype.l1","Sample"),
             anno_var = c("predicted.celltype.l1","Sample"),
             anno_colors = list('Set1',hue_pal()(5)),# These are the default colors of Seurat
             hm_limit = c(-2,0,2),
             hm_colors = c("green","black","red"))

