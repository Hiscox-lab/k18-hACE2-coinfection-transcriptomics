#coinfection analysis
library(edgeR)
library(limma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)
library(edgeR)
library(dplyr)
library(reshape2)
library(biomaRt)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GOSemSim)
library(cowplot)
library(statmod)
library(biomaRt)
library(DESeq2)
library(ggpubr)
library(tximport)
library(stringr)

#setwd
setwd("/home/rebee/projects/mice-transcriptomics/bin")
# read in coldata``
coldata <- read.csv("illumina.coldata.csv")
coldata<-coldata[!coldata$Tissue=='Blood',]

# change class of columns to factors as opposed to integers
sapply(coldata,class)
coldata$Day<-as.factor(coldata$Day)

#import GTF file
gtf_file <- ('/home/rebee/references/mouse_2021/gencode.vM27.annotation.gtf')
cat(readLines(gtf_file, n=10), ssep = "\n")


#read in GTF
gtf_content <- rtracklayer::import(gtf_file, feature.type = 'gene')
gtf_content
annotation <- data.frame(elementMetadata(gtf_content), stringsAsFactors = FALSE)
row.names(annotation) <- annotation$gene_id
annotation$gene_id

# import counts from salmon
annotation_transcript <- elementMetadata(rtracklayer::import(gtf_file, feature.type = 'transcript'))
tx2gene <- annotation_transcript[,c("transcript_id","gene_id")]
head(tx2gene)
quant <- list.files(list.dirs(path = "../data/illumina/quants/gencode_lungs/", full.names = TRUE, recursive = FALSE), pattern = "quant.sf", full.names = TRUE)
all(file.exists(quant))
quant
id <- list.dirs(path = "../data/illumina/quants/gencode_lungs/", recursive = FALSE, full.names = FALSE)
id <- as.data.frame(id)
id
coldata<-coldata[order(match(coldata$Sample, id$id)),]
coldata2<-coldata[!coldata$Virus=='FluMIST',]
coldata2<-coldata2[!coldata2$Virus=='IAV vax',]
coldata2<-coldata2[!coldata2$Virus=='FluMIST + SARS-CoV-2',]

quant2<-quant[!str_detect(quant,pattern="FluMist")]
quant2
id2<-id$id[!str_detect(quant,pattern="FluMist")]

txi_gene <- tximport(quant2, type = 'salmon', tx2gene = tx2gene, ignoreTxVersion = FALSE, txOut = FALSE)
y <- DGEList(counts=txi_gene$counts)
y$samples

barplot(y$samples$lib.size,names=colnames(y),las=2)
logcounts <- cpm(y,log=TRUE)
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")
plotMDS(y)


#add row and column names to quant files
salmon_counts<-salmon_counts<-matrix(as.integer(txi_gene$counts),nrow = length(txi_gene$counts[,1]))
head(salmon_counts)
rownames(salmon_counts)<-rownames(txi_gene$counts)
use_genes<-rownames(annotation)[rownames(annotation) %in% rownames(salmon_counts)]
annontation_salmon<-annotation[use_genes,]
filtered_salmon_counts<-salmon_counts[use_genes,]
colnames(filtered_salmon_counts)<-id2
head(filtered_salmon_counts)

coldata2$Sample %in% colnames(filtered_salmon_counts)

# Differential gene expression with EdgeR ####
library(edgeR)
dgList<- DGEList(
  counts = filtered_salmon_counts,
  samples = coldata2,
  genes = annontation_salmon,
  group = coldata2$Contrasts
)
dgList
dgList$samples
dgList$samples$lib.size

#normalisation function 
keep <- edgeR::filterByExpr(dgList, dgList[["samples"]]$group)
sum(keep)
dgList <- dgList[keep, , keep.lib.sizes = FALSE]
dgList <- calcNormFactors(dgList)

nrow(dgList)
dgList <- dgList[rowSums(cpm(dgList)>=0.1) >= 4, , keep.lib.sizes=FALSE]
dgList <- calcNormFactors(dgList)


#assess library size
barplot(dgList$samples$lib.size,names=colnames(dgList),las=2)
logcounts <- cpm(dgList,log=TRUE)
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (normalised)")
plotMDS(dgList)


#remove samples with low library sizes
dgList <- dgList[,which(!dgList$samples$Sample == "Sample_16-1723IAV_alone_hACE2_D7_RL_4")]
dgList <- dgList[,which(!dgList$samples$Sample == "Sample_27-1734Ctrl_hACE2_D7_RL_2")]

#plot PCA plot
edgeR.DDS <- DESeqDataSetFromMatrix(countData = round(dgList$counts), colData = dgList$samples, design = ~0+ group)
transform.edgeR.DDS <- rlog(edgeR.DDS, blind = TRUE)
pcaData <- plotPCA(transform.edgeR.DDS, intgroup = c("Virus","Day"), returnData = TRUE, ntop = 1000)  
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = Virus, shape = Day)) +
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  coord_fixed()+
  theme_pubr()+
  scale_color_brewer(palette = "Set1")

ggsave("../plots/illumina/Figure10B.tiff", dpi=300, device="tiff",width=7, height=5, units = "in")

#model matrix
design <- model.matrix(~0+ group, data = dgList$samples)
levels(dgList$samples$group)
dgGlm <- estimateDisp(dgList, design, robust = TRUE)


#Fitting data to the model
fit <- glmQLFit(dgGlm, design, robust = TRUE)

#make contrasts
my.contrasts <- makeContrasts(
  IAVD3VMock = groupIAV_d6_Lung-groupMock_d10_Lung,         
  IAVD7VMock = groupIAV_d10_Lung-groupMock_d10_Lung,         
  SARS2D3VMock = groupSARS_Cov_2_d6_Lung-groupMock_d10_Lung,         
  SARS2D7VMock = groupSARS_Cov_2_d10_Lung-groupMock_d10_Lung, 
  CoinfD3vmock = groupCoinfection_d6_Lung-groupMock_d10_Lung,
  CoinfD7vmock = groupCoinfection_d10_Lung-groupMock_d10_Lung,
  CoinfD3vIAVD3 = (groupCoinfection_d6_Lung-groupMock_d10_Lung)-(groupIAV_d6_Lung-groupMock_d10_Lung),
  CoinfD7vIAVD7 = (groupCoinfection_d10_Lung-groupMock_d10_Lung)-(groupIAV_d10_Lung-groupMock_d10_Lung),
  CoinfD3vSARSD3 = (groupCoinfection_d6_Lung-groupMock_d10_Lung)-(groupSARS_Cov_2_d6_Lung-groupMock_d10_Lung),
  CoinfD7vSARSD7 = (groupCoinfection_d10_Lung-groupMock_d10_Lung)-(groupSARS_Cov_2_d10_Lung-groupMock_d10_Lung),
  levels=design)
my.contrasts


fit <- glmQLFit(dgGlm, design, robust = TRUE)

de <- glmQLFTest(fit, contrast=my.contrasts)
de.IAVD3VMock     <-    glmQLFTest(fit, contrast=my.contrasts[,"IAVD3VMock"])
de.IAVD7VMock     <-    glmQLFTest(fit, contrast=my.contrasts[,"IAVD7VMock"])
de.SARS2D3VMock   <-  glmQLFTest(fit, contrast=my.contrasts[,"SARS2D3VMock"])
de.SARS2D7VMock   <-  glmQLFTest(fit, contrast=my.contrasts[,"SARS2D7VMock"])
de.CoinfD3vmock   <-  glmQLFTest(fit, contrast=my.contrasts[,"CoinfD3vmock"])
de.CoinfD7vmock   <-  glmQLFTest(fit, contrast=my.contrasts[,"CoinfD7vmock"])
de.CoinfD3vIAVD3  <- glmQLFTest(fit, contrast=my.contrasts[,"CoinfD3vIAVD3"])
de.CoinfD7vIAVD7  <- glmQLFTest(fit, contrast=my.contrasts[,"CoinfD7vIAVD7"])
de.CoinfD3vSARSD3 <- glmQLFTest(fit, contrast=my.contrasts[,"CoinfD3vSARSD3"])
de.CoinfD7vSARSD7 <- glmQLFTest(fit, contrast=my.contrasts[,"CoinfD7vSARSD7"])

#plot heatmaps

plot_heatmap <- function(plot_genes, title='', inputDgList,sample_annotation,gene_annotation){
  
  library(pheatmap)
  
  
  # Pick sample variables to plot
  annotation_col<-inputDgList$samples[sample_annotation]
  
  # Take out any genes not in the matrix
  plot_genes <- plot_genes[plot_genes %in% rownames(inputDgList)]
  
  # Make an expression mratix by calculating cpm and pulling out the selected genes, 
  #ordering by the sample annotation
  
  expression <- cpm(inputDgList)[plot_genes,rownames(annotation_col)]
  plotmatrix <- log2(expression + 0.1)
  
  # replace Ensembl IDs with gene names
  rownames(plotmatrix) <- gene_annotation[rownames(plotmatrix), 'gene_name'] 
  
  grid::grid.newpage()
  pheatmap(
    plotmatrix,
    show_rownames = T,
    annotation_col = annotation_col,
    border_color = NA,
    legend = FALSE,
    cluster_cols = TRUE, # change o TRUE to get a dendrogram of samples
    cluster_rows = TRUE,
    show_colnames = FALSE,
    scale = 'row',
    #color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
    main = title,
    cutree_cols = 6,
    cutree_rows = 3,
    cex = 1
  )
}

top_genes <- topTags(de, n = 50)

tiff('../plots/illumina/Figure10A.tiff', res = 300, width = 8, height = 12, units = "in")
plot_heatmap(rownames(top_genes), inputDgList = dgList,
             sample_annotation = c("Day","Virus"),gene_annotation=annotation)
dev.off()


results <- topTags(de, n=nrow(dgList), sort.by='none')$table

#make biomart object to annotate results with gene descriptions and entrezID for downstream analysis
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
genes<-results$gene_id

ConvertedGenes <- getBM(filters = "ensembl_gene_id_version", 
                        attributes = c("ensembl_gene_id_version", "entrezgene_id", "description"),
                        values = genes, 
                        mart = ensembl)
ConvertedGenes
colnames(ConvertedGenes)[1] <- "gene_id"


#left join results dataframes with annotations
results <- left_join(results, ConvertedGenes, by = "gene_id")
results.IAVD3VMock     <- topTags(de.IAVD3VMock    , n=nrow(dgList), sort.by = 'none')$table
results.IAVD7VMock     <- topTags(de.IAVD7VMock    , n=nrow(dgList), sort.by = 'none')$table
results.SARS2D3VMock   <- topTags(de.SARS2D3VMock  , n=nrow(dgList), sort.by = 'none')$table
results.SARS2D7VMock   <- topTags(de.SARS2D7VMock  , n=nrow(dgList), sort.by = 'none')$table
results.CoinfD3vmock   <- topTags(de.CoinfD3vmock  , n=nrow(dgList), sort.by = 'none')$table
results.CoinfD7vmock   <- topTags(de.CoinfD7vmock  , n=nrow(dgList), sort.by = 'none')$table
results.CoinfD3vIAVD3  <- topTags(de.CoinfD3vIAVD3 , n=nrow(dgList), sort.by = 'none')$table
results.CoinfD7vIAVD7  <- topTags(de.CoinfD7vIAVD7 , n=nrow(dgList), sort.by = 'none')$table
results.CoinfD3vSARSD3 <- topTags(de.CoinfD3vSARSD3, n=nrow(dgList), sort.by = 'none')$table
results.CoinfD7vSARSD7 <- topTags(de.CoinfD7vSARSD7, n=nrow(dgList), sort.by = 'none')$table
results.IAVD3VMock     <- left_join(results.IAVD3VMock    , ConvertedGenes, by = "gene_id")
results.IAVD7VMock     <- left_join(results.IAVD7VMock    , ConvertedGenes, by = "gene_id")
results.SARS2D3VMock   <- left_join(results.SARS2D3VMock  , ConvertedGenes, by = "gene_id")
results.SARS2D7VMock   <- left_join(results.SARS2D7VMock  , ConvertedGenes, by = "gene_id")
results.CoinfD3vmock   <- left_join(results.CoinfD3vmock  , ConvertedGenes, by = "gene_id")
results.CoinfD7vmock   <- left_join(results.CoinfD7vmock  , ConvertedGenes, by = "gene_id")
results.CoinfD3vIAVD3  <- left_join(results.CoinfD3vIAVD3 , ConvertedGenes, by = "gene_id")
results.CoinfD7vIAVD7  <- left_join(results.CoinfD7vIAVD7 , ConvertedGenes, by = "gene_id")
results.CoinfD3vSARSD3 <- left_join(results.CoinfD3vSARSD3, ConvertedGenes, by = "gene_id")
results.CoinfD7vSARSD7 <- left_join(results.CoinfD7vSARSD7, ConvertedGenes, by = "gene_id")     



write.csv(results.IAVD3VMock     , "../results/illumina/results.IAVD3VMock.csv")
write.csv(results.IAVD7VMock     , "../results/illumina/results.IAVD7VMock.csv")
write.csv(results.SARS2D3VMock   , "../results/illumina/results.SARS2D3VMock.csv")
write.csv(results.SARS2D7VMock   , "../results/illumina/results.SARS2D7VMock.csv")
write.csv(results.CoinfD3vmock   , "../results/illumina/results.CoinfD3vmock.csv")
write.csv(results.CoinfD7vmock   , "../results/illumina/results.CoinfD7vmock.csv")
write.csv(results.CoinfD3vIAVD3  , "../results/illumina/results.CoinfD3vIAVD3.csv")
write.csv(results.CoinfD7vIAVD7  , "../results/illumina/results.CoinfD7vIAVD7.csv")
write.csv(results.CoinfD3vSARSD3 , "../results/illumina/results.CoinfD3vSARSD3.csv")
write.csv(results.CoinfD7vSARSD7 , "../results/illumina/results.CoinfD7vSARSD7.csv")


#simple get de gene functions
get.de.genes <- function(results) {
  de.genes <- results[ which(results$logFC >= 2 | results$logFC <= -2),]
  de.genes <-de.genes[de.genes$FDR <=0.05,]
}

de.genes.IAVD3 <- get.de.genes(results.IAVD3VMock)
write.csv(de.genes.IAVD3, file = '../results/illumina/IAVD3.sig.csv')

de.genes.IAVD10 <- get.de.genes(results.IAVD7VMock)
write.csv(de.genes.IAVD10, file = '../results/illumina/IAVD7.sig.csv')

de.genes.SARSD3 <- get.de.genes(results.SARS2D3VMock)
write.csv(de.genes.SARSD3, file = '../results/illumina/SARSD3.sig.csv')

de.genes.SARSD7 <- get.de.genes(results.SARS2D7VMock)
write.csv(de.genes.SARSD7, file = '../results/illumina/SARSD7.sig.csv')

de.genes.CoinfD3 <- get.de.genes(results.CoinfD3vmock)
write.csv(de.genes.CoinfD3, file = '../results/illumina/CoinfD3.sig.csv')

de.genes.CoinfD7 <- get.de.genes(results.CoinfD7vmock)
write.csv(de.genes.CoinfD7, file = '../results/illumina/CoinfD7.sig.csv')

de.genes.CoinfD3vIAVD3 <- get.de.genes(results.CoinfD3vIAVD3)
write.csv(de.genes.CoinfD3vIAVD3, file = '../results/illumina/CoinfD3vIAVD3.sig.csv')

de.genes.CoinfD7vIAVD7 <- get.de.genes(results.CoinfD7vIAVD7)
write.csv(de.genes.CoinfD7vIAVD7, file = '../results/illumina/CoinfD7vIAVD7.sig.csv')

de.genes.CoinfD3vSARSD3 <- get.de.genes(results.CoinfD3vSARSD3)
write.csv(de.genes.CoinfD3vSARSD3, file = '../results/illumina/CoinfD3vSARSD3.sig.csv')

de.genes.CoinfD7vSARSD7 <- get.de.genes(results.CoinfD7vSARSD7)
write.csv(de.genes.CoinfD7vSARSD7, file = '../results/illumina/CoinfD7vSARSD7.sig.csv')

#volcano plot function
fdr_threshold <- 0.05
fc_threshold <- 2
fc_threshold2 <-2

volcano_plot<- function(results_table, fc_threshold, fdr_threshold,log=FALSE, label, ylim, xlim){
  results_table <- results_table %>% filter(!grepl('Gm|Rik|ENS', gene_name)) #filter out unreadable geneIDs for plot labelling
  results_table <- results_table %>% distinct(gene_name, .keep_all = TRUE)
  results_table$significant <- 'no'
  results_table$significant[ abs(results_table$logFC) >= log2(fc_threshold) &
                               results_table$FDR <= fdr_threshold  ] <- 'yes'
  results_table$plot_gene<- ifelse(results_table$logFC >= results_table$logFC[order(results_table$logFC, decreasing = T)[7]] & results_table$FDR < 0.05 |
                                     results_table$logFC <= results_table$logFC[order(results_table$logFC, decreasing = F)][7] & results_table$FDR < 0.05, results_table$gene_name, NA)
  # Un-log the values if necessary
  xlab<-"log fold change"
  if(! log){
    results_table$logFC <- sign(results_table$logFC)*2**abs(results_table$logFC) 
    xlab<-"fold change"
  }
  ggplot(results_table, aes(logFC, -log10(FDR), color=significant, label = plot_gene )) + 
    geom_point(alpha = 0.2) +  # This alters the transparency of the points
    scale_colour_manual(name = 'significant', values = setNames(c('red','grey'),c('yes', 'no'))) +
    xlab(xlab) +
    ylim(0,ylim)+
    xlim(-xlim,xlim)+
    geom_text_repel(color = "black",
                    max.overlaps = 20,
                    #nudge_x = 0.1,
                    #box.padding = 0.5,
                    #nudge_y = 2,
                    #segment.curvature = -0.1,
                    #segment.ncp = 3,
                    #segment.angle = 20
    )+
    geom_hline(yintercept = 1.3, linetype = "dashed")+
    geom_vline(xintercept = 1, linetype = "dashed")+
    geom_vline(xintercept = -1, linetype = "dashed")+
    ggtitle(label)+
    theme_classic()
}

IAVD3plot <- volcano_plot(results.IAVD3VMock, fc_threshold, fdr_threshold,log=TRUE, "IAV Day 6", 10,15) # Arguments are results table, fold change threshold,
IAVD7plot  <- volcano_plot(results.IAVD7VMock, fc_threshold, fdr_threshold,log=TRUE, "IAV Day 10",10,15) # Arguments are results table, fold change threshold,
SARSD3plot  <- volcano_plot(results.SARS2D3VMock, fc_threshold, fdr_threshold,log=TRUE, "SARS-CoV-2 Day 6",10,15) # Arguments are results table, fold change threshold,
SARSD7plot  <- volcano_plot(results.SARS2D7VMock, fc_threshold, fdr_threshold,log=TRUE, "SARS-CoV-2 Day 10",10,15) # Arguments are results table, fold change threshold,
CoinfD3plot <- volcano_plot(results.CoinfD3vmock, fc_threshold, fdr_threshold,log=TRUE, "Coinfection Day 6",10,15) # Arguments are results table, fold change threshold,
CoinfD7plot <- volcano_plot(results.CoinfD7vmock, fc_threshold, fdr_threshold,log=TRUE, "Coinfection Day 10",10,15) # Arguments are results table, fold change threshold,

plot_grid(IAVD3plot + theme(legend.position="none"),
          IAVD7plot + theme(legend.position="none"),
          SARSD3plot + theme(legend.position="none"),
          SARSD7plot + theme(legend.position="none"),
          CoinfD3plot + theme(legend.position="none"),
          CoinfD7plot + theme(legend.position="none"), ncol = 2)

ggsave("../plots/illumina/Figure10C.tiff", device = "tiff", dpi = 300, width = 8, height = 8)

IAVD3vCoinfD3 <- volcano_plot(results.CoinfD3vIAVD3, fc_threshold, fdr_threshold,log=TRUE, "IAV Day 6 vs Coinfection Day 6",10,15) # 
IAVD7vCoinfD7 <- volcano_plot(results.CoinfD7vIAVD7, fc_threshold, fdr_threshold,log=TRUE, "IAV Day 10 vs Coinfection Day 10",10,15) # 
SARSD3vCoinfD3 <- volcano_plot(results.CoinfD3vSARSD3, fc_threshold, fdr_threshold,log=TRUE, "SARS-CoV-2 Day 6 vs Coinfection Day 6",10,15) # 
SARSD7vCoinfD <- volcano_plot(results.CoinfD7vSARSD7, fc_threshold, fdr_threshold,log=TRUE, "SARS-CoV-2 Day 10 vs Coinfection Day 10",10,15) # 

plot_grid(IAVD3vCoinfD3 + theme(legend.position="none"), 
          IAVD7vCoinfD7 + theme(legend.position="none"), 
          SARSD3vCoinfD3 + theme(legend.position="none"),
          SARSD7vCoinfD + theme(legend.position="none"), ncol =2 )
ggsave("../plots/illumina/Figure11.tiff", device = "tiff", dpi = 300, width = 9, height = 6)

#enrichment plots
assayed.genes <- results

#subset up and down de genes
get.up.de.genes <- function(results) {
  de.genes <- results[ which(results$logFC >= 2),]
}

get.down.de.genes <- function(results) {
  de.genes <- results[ which(results$logFC <= 2),]
}

de.genes.CoinfD3.up<-   get.up.de.genes(de.genes.CoinfD3)
de.genes.CoinfD3.down<- get.down.de.genes(de.genes.CoinfD3)

de.genes.CoinfD7.up   <- get.up.de.genes(de.genes.CoinfD7)
de.genes.CoinfD7.down  <- get.down.de.genes(de.genes.CoinfD7)

de.genes.IAVD6.up     <-  get.up.de.genes(de.genes.IAVD3)
de.genes.IAVD6.down   <- get.down.de.genes(de.genes.IAVD3)

de.genes.IAVD10.up    <- get.up.de.genes(de.genes.IAVD10)
de.genes.IAVD10.down  <-get.down.de.genes(de.genes.IAVD10)

de.genes.SARSD3.up    <- get.up.de.genes(de.genes.SARSD3)
de.genes.SARSD3.down  <-get.down.de.genes(de.genes.SARSD3)

de.genes.SARSD7.up    <-get.up.de.genes(de.genes.SARSD7)
de.genes.SARSD7.down  <-get.down.de.genes(de.genes.SARSD7)

de.genes.CoinfD7_SARSD7.up    <-get.up.de.genes(de.genes.CoinfD7vSARSD7)
de.genes.CoinfD7_SARSD7.down  <-get.down.de.genes(de.genes.CoinfD7vSARSD7)

de.genes.CoinfD3_SARSD3.up    <-get.up.de.genes(de.genes.CoinfD3vSARSD3)
de.genes.CoinfD3_SARSD3.down  <-get.down.de.genes(de.genes.CoinfD3vSARSD3)

de.genes.CoinfD7_fluD7.up    <-get.up.de.genes(de.genes.CoinfD7vIAVD7)
de.genes.CoinfD7_fluD7.down  <-get.down.de.genes(de.genes.CoinfD7vIAVD7)

de.genes.CoinfD3_fluD3.up    <-get.up.de.genes(de.genes.CoinfD3vIAVD3)
de.genes.CoinfD3_fluD3.down  <-get.down.de.genes(de.genes.CoinfD3vIAVD3)

colnam<-c("Group","Up/Down","IAV Day 6", "IAV Day 10", "SARS-CoV-2 Day 6","SARS-CoV-2 Day 10", "Coinfection Day 6", "Coinfection Day 10")
mock_up<-c("Mock","Up",nrow(de.genes.IAVD6.up),nrow(de.genes.IAVD10.up),nrow(de.genes.SARSD3.up),nrow(de.genes.SARSD7.up), nrow(de.genes.CoinfD3.up),nrow(de.genes.CoinfD7.up))
mock_down<-c("","Down",nrow(de.genes.IAVD6.down),nrow(de.genes.IAVD10.down),nrow(de.genes.SARSD3.down),nrow(de.genes.SARSD7.down), nrow(de.genes.CoinfD3.down),nrow(de.genes.CoinfD7.down))
coinfD6_up<-c("Coinfection Day 6","Up",nrow(de.genes.CoinfD3_fluD3.up),"-",nrow(de.genes.CoinfD3_SARSD3.up),"-","-","-")
coinfD6_down<-c("","Down",nrow(de.genes.CoinfD3_fluD3.down),"-",nrow(de.genes.CoinfD3_SARSD3.down),"-","-","-")
coinfD10_up<-c("Coinfection Day 10","Up","-",nrow(de.genes.CoinfD7_fluD7.up),"-",nrow(de.genes.CoinfD7_SARSD7.up),"-","-")
coinfD10_down<-c("","Down","-",nrow(de.genes.CoinfD7_fluD7.down),"-",nrow(de.genes.CoinfD7_SARSD7.down),"-","-")
                                                                                                                                              
table1<-as.data.frame(rbind(mock_up,mock_down,coinfD6_up,coinfD6_down,coinfD10_up,coinfD10_down))
colnames(table1)<-colnam

table1<-gt(table1)%>%
  tab_header(
    title = md("Count of differentially expressed genes"))

table1 %>% gtsave("../results/illumina/table_1.docx")

#add groups for labelling
de.genes.CoinfD3.up$group <- "Up"
de.genes.CoinfD3.down$group <- "Down"
de.genes.CoinfD7.up$group <-"Up"
de.genes.CoinfD7.down$group <-"Down"
de.genes.IAVD6.up$group <-  "Up"
de.genes.IAVD6.down$group <-"Down"
de.genes.IAVD10.up$group <- "Up"
de.genes.IAVD10.down$group <-"Down"
de.genes.SARSD3.up$group <- "Up"
de.genes.SARSD3.down$group <- "Down"
de.genes.SARSD7.up$group <- "Up"
de.genes.SARSD7.down$group <- "Down"

de.genes.CoinfD3.up$othergroup <- "Coinfection Day 6"
de.genes.CoinfD3.down$othergroup <- "Coinfection Day 6"
de.genes.CoinfD7.up$othergroup <-"Coinfection Day 10"
de.genes.CoinfD7.down$othergroup <-"Coinfection Day 10"
de.genes.IAVD6.up$othergroup <-  "IAV Day 6"
de.genes.IAVD6.down$othergroup <-"IAV Day 6"
de.genes.IAVD10.up$othergroup <- "IAV Day 10"
de.genes.IAVD10.down$othergroup <-"IAV Day 10"
de.genes.SARSD3.up$othergroup <- "SARS-CoV-2 Day 6"
de.genes.SARSD3.down$othergroup <- "SARS-CoV-2 Day 6"
de.genes.SARSD7.up$othergroup <- "SARS-CoV-2 Day 10"
de.genes.SARSD7.down$othergroup <- "SARS-CoV-2 Day 10"

#rbind dataframes for input into clusterProfiler
GO.List <- rbind(de.genes.CoinfD3.up, de.genes.CoinfD3.down,de.genes.CoinfD7.up,de.genes.CoinfD7.down,de.genes.IAVD6.up,de.genes.IAVD6.down,
                 de.genes.IAVD10.up, de.genes.IAVD10.down, de.genes.SARSD3.up, de.genes.SARSD3.down, de.genes.SARSD7.up, de.genes.SARSD7.down)

#run clusterProfiler
compareCC<-compareCluster(entrezgene_id~group+othergroup, data=GO.List, fun = "enrichGO", OrgDb = "org.Mm.eg.db", ont="CC", pvalueCutoff = 0.01, qvalueCutoff  = 0.05, pAdjustMethod = "BH", universe = results$entrezgene_id, readable=TRUE)
compareMF<-compareCluster(entrezgene_id~group+othergroup, data=GO.List, fun = "enrichGO", OrgDb = "org.Mm.eg.db", ont="MF", pvalueCutoff = 0.01, qvalueCutoff  = 0.05, pAdjustMethod = "BH", universe = results$entrezgene_id, readable=TRUE)
compareBP<-compareCluster(entrezgene_id~group+othergroup, data=GO.List, fun = "enrichGO", OrgDb = "org.Mm.eg.db", ont="BP", pvalueCutoff = 0.01, qvalueCutoff  = 0.05, pAdjustMethod = "BH", universe = results$entrezgene_id, readable=TRUE)

#simplify terms
simplifyCC <- clusterProfiler::simplify(compareCC, cutoff=0.7, by="qvalue", select_fun=min)
simplifyMF <- clusterProfiler::simplify(compareMF, cutoff=0.7, by="qvalue", select_fun=min)
simplifyBP <- clusterProfiler::simplify(compareBP, cutoff=0.7, by="qvalue", select_fun=min)

#create a vector of terms for visualisation
terms<-c("regulation of interleukin-1 production","chemokine production","acute inflammatory response", "cellular response to interferon-gamma", "regulation of innate immune response", "cellular response to interferon-beta",
         "response to virus","activation of immune response","response to tumor necrosis factor","leukocyte chemotaxis",
         "myeloid leukocyte migration","cellular response to chemokine","cytokine-mediated signaling pathway","leukocyte migration", "mitotic nuclear division","positive regulation of cell cycle process")



#plot biological processes
dotplot(simplifyBP, x = "othergroup", showCategory = terms, color="qvalue", includeAll = TRUE)+
  scale_x_discrete(limits=c("Coinfection Day 6",
                            "Coinfection Day 10",
                            "IAV Day 6",
                            "IAV Day 10",
                            "SARS-CoV-2 Day 6",
                            "SARS-CoV-2 Day 10"))+
  theme_pubr()+
  xlab("")+
  guides(color = guide_colourbar(barwidth = 12, barheight = 1))

ggsave("../plots/illumina/Figure12.tiff", device = "tiff", dpi = 300, width = 12, height = 7)


#plot CC and MF terms for supplementary
dotplot(simplifyCC, x = "group", showCategory = 5, color="qvalue", includeAll = TRUE)+
  facet_grid(~factor(othergroup,levels = c("Coinfection Day 6",
                                "Coinfection Day 10",
                                "IAV Day 6",
                                "IAV Day 10",
                                "SARS-CoV-2 Day 6",
                                "SARS-CoV-2 Day 10")))+
  theme_pubr()+
  xlab("")+
  guides(color = guide_colourbar(barwidth = 12, barheight = 1))

  ggsave("../plots/illumina/Supplementary_Figure5.tiff", device = "tiff", dpi = 300, width = 12, height = 14)

dotplot(simplifyMF, x = "group", showCategory = 5, color="qvalue", includeAll = TRUE)+
  facet_grid(~factor(othergroup,levels = c("Coinfection Day 6",
                                           "Coinfection Day 10",
                                           "IAV Day 6",
                                           "IAV Day 10",
                                           "SARS-CoV-2 Day 6",
                                           "SARS-CoV-2 Day 10")))+
  theme_pubr()+
  xlab("")+
  guides(color = guide_colourbar(barwidth = 12, barheight = 1))

ggsave("../plots/illumina/Supplementary_Figure6.tiff", device = "tiff", dpi = 300, width = 12, height = 14)

write.csv(as.data.frame(simplifyBP@compareClusterResult), "../results/illumina/gene-ontology-BP.csv")
write.csv(as.data.frame(simplifyCC@compareClusterResult), "../results/illumina/gene-ontology-CC.csv")
write.csv(as.data.frame(simplifyMF@compareClusterResult), "../results/illumina/gene-ontology-MF.csv")









