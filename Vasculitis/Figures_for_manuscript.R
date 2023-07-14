
setwd("/home/rebee/projects/mice-transcriptomics/bin")
# read in coldata
coldata <- read.csv("illumina.coldata.csv")
coldata<-coldata[!coldata$Tissue=='Blood',]

# change class of columns to factors as opposed to integers
sapply(coldata,class)
coldata$Day<-as.factor(coldata$Day)
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
library(stringr)

#import GTF file
gtf_file <- ('/home/rebee/references/mouse_2021/gencode.vM27.annotation.gtf')
cat(readLines(gtf_file, n=10), ssep = "\n")


#read in GTF
library(rtracklayer)
gtf_content <- rtracklayer::import(gtf_file, feature.type = 'gene')
gtf_content
annotation <- data.frame(elementMetadata(gtf_content), stringsAsFactors = FALSE)
row.names(annotation) <- annotation$gene_id
annotation$gene_id
# import counts from salmon
library(tximport)
annotation_transcript <- elementMetadata(import(gtf_file, feature.type = 'transcript'))
tx2gene <- annotation_transcript[,c("transcript_id","gene_id")]
head(tx2gene)
quant <- list.files(list.dirs(path = "../data/illumina/quants/gencode_lungs/", full.names = TRUE, recursive = FALSE), pattern = "quant.sf", full.names = TRUE)
all(file.exists(quant))
quant
id <- list.dirs(path = "../data/illumina/quants/gencode_lungs/", recursive = FALSE, full.names = FALSE)
id <- as.data.frame(id)
id
coldata<-coldata[order(match(coldata$Sample, id$id)),]



volcano_plot<- function(results_table, fc_threshold, fdr_threshold,log=FALSE, label, ylim, xlim){
  results_table <- results_table %>% filter(!grepl('Gm|Rik|ENS', gene_name))
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


##### same comparisons made for illumina data!
#remove vax groups



coldata2<-coldata[!coldata$Virus=='FluMIST',]
coldata2<-coldata2[!coldata2$Virus=='IAV vax',]
coldata2<-coldata2[!coldata2$Virus=='FluMIST + SARS-CoV-2',]
coldata2<-coldata2[!coldata2$Virus=='Coinfection',]
coldata2<-coldata2[!coldata2$Contrasts=='IAV_d10_Lung',]


coldata2$Sample
coldata2$Day<-c("6","7","7","6","7","7","7","7","7","7","6","6","3","3","3","3")

quant2<-quant[!str_detect(quant,pattern="FluMist")]
quant2<-quant2[!str_detect(quant2,pattern="IAV_SC2")]
quant2<-quant2[!str_detect(quant2,pattern = "IAV_alone_hACE2_D7")]
quant2

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
head(filtered_salmon_counts)

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

nrow(dgList)

dgList <- dgList[,which(!dgList$samples$Sample == "Sample_31-2256FluMist_21_12_2020_FluMist_Alone_D3_RL_2")]
dgList <- dgList[,which(!dgList$samples$Sample == "Sample_16-1723IAV_alone_hACE2_D7_RL_4")]


barplot(dgList$samples$lib.size,names=colnames(dgList),las=2)
logcounts <- cpm(dgList,log=TRUE)
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (normalised)")
plotMDS(dgList)


plotMDS(
  dgList,
  gene.selection = 'pairwise',
  col = as.integer(dgList$samples$group),
  labels = dgList$samples$Label)

plotMDS(
  dgList,
  gene.selection = 'common',
  col = as.integer(dgList$samples$group),
  labels = dgList$samples$Sample)



library(ggplot2)
library(ggrepel)




design <- model.matrix(~0+ group, data = dgList$samples)

levels(dgList$samples$group)

design
dgGlm <- estimateDisp(dgList, design, robust = TRUE)
plotBCV(dgGlm)

#Fitting data to the model
fit <- glmQLFit(dgGlm, design, robust = TRUE)

design



my.contrasts <- makeContrasts(
  IAVD3VMock = groupIAV_d6_Lung-groupMock_d10_Lung,         
  SARS2D3VMock = groupSARS_Cov_2_d6_Lung-groupMock_d10_Lung,         
  SARS2D7VMock = groupSARS_Cov_2_d10_Lung-groupMock_d10_Lung, 
  levels=design)
my.contrasts


fit <- glmQLFit(dgGlm, design, robust = TRUE)

de <- glmQLFTest(fit, contrast=my.contrasts)
de.IAVD3VMock     <-    glmQLFTest(fit, contrast=my.contrasts[,"IAVD3VMock"])
de.SARS2D3VMock   <-  glmQLFTest(fit, contrast=my.contrasts[,"SARS2D3VMock"])
de.SARS2D7VMock   <-  glmQLFTest(fit, contrast=my.contrasts[,"SARS2D7VMock"])



summary(de)

de$dispersion

top_genes <- topTags(de, n = 50)
top_genes


top_genes.IAVD3VMock    <- topTags(de.IAVD3VMock, n = 75)
top_genes.SARS2D3VMock  <- topTags(de.SARS2D3VMock, n = 75)
top_genes.SARS2D7VMock  <- topTags(de.SARS2D7VMock, n = 75)


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
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
    main = title,
    cutree_cols = 5,
    cutree_rows = 2
  )
}


results <- topTags(de, n=nrow(dgList), sort.by='none')$table

library(biomaRt)
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
genes<-results$gene_id

ConvertedGenes <- getBM(filters = "ensembl_gene_id_version", 
                        attributes = c("ensembl_gene_id_version", "entrezgene_id", "description"),
                        values = genes, 
                        mart = ensembl)
ConvertedGenes
colnames(ConvertedGenes)[1] <- "gene_id"

colnames(ConvertedGenes)



results
results <- left_join(results, ConvertedGenes, by = "gene_id")
results


results.IAVD3VMock     <- topTags(de.IAVD3VMock    , n=nrow(dgList), sort.by = 'none')$table
results.SARS2D3VMock   <- topTags(de.SARS2D3VMock  , n=nrow(dgList), sort.by = 'none')$table
results.SARS2D7VMock   <- topTags(de.SARS2D7VMock  , n=nrow(dgList), sort.by = 'none')$table

results.IAVD3VMock    <- left_join(results.IAVD3VMock, ConvertedGenes, by = "gene_id")
results.SARS2D3VMock  <- left_join(results.SARS2D3VMock, ConvertedGenes, by = "gene_id")
results.SARS2D7VMock  <- left_join(results.SARS2D7VMock, ConvertedGenes, by = "gene_id")


write.csv(results.IAVD3VMock, "../vasculitis/IAV_DGE.csv")
write.csv(results.SARS2D3VMock, "../vasculitis/SC2_D3_DGE.csv")
write.csv(results.SARS2D7VMock, "../vasculitis/SC2_D7.DGE.csv")



get.de.genes <- function(results) {
  de.genes <- results[ which(results$logFC >= 2 | results$logFC <= -2),]
  de.genes <-de.genes[de.genes$FDR <=0.05,]
}

de.genes.IAVD3 <- get.de.genes(results.IAVD3VMock)

de.genes.SARSD3 <- get.de.genes(results.SARS2D3VMock)


de.genes.SARSD7 <- get.de.genes(results.SARS2D7VMock)



fdr_threshold <- 0.05
fc_threshold <- 2
fc_threshold2 <-2




IAVD3plot <- volcano_plot(results.IAVD3VMock, fc_threshold, fdr_threshold,log=TRUE, "IAV Day 6", 10,15) # Arguments are results table, fold change threshold,
SARSD3plot  <- volcano_plot(results.SARS2D3VMock, fc_threshold, fdr_threshold,log=TRUE, "SARS-CoV-2 Day 3",10,15) # Arguments are results table, fold change threshold,
SARSD7plot  <- volcano_plot(results.SARS2D7VMock, fc_threshold, fdr_threshold,log=TRUE, "SARS-CoV-2 Day 7",10,15) # Arguments are results table, fold change threshold,

library(cowplot)
plot_grid(IAVD3plot + theme(legend.position="none"),
          SARSD3plot + theme(legend.position="none"),
          SARSD7plot + theme(legend.position="none"), ncol = 2)


#enrichment plots

library(clusterProfiler)
library(org.Mm.eg.db)
library(GOSemSim)


assayed.genes <- results

get.up.de.genes <- function(results) {
  de.genes <- results[ which(results$logFC >= 2),]
}

get.down.de.genes <- function(results) {
  de.genes <- results[ which(results$logFC <= 2),]
}


de.genes.IAVD6.up     <-  get.up.de.genes(de.genes.IAVD3)
de.genes.IAVD6.down   <- get.down.de.genes(de.genes.IAVD3)


de.genes.SARSD3.up    <- get.up.de.genes(de.genes.SARSD3)
de.genes.SARSD3.down  <-get.down.de.genes(de.genes.SARSD3)

de.genes.SARSD7.up    <-get.up.de.genes(de.genes.SARSD7)
de.genes.SARSD7.down  <-get.down.de.genes(de.genes.SARSD7)


de.genes.IAVD6.up$group <-  "Up"
de.genes.IAVD6.down$group <-"Down"
de.genes.SARSD3.up$group <- "Up"
de.genes.SARSD3.down$group <- "Down"
de.genes.SARSD7.up$group <- "Up"
de.genes.SARSD7.down$group <- "Down"

de.genes.IAVD6.up$othergroup <-  "IAV Day 6"
de.genes.IAVD6.down$othergroup <-"IAV Day 6"
de.genes.SARSD3.up$othergroup <- "SARS-CoV-2 Day 3"
de.genes.SARSD3.down$othergroup <- "SARS-CoV-2 Day 3"
de.genes.SARSD7.up$othergroup <- "SARS-CoV-2 Day 7"
de.genes.SARSD7.down$othergroup <- "SARS-CoV-2 Day 7"

GO.List <- rbind(de.genes.IAVD6.up,de.genes.IAVD6.down,de.genes.SARSD3.up, de.genes.SARSD3.down, 
                 de.genes.SARSD7.up, de.genes.SARSD7.down)


compareBP<-compareCluster(entrezgene_id~group+othergroup, data=GO.List, fun = "enrichGO", 
                          OrgDb = "org.Mm.eg.db", ont="BP", pvalueCutoff = 0.01, qvalueCutoff  = 0.05, 
                          pAdjustMethod = "BH", universe = results$entrezgene_id, readable = TRUE)


simplifyBP <- clusterProfiler::simplify(compareBP, cutoff=0.7, by="qvalue", select_fun=min)

write.csv(simplifyBP@compareClusterResult, "../vasculitis/compareCluster_BP_output.csv")

brain_terms<-c("transmission of nerve impulse",
"cognition",
"regulation of sensory perception of pain",
"serotonin secretion",
"glial cell apoptotic process",
"glial cell activation",
"astrocyte development",
"neuron death",
"regulation of neuron death",
"synapse pruning",
"regulation of neuron apoptotic process",
"neuron apoptotic process",
"glial cell migration",
"negative regulation of glial cell apoptotic process",
"response to axon injury",
"immunological synapse formation",
"regulation of microglial cell activation",
"maintenance of location")


epithelial_terms<-c("cilium movement",
"cilium organization",
"microtubule bundle formation",
"cilium or flagellum-dependent cell motility",
"cilium-dependent cell motility",
"epithelial cilium movement involved in extracellular fluid movement",
"regulation of syncytium formation by plasma membrane fusion",
"epithelial cell migration",
"epithelium migration",
"syncytium formation by plasma membrane fusion",
"cell-cell fusion",
"cell killing")


epithelial_terms2<-c("cilium movement","cilium organization",
"cilium or flagellum-dependent cell motility","cilium-dependent cell motility",
"epithelial cilium movement involved in extracellular fluid movement","regulation of cilium beat frequency",
"respiratory system process","respiratory gaseous exchange by respiratory system",
"syncytium formation by plasma membrane fusion","cell-cell fusion",
"regulation of syncytium formation by plasma membrane fusion","syncytium formation by plasma membrane fusion",
"epithelial cell migration","epithelium migration",
"regulation of epithelial cell migration")


immune_terms<-c("cytokine-mediated signaling pathway",
"response to chemokine","tumor necrosis factor superfamily cytokine production",
"regulation of tumor necrosis factor superfamily cytokine production","response to interleukin-1",
"lymphocyte proliferation","negative regulation of T cell activation",
"leukocyte mediated cytotoxicity","production of molecular mediator of immune response",
"leukocyte activation involved in immune response","natural killer cell mediated cytotoxicity",
"viral entry into host cell","modulation by host of viral process",
"response to interferon-gamma","activation of immune response",
"lymphocyte mediated immunity","regulation of lymphocyte proliferation",
"regulation of leukocyte proliferation","leukocyte proliferation",
"T cell proliferation","natural killer cell mediated immunity",
"regulation of leukocyte mediated cytotoxicity","adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains",
"antigen processing and presentation","mononuclear cell differentiation",
"positive regulation of lymphocyte activation","inflammatory response to antigenic stimulus",
"cytokine production involved in immune response","regulation of antigen processing and presentation",
"cytokine production involved in inflammatory response","alpha-beta T cell activation",
"lymphocyte activation involved in immune response","regulation of cytokine production involved in immune response",
"T cell mediated immunity",
"positive regulation of B cell mediated immunity",
"positive regulation of immunoglobulin mediated immune response","negative regulation of viral genome replication",
"regulation of humoral immune response","response to virus",
"defense response to virus","response to interferon-beta",
"cellular response to interferon-beta","response to interferon-alpha",
"regulation of innate immune response","response to type I interferon",
"regulation of tumor necrosis factor production","tumor necrosis factor production",
"cellular response to interferon-alpha","regulation of response to cytokine stimulus",
"acute-phase response","negative regulation of immune system process",
"negative regulation of immune response","response to tumor necrosis factor",
"humoral immune response","organ or tissue specific immune response",
"regulation of T-helper 2 cell cytokine production","T-helper 2 cell cytokine production",
"regulation of immune effector process","regulation of adaptive immune response",
"regulation of type 2 immune response","CD4-positive, alpha-beta T cell activation",
"type 2 immune response","regulation of T cell activation",
"positive regulation of T cell activation")

vasculitis_terms<-c("leukocyte chemotaxis","neutrophil chemotaxis",
"regulation of phagocytosis","cytokine production involved in inflammatory response",
"positive regulation of apoptotic signaling pathway","negative regulation of myeloid cell apoptotic process",
"cell chemotaxis","leukocyte migration",
"regulation of cell-cell adhesion","cytokine-mediated signaling pathway",
"regulation of T cell activation","regulation of mononuclear cell migration",
"tumor necrosis factor production","negative regulation of cytokine production",
"regulation of tumor necrosis factor production","positive regulation of chemotaxis",
"response to tumor necrosis factor","acute inflammatory response",
"positive regulation of T cell activation","myeloid leukocyte activation",
"positive regulation of phagocytosis","leukocyte activation involved in inflammatory response",
"production of molecular mediator involved in inflammatory response",
"phagocytosis",
"negative regulation of cell motility","cell-substrate adhesion",
"positive regulation of macrophage activation","inflammatory cell apoptotic process",
"heterotypic cell-cell adhesion","homotypic cell-cell adhesion",
"lymphocyte chemotaxis","lymphocyte migration",
"response to chemokine","positive regulation of leukocyte migration",
"cellular response to interleukin-1","response to interleukin-1",
"cellular response to tumor necrosis factor","regulation of chemotaxis",
"positive regulation of cell-cell adhesion",
"regulation of inflammatory response","leukocyte cell-cell adhesion",
"macrophage activation")


vasculitis_terms2<-c("cytokine-mediated signaling pathway","response to interferon-gamma",
"leukocyte migration","leukocyte chemotaxis",
"regulation of innate immune response","tumor necrosis factor superfamily cytokine production","regulation of tumor necrosis factor superfamily cytokine production",
"response to chemokine","regulation of cell-cell adhesion",
"leukocyte cell-cell adhesion","negative regulation of immune response",
"myeloid leukocyte activation","positive regulation of leukocyte migration",
"leukocyte mediated cytotoxicity","regulation of leukocyte mediated cytotoxicity",
"negative regulation of cell activation","response to interleukin-1",
"cytokine production involved in immune response","regulation of cytokine production involved in immune response",
"humoral immune response","acute inflammatory response",
"response to tumor necrosis factor","mononuclear cell differentiation",
"positive regulation of chemotaxis","regulation of response to cytokine stimulus",
"production of molecular mediator involved in inflammatory response","cellular response to tumor necrosis factor",
"regulation of phagocytosis","leukocyte activation involved in inflammatory response",
"leukocyte activation involved in immune response","natural killer cell mediated cytotoxicity",
"positive regulation of lymphocyte activation","cytokine production involved in inflammatory response",
"regulation of leukocyte apoptotic process","positive regulation of macrophage activation",
"regulation of leukocyte degranulation","cell-substrate adhesion",
"cell adhesion mediated by integrin","heterotypic cell-cell adhesion",
"leukocyte degranulation","homotypic cell-cell adhesion",
"granulocyte chemotaxis","acute-phase response",
"cellular response to interleukin-1","regulation of leukocyte migration",
"negative regulation of cytokine production","regulation of T cell activation",
"negative regulation of T cell activation","cell-cell recognition",
"regulation of leukocyte tethering or rolling","cell killing",
"tumor necrosis factor production","regulation of tumor necrosis factor production",
"myeloid leukocyte differentiation","positive regulation of phagocytosis",
"phagocytosis","negative regulation of cell motility",
"regulation of macrophage activation")


apoptosis_terms<-c("positive regulation of lymphocyte apoptotic process","leukocyte apoptotic process",
                   "extrinsic apoptotic signaling pathway",
                   "regulation of leukocyte apoptotic process","myeloid cell apoptotic process",
                   "regulation of myeloid cell apoptotic process","regulation of endothelial cell apoptotic process",
                   "endothelial cell apoptotic process",
                   "apoptotic cell clearance")

simplifyBP@compareClusterResult$othergroup



library(ggpubr)

top5CP<-clusterProfiler::dotplot(simplifyBP, x = "group", showCategory = 3, color="qvalue", includeAll = TRUE) +
  xlab("")+
  facet_grid(~factor(othergroup, levels=c("IAV Day 6",
                                        
                                          "SARS-CoV-2 Day 3",
                                          "SARS-CoV-2 Day 7")))+
  theme_pubr()+
   guides(color = guide_colourbar(barwidth = 10, barheight = 1))


ggsave("../vasculitis/brain_terms_BP.tiff", device = "tiff", dpi=300, height=8,width = 11,units="in")

clusterProfiler::dotplot(simplifyBP, x = "group", showCategory = epithelial_terms, color="qvalue", includeAll = TRUE) +
  xlab("")+
  facet_grid(~factor(othergroup, levels=c("IAV Day 6",
                                         
                                          "SARS-CoV-2 Day 3",
                                          "SARS-CoV-2 Day 7")))+
  theme_pubr()+
  guides(color = guide_colourbar(barwidth = 10, barheight = 1))

ggsave("../vasculitis/epithelial_terms_BP.tiff", device = "tiff", dpi=300, height=13,width = 11,units="in")

EBP<-clusterProfiler::dotplot(simplifyBP, x = "group", showCategory = epithelial_terms2, color="qvalue", includeAll = TRUE) +
  xlab("")+
  facet_grid(~factor(othergroup, levels=c("IAV Day 6",
                                          
                                          "SARS-CoV-2 Day 3",
                                          "SARS-CoV-2 Day 7")))+
  theme_pubr()+
  guides(color = guide_colourbar(barwidth = 10, barheight = 1))

ggsave("../vasculitis/epithelial_terms2_BP.tiff", device = "tiff", dpi=300, height=13,width = 11,units="in")

clusterProfiler::dotplot(simplifyBP, x = "group", showCategory = immune_terms, color="qvalue", includeAll = TRUE) +
  xlab("")+
  facet_grid(~factor(othergroup, levels=c("IAV Day 6",
                                        
                                          "SARS-CoV-2 Day 3",
                                          "SARS-CoV-2 Day 7")))+
  theme_pubr()+
  guides(color = guide_colourbar(barwidth = 10, barheight = 1))

ggsave("../vasculitis/immune_terms_BP.tiff", device = "tiff", dpi=300, height=35,width = 15,units="in")
  
clusterProfiler::dotplot(simplifyBP, x = "group", showCategory = vasculitis_terms, color="qvalue", includeAll = TRUE) +
  xlab("")+
  facet_grid(~factor(othergroup, levels=c("IAV Day 6",
                                       
                                          "SARS-CoV-2 Day 3",
                                          "SARS-CoV-2 Day 7")))+
  theme_pubr()+
  guides(color = guide_colourbar(barwidth = 10, barheight = 1))

ggsave("../vasculitis/vasculitis_terms_BP.tiff", device = "tiff", dpi=300, height=25,width = 11,units="in")



clusterProfiler::dotplot(simplifyBP, x = "group", showCategory = vasculitis_terms2, color="qvalue", includeAll = TRUE) +
  xlab("")+
  facet_grid(~factor(othergroup, levels=c("IAV Day 6",
                                          
                                          "SARS-CoV-2 Day 3",
                                          "SARS-CoV-2 Day 7")))+
  theme_pubr()+
  guides(color = guide_colourbar(barwidth = 10, barheight = 1))

ggsave("../vasculitis/vasculitis_terms_BP2.tiff", device = "tiff", dpi=300, height=30,width = 11,units="in")


cnet1a<-cnetplot(simplifyBP,
         categorySize="pvalue",
         showCategory = epithelial_terms2,
         node_label="category",
         legend_n=2)
cnet1a<-cnet1a + scale_fill_brewer(palette = "Set1")

cnet1b<-cnetplot(simplifyBP,
                 categorySize="pvalue",
                 showCategory = epithelial_terms2,
                 node_label="gene",
                 legend_n=2)
cnet1b<-cnet1b + scale_fill_brewer(palette = "Set1")

legend_b<-get_legend(
  cnet1a + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)


net<-plot_grid(cnet1a + theme(legend.position = "none"),cnet1b+theme(legend.position = "none"), labels = c("B","C"))
cnet<-plot_grid(net,legend_b, ncol=1,rel_heights = c(1.5, .05))




ggsave("../vasculitis/epithelial_network.tiff", device="tiff", dpi=300,width=10, height=10)


cnet2<-cnetplot(simplifyBP,
               categorySize="pvalue",
               showCategory = vasculitis_terms2,
               node_label="category",
               legend_n=2)
cnet2<-cnet2 + scale_fill_brewer(palette = "Set1")

ggsave("../vasculitis/vasculitis_network.tiff", device="tiff", dpi=300,width=10, height=10)


cnet2b<-cnetplot(simplifyBP,
                categorySize="pvalue",
                showCategory = vasculitis_terms2,
                node_label="gene",
                legend_n=2)
cnet2b<-cnet2b + scale_fill_brewer(palette = "Set1")

legend_b<-get_legend(
  cnet3 + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

net<-plot_grid(cnet2+ theme(legend.position="none"),cnet2b+
                 theme(legend.position = "none"), ncol=2, labels="AUTO")

plot_grid(net,legend_b, ncol=1,rel_heights = c(1.5, .05))

ggsave("../vasculitis/vasculitis.networks.tiff", device="tiff", dpi =300)

library(cowplot)

plot_grid(EBP,cnet, labels = c("A"), ncol = 1, rel_heights = c(1,1))


ggsave("../vasculitis/Figure9.tiff",device="tiff",dpi=300,height=370,width = 210, units="mm")


cnet3<-cnetplot(simplifyBP,
                categorySize="pvalue",
                legend_n=2,              
                node_label = "category")
cnet3<-cnet3 + scale_fill_brewer(palette = "Set1")

cnet4<-cnetplot(simplifyBP,
                categorySize="pvalue",
                legend_n=2,
                node_label = "gene")
cnet4<-cnet4 + scale_fill_brewer(palette = "Set1")

legend_b <- get_legend(
  cnet3 + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

net<-plot_grid(cnet3+ theme(legend.position="none"),cnet4+
            theme(legend.position = "none"), ncol=2, labels=c("B","C"))

nets<-plot_grid(net,legend_b, ncol=1,rel_heights = c(1.8, .05))

nets

plot_grid(top5CP,nets, ncol = 1, rel_heights = c(1,1.2),labels = c("A"))

ggsave("../vasculitis/total.network.tiff", device="tiff", dpi =300, height=400,width = 350, units="mm")


# extra figure

terms <- c("- leukocyte migration", "leukocyte chemotaxis"
, "regulation of cell-cell adhesion"
, "leukocyte cell adhesion"
, "positive regulation of leukocyte migration", "positive regulation of chemotaxis"
, "cell-substrate adhesion", "cell adhesion mediated by integrin"
, "heterotypic cell-cell adhesion", "homotypic cell-cell adhesion"
, "granulocyte chemotaxis", "regulation of leukocyte migration"
, "regulation of leukocyte tethering and rolling"
, "negative regulation of cell motility")

fig10A<-clusterProfiler::dotplot(simplifyBP, x = "othergroup", showCategory = terms, color="qvalue", includeAll = TRUE) +
  xlab("")+

  theme_pubr()+
  guides(color = guide_colourbar(barwidth = 8, barheight = 1))

cnet10a<-cnetplot(simplifyBP,
                 categorySize="pvalue",
                 showCategory = terms,
                 node_label="category",
                 legend_n=2)
cnet10a<-cnet10a + scale_fill_brewer(palette = "Set1")

cnet10b<-cnetplot(simplifyBP,
                 categorySize="pvalue",
                 showCategory = terms,
                 node_label="gene",
                 legend_n=2)
cnet10b<-cnet10b + scale_fill_brewer(palette = "Set1")

legend_b<-get_legend(
  cnet10a + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)


net2<-plot_grid(cnet10a + theme(legend.position = "none"),cnet10b+theme(legend.position = "none"), labels = c("B","C"))
cnet2<-plot_grid(net2,legend_b, ncol=1,rel_heights = c(1.5, .05))

cnet2

plot_grid(fig10A,cnet2, labels = c("A"), ncol = 1, rel_heights = c(0.8,1))

ggsave("../vasculitis/Figure10.tiff",device="tiff",dpi=300,height=370,width = 210, units="mm")




# plot VCAM and ICAM

results %>% filter(gene_name=="Vcam1")
results %>% filter(gene_name=="Icam1")

library(ggplot2)
library(ggpubr)


plot_gene <- function(gene_id, inputDgList){
  
  expr <- log2(cpm(inputDgList))
  
  plot_data <- cbind(inputDgList$samples, expression = expr[gene_id,])
  plot_data <- plot_data[order(plot_data$group),]
  plot_data$sample <- factor(rownames(plot_data), levels=rownames(plot_data))
  
  p <- ggplot(plot_data, aes(x=group, y=expression+0.1, fill=group)) + 
    geom_boxplot()+
    geom_point(size=0.2)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + 
    scale_y_continuous(limits = c(2.5,12))+
    scale_x_discrete(limits=c("Mock_d10_Lung","IAV_d6_Lung","SARS_Cov_2_d6_Lung","SARS_Cov_2_d10_Lung"), labels=c("Mock", "IAV D6", "SARS-CoV-2 D3", "SARS-CoV-2 D7"))+
    scale_fill_grey()+
    xlab("")+ ylab("CPM (log2)")+
    ggtitle(paste(toupper(annotation$gene_name[match(gene_id, annotation$gene_id)])))+
      theme_classic() +
    theme(legend.position = "none")
  
  p_with_stats <- p + geom_signif(comparisons = list(c("Mock_d10_Lung", "IAV_d6_Lung"), 
                                                                   c("Mock_d10_Lung", "SARS_Cov_2_d6_Lung"), 
                                                                   c("Mock_d10_Lung", "SARS_Cov_2_d10_Lung"),
                                                                   c("IAV_d6_Lung","SARS_Cov_2_d6_Lung"),
                                                                   c("IAV_d6_Lung","SARS_Cov_2_d10_Lung")), 
                                  test="wilcox.test",
                                  map_signif_level=TRUE, 
                                  textsize = 3,
                                  tip_length = 0.01,
                                  vjust = 0,
                                  step_increase = 0.2)
  
  print(p_with_stats)
}


p1<-plot_gene("ENSMUSG00000027962.15", dgList)
p2<-plot_gene("ENSMUSG00000037405.9",dgList)

plot_grid(p1,p2)

ggsave("../vasculitis/icam-vcam-plot.tiff", device="tiff",dpi=300,width=10,height=5)

