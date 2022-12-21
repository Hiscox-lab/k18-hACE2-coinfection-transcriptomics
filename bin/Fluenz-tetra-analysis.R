library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)
library(edgeR)
library(tidyverse)
library(reshape2)
library(biomaRt)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GOSemSim)
library(cowplot)
library(statmod)
library(biomaRt)
library(gt)
library(tidyverse)
library(glue)


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

txi_gene <- tximport(quant, type = 'salmon', tx2gene = tx2gene, ignoreTxVersion = FALSE, txOut = FALSE)
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
  samples = coldata,
  genes = annontation_salmon,
  group = coldata$Contrasts
)
dgList
dgList$samples
dgList$samples$lib.size

#normalisation function 
keep <- edgeR::filterByExpr(dgList, dgList[["samples"]]$group)
sum(keep)
dgList <- dgList[keep, , keep.lib.sizes = FALSE]
dgList <- calcNormFactors(dgList)
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
dgList <- dgList[,which(!dgList$samples$Sample == "Sample_31-2256FluMist_21_12_2020_FluMist_Alone_D3_RL_2")]
dgList <- dgList[,which(!dgList$samples$Sample == "Sample_16-1723IAV_alone_hACE2_D7_RL_4")]
dgList <- dgList[,which(!dgList$samples$Sample == "Sample_27-1734Ctrl_hACE2_D7_RL_2")]

#model matrix
design <- model.matrix(~0+ group, data = dgList$samples)
design
dgGlm <- estimateDisp(dgList, design, robust = TRUE)


#Fitting data to the model
fit <- glmQLFit(dgGlm, design, robust = TRUE)


#make contrasts
my.contrasts <- makeContrasts(
  d6fluenz_SARS_v_SARS = (groupSARS_Cov_2_d6_Lung-groupMock_d10_Lung)-(groupIAV_Vax_coinfection_d6_Lung-groupMock_d10_Lung),
  d6fluenz_SARS_v_Coinf= (groupCoinfection_d6_Lung-groupMock_d10_Lung)-(groupIAV_Vax_coinfection_d6_Lung-groupMock_d10_Lung),
  d10fluenz_SARS_v_SARS = (groupSARS_Cov_2_d10_Lung-groupMock_d10_Lung)-(groupIAV_Vax_coinfection_d10_Lung-groupMock_d10_Lung),
  d10fluenz_SARS_v_coinf=(groupCoinfection_d10_Lung-groupMock_d10_Lung)-(groupIAV_Vax_coinfection_d10_Lung-groupMock_d10_Lung),
  levels=design)

#de genes in one object
de <- glmQLFTest(fit, contrast=my.contrasts)

#de genes for each contrast
de.d6fluenz_SARS_v_SARS            <-    glmQLFTest(fit, contrast=my.contrasts[,"d6fluenz_SARS_v_SARS"])
de.d10fluenz_SARS_v_SARS            <-    glmQLFTest(fit, contrast=my.contrasts[,"d10fluenz_SARS_v_SARS"])

 de.d6fluenz_SARS_v_coinf            <-     glmQLFTest(fit, contrast=my.contrasts[,"d6fluenz_SARS_v_Coinf"])
de.d10fluenz_SARS_v_coinf            <-    glmQLFTest(fit, contrast=my.contrasts[,"d10fluenz_SARS_v_coinf"])


#make biomark object for annotation and left join to results dataframes 
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
genes<-results$gene_id

ConvertedGenes <- getBM(filters = "ensembl_gene_id_version", 
                        attributes = c("ensembl_gene_id_version", "entrezgene_id", "description"),
                        values = genes, 
                        mart = ensembl)
ConvertedGenes
colnames(ConvertedGenes)[1] <- "gene_id"

results <- topTags(de, n=nrow(dgList), sort.by='none')$table
results <- left_join(results, ConvertedGenes, by = "gene_id")


results.d6fluenz_SARS_v_SARS   <- topTags(de.d6fluenz_SARS_v_SARS, n=nrow(dgList), sort.by = 'none')$table
results.d10fluenz_SARS_v_SARS   <- topTags(de.d10fluenz_SARS_v_SARS, n=nrow(dgList), sort.by = 'none')$table

results.d6fluenz_SARS_v_coinf   <-  topTags(de.d6fluenz_SARS_v_coinf, n=nrow(dgList), sort.by = 'none')$table
results.d10fluenz_SARS_v_coinf   <- topTags(de.d10fluenz_SARS_v_coinf, n=nrow(dgList), sort.by = 'none')$table


results.d6fluenz_SARS_v_SARS      <-left_join(results.d6fluenz_SARS_v_SARS, ConvertedGenes, by = "gene_id")
results.d10fluenz_SARS_v_SARS       <-left_join(results.d10fluenz_SARS_v_SARS, ConvertedGenes, by = "gene_id")

results.d6fluenz_SARS_v_coinf       <-left_join(results.d6fluenz_SARS_v_coinf, ConvertedGenes, by = "gene_id")
results.d10fluenz_SARS_v_coinf      <-left_join(results.d10fluenz_SARS_v_coinf, ConvertedGenes, by = "gene_id")

#write de genes to csv
write.csv(results.d6fluenz_SARS_v_SARS     , file='../results/illumina/differential_results.d6fluenz_SARS_v_SARS.csv   ')
write.csv(results.d10fluenz_SARS_v_SARS      , file='../results/illumina/differential_results.d10fluenz_SARS_v_SARS.csv    ')

 write.csv(results.d6fluenz_SARS_v_coinf     ,   file='../results/illumina/differential_results.d6fluenz_SARS_v_coinf.csv   ')
write.csv(results.d10fluenz_SARS_v_coinf      , file='../results/illumina/differential_results.d10fluenz_SARS_v_coinf.csv    ')

#basic get de genes function to subset significant genes only and then write to csv
get.de.genes <- function(results) {
  de.genes <- results[ which(results$logFC >= 2 | results$logFC <= -2),]
  de.genes <-de.genes[de.genes$FDR <=0.05,]
}

de.genes.d6fluenz_SARS_v_SARS <- get.de.genes(results.d6fluenz_SARS_v_SARS)
write.csv(de.genes.d6fluenz_SARS_v_SARS, file = '../results/illumina/d6fluenz_SARS_v_SARS.sig.csv')

de.genes.d10fluenz_SARS_v_SARS <- get.de.genes(results.d10fluenz_SARS_v_SARS)
write.csv(de.genes.d10fluenz_SARS_v_SARS, file = '../results/illumina/d10fluenz_SARS_v_SARS.sig.csv')


de.genes.d6fluenz_SARS_v_coinf <- get.de.genes(results.d6fluenz_SARS_v_coinf)
write.csv(de.genes.d6fluenz_SARS_v_coinf, file = '../results/illumina/d6fluenz_SARS_v_coinf.sig.csv')

de.genes.d10fluenz_SARS_v_coinf <- get.de.genes(results.d10fluenz_SARS_v_coinf)
write.csv(de.genes.d10fluenz_SARS_v_coinf, file = '../results/illumina/d10fluenz_SARS_v_coinf.sig.csv')


# set thresholds for volcano plots
fdr_threshold <- 0.05
fc_threshold <- 2
fc_threshold2 <-2

#volcano plot function
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


p9<-volcano_plot(results.d6fluenz_SARS_v_SARS, fc_threshold, fdr_threshold,log=TRUE, "Fluenz Tetra + SARS-CoV-2 vs SARS-CoV-2 Day 6", 5,15) # Arguments are results table, fold change threshold,
p10<-volcano_plot(results.d10fluenz_SARS_v_SARS, fc_threshold, fdr_threshold,log=TRUE, "Fluenz Tetra + SARS-CoV-2 vs SARS-CoV-2 Day 10",5,15) # Arguments are results table, fold change threshold,

p11<-volcano_plot(results.d6fluenz_SARS_v_coinf, fc_threshold, fdr_threshold,log=TRUE, "Fluenz Tetra + SARS-CoV-2 vs Coinfection Day 6", 10,15) # Arguments are results table, fold change threshold,
p12<-volcano_plot(results.d10fluenz_SARS_v_coinf, fc_threshold, fdr_threshold,log=TRUE, "Fluenz Tetra + SARS-CoV-2 vs Coinfection Day 10",10,15) # Arguments are results table, fold change threshold,


plot_grid(p9+theme(legend.position="none"),
          p10+theme(legend.position="none"),
          p11 + theme(legend.position="none"),
          p12 + theme(legend.position="none"))

ggsave("../plots/illumina/Figure13.tiff", device="tiff", dpi=300,width = 10,height = 8, units="in")

#make a df with universe genes for later using in clusterProfiler as universe
assayed.genes <- results

#get up and down regulated gene functions
get.up.de.genes <- function(results) {
  de.genes <- results[ which(results$logFC >= 2),]
}

get.down.de.genes <- function(results) {
  de.genes <- results[ which(results$logFC <= 2),]
}


de.genes.d6fluenz_SARS_v_SARS.up          <- get.up.de.genes(de.genes.d6fluenz_SARS_v_SARS)
de.genes.d6fluenz_SARS_v_SARS.down        <-get.down.de.genes(de.genes.d6fluenz_SARS_v_SARS)
de.genes.d10fluenz_SARS_v_SARS.up          <- get.up.de.genes(de.genes.d10fluenz_SARS_v_SARS)
de.genes.d10fluenz_SARS_v_SARS.down        <-get.down.de.genes(de.genes.d10fluenz_SARS_v_SARS)


de.genes.d6fluenz_SARS_v_coinf.up            <- get.up.de.genes(de.genes.d6fluenz_SARS_v_coinf)
de.genes.d6fluenz_SARS_v_coinf.down         <-get.down.de.genes(de.genes.d6fluenz_SARS_v_coinf)
de.genes.d10fluenz_SARS_v_coinf.up          <- get.up.de.genes(de.genes.d10fluenz_SARS_v_coinf)
de.genes.d10fluenz_SARS_v_coinf.down       <-get.down.de.genes(de.genes.d10fluenz_SARS_v_coinf)




colnam<-c("Group","Up/Down","SARS-CoV-2 Day 6 vs Mock", "SARS-CoV-2 Day 10 vs Mock","Coinfection Day 6 vs Mock","Coinfection Day 10 vs Mock")
D6_up<-c("Fluenz Tetra + SARS-CoV-2 Day 6 vs Mock","Up",nrow(de.genes.d6fluenz_SARS_v_SARS.up),"-",nrow(de.genes.d6fluenz_SARS_v_coinf.up),"-")
D6_down<-c(" ","Down",nrow(de.genes.d6fluenz_SARS_v_SARS.down),"-",nrow(de.genes.d6fluenz_SARS_v_coinf.down),"-")
D10_up<-c("Fluenz Tetra + SARS-CoV-2 Day 10 vs Mock","Up","-",nrow(de.genes.d10fluenz_SARS_v_SARS.up),"-",nrow(de.genes.d10fluenz_SARS_v_coinf.up))
D10_down<-c(" ","Down","-",nrow(de.genes.d10fluenz_SARS_v_SARS.down),"-",nrow(de.genes.d10fluenz_SARS_v_coinf.down))


table2<-as.data.frame(rbind(D6_up,D6_down,D10_up,D10_down))
colnames(table2)<-colnam

table2<-gt(table2)%>%
  tab_header(
    title = md("Count of differentially expressed genes"))

table2

table2 %>% gtsave("../results/illumina/table_2.docx")





#add in new column with descriptions before rbinding into one dataframe
de.genes.d6fluenz_SARS_v_SARS.up$group   <- "Up"
de.genes.d6fluenz_SARS_v_SARS.down$group <- "Down"
de.genes.d10fluenz_SARS_v_SARS.up$group  <- "Up"
de.genes.d10fluenz_SARS_v_SARS.down$group <- "Down"


 de.genes.d6fluenz_SARS_v_coinf.up$group   <- "Up"
 de.genes.d6fluenz_SARS_v_coinf.down$group <- "Down"
de.genes.d10fluenz_SARS_v_coinf.up$group  <- "Up"
de.genes.d10fluenz_SARS_v_coinf.down$group <- "Down"



de.genes.d6fluenz_SARS_v_SARS.up$othergroup   <- "Fluenz Tetra + SARS-CoV-2 vs SARS-CoV-2 only"
de.genes.d6fluenz_SARS_v_SARS.down$othergroup <- "Fluenz Tetra + SARS-CoV-2 vs SARS-CoV-2 only"
de.genes.d10fluenz_SARS_v_SARS.up$othergroup  <- "Fluenz Tetra + SARS-CoV-2 vs SARS-CoV-2 only"
de.genes.d10fluenz_SARS_v_SARS.down$othergroup <- "Fluenz Tetra + SARS-CoV-2 vs SARS-CoV-2 only"

 de.genes.d6fluenz_SARS_v_coinf.up$othergroup   <- "Fluenz Tetra + SARS-CoV-2 vs Coinfection"
 de.genes.d6fluenz_SARS_v_coinf.down$othergroup <- "Fluenz Tetra + SARS-CoV-2 vs Coinfection"
de.genes.d10fluenz_SARS_v_coinf.up$othergroup  <-  "Fluenz Tetra + SARS-CoV-2 vs Coinfection"
de.genes.d10fluenz_SARS_v_coinf.down$othergroup <- "Fluenz Tetra + SARS-CoV-2 vs Coinfection"



de.genes.d6fluenz_SARS_v_SARS.up$day    <- "Day 6"
de.genes.d6fluenz_SARS_v_SARS.down$day   <- "Day 6"
de.genes.d10fluenz_SARS_v_SARS.up$day   <- "Day 10"
de.genes.d10fluenz_SARS_v_SARS.down$day <- "Day 10"

 de.genes.d6fluenz_SARS_v_coinf.up$day    <- "Day 6"
 de.genes.d6fluenz_SARS_v_coinf.down$day   <- "Day 6"
de.genes.d10fluenz_SARS_v_coinf.up$day   <- "Day 10"
de.genes.d10fluenz_SARS_v_coinf.down$day <- "Day 10"


#rbind de genes with identifiers to prepare input for clusterProfiler
GO.List2<-rbind(de.genes.d6fluenz_SARS_v_SARS.up,
de.genes.d6fluenz_SARS_v_SARS.down,
de.genes.d10fluenz_SARS_v_SARS.up,
de.genes.d10fluenz_SARS_v_SARS.down, 
de.genes.d6fluenz_SARS_v_coinf.up,  
de.genes.d6fluenz_SARS_v_coinf.down,
de.genes.d10fluenz_SARS_v_coinf.up,  
de.genes.d10fluenz_SARS_v_coinf.down)

#run enrichment analysis
compareCC<-compareCluster(entrezgene_id~group+othergroup+day, data=GO.List2, fun = "enrichGO", OrgDb = "org.Mm.eg.db", ont="CC", pvalueCutoff = 0.01, qvalueCutoff  = 0.05, pAdjustMethod = "BH", universe = results$entrezgene_id, readable=TRUE)
compareMF<-compareCluster(entrezgene_id~group+othergroup+day, data=GO.List2, fun = "enrichGO", OrgDb = "org.Mm.eg.db", ont="MF", pvalueCutoff = 0.01, qvalueCutoff  = 0.05, pAdjustMethod = "BH", universe = results$entrezgene_id, readable=TRUE)
compareBP<-compareCluster(entrezgene_id~group+othergroup+day, data=GO.List2, fun = "enrichGO", OrgDb = "org.Mm.eg.db", ont="BP", pvalueCutoff = 0.01, qvalueCutoff  = 0.05, pAdjustMethod = "BH", universe = results$entrezgene_id, readable=TRUE)

#remove similar terms 
simplifyCC <- clusterProfiler::simplify(compareCC, cutoff=0.7, by="qvalue", select_fun=min)
simplifyMF <- clusterProfiler::simplify(compareMF, cutoff=0.7, by="qvalue", select_fun=min)
simplifyBP <- clusterProfiler::simplify(compareBP, cutoff=0.7, by="qvalue", select_fun=min)


#create vector of terms to visualise
FTterms<-c("response to virus","regulation of T cell activation", "leukocyte aggregation", "leukocyte migration involved in inflammatory response",
           "leukocyte chemotaxis", "cilium movement", "cilium organization","cytokine-mediated signaling pathway","regulation of innate immune response","skin development","intermediate filament organization", "cell chemotaxis",
           "microtubule bundle formation", "respiratory system process","homeostasis of number of cells","lymphocyte proliferation")



#plot biological processes
dotplot(simplifyBP, x = "group", showCategory = FTterms, color="qvalue", includeAll = TRUE) + ggplot2::facet_grid(~factor(othergroup)+factor(day, levels=c('Day 6', 'Day 10')), 
                                                                                                                 labeller = label_wrap_gen())+
  theme_pubr()+
  xlab("")

ggsave("../plots/illumina/Figure15.tiff", device = "tiff", dpi=300,height=9,width=11,units="in")


# generate table for publication to show dge genes for fluenztetra

suptab1<-de.genes.d6fluenz_SARS_v_SARS %>%
  select(gene_name, logFC,FDR, description) %>%
  mutate(across(logFC, round, 2)) %>%
  arrange(desc(logFC))%>%
  rename("gene name"=gene_name)%>%
  mutate_at(c('description'), ~replace_na(.,""))%>%
  gt()%>%
  tab_header(
    title = md("Differentially expressed genes"),
    subtitle = "Fluenz tetra and SARS-CoV-2 infection vs SARS-CoV-2 only infection at day 6"
  )

suptab1
suptab2<-de.genes.d10fluenz_SARS_v_SARS %>%
  select(gene_name, logFC,FDR, description) %>%
  mutate(across(logFC, round, 2)) %>%
  arrange(desc(logFC))%>%
  rename("gene name"=gene_name)%>%
  mutate_at(c('description'), ~replace_na(.,""))%>%
  gt()%>%
  tab_header(
    title = md("Differentially expressed genes"),
    subtitle = "Fluenz tetra and SARS-CoV-2 infection vs SARS-CoV-2 only infection at day 10"
  )

suptab1 %>% gtsave("../results/illumina/supplementary_table_1.docx")
suptab2 %>% gtsave("../results/illumina/supplementary_table_2.docx")



