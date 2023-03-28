setwd("/home/rebee/projects/mice-transcriptomics/bin")
# read in coldata and remove groups not included in analysis
coldata <- read.csv("illumina.coldata.csv")
coldata<-coldata[!coldata$Tissue=='Blood',]
coldata2<-coldata[!coldata$Virus=='FluMIST',]
coldata2<-coldata2[!coldata2$Virus=='IAV vax',]
coldata2<-coldata2[!coldata2$Virus=='FluMIST + SARS-CoV-2',]
coldata2<-coldata2[!coldata2$Virus=='Coinfection',]
coldata2<-coldata2[!coldata2$Virus=='IAV',]

# change class of columns to factors as opposed to integers
sapply(coldata2,class)
coldata2$Day<-as.factor(coldata2$Day)
coldata2$Virus<-as.factor(coldata2$Virus)
coldata2$Contrasts<-as.factor(coldata2$Contrasts)


# load libraries
library(edgeR)
library(limma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)
library(edgeR)
library(tidyverse)
library(reshape2)
library(biomaRt)
library(ggplot2)
library(EnrichmentBrowser)
library("AnnotationDbi")
library("org.Mm.eg.db")
library("RColorBrewer")
library(DESeq2)


#read in annotation file for edgeR analysis later on in script

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

#read in featureCounts output

countfiles <- list.files(list.dirs(path = "../../cambridge-upr-project/data/featureCounts/", full.names = TRUE), pattern = ".txt", full.names = TRUE)
countfiles<-countfiles[!str_detect(countfiles,pattern="summary")]
counts <- lapply(countfiles, function(cf){
  read.delim(cf, header = TRUE, skip = 1, row.names = 1)[,6,drop=FALSE]
})

countfiles
head(counts)

#make a dataframe of sample IDs and then match to order in coldata2
id <- list.files(path = "../../cambridge-upr-project/data/featureCounts/", full.names = FALSE)
id<-id[!str_detect(id,pattern="summary")]
id<-str_remove_all(id, ".txt")
id <- as.data.frame(id)
id
coldata2<-coldata2[order(match(coldata2$Sample, id$id)),]

# get rownames for counts
names <- sort(Reduce(union, lapply(counts, rownames)))
#convert list of counts into dataframe
raw_counts <- do.call(cbind, lapply(counts, function(x) x[names,,drop=FALSE]))
head(raw_counts)
colnames(raw_counts)
coldata2$Sample
colnames(raw_counts)<-c("Sample_18-1725SC2_alone_hACE2_D7_RL_1","Sample_19-1726SC2_alone_hACE2_D7_RL_2",
"Sample_20-1727SC2_alone_hACE2_D7_RL_3", "Sample_21-1728SC2_alone_hACE2_D7_RL_4", "Sample_26-1733Ctrl_hACE2_D7_RL_1", 
                        "Sample_27-1734Ctrl_hACE2_D7_RL_2", "Sample_28-1735Ctrl_hACE2_D7_RL_3", 
                        "Sample_29-1736Ctrl_hACE2_D7_RL_4" , "Sample_5-1712SC2_alone_hACE2_D3_RL_1", 
                        "Sample_6-1713SC2_alone_hACE2_D3_RL_2", "Sample_7-1714SC2_alone_hACE2_D3_RL_3",
                        "Sample_8-1715SC2_alone_hACE2_D3_RL_4")

#remove the extension of the gene ID to allow for downstream analysis                        
rownames(raw_counts) <- gsub("\\..*","",rownames(raw_counts))


#download geneset for enrichment analysis

go.gs <- getGenesets(org = "mmu", db = "go", onto = "BP", mode = "biomart")

#drop D6 samples from counts and coldata and check colnames
raw_countsD3<-raw_counts[-c(1:4)]
colnames(raw_countsD3)
coldataD3<-coldata2 %>% filter(Contrasts!="SARS_Cov_2_d10_Lung")


#make DESeq object
dds_D3 <- DESeqDataSetFromMatrix(raw_countsD3, 
                              colData = coldataD3,
                              design = ~Contrasts)

#filter low counts
keep <- rowSums(counts(dds_D3)) >= 10  # >=10 because recommended in the vignette
dds_D3 <- dds_D3[keep,]

dds_D3$Virus 

resultsNames(dds_D3)

#Do DESeq and make results object
dds_D3 <- DESeq(dds_D3)
res_D3 <- results(dds_D3)

res_D3

#import for enrichmentBrowser
se_D3 <- import(dds_D3, res_D3)
se_D3 <- normalize(se_D3, norm.method = "vst")

boxplot(assay(se_D3, "raw"))
boxplot(assay(se_D3, "norm"))

#convert ENSEMBL ID to ENTREZ for camera analysis
SE_D3 <- idMap(se_D3, org = "mmu", from = "ENSEMBL", to = "ENTREZID")

#conduct camera analysis and extract significant results
camera_D3<-sbea(method = "camera", se = SE_D3, gs = go.gs)
camera.sig_D3<-gsRanking(camera_D3, signif.only = TRUE)

#convert to dataframe and filter terms associated with protein folding responses
BPdf_D3<-as.data.frame(camera.sig_D3)


UPR_GO_D3<-BPdf_D3 %>% filter(str_detect(GENE.SET, 'chaperone|glycosylation|protein folding|unfolded|response to stress'))


#edgeR analysis as well as DESeq

# (1) import from edgeR (RNA-seq count data)
# (1a) create the expression data object
library(edgeR)
d <- DGEList(counts = raw_countsD3,
             samples = coldataD3,
             genes = annotation,
             group = coldataD3$Contrasts)
d <- calcNormFactors(d)
rdesign<- model.matrix(~0+ group, data = d$samples)
d <- estimateDisp(d, rdesign)

# (1b) obtain differential expression results 
fit <- glmQLFit(d, rdesign)
qlf <- glmQLFTest(fit)
res <- topTags(qlf, n = nrow(d), sort.by = "none")

# (1c) import
se <- import(d, res)

se <- normalize(se, norm.method = "vst")

boxplot(assay(se, "raw"))
boxplot(assay(se, "norm"))


SE <- idMap(se, org = "mmu", from = "ENSEMBL", to = "ENTREZID")

camera<-sbea(method = "camera", se = SE, gs = go.gs)
camera.sig_D3_edgeR<-gsRanking(camera, signif.only = TRUE)


BPdf_D3<-as.data.frame(camera.sig_D3_edgeR)

UPR_GO_D3_edgeR<-BPdf_D3 %>% filter(str_detect(GENE.SET, 'chaperone|glycosylation|protein folding|unfolded|response to stress'))


#filter data for second timepoint
coldata2$Day <- as.factor(coldata2$Day)
head(coldata2)
head(raw_counts)

raw_countsD7<-raw_counts[-c(9:12)]
coldataD7<-coldata2 %>% filter(Contrasts!="SARS_Cov_2_d6_Lung")

dds_D7 <- DESeqDataSetFromMatrix(raw_countsD7, 
                              colData = coldataD7,
                              design = ~Contrasts)
dds_D7 <- DESeq(dds_D7)
res_D7 <- results(dds_D7)

se_D7 <- import(dds_D7, res_D7)
se_D7 <- normalize(se_D7, norm.method = "vst")

boxplot(assay(se_D7, "raw"))
boxplot(assay(se_D7, "norm"))


SE_D7 <- idMap(se_D7, org = "mmu", from = "ENSEMBL", to = "ENTREZID")

camera_D7<-sbea(method = "camera", se = SE_D7, gs = go.gs)
camera.sig_D7<-gsRanking(camera_D7, signif.only = TRUE)


BPdf_D7<-as.data.frame(camera.sig_D7)

UPR_GO_D7<-BPdf_D7 %>% filter(str_detect(GENE.SET, 'chaperone|glycosylation|protein folding|unfolded|response to stress'))

write.csv(UPR_GO_D3, "../../cambridge-upr-project/data/UPR-camera-D3.csv")
write.csv(UPR_GO_D7, "../../cambridge-upr-project/data/UPR-camera-D7.csv")



#edgeR

# (1) import from edgeR (RNA-seq count data)
# (1a) create the expression data object
library(edgeR)
d <- DGEList(counts = raw_countsD7,
             samples = coldataD7,
             genes = annotation,
             group = coldataD7$Contrasts)
d <- calcNormFactors(d)
rdesign<- model.matrix(~0+ group, data = d$samples)
d <- estimateDisp(d, rdesign)

# (1b) obtain differential expression results 
fit <- glmQLFit(d, rdesign)
qlf <- glmQLFTest(fit)
res <- topTags(qlf, n = nrow(d), sort.by = "none")

# (1c) import
se <- import(d, res)

se <- normalize(se, norm.method = "vst")

boxplot(assay(se, "raw"))
boxplot(assay(se, "norm"))


SE <- idMap(se, org = "mmu", from = "ENSEMBL", to = "ENTREZID")

camera<-sbea(method = "camera", se = SE, gs = go.gs, browse = FALSE)
camera.sig_D7_edgeR<-gsRanking(camera, signif.only = TRUE)


BPdf_D7<-as.data.frame(camera.sig_D7_edgeR)

UPR_GO_D7_edgeR<-BPdf_D7 %>% filter(str_detect(GENE.SET, 'chaperone|glycosylation|protein folding|unfolded|response to stress'))


write.csv(UPR_GO_D3_edgeR, "../../cambridge-upr-project/data/UPR-camera-D3-edgeR.csv")
write.csv(UPR_GO_D7_edgeR, "../../cambridge-upr-project/data/UPR-camera-D7-edgeR.csv")


#DGE analysis in DESeq2

dds <- DESeqDataSetFromMatrix(raw_counts, 
                              colData = coldata2,
                             design = ~Contrasts)



keep <- rowSums(counts(dds)) >= 10  # >=10 because recommended in the vignette
dds <- dds[keep,]
dds 


dds$Contrasts <- factor(dds$Contrasts, levels = c("Mock_d10_Lung","SARS_Cov_2_d6_Lung","SARS_Cov_2_d10_Lung"))

dds$Contrasts

dds <- DESeq(dds)

res_D3_a <- results(dds, contrast = c("Contrasts", "Mock_d10_Lung", "SARS_Cov_2_d10_Lung"))
summary(res_D3_a)
res_D7_a <- results(dds, contrast = c("Contrasts", "Mock_d10_Lung",  "SARS_Cov_2_d6_Lung"))
summary(res_D7_a)


#add gene symbols
ens.str <- rownames(res_D3_a)
res_D3_a$symbol <- mapIds(org.Mm.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res_D3_a$entrez <- mapIds(org.Mm.eg.db,
                     keys=ens.str,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
res_D7_a$symbol <- mapIds(org.Mm.eg.db,
                          keys=ens.str,
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")

ens.str <- rownames(dds)
dds$symbol <- mapIds(org.Mm.eg.db,
                          keys=ens.str,
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")
dds$entrez <- mapIds(org.Mm.eg.db,
                          keys=ens.str,
                          column="ENTREZID",
                          keytype="ENSEMBL",
                          multiVals="first")


library(apeglm)
resultsNames(dds)

#shrink lfc
resLFC_D7 <- lfcShrink(dds, coef="Contrasts_SARS_Cov_2_d10_Lung_vs_Mock_d10_Lung", type="apeglm")
resLFC_D7 


resLFC_D7$symbol <- mapIds(org.Mm.eg.db,
                     keys=rownames(resLFC_D7),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

LFCD7<-as.data.frame(resLFC_D7)

write.csv(LFCD7,"../../cambridge-upr-project/data/LFC_D7.csv")

resLFC_D3 <- lfcShrink(dds, coef="Contrasts_SARS_Cov_2_d6_Lung_vs_Mock_d10_Lung", type="apeglm")
resLFC_D3

resLFC_D3$symbol <- mapIds(org.Mm.eg.db,
                           keys=rownames(resLFC_D3),
                           column="SYMBOL",
                           keytype="ENSEMBL",
                           multiVals="first")


LFCD3<-as.data.frame(resLFC_D3)

write.csv(LFCD3, "../../cambridge-upr-project/data/LFC_D3.csv")


library(EnhancedVolcano)

EnhancedVolcano(resLFC_D3,
                lab = resLFC_D3$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                title = "SARS-CoV-2 infection in k18-hACE2",
                subtitle = "Day 3 vs Mock")

ggsave("../../cambridge-upr-project/plots/enhancedvolcano_padj_D3.png", device = "png", width=8,height=5,units = "in", dpi=300)

EnhancedVolcano(resLFC_D7,
                lab = resLFC_D7$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                title = "SARS-CoV-2 infection in k18-hACE2",
                subtitle = "Day 7 vs Mock")


ggsave("../../cambridge-upr-project/plots/enhancedvolcano_padj_D7.png", device = "png", width=8,height=5,units = "in", dpi=300)

#rlog tranformations

rld <- rlog(dds, blind=FALSE)
rld

pcaData <- plotPCA(rld, intgroup=c("Contrasts"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Contrasts)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  theme_pubr()+
  scale_color_brewer(palette="Set1",labels=c('Mock', 'SARS-CoV-2 Day 3', 'SARS-CoV-2 Day 7'))
  
ggsave("../../cambridge-upr-project/plots/PCA_plot.tiff", device = "tiff", dpi=300,width=7,height=6,units="in")

