setwd("/home/rebee/projects/mice-transcriptomics/bin")
# read in coldata
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


library(DESeq2)

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

#eaBrowse(camera.sig_D7_edgeR, out.dir="../../cambridge-upr-project/Camera/D7", report.name="camera-edgeR-D7")

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


res_D7_a
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
LFCD3 %>% filter(rownames(LFCD3) %in% genes2)

res_D3_ordered <- res_D3_a[order(res_D3_a$pvalue),]
res_D3_ordered

res_D7_ordered <- res_D7_a[order(res_D7_a$pvalue),]
res_D7_ordered

library(EnhancedVolcano)


genes<-c("Xbp1","Hspa5","Hspa8","Calr","Pdia4",
         "Tor1b","Sec61b","Hsp90ab1","Asns",
         "Atf4","Chop","Gadd34","Ppp1r15a","Ddit3")

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


summary(res_D3_a)

sum(res_D3_a$padj < 0.1, na.rm=TRUE)
sum(res_D7_a$padj < 0.1, na.rm=TRUE)

res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)
res05

#MA plots

plotMA(res_D3_a, ylim=c(-2,2))
plotMA(res_D7_a, ylim=c(-2,2))

plotMA(resLFC_D3, ylim=c(-2,2))
plotMA(resLFC_D7, ylim=c(-2,2))



#rlog tranformations


rld <- rlog(dds, blind=FALSE)

meanSdPlot(assay(rld))

library(pheatmap)


select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Virus","Day")])
rld

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, show_colnames = FALSE, cutree_cols = 2)


#make heatmap of genes
D3LFC<-data.frame()
D3LFC<-data.frame(resLFC_D3[,'log2FoldChange'])
D3LFC
colnames(D3LFC)<-"3"
rownames(D3LFC)<-rownames(res_D3_a)

D6LFC<-data.frame()
D6LFC<-data.frame(resLFC_D7[,'log2FoldChange'])
D6LFC
colnames(D6LFC)<-"7"
rownames(D6LFC)<-rownames(res_D7_a)

genes2<-mapIds(org.Mm.eg.db,
               keys=genes,
               column="ENSEMBL",
               keytype="SYMBOL",
               multiVals="first")

D6UPRLFC<-D6LFC %>% filter(rownames(D6LFC) %in% genes2)
D3UPRLFC<-D3LFC %>% filter(rownames(D3LFC) %in% genes2)

D3UPRLFC1 <- D3UPRLFC %>% add_rownames()
D6UPRLFC1 <- D6UPRLFC %>% add_rownames()
LFC <- full_join(D3UPRLFC1,D6UPRLFC1)

LFC

LFC$symbol <- mapIds(org.Mm.eg.db,
                 keys=LFC$rowname,
                 column="SYMBOL",
                 keytype="ENSEMBL",
                 multiVals="first")

LFC

plot<-gather(LFC, day, Log2FC, 2:3)

library(ggpubr)
ggplot(plot)+
  geom_tile(aes(x=day,y=symbol,fill=Log2FC))+
  xlab("dpi")+
  ylab("")+
  scale_fill_gradient(low = "black", high = "red", limits =c(-2,2))+
  theme_pubclean()

ggsave("../../cambridge-upr-project/plots/logFC_plot.png",device = "png", dpi=300, height = 7,width=5,units = "in")

sampleDists <- dist(t(assay(vsd)))


library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


pcaData <- plotPCA(rld, intgroup=c("Contrasts"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Contrasts)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  theme_pubr()+
  scale_color_brewer(palette="Set1",labels=c('Mock', 'SARS-CoV-2 Day 3', 'SARS-CoV-2 Day 7'))
  
ggsave("../../cambridge-upr-project/plots/PCA_plot.png", device = "png", dpi=300,width=7,height=6,units="in")





# Differential gene expression with EdgeR ####
library(edgeR)
dgList<- DGEList(
  counts = raw_counts,
  samples = coldata2,
  genes = annotation,
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
dgList$samples
dgList$samples$lib.size


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
  labels = dgList$samples$Label)

library(ggplot2)
library(ggrepel)


pca_plot<-function(inputDgList,labels,color,shape,PC1,PC2){
  
  
  ## PCA analysis
  dgList<-inputDgList
  pcavals <- log2(cpm(dgList$counts) + 0.1)   # log transforming the counts makes it easy to compare very high and very
  # low levels of expression. Adding 0.1 is necessary to avoid log transforming zeros.
  pcavals_var <- pcavals[apply(pcavals, 1, function(x) length(unique(x))) > 1, ]  #  remove genes that have the same expression in every sample
  pcavals_t <- t(pcavals_var)                          # transpose  (swap row and columns) the data frame
  pcavals_m<- as.matrix(pcavals_t, scale =T)        # conver the data from to a matrix
  pca <- prcomp(pcavals_m)               #  Perform the principle components analysis
  fraction_explained <- round((pca$sdev)^2/sum(pca$sdev^2), 3) * 100 # Get the fraction of variance explained by each princple component 
  
  #guides<-guides(shape=FALSE)
  #if(PC2=="PCA2"){guides<-guides(shape=FALSE,color=FALSE)}
  plotdata<-as.data.frame(pca$x)
  plotdata$label = dgList$samples[[labels]]         # Make a list of labels for each point
  plotdata[[color]] <-as.factor(dgList$samples[[color]])     # Make a list of factors to determine the color of each point
  plotdata[[shape]] <-as.factor(dgList$samples[[shape]])  # Make a list of factors to determine the shape of each point
  pcaplot<-ggplot(plotdata, aes_string(PC1, PC2, color = color, shape = shape, label = 'label')) + 
    geom_point(size = 3) + # change this numnber to alter the size of the points
    stat_ellipse()+
    geom_hline(yintercept = 0, linetype="dashed")+
    geom_vline(xintercept = 0, linetype="dashed")+
    geom_text_repel(size = 3, segment.size = 0.1, nudge_y = 0.05, nudge_x = 0.05, point.padding = unit(0.5, 'lines'), force = 2, max.iter = 10000, max.overlaps = Inf)+
    theme_classic()+
    ylim(-450,450)+
    xlim(-500,500)+
    theme(
      axis.text.x=element_text(size=14),
      axis.text.y=element_text(size=14),
      axis.title.x=element_text(size=16),
      axis.title.y=element_text(size=16),
      legend.text = element_text(size=12),
      legend.title = element_text(size=14),
      legend.position="bottom"
    )
  #By addgng `geom_text_repel` the labels will be added in sensible places, join to the point with a line if necessary.
  
  
  ## make a scree plot
  variance<- pca$sdev^2                 # get the variance for each principle component.
  variance <- variance[1:(length(variance)-1)]     # Miss the last one because it's 0 and makes the plot look weird
  pev<-(cumsum(variance)/sum(variance))*100       # Calculate the culmative variance explained by increasing numbers of principle components.
  plotdata <- data.frame(
    PC = c(1:length(pev)),
    pev = pev,
    type = 'scree'
  )
  
  origin_line <- data.frame(PC = c(0, 1), pev = c(0, pev[1]), type = 'orign')
  screeplot<-  ggplot(rbind(plotdata, origin_line), aes(PC, pev,color = type, linetype = type)) +
    geom_point(size = 3) + geom_path()  + scale_color_manual(values = c('black', 'grey'))+
    ylab('Cumulative variance explained (%)') + xlab('Principal component') + theme_bw() + theme(legend.position="none") + 
    theme_classic()
  
  list(
    pca = pca,
    pcaplot = pcaplot,
    screeplot  = screeplot
  )
  
}
pca_output<-pca_plot(dgList, labels = 'Label', color = 'Label', shape = 'Day', 'PC1','PC2')
print(pca_output$pcaplot)




design <- model.matrix(~0+group, data = dgList$samples)

levels(dgList$samples$group)

#Estimating dispersion
library(statmod)

design
dgGlm <- estimateDisp(dgList, design, robust = TRUE)
plotBCV(dgGlm)

#Fitting data to the model
fit <- glmQLFit(dgGlm, design, robust = TRUE)

design

my.contrasts <- makeContrasts(SARSvmock_d6 = groupSARS_Cov_2_d6_Lung-groupMock_d10_Lung,
                               SARSvmock_d10 = groupSARS_Cov_2_d10_Lung-groupMock_d10_Lung,
                               levels=design)
my.contrasts


fit <- glmQLFit(dgGlm, design, robust = TRUE)

de <- glmQLFTest(fit, contrast=my.contrasts)

de.SARSvmock_d10           <-    glmQLFTest(fit, contrast=my.contrasts[,"SARSvmock_d10"])
de.SARSvmock_d6           <-    glmQLFTest(fit, contrast=my.contrasts[,"SARSvmock_d6"])



summary(de)

de$dispersion

top_genes <- topTags(de, n = 75)
top_genes




results <- topTags(de, n=nrow(dgList), sort.by='none')$table
results


results.SARSvmock_d10           <- topTags(de.SARSvmock_d10       , n=nrow(dgList), sort.by = 'none')$table
results.SARSvmock_d6           <- topTags(de.SARSvmock_d6       , n=nrow(dgList), sort.by = 'none')$table





get.de.genes <- function(results) {
  de.genes <- results[ which(results$logFC >= 2 | results$logFC <= -2),]
  de.genes <-de.genes[de.genes$FDR <=0.05,]
}

de.genes.SARS_d10 <- get.de.genes(results.SARSvmock_d10)


de.genes.SARS_d6 <- get.de.genes(results.SARSvmock_d6)

fdr_threshold <- 0.05
fc_threshold <- 2
fc_threshold2 <-2


volcano_plot<- function(results_table, fc_threshold, fdr_threshold,log=FALSE, label, ylim, xlim){
  results_table <- results_table %>% filter(!grepl('Gm|Rik|ENS', gene_name))
  results_table <- results_table %>% distinct(gene_name, .keep_all = TRUE)
  results_table$significant <- 'no'
  results_table$significant[ abs(results_table$logFC) >= log2(fc_threshold) &
                               results_table$FDR <= fdr_threshold  ] <- 'yes'
  results_table$plot_gene<- ifelse(results_table$logFC >= results_table$logFC[order(results_table$logFC, decreasing = T)[10]] & results_table$FDR < 0.05 |
                                     results_table$logFC <= results_table$logFC[order(results_table$logFC, decreasing = F)][10] & results_table$FDR < 0.05, results_table$gene_name, NA)
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

volcano_plot(results.SARSvmock_d6,fc_threshold, fdr_threshold,log=TRUE, "Day 3",5,8)
volcano_plot(results.SARSvmock_d10, fc_threshold, fdr_threshold,log=TRUE, "Day 10",5,8)


assayed.genes <- results

get.up.de.genes <- function(results) {
  de.genes <- results[ which(results$logFC >= 2),]
}

get.down.de.genes <- function(results) {
  de.genes <- results[ which(results$logFC <= 2),]
}

de.genes.SARSD3.up<-   get.up.de.genes(de.genes.SARS_d6)
de.genes.SARSD3.down<- get.down.de.genes(de.genes.SARS_d6)

de.genes.SARSD7.up<-   get.up.de.genes(de.genes.SARS_d10)
de.genes.SARSD7.down<- get.down.de.genes(de.genes.SARS_d10)



GO.genes <- function(de.genes){
  go.genes<-bitr(geneID = de.genes$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db", drop = TRUE)
}

de.genes.SARSD3.up

#bitr(geneID, fromType, toType, OrgDb, drop = TRUE)
GO.genes.SARSD3.up    <- bitr(geneID = rownames(de.genes.SARSD3.up), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db", drop = TRUE)
GO.genes.SARSD3.down  <- bitr(geneID = rownames(de.genes.SARSD3.down), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db", drop = TRUE)
GO.genes.SARSD7.up    <- bitr(geneID = rownames(de.genes.SARSD7.up), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db", drop = TRUE)
GO.genes.SARSD7.down  <- bitr(geneID = rownames(de.genes.SARSD7.down), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db", drop = TRUE)

assayed.genes2 <- bitr(geneID = rownames(assayed.genes), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db", drop = TRUE)
GO.genes.SARSD3.up$group <- "Up"
GO.genes.SARSD3.down$group <- "Down"
GO.genes.SARSD7.up$group <- "Up"
GO.genes.SARSD7.down$group <- "Down"


GO.genes.SARSD3.up$othergroup <- "SARS-CoV-2 Day 3"
GO.genes.SARSD3.down$othergroup <- "SARS-CoV-2 Day 3"
GO.genes.SARSD7.up$othergroup <- "SARS-CoV-2 Day 7"
GO.genes.SARSD7.down$othergroup <- "SARS-CoV-2 Day 7"

GO.List <- rbind(GO.genes.SARSD3.up, GO.genes.SARSD3.down, GO.genes.SARSD7.up, GO.genes.SARSD7.down)



GO.List$othergroup_f = factor(GO.List$othergroup, levels=c("SARS-CoV-2 Day 3",
                                                           "SARS-CoV-2 Day 7"))

compareCC<-compareCluster(ENTREZID~group+othergroup_f, data=GO.List, fun = "enrichGO", OrgDb = "org.Mm.eg.db", ont="CC", pvalueCutoff = 0.01, qvalueCutoff  = 0.05, pAdjustMethod = "BH", universe = assayed.genes2$ENTREZID, readable = TRUE)
compareMF<-compareCluster(ENTREZID~group+othergroup_f, data=GO.List, fun = "enrichGO", OrgDb = "org.Mm.eg.db", ont="MF", pvalueCutoff = 0.01, qvalueCutoff  = 0.05, pAdjustMethod = "BH", universe = assayed.genes2$ENTREZID, readable = TRUE)
compareBP<-compareCluster(ENTREZID~group+othergroup_f, data=GO.List, fun = "enrichGO", OrgDb = "org.Mm.eg.db", ont="BP", pvalueCutoff = 0.01, qvalueCutoff  = 0.05, pAdjustMethod = "BH", universe = assayed.genes2$ENTREZID, readable = TRUE)


simplifyCC <- simplify(compareCC, cutoff=0.7, by="qvalue", select_fun=min)
simplifyMF <- simplify(compareMF, cutoff=0.7, by="qvalue", select_fun=min)
simplifyBP <- simplify(compareBP, cutoff=0.7, by="qvalue", select_fun=min)



CC <-clusterProfiler::dotplot(simplifyCC, x = "group", showCategory = 25, color="qvalue", includeAll = TRUE) + ggplot2::facet_grid(~factor(othergroup_f, levels=c("SARS-CoV-2 Day 3",
                                                                                                                                                                  "SARS-CoV-2 Day 7"))) + ggtitle("Molecular Function GO Terms",
                                                                                                                                                                                                  xlab(""))


MF <-clusterProfiler::dotplot(simplifyMF, x = "group", showCategory = 25, color="qvalue", includeAll = TRUE) + ggplot2::facet_grid(~factor(othergroup_f, levels=c("SARS-CoV-2 Day 3",
                                                                                                                                                 "SARS-CoV-2 Day 7"))) + ggtitle("Molecular Function GO Terms", xlab(""))

BP<-clusterProfiler::dotplot(simplifyBP, x = "group", showCategory = 10, color="qvalue", includeAll = TRUE) + ggplot2::facet_grid(~factor(othergroup_f, levels=c("SARS-CoV-2 Day 3",
    
                                                                                                                                                                                                                                                                                                             "SARS-CoV-2 Day 7"))) + ggtitle("Biological Process GO Terms", xlab(""))

CC
ggsave("../../cambridge-upr-project//plots/GO_CC.tiff", device = "tiff", dpi = 300, width = 12, height = 8)
BP
ggsave("../../cambridge-upr-project/plots/GO_BP.tiff", device = "tiff", dpi = 300, width = 14, height = 10)
MF
ggsave("../../cambridge-upr-project/plots/GO_MF.tiff", device = "tiff", dpi = 300, width = 15, height = 8)



MFres<-simplifyMF@compareClusterResult
BPres<-simplifyBP@compareClusterResult


LFCplot<- results %>% filter(gene_name %in% genes)

plot2<-gather(LFCplot, group, LFC, 12:13 )

ggplot(plot2)+
  geom_tile(aes(x=group,y=gene_name,fill=LFC))+
  xlab("dpi")+
  ylab("")+
  scale_fill_gradient(low = "black", high = "red", limits =c(-2,2))+
  theme_pubclean()
