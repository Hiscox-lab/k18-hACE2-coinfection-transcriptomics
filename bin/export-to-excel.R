setwd("/home/rebee/projects/mice-transcriptomics/results/illumina/")

library(openxlsx)

df1<-read.csv("CoinfD3.sig.csv") 
df2<-read.csv("CoinfD3vIAVD3.sig.csv")
df3<-read.csv("CoinfD3vSARSD3.sig.csv")
df4<-read.csv("CoinfD7.sig.csv")
df5<-read.csv("CoinfD7vIAVD7.sig.csv")
df6<-read.csv("CoinfD7vSARSD7.sig.csv")
df7<-read.csv("SARSD3.sig.csv")
df8<-read.csv("SARSD7.sig.csv")
df9<-read.csv("IAVD3.sig.csv")
df10<-read.csv("IAVD7.sig.csv")
df11<-read.csv("d6vax_SARS_v_Coinf.sig.csv")
df12<-read.csv("d6vax_SARS_v_SARS.sig.csv")
df13<-read.csv("d10vax_SARS_v_Coinf.sig.csv")
df14<-read.csv("d10vax_SARS_v_SARS.sig.csv")



dataset_names <- list('Coinf-v-Mock-D6'= df1,
     'Coinf-v-IAV-D6'= df2,
     'Coinf-v-SC2-D6'= df3,
     'Coinf-v-Mock-D10'= df4,
     'Coinf-v-IAV-D10'= df5,
     'Coinf-v-SC2-D10'= df6,
     'SC2-v-Mock-D6'= df7,
     'SC2-v-Mock-D10'= df8,
     'IAV-v-Mock-D6'= df9,
     'IAV-v-Mock-D10'=df10,
     'Fluenz-Coinf-v-Coinf-D6'=df11,
     'Fluenz-Coinf-v-SCV2-D6' =df12,
     'Fluenz-Coinf-v-Coinf-D10' =df13,
     'Fluenz-Coinf-v-SCV2-D10'=df14)


write.xlsx(dataset_names, file = 'Supplementary-data-DGE.xlsx')

GO1<-"../results/illumina/gene-ontology-BP.csv"
GO2<-"../results/illumina/gene-ontology-CC.csv"
GO3<-"../results/illumina/gene-ontology-MF.csv"
GO4<-"../results/illumina/gene-ontology-fluenz-BP.csv"
GO5<-"../results/illumina/gene-ontology-fluenz-CC.csv"
GO6<-"../results/illumina/gene-ontology-fluenz-MF.csv"

dataset_names<-list('Coinfection BP'=GO1,
                    'Coinfection CC'=GO2,
                    'Coinfection MF'=GO3,
                    'Fluenz BP'= GO4,
                    'Fluenz CC'= GO5,
                    'FLuenz MF'= GO6)

write.xlsx(dataset_names, file = 'Supplementary-data-GO.xlsx')
