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

GO1<-read.csv("gene-ontology-BP.csv")
GO2<-read.csv("gene-ontology-CC.csv")
GO3<-read.csv("gene-ontology-MF.csv")
GO4<-read.csv("gene-ontology-fluenz-BP.csv")
GO5<-read.csv("gene-ontology-fluenz-CC.csv")
GO6<-read.csv("gene-ontology-fluenz-MF.csv")

dataset_names<-list('Coinfection BP'=GO1,
                    'Coinfection CC'=GO2,
                    'Coinfection MF'=GO3,
                    'Fluenz BP'= GO4,
                    'Fluenz CC'= GO5,
                    'FLuenz MF'= GO6)

write.xlsx(dataset_names, file = 'Supplementary-data-GO.xlsx')






df1<-read.csv("results.IAVD3VMock.csv")
df2<-read.csv("results.IAVD7VMock.csv")
df3<-read.csv("results.SARS2D3VMock.csv")
df4<-read.csv("results.SARS2D7VMock.csv")
df5<-read.csv("results.CoinfD3vmock.csv")
df6<-read.csv("results.CoinfD7vmock.csv")
df7<-read.csv("results.CoinfD3vIAVD3.csv")
df8<-read.csv("results.CoinfD7vIAVD7.csv")
df9<-read.csv("results.CoinfD3vSARSD3.csv")
df10<-read.csv("results.CoinfD7vSARSD7.csv")
df11<-read.csv("differential_results.d6fluenz_SARS_v_coinf.csv   ")
df12<-read.csv("differential_results.d6fluenz_SARS_v_SARS.csv   ")
df13<-read.csv("differential_results.d10fluenz_SARS_v_coinf.csv    ")
df14<-read.csv("differential_results.d10fluenz_SARS_v_SARS.csv    ")



dataset_names <- list('IAV-v-Mock-D6'= df1,
                      'IAV-v-Mock-D10'= df2,
                      'SC2-v-Mock-D6'= df3,
                      'SC2-v-Mock-D10'= df4,
                      'Coinf-v-Mock-D6'= df5,
                      'Coinf-v-Mock-D10'= df6,
                      'Coinf-v-IAV-D6'= df7,
                      'Coinf-v-IAV-D10'= df8,
                      'Coinf-v-SC2-D6'= df9,
                      'Coinf-v-SC2-D10'=df10,
                      'Fluenz-Coinf-v-Coinf-D6'=df11,
                      'Fluenz-Coinf-v-SCV2-D6' =df12,
                      'Fluenz-Coinf-v-Coinf-D10' =df13,
                      'Fluenz-Coinf-v-SCV2-D10'=df14)

write.xlsx(dataset_names, file = 'Supplementary-data-DGE-all.xlsx')
