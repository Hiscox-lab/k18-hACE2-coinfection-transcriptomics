setwd("/home/rebee/projects/mice-transcriptomics/vasculitis/")

library(openxlsx)

df1<-read.csv("IAV_DGE.csv") 
df2<-read.csv("SC2_D3_DGE.csv")
df3<-read.csv("SC2_D7.DGE.csv")



dataset_names <- list('IAV_D6_DGE'= df1,
     'SC2_D3_DGE'= df2,
     'SC2_D7_DGE'= df3)


write.xlsx(dataset_names, file = 'Supplementary-data-DGE.xlsx')

GO1<-"compareCluster_BP_output.csv"


dataset_names<-list('BP'=GO1)

write.xlsx(dataset_names, file = 'Supplementary-data-GO.xlsx')
