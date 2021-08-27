BiocManager::install("biomaRt")
install.packages("dplyr")
install.packages("ellipsis")

library("ellipsis")
library("biomaRt")
library("dplyr")

setwd("/Users/chixiangchen/Desktop/phd_in_psu/research/Prof._Wu/Network/tissue/blood_vessel/new_analysis_07082021/")
gene_exact_names<-read.csv(file = "85_gene_information.csv")
head(gene_exact_names)
dim(gene_exact_names)

mart<-useDataset("hsapiens_gene_ensembl",useMart("ensembl"))
#genes<-gene_exact_names$x
genes<-sub("\\..*", "", gene_exact_names$x)
gene_IDs<-getBM(filters = "ensembl_gene_id",attributes=c("ensembl_gene_id","hgnc_symbol"),
                values=genes, mart=mart)
gene_IDs$hgnc_symbol
which(genes %in% gene_IDs$ensembl_gene_id)

doc.gene<-data.frame(ensg=genes)
head(doc.gene)
doc.gene$gene.name<-0
rownames(doc.gene)<-doc.gene$ensg
doc.gene[gene_IDs$ensembl_gene_id,2]<-gene_IDs$hgnc_symbol
rownames(doc.gene)<-1:nrow(doc.gene)
doc.gene

write.csv(doc.gene,"85_gene_name_information.csv")





