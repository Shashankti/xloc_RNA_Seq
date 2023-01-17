# DESeq2 analysis from the ZARP counts dta

library (DESeq2)
library (tidyverse)
library (airway)
library(biomaRt)

# Step 1: preparing count data ----------------
# read in counts data
counts_data <- read.table ('./Data/genes_numreads.tsv',header = T)
head(counts_data)
counts_data <- counts_data %>% column_to_rownames("Name")

# read in sample info
colData_2f<-read.table ('./Data/sample_info.tsv')
#view (colData_2f)



#making sure the row names in colData matches the column names in counts_data 

all(colnames(counts_data) %in% rownames(colData_2f))

#check the order
all(colnames(counts_data) == rownames(colData_2f))

# Step 2: construct a DESeqDataSet object
dds_2f <- DESeqDataSetFromMatrix (countData = round(counts_data), colData = colData_2f , design = ~ condition )

keep_2f <- rowSums(counts(dds_2f)) >= 10
dds_2f <- dds_2f [keep_2f,]
dds_2f

#get normalised counts information
dds_2f <- estimateSizeFactors(dds_2f)
sizeFactors(dds_2f)
dds_2f$condition <- relevel(dds_2f$condition,ref="control")

# Step 3: Run DESeq ----------------------

dds_2f <- DESeq(dds_2f)
resultnames_2f <- resultsNames(dds_2f)
resultnames_2f
#day3 vs day5 cannot merge replicattes because of 2 replicates minimum requirement
res_D3_vs_D5 <- results(dds_2f,contrast = c("condition","sgRNA_d5","sgRNA_d3"))
normalized_counts_2 <- counts(dds_2f,normalized=TRUE)
normalized_counts_2 <- as.data.frame(normalized_counts_2)
normalized_counts_2$GeneNames <- normalized_counts$GeneName
#merge technical replicates
dds_2f <- collapseReplicates(dds_2f,dds_2f$sample)
normalized_counts <- counts(dds_2f,normalized=TRUE)
normalized_counts <- as.data.frame(normalized_counts) %>% rownames_to_column("GeneID")
write.table(normalized_counts,file = "./Analysis/Normalized_counts_merged.txt",sep = "\t",quote = F,col.names = NA)

res_control_vs_D3 <- results(dds_2f,contrast = c("condition","sgRNA_d3","control"))
res_control_vs_D5 <- results(dds_2f,contrast = c("condition","sgRNA_d5","control"))


#summary of results
summary(res_control_vs_D3)
summary(res_control_vs_D5)
summary(res_D3_vs_D5)

#convert to dataframe
d3_df <- as.data.frame(res_control_vs_D3)
d5_df <- as.data.frame(res_control_vs_D5)
d3_vs_d5_df <- as.data.frame(res_D3_vs_D5)

# LFC shrinkage for visualisation

#annotate geneIDs

## control vs d3
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = rownames(d3_df), mart= mart)

d3_df$GeneName <- gene_IDs[match(rownames(d3_df),gene_IDs[,1]),2]
d3_df <- arrange(d3_df,padj)

## control vs d5
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = rownames(d5_df), mart= mart)

d5_df$GeneName <- gene_IDs[match(rownames(d5_df),gene_IDs[,1]),2]
d5_df <- arrange(d5_df,padj)

## d3 vs d5
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = rownames(d3_vs_d5_df), mart= mart)

d3_vs_d5_df$GeneName <- gene_IDs[match(rownames(d3_vs_d5_df),gene_IDs[,1]),2]
d3_vs_d5_df <- arrange(d3_vs_d5_df,padj)

# Add gene names to normalised counts
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = normalized_counts[,1], mart= mart)

#!normalized_counts <- merge(normalized_counts,gene_IDs, by.x="GeneID",by.y="ensembl_gene_id")
normalized_counts$GeneName <- gene_IDs[match(normalized_counts[,1],gene_IDs[,1]),2]
?corr.test

