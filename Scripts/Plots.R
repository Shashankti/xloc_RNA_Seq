# make plots for EDA
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(DEGreport)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(biomaRt)
library(viridis)
library(genefilter)
library(ComplexHeatmap)
library(data.table)
library(cowplot)
library(InteractiveComplexHeatmap)
library("glmpca")

# perform lfc shrinkage for better visualisation
res_control_vs_D3_shrunk <- lfcShrink(dds_2f,contrast = c("condition","sgRNA_d3","control"),type = "ashr")
res_control_vs_D5_shrunk <- lfcShrink(dds_2f,contrast = c("condition","sgRNA_d5","control"),type="ashr")
res_D3_vs_D5_shrunk <- lfcShrink(dds_2f,contrast = c("condition","sgRNA_d5","sgRNA_d3"),type="ashr")

# generate pca plots 
## create transformed values
vsd <- vst(dds_2f,blind=FALSE)
rld <- rlog(dds_2f,blind=F)
head(assay(vsd),3)
select <- order(rowMeans(counts(dds_2f,normalized=TRUE)),
                decreasing=TRUE)[1:40]

#$PCA Plots

pcaData <- plotPCA(vsd,intgroup=c("sample"),returnData=F)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
vst_pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")+
  theme_cowplot()
# GLM-PCA plots
tiff("Plots/GPCA_plot_vst.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(vst_pca_plot)
dev.off()

gpca <- glmpca(counts(dds_2f), L=2)
gpca.dat <- gpca$factors
gpca.dat$condition <- dds_2f$condition

vst_pca_plot <- ggplot(gpca.dat, aes(x = dim1, y = dim2, color = condition)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")+
  theme_cowplot()

#adjusting normalized counts
rownames(normalized_counts) <- normalized_counts$GeneID

# Making heatmaps for D3 vs control

padj.cutoff <- 0.05
lfc.cutoff <- 1

threshold <- d3_df$padj < padj.cutoff & abs(d3_df$log2FoldChange) > lfc.cutoff
length(which(threshold))
d3_df$threshold <- threshold
d3_sig <- data.frame(subset(d3_df,threshold==TRUE))
d3_sig <- d3_sig[order(d3_sig$padj),]
only_d3 <- normalized_counts[rownames(d3_sig),c(1,2,3,5,7)]
only_d3$padj <- d3_sig[match(rownames(only_d3),rownames(d3_sig)),6]
#melted_norm <- data.frame(reshape::melt(sigOE_ordered_qValue)) 
only_low_mat <- as.matrix(only_low[,c(2,3,4,5)])
only_d3_mat_scaled <- t(scale(t(only_d3[,c(2,3,4)])))
rownames(only_d3_mat_scaled) <- only_d3$GeneName
only_d3_mat_scaled <- cbind(only_d3_mat_scaled,only_d3[match(rownames(only_d3_mat_scaled),only_d3$GeneName),6])
ht_d3 <-ComplexHeatmap::Heatmap(only_d3_mat_scaled[,c(1,2,3)],
                                 column_names_side = "bottom",
                                 cluster_rows = T,
                                 cluster_columns = F,
                                 column_names_rot = 45,
                                 column_names_gp = gpar(fontsize=7),
                                 row_names_gp = gpar(fontsize=7),
                                row_names_rot = 0,
                                name = "Z-score",width = 4*unit(25,'mm'),
                                column_title = "Control vs Day3-sgRNA")
ht_d3
tiff("Plots/D3_vs_ct.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(ht_d3)
dev.off()

#making a volcano plot
d3_df$GeneName <- normalized_counts[match(rownames(d3_df),normalized_counts$GeneID),7]
E1 <- EnhancedVolcano(d5_df,
                      lab = d5_df$GeneName,
                      x='log2FoldChange',
                      y='padj',
                      title = 'Day5 vs control',
                      subtitle = 'Differential Expression',
                      caption = bquote(~Log[2]~ 'foldchange cutoff 1; p-value cutoff,0.05'),
                      pCutoff = 0.05,
                      FCcutoff = 1,
                      col = c('black','black','green','red'),
                      hline = 10e-03,
                      labSize = 3,
                      pointSize = 1,
                      shape = c(1,1,1,25),
                      colAlpha = 1,
                      legendPosition = 'top',
                      legendIconSize = 3,
                      drawConnectors = F,
                      widthConnectors = 0.5)+
  coord_cartesian(xlim = c(-10,15))
tiff("Plots/Day5_vs_control_volcano.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(E1)
dev.off()



#Making heatmaps for D5 vs control
padj.cutoff <- 0.05
lfc.cutoff <- 1

threshold <- d5_df$padj < padj.cutoff & abs(d5_df$log2FoldChange) > lfc.cutoff
length(which(threshold))
d5_df$threshold <- threshold
d5_sig <- data.frame(subset(d5_df,threshold==TRUE))
d5_sig <- d5_sig[order(d5_sig$padj),]
only_d5 <- normalized_counts[rownames(d5_sig),c(1,2,4,6,7)]
only_d5$padj <- d5_sig[match(rownames(only_d5),rownames(d5_sig)),6]
#melted_norm <- data.frame(reshape::melt(sigOE_ordered_qValue)) 
#only_d5_mat <- as.matrix(only_low[,c(2,3,4,5)])
only_d5_mat_scaled <- t(scale(t(only_d5[,c(2:4)])))
rownames(only_d5_mat_scaled) <- only_d5$GeneName

d5_sig$GeneName <- only_d5[match(rownames(d5_sig),only_d5$GeneID),6]

#save output as csv file
write.csv(d5_sig,file = "Analysis/res_ct_vs_d5.csv")

only_d5_mat_scaled <- cbind(only_d5_mat_scaled,only_d5[match(rownames(only_d5_mat_scaled),only_d5$GeneName),6])
ht_d5 <-ComplexHeatmap::Heatmap(na.omit(only_d5_mat_scaled[,c(1,2,3)]),
                                column_names_side = "bottom",
                                cluster_rows = T,
                                cluster_columns = F,
                                column_names_rot = 45,
                                column_names_gp = gpar(fontsize=7),
                                row_names_gp = gpar(fontsize=7),
                                row_names_rot = 0,
                                name = "Z-score",width = 4*unit(25,'mm'),
                                column_title = "Control vs Day5-sgRNA")
ht_d5



#day 3 vs day5
padj.cutoff <- 0.05
lfc.cutoff <- 1

threshold <- d3_vs_d5_df$padj < padj.cutoff & abs(d3_vs_d5_df$log2FoldChange) > lfc.cutoff
length(which(threshold))
d3_vs_d5_df$threshold <- threshold
d3_vs_d5_sig <- data.frame(subset(d3_vs_d5_df,threshold==TRUE))
d3_vs_d5_sig <- d3_vs_d5_sig[order(d3_vs_d5_sig$padj),]
only_d3_d5 <- normalized_counts_2[rownames(d3_vs_d5_sig),c(3:11)]
#melted_norm <- data.frame(reshape::melt(sigOE_ordered_qValue)) 
only_d3_d5_mat_scaled <- t(scale(t(only_d3_d5[,c(1:8)])))
rownames(only_d3_d5_mat_scaled) <- only_d3_d5$GeneName
#only_d3_d5_mat_scaled <- cbind(only_d3_d5_mat_scaled,only_d3_d5[match(rownames(only_d3_d5_mat_scaled),only_d3_d5$GeneName),6])
ht_d3_d5 <-ComplexHeatmap::Heatmap(only_d3_d5_mat_scaled,
                                column_names_side = "bottom",
                                cluster_rows = T,
                                cluster_columns = F,
                                column_names_rot = 45,
                                column_names_gp = gpar(fontsize=7),
                                row_names_gp = gpar(fontsize=7),
                                row_names_rot = 0,
                                name = "Z-score",width = 6*unit(25,'mm'),
                                column_title = "Day5 vs Day3-sgRNA")
ht_d3_d5
tiff("Plots/D3_vs_D5.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(ht_d3_d5)
dev.off()

#making a volcano plot
d3_df$GeneName <- normalized_counts[match(rownames(d3_df),normalized_counts$GeneID),7]
E1 <- EnhancedVolcano(d3_vs_d5_df,
                      lab = d3_vs_d5_df$GeneName,
                      x='log2FoldChange',
                      y='padj',
                      title = 'Day5 vs Day3',
                      subtitle = 'Differential Expression',
                      caption = bquote(~Log[2]~ 'foldchange cutoff 1; p-value cutoff,0.05'),
                      pCutoff = 0.05,
                      FCcutoff = 1,
                      col = c('black','black','green','red'),
                      hline = 10e-03,
                      labSize = 3,
                      pointSize = 1,
                      shape = c(1,1,1,25),
                      colAlpha = 1,
                      legendPosition = 'top',
                      legendIconSize = 3,
                      drawConnectors = F,
                      widthConnectors = 0.5)+
  coord_cartesian(xlim = c(-10,15))
tiff("Plots/Day5_vs_Day3_volcano.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(E1)
dev.off()


#making heatmaps with all the columns/samples
all_d5 <- normalized_counts[rownames(d5_sig),]
all_d5$padj <- d5_sig[match(rownames(all_d5),rownames(d5_sig)),6]
#melted_norm <- data.frame(reshape::melt(sigOE_ordered_qValue)) 
#only_d5_mat <- as.matrix(only_low[,c(2,3,4,5)])
all_d5_mat_scaled <- t(scale(t(all_d5[,c(2:6)])))
rownames(all_d5_mat_scaled) <- all_d5$GeneName

all_d5_mat_scaled <- cbind(all_d5_mat_scaled,all_d5[match(rownames(all_d5_mat_scaled),all_d5$GeneName),8])
ht_d5_all <-ComplexHeatmap::Heatmap(na.omit(all_d5_mat_scaled[c(1:30),c(1,2,4,3,5)]),
                                column_names_side = "bottom",
                                cluster_rows = T,
                                cluster_columns = F,
                                column_names_rot = 45,
                                column_names_gp = gpar(fontsize=7),
                                row_names_gp = gpar(fontsize=7),
                                row_names_rot = 0,
                                name = "Z-score",width = 6*unit(25,'mm'),
                                column_title = "Control vs Day5-sgRNA")
ht_d5_all
tiff("Plots/Day5_vs_control_wt_d3_readable.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(ht_d5_all)
dev.off()



##saving output of sig genes in csv files
write.csv(d3_sig,"Analysis/res_d3_vs_ct_sig.csv")
write.csv(d5_sig,"Analysis/res_d5_vs_ct_sig.csv")
write.csv(d3_vs_d5_sig,"Analysis/res_d5_vs_d3_sig.csv")







#make list of neighboring genes
mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart)
attributes <- c("ensembl_gene_id","start_position","end_position","strand","hgnc_symbol","chromosome_name","ucsc","band")
filters <- c("chromosome_name","start","end")
values <- list(chromosome="18",start=" 15169976",end="16169976")
all.genes_up <- getBM(attributes=attributes, filters=filters, values=values, mart=mart)
all.genes_up
values <- list(chromosome="18",start=" 14165346",end="15165346")
all.genes_d <- getBM(attributes=attributes, filters=filters, values=values, mart=mart)
all.genes_d

all.genes <- rbind(all.genes_d,all.genes_up)


# overlap the neighboring genes to find expression changes

##overlap with d3
d3_sig[rownames(d3_sig) %in% all.genes_d$ensembl_gene_id,]
d3_df[rownames(d3_df) %in% all.genes_d$ensembl_gene_id,]