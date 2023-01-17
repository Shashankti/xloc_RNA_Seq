# GO and GSEA for RNA-Seq data

library(clusterProfiler)
library(goseq)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggridges)
library(pathview)
library(ggplot2)
set.seed(2022)

#Select up regulated gene
up_genes <- rownames(res_control_vs_D5_shrunk[res_control_vs_D5_shrunk$log2FoldChange>0,])
#make the over representation test
GO_res <- enrichGO(gene = up_genes,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENSEMBL",
                   ont = "BP")
as.data.frame(GO_res)
# plot the results
fit1 <- plot(barplot(GO_res,showCategory = 30,font.size = 7,title = "Significantly up-regulated genes control vs Day5"))
fit1

# same process for the downregulated genes
down_genes <- rownames(res_control_vs_D3_shrunk[res_control_vs_D3_shrunk$log2FoldChange<0,])

GO_res_d <- enrichGO(gene = down_genes,
                     OrgDb = "org.Hs.eg.db",
                     keyType = "ENSEMBL",
                     ont = "BP")
as.data.frame(GO_res_d)
fit2 <- plot(barplot(GO_res_d,showCategory = 30,font.size = 7,title = "Significantly down-regulated genes in control vs Day5"))
fit2

# make a dot plot
dot_go <- dotplot(GO_res,showCategory=20,title="Up-regulated GO-enrichment for control vs Day5",font.size=7)
dot_go <- dotplot(GO_res_d,showCategory=20,title="Down-regulated GO-enrichment for controld vs Day5",font.size=9)
tiff("Plots/GO_bar_Day5_vs_control_up.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(fit1)
dev.off()
#GSEA

# start the GSEA with DESEQ2 output
d3_df_GSEA <- d5_df[d5_df$baseMean>50,]
d3_df_GSEA <- d3_df_GSEA[order(-d3_df_GSEA$stat),]
gene_list_up <- d3_df_GSEA$stat
names(gene_list_up) <- rownames(d3_df_GSEA)
gse <- gseGO(gene_list_up,
             ont = "BP",
             OrgDb = "org.Hs.eg.db",
             keyType = "ENSEMBL",
             eps = 1e-300)
gse_df <- as.data.frame(gse)

#make a pubmed trends plot
terms <- gse$Description[1:3]
pmcplot(terms,2010:2020,proportion = F)+theme_cowplot()
#make a dotplot
gsea_dot <- dotplot(gse,showCategory=20,split=".sign",font.size=5) +
  facet_grid(.~.sign)
tiff("Plots/GSEA_dot_Day5_vs_control.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(gsea_dot)
dev.off()

#save gsea plots
g1 <- gseaplot(gse,geneSetID = 15,by='all',title = gse_df$Description[15])+theme_cowplot()
tiff("Plots/ER to Golgi vesicle mediated transport.tiff",width = 35,height = 21,units = 'cm',res = 300)
print(g1)
dev.off()


#make an enrichment map
x2 <- pairwise_termsim(gse)
emapplot(x2,showCategory = 40)




# load annotate library
library(annotate)
findNeighbors("humanCHRLOC", chromosome = 18, upBase = 15167128, downBase = 15168337,organism)

library(biomaRt)


