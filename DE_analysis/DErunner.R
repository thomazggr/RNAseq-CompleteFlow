# Main script was made by:
# @ Thomas W. Battaglia
# tb1280@nyu.edu
# Changes made by:
# @ Thomaz G. Ramalheira
# thomaz@vivaldi.net

# Install R base and R dev
# Install BiocManager first
# sudo apt install r-cran-xml r-cran-xml2 libxml2-dev libcurl4-openssl-dev libssl-dev 
# install.packages(c("httr", "curl", "RCurl", "openssl", "XML"))
# BiocManager::install(c("DESeq2","ggplot2","clusterProfiler","biomaRt","ReactomePA","ggsci","gage","dplyr","topGO","DOSE","org.Hs.eg.db","org.Mm.eg.db","pheatmap","genefilter","GO.db","KEGG.db","RColorBrewer"))
# Load required libraries
library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(ReactomePA)
library(DOSE)
library(KEGG.db)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(pheatmap)
library(genefilter)
library(RColorBrewer)
library(GO.db)
library(topGO)
library(dplyr)
library(gage)
library(ggsci)


# - - - - - - - - - - - - - 
# Import gene counts table
# - - - - - - - - - - - - - 
# - skip first row (general command info)
# - make row names the gene identifiers
countdata <- read.table("final_counts.txt", header = TRUE, skip = 1, row.names = 1)

# Remove .bam + '..' from column identifiers
colnames(countdata) <- gsub(".sam", "", colnames(countdata), fixed = T)
colnames(countdata) <- gsub("X", "", colnames(countdata), fixed = T)
colnames(countdata) <- gsub("..", "", colnames(countdata), fixed = T)

# Remove length/char columns
countdata <- countdata[ ,c(-1:-5)]

# Make sure ID's are correct
# head(countdata)

# - - - - - - - - - - - - - 
# Import metadata file
# - - - - - - - - - - - - - 

# Import and make row names the matching sampleID's from the countdata
metadata <- read.delim("metadata.txt", sep=",",row.names = 1)

# Add sampleID's to the mapping file
metadata$sampleid <- row.names(metadata)

# Reorder sampleID's to match featureCounts column order. 
metadata <- metadata[match(colnames(countdata), metadata$sampleid), ]

# Make sure ID's are correct
# head(metadata)


# - - - - - - - - - - - - - 
# Run Deseq2
# - - - - - - - - - - - - - 
# - countData : count dataframe
# - colData : sample metadata in the dataframe with row names as sampleID's
# - design : The design of the comparisons to use. 
#            Use (~) before the name of the column variable to compare
print("Run DESeq Matrix")
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = metadata,
                                 design = ~Group)


# Find differential expressed genes
# Run DESEq2
print("Run DESeq")
ddsMat <- DESeq(ddsMat)

# Compute pairwise comparison for both 3 possibilities
res.CT.T1 <- results(ddsMat, contrast = c("Group", "CT", "T1"), pAdjustMethod = "fdr", alpha = 0.05)
res.CT.T2 <- results(ddsMat, contrast = c("Group", "CT", "T2"), pAdjustMethod = "fdr", alpha = 0.05)
res.T1.T2 <- results(ddsMat, contrast = c("Group", "T1", "T2"), pAdjustMethod = "fdr", alpha = 0.05)

# Get results from testing with Likelihood Ratio Test for 3 groups
ddsLRT <- DESeq(ddsMat, test="LRT", reduced = ~1)
resLRT <- results(ddsLRT)

# Get results from testing with FDR adjust pvalues
results <- results(ddsMat, pAdjustMethod = "fdr", alpha = 0.05)

# Generate summary of testing. 
# summary(results)

# Check directionality of fold change
# mcols(results, use.names = T)


# - - - - - - - - - - - - - 
# Gene annotation
# - - - - - - - - - - - - - 
library(AnnotationDbi)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

# Add gene full name
results$description <- mapIds(x = org.Mm.eg.db,
                              keys = row.names(results),
                              column = "GENENAME",
                              keytype = "SYMBOL",
                              multiVals = "first")

res.CT.T1$description <- mapIds(x = org.Mm.eg.db,
                              keys = row.names(res.CT.T1),
                              column = "GENENAME",
                              keytype = "SYMBOL",
                              multiVals = "first")

res.CT.T2$description <- mapIds(x = org.Mm.eg.db,
                              keys = row.names(res.CT.T2),
                              column = "GENENAME",
                              keytype = "SYMBOL",
                              multiVals = "first")

res.T1.T2$description <- mapIds(x = org.Mm.eg.db,
                              keys = row.names(res.T1.T2),
                              column = "GENENAME",
                              keytype = "SYMBOL",
                              multiVals = "first")

resLRT$description <- mapIds(x = org.Mm.eg.db,
                              keys = row.names(resLRT),
                              column = "GENENAME",
                              keytype = "SYMBOL",
                              multiVals = "first")

# Add gene symbol
results$symbol <- row.names(results)

res.CT.T1$symbol <- row.names(res.CT.T1)

res.CT.T2$symbol <- row.names(res.CT.T2)

res.T1.T2$symbol <- row.names(res.T1.T2)

resLRT$symbol <- row.names(resLRT)

# Add ENTREZ ID
results$entrez <- mapIds(x = org.Mm.eg.db,
                         keys = row.names(results),
                         column = "ENTREZID",
                         keytype = "SYMBOL",
                         multiVals = "first")

res.CT.T1$entrez <- mapIds(x = org.Mm.eg.db,
                         keys = row.names(res.CT.T1),
                         column = "ENTREZID",
                         keytype = "SYMBOL",
                         multiVals = "first")

res.CT.T2$entrez <- mapIds(x = org.Mm.eg.db,
                         keys = row.names(res.CT.T2),
                         column = "ENTREZID",
                         keytype = "SYMBOL",
                         multiVals = "first")

res.T1.T2$entrez <- mapIds(x = org.Mm.eg.db,
                         keys = row.names(res.T1.T2),
                         column = "ENTREZID",
                         keytype = "SYMBOL",
                         multiVals = "first")

resLRT$entrez <- mapIds(x = org.Mm.eg.db,
                         keys = row.names(resLRT),
                         column = "ENTREZID",
                         keytype = "SYMBOL",
                         multiVals = "first")

# Subset for only significant genes (q < 0.05)
results_sig <- subset(results, pvalue < 0.01)

res.CT.T1_sig <- subset(res.CT.T1, pvalue < 0.01)

res.CT.T2_sig <- subset(res.CT.T2, pvalue < 0.01)

res.T1.T2_sig <- subset(res.T1.T2, pvalue < 0.01)

resLRT_sig <- subset(resLRT, pvalue < 0.01)

# Extracting Rik and Gm genes that has been annotated but have no canonical name
results_sig <- results_sig[!grepl("Gm[0-9]{3}", results_sig$symbol),]
results_sig <- results_sig[!grepl("[0-9]Rik", results_sig$symbol),]

res.CT.T1_sig <- res.CT.T1_sig[!grepl("Gm[0-9]{3}", res.CT.T1_sig$symbol),]
res.CT.T1_sig <- res.CT.T1_sig[!grepl("[0-9]Rik", res.CT.T1_sig$symbol),]

res.CT.T2_sig <- res.CT.T2_sig[!grepl("Gm[0-9]{3}", res.CT.T2_sig$symbol),]
res.CT.T2_sig <- res.CT.T2_sig[!grepl("[0-9]Rik", res.CT.T2_sig$symbol),]

res.T1.T2_sig <- res.T1.T2_sig[!grepl("Gm[0-9]{3}", res.T1.T2_sig$symbol),]
res.T1.T2_sig <- res.T1.T2_sig[!grepl("[0-9]Rik", res.T1.T2_sig$symbol),]

resLRT_sig <- resLRT_sig[!grepl("Gm[0-9]{3}", resLRT_sig$symbol),]
resLRT_sig <- resLRT_sig[!grepl("[0-9]Rik", resLRT_sig$symbol),]

# head(results_sig)
# print("Generate normalized counts")
# Write normalized gene counts to a .txt file
# write.table(x = as.data.frame(counts(ddsMat), normalized = T), 
#             file = 'normalized_counts.txt', 
#             sep = '\t', 
#             quote = F,
#             col.names = NA)
# print("Generate normalized counts sign")
# Write significant normalized gene counts to a .txt file
# write.table(x = counts(ddsMat[row.names(results_sig)], normalized = T), 
#             file = 'normalized_counts_significant.txt', 
#             sep = '\t', 
#             quote = F, 
#             col.names = NA)
# print("Generate gene annotated")
# Write the annotated results table to a .txt file
# write.table(x = as.data.frame(results), 
#             file = "results_gene_annotated.txt", 
#             sep = '\t', 
#             quote = F,
#             col.names = NA)
print("Generate gene annotated sign")
# Write significant annotated results table to a .txt file
write.table(x = as.data.frame(results_sig), 
            file = "DESeq3G_gene_annotated_significant.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)

write.table(x = as.data.frame(res.CT.T1_sig), 
            file = "CTvT1_gene_annotated_significant.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)

write.table(x = as.data.frame(res.CT.T2_sig), 
            file = "CTvT2_gene_annotated_significant.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)

write.table(x = as.data.frame(res.T1.T2_sig), 
            file = "T1vT2_gene_annotated_significant.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)

write.table(x = as.data.frame(resLRT_sig), 
            file = "resLRT3G_gene_annotated_significant.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)

# - - - - - - - - - - - - -
# Venn Diagram
# - - - - - - - - - - - - -
# Write tables used for Venn diagram
write.table(row.names(res.T1.T2_sig), "T1vsT2_genes_VENN.txt", sep="\n",row.names = F, quote = F)
write.table(row.names(res.CT.T2_sig), "CTvsT2_gene_VENN.txt", sep="\n",row.names = F, quote = F)
write.table(row.names(res.CT.T1_sig), "CTvsT1_gene_VENN.txt", sep="\n",row.names = F, quote = F)

# - - - - - - - - - - - - - 
# Heatmap plot
# - - - - - - - - - - - - - 
# Load libraries
library(pheatmap) 
library(RColorBrewer) 

# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)

# Gather 30 significant genes and make matrix
mat <- assay(ddsMat_rlog[row.names(results_sig)])[1:30, ]

# Choose which column variables you want to annotate the columns by.
annotation_col = data.frame(
  Group = factor(colData(ddsMat_rlog)$Group),
  Replicate = factor(colData(ddsMat_rlog)$Replicate),
  row.names = colData(ddsMat_rlog)$sampleid
)

# Specify colors you want to annotate the columns by.
ann_colors = list(
  Group = c("CT" = "blue", "T1" = "orange", "T2" = "black"),
  Replicate = c(R1 = "red", R2 = "green")
)
print("Generate heatmap")
# Make Heatmap with pheatmap function.
# See more in documentation for customization
pheatmap(mat = mat, 
         color = colorRampPalette(brewer.pal(9, "YlOrBr"))(255), 
         scale = "row", 
         annotation_col = annotation_col, 
         annotation_colors = ann_colors, 
         show_colnames = F,
         filename="heatmap.png",
         silent=TRUE)


# - - - - - - - - - - - - - 
# Volcano plot
# - - - - - - - - - - - - - 
# Load libraries
library(ggplot2)
library(RColorBrewer)

# Gather Log-fold change and FDR-corrected pvalues from deseq2 results
data <- data.frame(pval = -log10(resLRT$padj), 
                   lfc = resLRT$log2FoldChange, 
                   row.names = row.names(resLRT))

# Remove any rows that have NA as an entry
data <- na.omit(data)

# Color the points which are up or down
## If fold-change > 0 and pvalue > 1.3 (Increased significant)
## If fold-change < 0 and pvalue > 1.3 (Decreased significant)
data <- mutate(data, color = case_when(data$lfc > 1 & data$pval > 1.3 ~ "Increased",
                                       data$lfc < -1 & data$pval > 1.3 ~ "Decreased",
                                       data$pval < 1.3 ~ "nonsignificant"))
# Make a basic ggplot2 object with x-y values
vol <- ggplot(data, aes(x = lfc, y = pval, color = color))

# Add ggplot2 layers
volcano_plot <- vol +   
  ggtitle(label = "Volcano Plot", subtitle = "Colored by fold-change direction") +
  geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
  scale_color_manual(name = "Directionality",
                     values = c(Increased = "#008B00", Decreased = "#CD4F39", nonsignificant = "darkgray")) +
  theme_bw(base_size = 14) + # change overall theme
  theme(legend.position = "right") + # change the legend
  xlab(expression(log[2]("case1" / "case2"))) + # Change X-Axis label
  ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
  geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
  scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values
print("Generate volcano")
ggsave("volcano.png", plot=volcano_plot, device="png")


# - - - - - - - - - - - - - - -
# Pathway analysis of DE genes
# - - - - - - - - - - - - - - -

# Load required libraries
library(clusterProfiler)
library(ReactomePA)
library(KEGG.db)
library(DOSE)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

# Remove any genes that do not have any entrez identifiers
results_sig_entrez <- subset(resLRT_sig, is.na(entrez) == FALSE)

# Create a matrix of gene log2 fold changes
gene_matrix <- results_sig_entrez$log2FoldChange

# Add the entrezID's as names for each logFC entry
names(gene_matrix) <- results_sig_entrez$entrez


# - - - - - - - - - - - - -
# Enrich with KEGG database
# - - - - - - - - - - - - -
print("Generate KEGG barplot")

kegg_enrich <- enrichKEGG(gene = names(gene_matrix),
                          organism = 'mmu',
                          pvalueCutoff = 0.05)

# Get table of results
kegg_table <- head(as.data.frame(kegg_enrich), n=10) %>% arrange(desc(-log10(pvalue)))

# KEGG plot
kegg_bar <- ggplot(kegg_table, aes(x=reorder(Description, -log10(pvalue)), y=-log10(pvalue))) + geom_bar(stat='identity', fill="#6B2525") + geom_col(width=0.7) + labs(title="KEGG Enrichment Pathways", x="Termos de KEGG") + coord_flip()

ggsave("kegg_bar.png", plot=kegg_bar, device="png")


# - - - - - - - - - - - - -
# Enrich with GO
# - - - - - - - - - - - - -
print("Generate GO barplot")

go_enrich <- enrichGO(gene = names(gene_matrix), 
                      OrgDb = "org.Mm.eg.db",
                      ont = "BP",
                      pvalueCutoff = 0.05)

# Get table of results
go_table <- head(as.data.frame(go_enrich), n=10) %>% arrange(desc(-log10(pvalue)))

# Plot results
go_bar <- ggplot(go_table, aes(x=reorder(Description, -log10(pvalue)), y=-log10(pvalue))) + geom_bar(stat='identity', fill="#157296") + geom_col(width=0.7) + labs(title="GO Biological Pathways", x="Termos de GO") + coord_flip()

ggsave("go_bar.png", plot=go_bar, device="png")