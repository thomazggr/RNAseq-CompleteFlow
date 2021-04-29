# Adapted script for DE analysis and data visualization

# Install R base and R dev
# Install BiocManager first
# sudo apt install r-cran-xml r-cran-xml2 libxml2-dev libcurl4-openssl-dev libssl-dev 
# install.packages(c("httr", "curl", "RCurl", "openssl", "XML"))
# BiocManager::install(c("DESeq2","ggplot2","clusterProfiler","AnnotationDbi","ReactomePA","gage","ggsci","dplyr","DOSE","org.Hs.eg.db","org.Mm.eg.db","pheatmap","KEGG.db","RColorBrewer"))

# Load general libraries
library(ggplot2)
library(dplyr)
library(gage)
library(ggsci)

# Differential expression analysis library
library(DESeq2)

# Gene annotation libraries
library(AnnotationDbi)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

# Heatmap libraries
library(pheatmap) 
library(RColorBrewer)

# Volcano libraries
library(ggplot2)
library(RColorBrewer)

# Pathway analysis libraries
library(clusterProfiler)
library(ReactomePA)
library(KEGG.db)
library(DOSE)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

# - - - - - - - - - - - - - 
# Import gene counts table
# - - - - - - - - - - - - - 
# - skip first row (general command info)
# - make row names the gene identifiers
countdata <- read.table("./DE_analysis/final_counts_m24.txt", header = TRUE, skip = 1, row.names = 1)

# Remove .sam + '..' from column identifiers
colnames(countdata) <- gsub(".sam", "", colnames(countdata), fixed = T)
colnames(countdata) <- gsub("X", "", colnames(countdata), fixed = T)

# Remove length/char columns
countdata <- countdata[ ,c(-1:-5)]

# Make sure ID's are correct
print(head(countdata))
idcd <- readline(prompt="Count data ID are correct? [y/n] ")
if(idcd == "y"){ } else if(idcd == "n"){ stop("Count data IDs are not correct. Finishing execution.") }

# - - - - - - - - - - - - - 
# Import metadata file
# - - - - - - - - - - - - - 

# Import and make row names the matching sampleID's from the countdata
metadata <- read.delim("./DE_analysis/metadata.txt", sep=",",row.names = 1)

# Add sampleID's to the mapping file
metadata$sampleid <- row.names(metadata)

# Reorder sampleID's to match featureCounts column order. 
metadata <- metadata[match(colnames(countdata), metadata$sampleid), ]

# Make sure ID's are correct
print(head(metadata))
idmd <- readline(prompt="Count data ID are correct? [y/n] ")
if(idmd == "y"){ } else if(idmd == "n"){ stop("Metadata IDs are not correct. Finishing execution.") }


# - - - - - - - - - - - - - 
# Run Deseq2
# - - - - - - - - - - - - - 
# - countData : count dataframe
# - colData : sample metadata in the dataframe with row names as sampleID's
# - design : The design of the comparisons to use. 
#            Use (~) before the name of the column variable to compare
deMatrix <- function(countdata, metadata) {
  print("Run DESeq Matrix")
  ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                  colData = metadata,
                                  design = ~Group)

  # Find differential expressed genes
  # Run DESEq2
  print("Run DESeq")
  ddsMat <- DESeq(ddsMat)
  return(ddsMat)
}

de2gp <- function(ddsMat) {
  print("Get results from Walden's test for 2 groups")
  results <- results(ddsMat, pAdjustMethod = "fdr", alpha = 0.05)
  return(results)
}

de3gpPWC_LRT <- function(ddsMat, groups){
  print("Get results from pairwise comparison for 3 groups")

  names <- c(paste0(groups[1],"vs",groups[2])
            , paste0(groups[1],"vs",groups[3])
            , paste0(groups[2],"vs",groups[3])
            , "resLRT")

  assign(names[1], results(ddsMat, contrast = c("Group", groups[1], groups[2]), pAdjustMethod = "fdr", alpha = 0.05))
  assign(names[2], results(ddsMat, contrast = c("Group", groups[1], groups[3]), pAdjustMethod = "fdr", alpha = 0.05))
  assign(names[3], results(ddsMat, contrast = c("Group", groups[2], groups[3]), pAdjustMethod = "fdr", alpha = 0.05))
  
  print("Get results from likelihood ratio test for 3 groups")

  ddsLRT <- DESeq(ddsMat, test="LRT", reduced= ~1)
  resLRT <- results(ddsLRT)

  results <- lapply(names, get)

  names(results) <- names

  return(results)
}

# # Compute pairwise comparison for both 3 possibilities
# res.CT.Prop <- results(ddsMat, contrast = c("Group", "CT", "Prop"), pAdjustMethod = "fdr", alpha = 0.05)
# res.CT.Res <- results(ddsMat, contrast = c("Group", "CT", "Res"), pAdjustMethod = "fdr", alpha = 0.05)
# res.Prop.Res <- results(ddsMat, contrast = c("Group", "Prop", "Res"), pAdjustMethod = "fdr", alpha = 0.05)

# # Get results from testing with Likelihood Ratio Test for 3 groups
# ddsLRT <- DESeq(ddsMat, test="LRT", reduced = ~1)
# resLRT <- results(ddsLRT)

# # Get results from testing with FDR adjust pvalues
# results <- results(ddsMat, pAdjustMethod = "fdr", alpha = 0.05)

# Generate summary of testing. 
# summary(results)

# Check directionality of fold change
# mcols(results, use.names = T)

# - - - - - - - - - - - - - 
# Gene annotation
# - - - - - - - - - - - - - 

geneAnnotation <- function(results, organism){
  if (organism == "mmu") { db <- org.Mm.eg.db } else if(organism == "hsa") { db <- org.Hs.eg.db }
  
  for (idx in 1:length(results)){
    results[[idx]]$description <- mapIds(x = db,
                              keys = row.names(results[[idx]]),
                              column = "GENENAME",
                              keytype = "SYMBOL",
                              multiVals = "first")

    results[[idx]]$symbol <- row.names(results[[idx]])

    results[[idx]]$entrez <- mapIds(x = org.Mm.eg.db,
                         keys = row.names(results[[idx]]),
                         column = "ENTREZID",
                         keytype = "SYMBOL",
                         multiVals = "first")

    signame <- paste0(names(results)[idx], "_sig")

    assign(paste0(signame), subset(results[[idx]], pvalue < 0.01))
    eval(as.character(signame)) <- eval(as.character(signame))[!grepl("Gm[0-9]{3}"
                                                              , eval(as.character(signame))$symbol),]
    eval(as.character(signame)) <- eval(as.character(signame))[!grepl("[0-9]Rik"
                                                              , eval(as.character(signame))$symbol),]
    
    names[idx] <- signame
  }
  results_sig <- lapply(names, get)
  names(results_sig) <- names
  return(results_sig)
}

stop()
# # Add gene full name
# results$description <- mapIds(x = org.Mm.eg.db,
#                               keys = row.names(results),
#                               column = "GENENAME",
#                               keytype = "SYMBOL",
#                               multiVals = "first")

# res.CT.Prop$description <- mapIds(x = org.Mm.eg.db,
#                               keys = row.names(res.CT.Prop),
#                               column = "GENENAME",
#                               keytype = "SYMBOL",
#                               multiVals = "first")

# res.CT.Res$description <- mapIds(x = org.Mm.eg.db,
#                               keys = row.names(res.CT.Res),
#                               column = "GENENAME",
#                               keytype = "SYMBOL",
#                               multiVals = "first")

# res.Prop.Res$description <- mapIds(x = org.Mm.eg.db,
#                               keys = row.names(res.Prop.Res),
#                               column = "GENENAME",
#                               keytype = "SYMBOL",
#                               multiVals = "first")

# resLRT$description <- mapIds(x = org.Mm.eg.db,
#                               keys = row.names(resLRT),
#                               column = "GENENAME",
#                               keytype = "SYMBOL",
#                               multiVals = "first")

# # Add gene symbol
# results$symbol <- row.names(results)

# res.CT.Prop$symbol <- row.names(res.CT.Prop)

# res.CT.Res$symbol <- row.names(res.CT.Res)

# res.Prop.Res$symbol <- row.names(res.Prop.Res)

# resLRT$symbol <- row.names(resLRT)

# # Add ENTREZ ID
# results$entrez <- mapIds(x = org.Mm.eg.db,
#                          keys = row.names(results),
#                          column = "ENTREZID",
#                          keytype = "SYMBOL",
#                          multiVals = "first")

# res.CT.Prop$entrez <- mapIds(x = org.Mm.eg.db,
#                          keys = row.names(res.CT.Prop),
#                          column = "ENTREZID",
#                          keytype = "SYMBOL",
#                          multiVals = "first")

# res.CT.Res$entrez <- mapIds(x = org.Mm.eg.db,
#                          keys = row.names(res.CT.Res),
#                          column = "ENTREZID",
#                          keytype = "SYMBOL",
#                          multiVals = "first")

# res.Prop.Res$entrez <- mapIds(x = org.Mm.eg.db,
#                          keys = row.names(res.Prop.Res),
#                          column = "ENTREZID",
#                          keytype = "SYMBOL",
#                          multiVals = "first")

# resLRT$entrez <- mapIds(x = org.Mm.eg.db,
#                          keys = row.names(resLRT),
#                          column = "ENTREZID",
#                          keytype = "SYMBOL",
#                          multiVals = "first")

# # Subset for only significant genes (q < 0.05)
# results_sig <- subset(results, pvalue < 0.01)

# res.CT.Prop_sig <- subset(res.CT.Prop, pvalue < 0.01)

# res.CT.Res_sig <- subset(res.CT.Res, pvalue < 0.01)

# res.Prop.Res_sig <- subset(res.Prop.Res, pvalue < 0.01)

# resLRT_sig <- subset(resLRT, pvalue < 0.01)

# # Extracting Rik and Gm genes that has been annotated but have no canonical name
# results_sig <- results_sig[!grepl("Gm[0-9]{3}", results_sig$symbol),]
# results_sig <- results_sig[!grepl("[0-9]Rik", results_sig$symbol),]

# res.CT.Prop_sig <- res.CT.Prop_sig[!grepl("Gm[0-9]{3}", res.CT.Prop_sig$symbol),]
# res.CT.Prop_sig <- res.CT.Prop_sig[!grepl("[0-9]Rik", res.CT.Prop_sig$symbol),]

# res.CT.Res_sig <- res.CT.Res_sig[!grepl("Gm[0-9]{3}", res.CT.Res_sig$symbol),]
# res.CT.Res_sig <- res.CT.Res_sig[!grepl("[0-9]Rik", res.CT.Res_sig$symbol),]

# res.Prop.Res_sig <- res.Prop.Res_sig[!grepl("Gm[0-9]{3}", res.Prop.Res_sig$symbol),]
# res.Prop.Res_sig <- res.Prop.Res_sig[!grepl("[0-9]Rik", res.Prop.Res_sig$symbol),]

# resLRT_sig <- resLRT_sig[!grepl("Gm[0-9]{3}", resLRT_sig$symbol),]
# resLRT_sig <- resLRT_sig[!grepl("[0-9]Rik", resLRT_sig$symbol),]

# - - - - - - - - - - - - - 
# Not needed tables for now
# - - - - - - - - - - - - - 
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
# - - - - - - - - - - - - - 

print("Generate gene annotated sign")
# Write significant annotated results table to a .txt file
write.table(x = as.data.frame(results_sig), 
            file = "DESeq3G_gene_annotated_significant_prop_res.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)

write.table(x = as.data.frame(res.CT.Prop_sig), 
            file = "CTvProp_gene_annotated_significant.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)

write.table(x = as.data.frame(res.CT.Res_sig), 
            file = "CTvRes_gene_annotated_significant.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)

write.table(x = as.data.frame(res.Prop.Res_sig), 
            file = "PropvRes_gene_annotated_significant.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)

write.table(x = as.data.frame(resLRT_sig), 
            file = "resLRT3G_gene_annotated_significant_prop_res.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)

# - - - - - - - - - - - - -
# Venn Diagram
# - - - - - - - - - - - - -
# > 3 tables that gives us the gene names from each pairwise comparison in case of 3 groups
vennTables <- function(results){
  for (idx in 1:length(results)){
    write.table(row.names(results[idx]), paste("",names(results[idx]),"_gene_names.txt"), sep="\n",row.names = F, quote = F)
  }
}
# Older version
write.table(row.names(res.CT.Res_sig), "PropvsRes_gene_names.txt", sep="\n",row.names = F, quote = F)
write.table(row.names(res.CT.Res_sig), "CTvsRes_gene_names.txt", sep="\n",row.names = F, quote = F)
write.table(row.names(res.CT.Prop_sig), "CTvsProp_gene_names.txt", sep="\n",row.names = F, quote = F)

# - - - - - - - - - - - - - 
# Heatmap plot
# - - - - - - - - - - - - - 
# Load libraries

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
  Group = c("CT" = "blue", "Prop" = "orange", "Res" = "black"),
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
         filename="heatmap_Prop_Res.png",
         silent=TRUE)


# - - - - - - - - - - - - - 
# Volcano plot
# - - - - - - - - - - - - - 


volcanos_data <- list(res.CT.Prop, res.CT.Res, res.Prop.Res)
volcanos_names <- list("CTvsProp", "CTvsRes", "PropvsRes")
idx <- 1
for (i in volcanos_data){
  # Gather Log-fold change and FDR-corrected pvalues from deseq2 results
  data <- data.frame(pval = -log10(i$padj), 
                    lfc = i$log2FoldChange, 
                    row.names = row.names(i))

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
    xlab(expression(log[2](volcanos_names[idx]))) + # Change X-Axis label
    ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
    geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
    scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values
  print(paste("Generate volcano for", volcanos_names[idx]))
  ggsave(paste0("volcano",volcanos_names[idx],".png"), plot=volcano_plot, device="png")
  idx <- idx + 1
}



# - - - - - - - - - - - - - - -
# Pathway analysis of DE genes
# - - - - - - - - - - - - - - -
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

ggsave("kegg_bar_Prop_Res.png", plot=kegg_bar, device="png")


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

ggsave("go_bar_Prop_Res.png", plot=go_bar, device="png")