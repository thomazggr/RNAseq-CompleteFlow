# This script has been firstly made by Thomas W. Battaglia (https://github.com/twbattaglia/RNAseq-workflow)
# It has been shortened and adapted for DE analysis and data visualization only
# Refactored creating functions and adaptation for ~ 3+ groups  analysis computing pairwise and LRT tests
# Adapted by Thomaz Guadagnini (https://github.com/ThomazGR | thomaz@vivaldi.net)

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
idcd <- as.character(readline(prompt="Count data ID are correct? [y/n] "))
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
idmd <- as.character(readline(prompt="Count data ID are correct? [y/n] "))
if(idmd == "y"){ } else if(idmd == "n"){ stop("Metadata IDs are not correct. Finishing execution.") }

runner <- function(){
  if(length(unique(metadata$Group)) > 2){
    ddsMat <- deMatrix(countdata, metadata)
    results <- de3gpPWC_LRT(ddsMat, unique(metadata$Group))
    results_sig <- geneAnnotation(results, "mmu")
    generateTables(results_sig)
    vennTables(results_sig)
    generateHeatmap(ddsMat, results_sig)
    generateVolcanos(results)
    pathwayAnalysis(results_sig, "mmu")
    print("Everything has ended!")
  } else if(length(unique(metadata$Group)) == 2){
    ddsMat <- deMatrix(countdata, metadata)
    results <- de2gp(ddsMat)
    results_sig <- geneAnnotation(results, "mmu")
    generateTables(results_sig)
    vennTables(results_sig)
    generateHeatmap(ddsMat, results_sig)
    generateVolcanos(results)
    pathwayAnalysis(results_sig, "mmu")
    print("Everything has ended!")
  }
}

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
  print("|-----DONE-----|")
  # Find differential expressed genes
  # Run DESEq2
  print("Run DESeq")
  ddsMat <- DESeq(ddsMat)
  print("|-----DONE-----|")
  return(ddsMat)
}

# - - - - - - - - - - - - - 
# Results for 2 groups only
# - - - - - - - - - - - - - 
# - ddsMat : Matrix from DESeq2 function
de2gp <- function(ddsMat) {
  print("Get results from Walden's test for 2 groups")
  results <- list(results(ddsMat, pAdjustMethod = "fdr", alpha = 0.05))
  names(results) <- paste(unique(metadata$Group), collapse = "vs")
  print("|-----DONE-----|")
  return(results)
}

# - - - - - - - - - - - - - 
# Gather results for ~+3groups
# - - - - - - - - - - - - - 
# - ddsMat : Matrix from DESeq2 function
# - groups : Vector with string as groups, e.g. c("CT", "G1", "G2") or pass unique(metadata$Group)
de3gpPWC_LRT <- function(ddsMat, groups){
  print("Get results from pairwise comparison for 3 groups")

  names <- c(paste0(groups[1],"vs",groups[2])
            , paste0(groups[1],"vs",groups[3])
            , paste0(groups[2],"vs",groups[3])
            , "resLRT")

  results <- list()
  
  results[1] <- results(ddsMat, contrast = c("Group", groups[1], groups[2]), pAdjustMethod = "fdr", alpha = 0.05)
  results[2] <- results(ddsMat, contrast = c("Group", groups[1], groups[3]), pAdjustMethod = "fdr", alpha = 0.05)
  results[3] <- results(ddsMat, contrast = c("Group", groups[2], groups[3]), pAdjustMethod = "fdr", alpha = 0.05)
  print("|-----DONE-----|")
  
  print("Compute LRT for 3+ groups and save results")
  ddsLRT <- DESeq(ddsMat, test="LRT", reduced= ~1)
  results[4] <- results(ddsLRT)
  names(results) <- names
  print("|-----DONE-----|")
  return(results)
}

# - - - - - - - - - - - - - 
# Gene annotation and retrieves significant genes
# - - - - - - - - - - - - - 
# - results : Results gathered from 2 groups only or +3 groups functions
# - organism : Passing organism to be used for annotation. "hsa" for homo sapiens or "mmu" mus musculus
geneAnnotation <- function(results, organism){
  # Check the organism passed to select a DB
  if (organism == "mmu") { db <- org.Mm.eg.db } else if(organism == "hsa") { db <- org.Hs.eg.db }
  
  # Start loop in every result passed
  print("Get gene description, symbol and ENTREZID for every significant result. Also filter Gm and Rik genes.")
  for (idx in 1:length(results)){
    results[[idx]] <- subset(results[[idx]], pvalue < 0.01)
    
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
    
    results[[idx]] <- results[[idx]][!grepl("Gm[0-9]{3}", results[[idx]]$symbol),]
    results[[idx]] <- results[[idx]][!grepl("[0-9]Rik", results[[idx]]$symbol),]
  }
  for(i in 1:length(names(results))){ names(results)[i] <- paste0(names(results)[i], "_sig") }
  print("|-----DONE-----|")
  return(results)
}

# - - - - - - - - - - - - - 
# Generate table(s) from significant genes
# - - - - - - - - - - - - - 
# - res_sig : Gathered results from significant genes
generateTables <- function(res_sig){
  print("Generate annotated significant genes")
  for (idx in 1:length(res_sig)){
    write.table(x = as.data.frame(res_sig[[idx]]), 
                file = paste0(names(res_sig)[idx],"_annotated.txt"), 
                sep = '\t', 
                quote = F,
                col.names = NA)
  }
  print("|-----DONE-----|")
}

# - - - - - - - - - - - - - 
# Generate table(s) with gene name for external venn diagram
# - - - - - - - - - - - - - 
# - res_sig : Gathered results from significant genes
vennTables <- function(res_sig){
  print("Writing gene names table for easy copy-paste")
  for (idx in 1:length(res_sig)){
    write.table(row.names(res_sig[[idx]]), paste0("Gene_names_",names(res_sig)[idx],".txt"), sep="\n",row.names = F, quote = F)
    
  }
  print("|-----DONE-----|")
}


# - - - - - - - - - - - - - 
# Heatmap plot
# - - - - - - - - - - - - - 
# - ddsMat : Matrix generated in the beginning
# - res_sig : Gathered results from significant genes
generateHeatmap <- function(ddsMat, res_sig){
  if(length(res_sig) > 2){
    results_sig <- res_sig[[length(res_sig)]]
    gps <- unique(metadata$Group)
    # Specify colors you want to annotate the columns by.
    
  }  else if(length(res_sig) == 1) {
      results_sig <- res_sig[[1]]
      gps <- unique(metadata$Group)
      # Specify colors you want to annotate the columns by.
      }
  
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
  ann_colors = list(
      Group = c("Control" = "blue", "Group1" = "orange"), # Add "Group2" = "black" if 3 groups are being compared
      Replicate = c(R1 = "red", R2 = "green")
    )
  nms <- paste(gps, collapse = "-")
  print("Generate heatmap")
  # Make Heatmap with pheatmap function.
  # See more in documentation for customization
  pheatmap(mat = mat, 
           color = colorRampPalette(brewer.pal(9, "YlOrBr"))(255), 
           scale = "row", 
           annotation_col = annotation_col, 
           annotation_colors = ann_colors, 
           show_colnames = F,
           filename=paste0("heatmap_",nms,".png"),
           silent=TRUE)
  print("|-----DONE-----|")
}

# - - - - - - - - - - - - - 
# Volcano plot
# - - - - - - - - - - - - - 
# - results : Results gathered from 2 groups only or +3 groups functions
generateVolcanos <- function(results){
  if(length(results) > 2){ res <- results[-4] } else{ res <- results }
  print("Starting generating volcano plot(s)")
  for (idx in 1:length(res)){
    # Gather Log-fold change and FDR-corrected pvalues from deseq2 results
    data <- data.frame(pval = -log10(res[[idx]]$padj), 
                      lfc = res[[idx]]$log2FoldChange, 
                      row.names = row.names(res[[idx]]))
  
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
      xlab(expression(log[2](paste0(names(res)[idx])))) + # Change X-Axis label
      ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
      geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
      scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values
    print(paste("Generate volcano for", names(res)[idx]))
    ggsave(paste0("volcano",names(res)[idx],".png"), plot=volcano_plot, device="png")
  }
  print("|-----DONE-----|")
}

# - - - - - - - - - - - - - 
# Pathway analysis of DE genes
# - - - - - - - - - - - - - 
# - res_sig : Gathered results from significant genes
# - organism : Passing organism to be used for annotation. "hsa" for homo sapiens or "mmu" mus musculus
pathwayAnalysis <- function(results_sig, organism){
  print("Preprocessing LRT results data if exists")
  if(length(results_sig) > 2){
    res_sig <- results_sig[["resLRT_sig"]]
  }  else if(length(results_sig) == 1) {
    res_sig <- results_sig[[1]]
  }
  # Remove any genes that do not have any entrez identifiers
  res_sig_entrez <- subset(res_sig, is.na(entrez) == FALSE)
  # Create a matrix of gene log2 fold changes
  gene_matrix <- res_sig_entrez$log2FoldChange
  # Add the entrezID's as names for each logFC entry
  names(gene_matrix) <- res_sig_entrez$entrez
  
  gps <- unique(metadata$Group)
  nms <- paste(gps, collapse = "-")
  
  # Check the organism passed to select a DB
  if (organism == "mmu") { 
    db <- "org.Mm.eg.db"
  } else if(organism == "hsa") { 
      db <- "org.Hs.eg.db"
  } else { 
        stop("Need to pass hsa or mmu as organism") 
    }
  # - - - - - - - - - - - - -
  # Enrich with KEGG database
  # - - - - - - - - - - - - -
  print("Generate KEGG pathway results")
  kegg_enrich <- enrichKEGG(gene = names(gene_matrix),
                            organism = paste0(organism),
                            pvalueCutoff = 0.05)
  
  # Get table of results
  kegg_table <- head(as.data.frame(kegg_enrich), n=10) %>% 
                arrange(desc(-log10(pvalue)))
  
  # KEGG plot
  kegg_bar <- ggplot(kegg_table, aes(x=reorder(Description, -log10(pvalue)), y=-log10(pvalue))) + 
              geom_bar(stat='identity', fill="#6B2525") + 
              geom_col(width=0.7) + 
              labs(title="KEGG Enrichment Pathways", x="Termos de KEGG") + 
              coord_flip()
  
  ggsave(paste0("kegg_bar_",nms,".png"), plot=kegg_bar, device="png")
  print("|-----DONE-----|")
  
  # - - - - - - - - - - - - -
  # Enrich with GO
  # - - - - - - - - - - - - -
  print("Generate GO barplot")
  
  go_enrich <- enrichGO(gene = names(gene_matrix), 
                        OrgDb = paste0(db),
                        ont = "BP",
                        pvalueCutoff = 0.05)
  
  # Get table of results
  go_table <- head(as.data.frame(go_enrich), n=10) %>% 
              arrange(desc(-log10(pvalue)))
  
  # Plot results
  go_bar <- ggplot(go_table, aes(x=reorder(Description, -log10(pvalue)), y=-log10(pvalue))) + 
            geom_bar(stat='identity', fill="#157296") + 
            geom_col(width=0.7) + 
            labs(title="GO Biological Pathways", x="Termos de GO") + 
            coord_flip()
  
  ggsave(paste0("go_bar_",nms,".png"), plot=go_bar, device="png")
  print("|-----DONE-----|")
}

runner()