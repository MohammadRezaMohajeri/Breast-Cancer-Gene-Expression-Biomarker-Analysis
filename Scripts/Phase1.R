# Phase 1
# R Code (1-Dataset Analysis - limma)

# Differential Expression Analysis of Breast Cancer Dataset (GSE25055) using limma in R: Comparing Molecular Subtypes and Grading System
#==========================================================================================================================================


# Load required packages
require(GEOquery)    # For downloading and processing GEO datasets
require(limma)       # For performing differential expression analysis
require(tidyverse)   # For data manipulation and visualization
require(plotly)      # For interactive plots


# Function to download a GEO dataset
DownloadGEO <- function(GSEname) {
  # Create a folder for the dataset
  dir.create(GSEname, showWarnings = FALSE)
  
  # Increase memory allocation
  Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10) 
  
  # Download the dataset
  gset <- getGEO(
    GSEname,
    GSEMatrix = TRUE,
    AnnotGPL = TRUE,
    destdir = paste0(getwd(), "/", GSEname)
  )
 
  if(!is.null(gset)) {
    # Extract matrix, phenotype, and annotation
    matrix <- gset[[paste0(GSEname, "_series_matrix.txt.gz")]]@assayData[["exprs"]]
    matrix <- cbind(rownames(matrix), matrix)
    colnames(matrix)[1] <- "ID"
    
    phenotype <- gset[[paste0(GSEname, "_series_matrix.txt.gz")]]@phenoData@data
    phenotype <- cbind(rownames(phenotype), phenotype)
    colnames(phenotype)[1] <- "ID"
    
    annot <- gset[[paste0(GSEname, "_series_matrix.txt.gz")]]@featureData@data
    
    # Write the extracted variables to files
    write.csv(matrix, paste0(GSEname, "/", "matrix.csv"), quote = FALSE, row.names = FALSE)
    write_tsv(phenotype, paste0(GSEname, "/", "phenotype.tsv")) 
    write_tsv(annot, paste0(GSEname, "/", "annot.tsv")) 
    
    # Remove variables from memory
    rm(gset, matrix, phenotype, annot)
    
    message("+++++ Dataset is ready to analyze. +++++")
  }
}

# # Function for performing differential expression analysis
DEGanalysis <- function(projname,
                        compare = "subtype",
                        groups,
                        adjust = "fdr") {
  # Load the matrix data from the CSV file
  matrix <- read_csv(paste0(projname, "/", "matrix.csv"))
  matrix <- data.frame(matrix)
  rownames(matrix) <- matrix[,1]
  matrix <- matrix[,-1]

  # Load the phenotype data from the TSV file
  phenotype <- read_tsv(paste0(projname, "/", "phenotype.tsv"))
  phenotype <- data.frame(phenotype)
  rownames(phenotype) <- phenotype[,1]
  phenotype <- phenotype[,-1]

  # Load the annotation data from the TSV file
  annot <- read_tsv(paste0(projname, "/", "annot.tsv"))
  annot <- data.frame(annot)
  rownames(annot) <- annot[,1]
  annot <- annot[,-1]
   
  (subtype) Extract the relevant column from the phenotype data and rename it as "subtype"
  phecoln <- grep("pam50", phenotype)
  pheno <- data.frame(phenotype[ ,phecoln])
  colnames(pheno) <- "subtype"
  pheno$subtype <- gsub("pam50.+:", "", pheno$subtype)
  rownames(pheno) <- rownames(phenotype)

  # (grade) Extract the relevant column from the phenotype data and rename it as "grade"
  gradecoln <- grep("grade:", phenotype)
  grade <- data.frame(phenotype[ ,gradecoln])
  colnames(grade) <- "grade"
  grade$grade <- gsub("grade:", "", grade$grade)
  grade$grade <- gsub("=.+", "", grade$grade)
  rownames(grade) <- rownames(phenotype)
  
  # Perform log2 transformation if the matrix is not already transformed
  qx <- as.numeric(quantile(matrix, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
  if(LogC){
    matrix[which(matrix <= 0)] <- NaN
    matrix <- log2(matrix)
  }

  # Create a list to store the sample IDs for each group
  if(compare == "subtype"){
    grouplis <- list() 
    for(i in seq_along(groups)){
      name <- paste0(groups[[i]], collapse = "|")
      rownum <- grep(name, pheno$subtype)
      grouplis[[names(groups)[i]]] <- rownames(pheno)[rownum]
    }
  } else{
    grouplis <- list()
    for(i in seq_along(groups)){
      name <- paste0(groups[[i]], collapse = "|")
      rownum <- grep(name, grade$grade)
      grouplis[[names(groups)[i]]] <- rownames(grade)[rownum]
    }
  }
  
  # Rename the matrix file according to the group names specified by the use 
  for(i in seq_along(grouplis)){
    matname <- paste0(grouplis[[i]], collapse = "|")
    matcoln <- grep(matname, colnames(matrix))
    colnames(matrix)[matcoln] <- paste0(names(grouplis)[i], 1:length(grouplis[[i]]))
  } 
  
  # Keep only the samples that are renamed and order them
  matnames <- grep(paste0(names(grouplis), collapse = "|"), colnames(matrix))
  matrix <- matrix[, matnames]
  matrix <- matrix[, order(colnames(matrix))]
  
  # Create a factor variable with 0 and 1 based on the length of groups specified by user
  gname <- list()
  for(i in seq_along(grouplis)){
    number <- length(grep(names(grouplis)[i], colnames(matrix)))
    gname[[i]] <- factor(rep(i - 1, number))
  }
  grouping <- unlist(gname)
  
  # Create the design matrix for the linear model
  design <- model.matrix(~ 0 + grouping)
  colnames(design) <- names(grouplis)
  
  # Fit the linear model
  fit1 <- lmFit(matrix, design)
  l <- length(names(grouplis))
  cts <- paste0(names(grouplis)[l:1], collapse = "-")
  cont.matrix <- makeContrasts(contrasts = cts, levels = design)
  fit2 <- contrasts.fit(fit1, cont.matrix)
  fit2 <- eBayes(fit2)
  
  # Create a top table with all genes
  DEGs <- topTable(fit2, adjust = adjust, number = Inf)
  DEGs <- cbind(rownames(DEGs), DEGs)
  colnames(DEGs)[1] <- "ID" 
  
  # Prepare the annotation variable
  annotation <- annot[, c("Gene.symbol", "Chromosome.location", "GO.Function")]
  annotation <- cbind(rownames(annotation), annotation)
  colnames(annotation)[1] <- "ID" 
  
  # Ensure that the row names of DEGs and annotation are the same
  annotation <- annotation[order(annotation$ID), ]
  DEGs <- DEGs[order(DEGs$ID), ]
  
  # Add the annotation to the top table
  DEGs <- cbind(DEGs,  annotation[, c(2, 3, 4)])
  
  # Removing missing gene symbols
  genemissing <- which(is.na(DEGs$Gene.symbol) == TRUE)
  DEGs <- DEGs[-genemissing, ]
  
  # Write the top table to a TSV file
  write_tsv(DEGs, paste0(projname, "/", "DEGs.tsv"))
  
  message("+++++ Analysis has completed. +++++")
}



# Function for generating a volcano plot
Makevolcano <- function(projname,
                        padjlevel = 0.05,
                        uplogFC = 1,
                        downlogFC = -1 ,
                        ntop = 10) {

  # Read the top table from the TSV file
  table <- read_tsv(paste0(projname, "/", "DEGs.tsv"))
  table <- data.frame(table)
	   
  # Define color palette for the plot
  my_pal <- c("#1B9E77","#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#9A7D0A")

  # Add a "DEGs" column based on the adjusted p-value and log fold change
  volcano <- table %>%
    mutate(AdjustedPvalue = -log10(adj.P.Val)) %>%
    mutate(DEGs = "Not") %>%
    mutate(DEGs = ifelse(AdjustedPvalue > -log10(padjlevel) & logFC > uplogFC, "Upregulated", DEGs)) %>%
    mutate(DEGs = ifelse(AdjustedPvalue > -log10(padjlevel) & logFC < downlogFC, "Downregulated", DEGs))

  # Create a ggplot object, add layers, and annotate the plot
  p <- ggplot(data = volcano, aes(x = logFC, y = AdjustedPvalue, color = DEGs, fill = DEGs, text = paste("ID: ", ID,
                                                                                                       "<br>Gene: ", Gene.symbol,
                                                                                                       "<br>Chromosome: ", Chromosome.location,
                                                                                                       "<br>GO function: ", GO.Function))) +
    labs(x= 'log2 (Fold Change)', y = "-log10(Adjusted P-value)") + 
    geom_point(size = 1, shape = 21) +
    scale_color_manual(values = c(my_pal)) +
    scale_fill_manual(values = c(paste(my_pal, "66" , sep = "" ))) +
    theme_classic() +
    theme(axis.text = element_text(family = "Times", size = 15 , colour = "black"),
          axis.text.x = element_text(family = "Times", colour = "black", size = 15),
          axis.text.y = element_text(family = "Times", colour = "black", size = 15),
          plot.subtitle = element_text(family = "Times", size = 20, colour = "black", hjust = 0.5),
          axis.title.y = element_text(family = "Times", size = rel(1.8), angle = 90),
          axis.title.x = element_text(family = "Times", size = rel(1.8), angle = 00)) +
    labs(subtitle = "Volcano plot")

  # Convert ggplot to plotly for interactivity
  p <- ggplotly(p, tooltip = "text")

  # Calculate the number of up-regulated and down-regulated genes
  up <- volcano[which(volcano$DEGs == "Upregulated"), ]
  up <- up[order(up$logFC, decreasing = T), ]
  down <- volcano[which(volcano$DEGs == "Downregulated"), ]
  down <- down[order(down$logFC, decreasing = F), ]

  # Creating the top table (up and down) based on the number specified (ntop)
  top <- rbind(up[1:ntop, ], down[1:ntop, ])
  top <- top[, c(1:6, 11, 7, 8, 9, 10, 12)]
  colnames(top)[7] <- "n.log10(adj.P.Val)"

  # Save the top table as a CSV file in the project folder
  write.csv(top, paste0(projname, "/", "top.csv"), quote = FALSE, row.names = FALSE)

  # Display a message with the number of up-regulated and down-regulated genes
  message(sprintf("Up-regulated genes: %s; Down-regulated genes: %s", nrow(up), nrow(down)))

# Return the interactive volcano plot
  return(p)
}




# Download GEO dataset with the accession number "GSE25055"
DownloadGEO ("GSE25055")

#============================ Q1 =========================#
# other adjust methods ---> "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"

# (Lum A, Normal-like subtypes) VS (Basal-like, HER2, Lum B)
DEGanalysis(projname = "GSE25055",
            compare = "subtype",
            groups = list(normal = c("LumA", "Normal"), tumor = c("Basal", "Her2", "LumB")),
            adjust = "fdr")

# alpha = 0.01 - logFC 0
volcano1 <- Makevolcano(projname = "GSE25055",
                        padjlevel = 0.01,
                        uplogFC = 1.5,
                        downlogFC = -1.5,
                        ntop = 10)
volcano1

# alpha = 0.01 - logFC 1
volcano2 <- Makevolcano(projname = "GSE25055",
                        padjlevel = 0.00001,
                        uplogFC = 1,
                        downlogFC = -1)
volcano2

# alpha = 0.05 - logFC 0
volcano3 <- Makevolcano(projname = "GSE25055",
                        padjlevel = 0.05,
                        uplogFC = 0,
                        downlogFC = 0)
volcano3

# alpha = 0.05 - logFC 1
volcano4 <- Makevolcano(projname = "GSE25055",
                        padjlevel = 0.05,
                        uplogFC = 1,
                        downlogFC = -1
volcano4

#============================ Q2 =========================#
# Grade 1 vs Grade 3
DEGanalysis(projname = "GSE25055",
            compare = "grade",
            groups = list(normal = c("1"), tumor = c("3")),
            adjust = "fdr")

# alpha = 0.01 - logFC 0
volcano5 <- Makevolcano(projname = "GSE25055",
                        padjlevel = 0.01,
                        uplogFC = 0,
                        downlogFC = 0)
volcano5

# alpha = 0.01 - logFC 1.5
volcano6 <- Makevolcano(projname = "GSE25055",
                        padjlevel = 0.01,
                        uplogFC = 1.5,
                        downlogFC = -1.5)
volcano6

# alpha = 0.05 - logFC 0
volcano7 <- Makevolcano(projname = "GSE25055",
                        padjlevel = 0.05,
                        uplogFC = 0,
                        downlogFC = 0)
volcano7

# alpha = 0.01 - logFC 1
volcano8 <- Makevolcano(projname = "GSE25055",
                        padjlevel = 0.01,
                        uplogFC = 1,
                        downlogFC = -1)
volcano8
