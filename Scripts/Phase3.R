# Phase 3
# R Code For Meta-analysis (4 datasets)

# Integration of Gene Expression Data through Meta-Analysis for Robust Biomarker Discovery in Breast Cancer (Grade 1 vs Grade 3) using ‘’GeneMeta’’ package in R
#================================================================================================================================================================

####################################################### "Package Installation" #######################################################span>
############# Uncomment the following lines if any of the packages are not already installed or need to be reinstalled. #############
# BiocManager packages 
# Install required Bioconductor packages
# BiocManager::install("sva", force = TRUE)
# BiocManager::install("GeneMeta", force = TRUE)
# BiocManager::install("ComplexHeatmap", force = TRUE)
# BiocManager::install("fgsea", force = TRUE)
# BiocManager::install("clusterProfiler", force = TRUE)
			
# Other packages
# Install required CRAN packages:
# install.packages("RColorBrewer")
# install.packages("msigdbr")
# install.packages("ggpubr")
# install.packages("reshape2")
# install.packages("caret")
# install.packages("rgl")

# sva package installation
# Install "sva" package using BiocManager
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("sva")

# VSN package installation
# Install "vsn" package using BiocManager
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("vsn")

# tidyverse package installation
# install.packages("tidyverse")

# GEOquery package installation
# Install "GEOquery" package using BiocManager
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("GEOquery")

# limma package installation
# install.packages("limma")

# plotly package installation
# install.packages("plotly")

# ggvenn package installation
# install.packages("ggvenn")

# ggpubr package installation
# install.packages("ggpubr")

# rgl package installation
# For rgl installation, require remotes package first
# install.packages("remotes")
# require(remotes)
# remotes::install_github("dmurdoch/rgl")


# Package Loading
# Load the required packages for data analysis and visualization.
			
require(tidyverse)       # Package for comprehensive data manipulation and visualization, including dplyr, ggplot2, and other useful packages
require(GEOquery)        # Package for accessing and retrieving gene expression data from the Gene Expression Omnibus (GEO) database
require(reshape2)        # Package for reshaping data from wide to long format and vice versa, useful for data preprocessing
require(caret)           # Package for machine learning and predictive modeling, providing a unified interface for various algorithms and performance evaluation
require(GeneMeta)        # Package for meta-analysis of gene expression data, allowing combining results from multiple studies
require(limma)           # Package for analyzing microarray data using linear models, including differential expression analysis
require(rgl)             # Package for creating interactive 3D visualizations, useful for exploring complex data structures
require(sva)             # Package for surrogate variable analysis, used to identify and adjust for unwanted variation in high-dimensional datasets
require(plotly)          # Package for creating interactive and dynamic plots, including scatter plots, bar plots, and more
require(ggvenn)          # Package for creating Venn diagrams using ggplot2, useful for visualizing set relationships
require(ggpubr)          # Package for creating publication-ready plots using ggplot2, providing additional customization options
require(ComplexHeatmap)  # Package for creating complex heatmaps, allowing visualization of multiple layers of data and annotations
require(RColorBrewer)    # Package for providing color palettes suitable for data visualization
require(msigdbr)         # Package for accessing and using the Molecular Signatures Database (MSigDB) for gene set enrichment analysis
require(fgsea)           # Package for performing gene set enrichment analysis (GSEA) using fast gene set testing algorithms
require(locfit)          # Package for local regression modeling, useful for fitting flexible curves to data
require(vsn)             # Package for variance stabilization and normalization, specifically for microarray data analysis


	
######################################################################################################################################
####################################################### 1- DownloadGEO Function#######################################################
# DownloadGEO: Downloads a Gene Expression Omnibus (GEO) dataset and prepares it for analysis 

# Input (Arguments): 
#   - GSEname: A character vector specifying the GEO dataset name (default: "GSE25065")

# Output: 
#   - Extracted variables saved as separate files (matrix.csv, phenotype.tsv, annot.tsv) 

# Function Description: 
# This function downloads a specified GEO dataset from Gene Expression Omnibus (GEO) and prepares it for analysis. 
# It creates a folder to store the downloaded dataset, increases memory allocation for efficient processing,  
# and downloads the specified GEO dataset. It then extracts the gene expression matrix, phenotype information, 
# and annotation data from the downloaded dataset. The extracted variables are saved as separate files. 
	
DownloadGEO <- function(GSEname = "GSE25065"){
  # Create a folder to store the downloaded dataset 
  dir.create(GSEname, showWarnings = FALSE)
	
  # Increase the memory allocation for efficient processing 
  Sys.setenv("VROOM_CONNECTION_SIZE"=131072 *131072)
  
  #Download dataset 
  gset <- getGEO(GSEname, GSEmatrix = TRUE, AnnotGPL = TRUE,
                 destdir = paste0(getwd(), "/", GSEname))
	
 if(!is.null(gset)){
    # Extract the gene expression matrix, phenotype information, and annotation data: 
    # Extract the gene expression matrix from the downloaded dataset 
    matrix <- gset[[paste0(GSEname, "_series_matrix.txt.gz")]]@assayData[["exprs"]]
	
    # Check if the matrix has undergone log2 transformation 
    qx <- as.numeric(quantile(matrix, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm= TRUE))
    LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
    if(!LogC){
      # Reverse the log2 transformation if applied 
      matrix <- 2^(matrix)
    }

    # Add row names to the matrix and set the column name for the identifier column 
    matrix <- cbind(rownames(matrix), matrix)
    colnames(matrix)[1] <- "ID"

    # Extract the phenotype information from the downloaded dataset 
    phenotype <- gset[[paste0(GSEname, "_series_matrix.txt.gz")]]@phenoData@data
    phenotype <- cbind(rownames(phenotype), phenotype)
    colnames(phenotype)[1] <- "ID"

    # Extract the annotation data from the downloaded dataset 
    annot<- gset[[paste0(GSEname, "_series_matrix.txt.gz")]]@featureData@data
    
	
    # Write the extracted variables to separate files: 
    # Write the gene expression matrix to a CSV file 
    write.csv(matrix, paste0(GSEname, "/", "matrix.csv"), quote = FALSE, row.names = FALSE)

    # Write the phenotype information to a tab-separated values (TSV) file 
    write_tsv(phenotype, paste0(GSEname, "/", "phenotype.tsv"))

    # Write the annotation data to a TSV file 
    write_tsv(annot, paste0(GSEname, "/", "annot.tsv"))
    
    # Remove unnecessary variables from memory 
    rm(gset, matrix, phenotype, annot)

    # Display a message indicating that the dataset is ready for analysis 
    message ("+++++ Datatset is ready to analyse. +++++")
  }
}


######################################################################################################################################
######################################################### 2- ReadGEO Function#########################################################   
# ReadGEO: Reads the gene expression matrix file for a given project 

# Input (Arguments): 
#   - projname: A character vector specifying the project name 

# Output: 
#   - A data frame containing the gene expression matrix 

# Function Description: 
# This function reads the gene expression matrix file for a specified project.  
# It takes the argument `projname`, which is a character vector specifying the project name. 
# The function reads the matrix file in CSV format using the `read_csv` function and stores it in the variable `matrix`.  
# The gene expression matrix is assumed to be stored in a file named "matrix.csv" in the project folder. 
# Finally, the matrix is converted to a data frame using the `data.frame` function before being returned as the output. 
# The function does not modify the original matrix file; it only reads the file and returns the matrix as a data frame. 

ReadGEO  <- function(projname){
  # Reading the gene expression matrix file 
  matrix <- read_csv (paste0(projname , "/", "matrix.csv"))

  # Converting the matrix to a data frame 	
  matrix <- data.frame (matrix)
}

	
	
###################################################################################################################################### 
###################################################### 3- MakePhenotype Function ##################################################### 
# MakePhenotype: Creates a phenotype file based on project-specific phenotype data 

# Input (Arguments): 
#   - projname: A character vector of project names 
#   - compare: A string indicating whether to compare "subtype" or "grade" 

# Output: 
#   - A phenotype file containing the selected comparison (subtype or grade) 

# Function Description: 
# This function creates a phenotype file based on project-specific phenotype data. It takes the argument `projname`,  
# which is a character vector of project names. The function reads the phenotype file for each project, converts  
# the phenotype data to a data frame, and sets row names. Depending on the value of the `compare` argument, the function  
# either compares subtype or grade information. For subtype comparison, it extracts subtype data, formats it, and stores 
# it in a new list. For grade comparison, it extracts grade data, handles special cases (such as Greek symbols), and  
# formats it. The formatted phenotype data is stored in a new list for each project. The function then combines the  
# phenotype data from all projects into a single data frame. Finally, it writes the phenotype data to a CSV file and  
# displays a message indicating that the phenotype file is created. 

Makephenotype  <- function(projname , compare= "grade"){
  # Create a folder to store the phenotype file
  dir.create("phenotype", showWarnings = FALSE)

  # Initialize an empty list to store the phenotype data for each project	
  pheno <-list()

  # Iterate over each project name	
  for(i in seq_along(projname)){
    # Read the phenotype file for the current project
    phenotype <- read_tsv(paste0(projname [i], "/", "phenotype.tsv"))

    # Convert the phenotype data to a data frame and set row names
    phenotype <- data.frame (phenotype)
    rownames(phenotype) <- phenotype[,1]
    phenotype <- phenotype[,-1]

    # Add the phenotype data to the list using the project name as the key
    pheno[[projname [i]]] <- phenotype
  }
	
 # Choose either subtype or grade for comparison
 if(compare  == "subtype"){
    pheno2 <-list()

    # Iterate over each project in the phenotype list
    for (i in seq_along(pheno)){
      # Find columns related to subtype information (e.g., pam50)
      phecoln <- grep("pam50", pheno[[i]])

     # Check if any subtype columns were found
     if(length(phecoln) > 0){
	# Extract the subtype data and format it
        pheno1 <- data.frame (pheno[[i]][ ,phecoln])
        colnames(pheno1) <- "subtype"
        pheno1$subtype <- gsub("pam50.+:", "", pheno1$subtype)

	# Prepare the final phenotype data by adding sample identifiers and column names
        pheno1 <- cbind(rownames(pheno[[i]]), pheno1)
        colnames(pheno1)[1] <- "Sample"

	# Store the formatted phenotype data in a new list using the project name as the key
        pheno2[[names(pheno)[i]]] <- pheno1
      }
    }

    # Combine the phenotype data from all projects into a single data frame
    pheno3 <- do.call(rbind, pheno2)
    rownames(pheno3) <- NULL 

	
  } elseif(compare  == "grade"){
    pheno2 <-list()

    # Iterate over each project in the phenotype list	
    for (i in seq_along(pheno)){
      # Find columns related to grade information
      gradecoln <- grep("grade:", pheno[[i]])

     # Check if any grade columns were found
     if(length(gradecoln) > 0){
	# Extract the grade data and format it
        pheno1 <- data.frame (pheno[[i]][ ,gradecoln])
        colnames(pheno1) <- "grade"
        pheno1$grade <- gsub("grade:|grade: |B-R grade: ", "", pheno1$grade)
        pheno1$grade <- gsub("=.+", "", pheno1$grade)

	# Handle special cases where grades are represented as Greek symbols
        greek <- grep("III", pheno1$grade)
        if(length(greek) > 0){
          one <- which(factor(pheno1$grade) == "I")
          two <- which(factor(pheno1$grade) == "II")
          three <- which(factor(pheno1$grade) == "III")
          pheno1[one, ] <- "1"
          pheno1[two, ] <- "2"
          pheno1[three, ] <- "3"
        }

	# Prepare the final phenotype data by adding sample identifiers and column names
        pheno1 <- cbind(rownames(pheno[[i]]), pheno1)
        colnames(pheno1)[1] <- "Sample"

	# Store the formatted phenotype data in a new list using the project name as the key
        pheno2[[names(pheno)[i]]] <- pheno1
      }
    }

    # Combine the phenotype data from all projects into a single data frame
    pheno3 <- do.call(rbind, pheno2)
    rownames(pheno3) <- NULL 
  }
	
  # Write the phenotype data to a CSV file
  write.csv(pheno3, paste0("phenotype", "/", "phenotype_",compare , ".csv"), row.names = FALSE, quote = FALSE)

  # Display a message indicating that the phenotype file is created	
  message ("phenotype file is created into phenotype folder")
}


	
###################################################################################################################################### 
####################################################### 4- MakeBoxPlot Function ###################################################### 
# MakeBoxPlot: Creates a box plot to visualize gene expression values from a gene expression matrix 

# Input (Arguments): 
#   - GEOmatrix: A gene expression matrix in data frame format 
#   - projname: A character vector specifying the project name 

# Output: 
#   - A box plot visualizing the gene expression values 

# Function Description: 
# This function creates a box plot to visualize the gene expression values from a gene expression matrix.  
# It takes two arguments: `GEOmatrix`, which is the gene expression matrix in data frame format, and `projname`,  
# which is a character vector specifying the project name. The function converts the gene expression matrix into  
# a long format data frame using the `melt` function. It then uses `ggplot2` to create the box plot, with the  
# x-axis representing the samples and the y-axis representing the gene expression values. The plot is customized  
# with color and fill options, axis labels, y-axis scale labels, and theme settings. The resulting box plot  
# is displayed, with the subtitle indicating the project name. 
	
MakeBoxPlot <- function(GEOmatrix,projname){
  # Convert the gene expression matrix into a long format data frame
  melted <- melt(GEOmatrix, id.vars= "ID")

  # Define a color palette for the plot	
  my_pal <- c("#1B9E77")

  # Plot the box plot using ggplot2 
  ggplot(data = melted, aes(x = variable, y = value)) +
    # Add the box plot layer with color and fill options
    geom_boxplot(aes(color = "#45B39D", fill= "#45B39D"), outlier.alpha = 0.1)+
    # Set the x and y axis labels
    labs(x= 'SampleID', y= 'Values') +
    # Customize the y-axis scale labels
    scale_y_continuous(labels = scales::comma) +
    # Apply a classic theme to the plot	
    theme_classic() +
    # Set the color scale manually
    scale_color_manual(values=c(my_pal)) +
    # Set the fill scale manually with a modified transparency value
    scale_fill_manual(values=c(paste(my_pal, "66", sep = ""))) +
    # Customize the appearance of the axes, text, and title
    theme(
	# Customize axis labels and text styles:
        # Customize axis text font, size, color, and angle
	axis.text = element_text(family = "Arial",size = 24 , colour = "black", angle = 90),
	# Customize x-axis text font, color, and size
        axis.text.x = element_text(family = "Arial",colour = "black", size = 6),
	# Customize y-axis text font, color, and size
        axis.text.y = element_text(family = "Arial",colour = "black", size = 12),
	# Customize plot subtitle font, size, color, and horizontal justification
        plot.subtitle = element_text(family = "Arial",size = 24, colour = "black", hjust = 0.5),
	# Customize y-axis title font, size, and angle
        axis.title.y = element_text(family = "Arial", size = 24, angle = 90),
	# Customize x-axis title font, size, and angle
        axis.title.x = element_text(family = "Arial", size = 24, angle = 00),
	# Remove the legend from the plot
        legend.position="none"
    ) +

    # Add a subtitle to the plot indicating the project name
    labs(subtitle = paste0("Box plot - ", projname))
}


	
###################################################################################################################################### 
####################################################### 5- the MakePCA Function ###################################################### 
# MakePCA: Performs principal component analysis (PCA) and generates a plot 
 
# Input (Arguments): 
#   - GEOmatrix: A gene expression matrix in data frame format 
#   - studies: A list of gene expression matrices from different studies (optional) 
#   - projname: A character vector specifying the project name(s) 
#   - compare: A character specifying the comparison type (e.g., "grade") 

# Output: 
#   - A PCA plot or a 3D plot (based on the number of studies) 

# Function Description: 
# This function performs principal component analysis (PCA) on the gene expression matrix provided.  
# It takes the gene expression matrix, `GEOmatrix`, as the primary input. Additionally, it can take a  
# list of gene expression matrices from different studies, `studies`, the project name(s), `projname`,  
# and the comparison type, `compare` (e.g., "grade"). The function first prepares the data by scaling  
# the gene expression matrix and performing PCA using the `preProcess` function from the `caret` package.  
# It then plots the PCA results either as a 2D plot for multiple studies or a 3D plot for a single study.  
# The plot is customized with colors, legends, axis labels, and theme settings. The resulting PCA plot  
# or 3D plot is displayed or saved as an HTML widget based on the number of studies provided. 
	
MakePCA<- function(GEOmatrix, studies =  NULL , projname , compare  = "grade"){

  # Create a directory to store the plots	
  dir.create("Plots", showWarnings = FALSE)
	
  # Read the gene expression matrix file
  set.seed(123)

  # Check if the first column of GEOmatrix is "ID"	
  if(colnames(GEOmatrix)[1] == "ID"){
    # Set the row names of GEOmatrix to the first column and remove the ID column
    rownames(GEOmatrix) <- GEOmatrix[,1]
    GEOmatrix <- GEOmatrix[,-1]
  }
	
  # Read the phenotype file
  pheno <- read_csv (paste0("phenotype/phenotype_", compare , ".csv"))
  pheno <- data.frame (pheno)

  # Subset the phenotype data for samples present in the gene expression matrix
  paired_pheno <- which(pheno$Sample %in% colnames(GEOmatrix))
  pheno <- pheno[paired_pheno, ]
	
  # Perform PCA analysis:
  # Convert the gene expression matrix to a data frame and transpose it
  GEOmatrix <- data.frame (t(GEOmatrix))

  # Scale the data for PCA	
  pca_df <- scale(GEOmatrix)

  # Create a PCA model
  pca <- preProcess(x = pca_df, method = 'pca', pcaComp = 3)

  # Apply the PCA transformation to the data	
  pca_df <- data.frame (predict(pca, pca_df))

  # Sort the PCA data frame and the phenotype data frame by row names
  pca_df <- pca_df[order(rownames(pca_df)), ]
  pheno <- pheno[order(pheno$Sample), ]
 
  # Handling missing values in phenotype data:
  # Identify rows with NA values in the "grade" column	
  na_row <- which(is.na(pheno$grade))

  # Execute the following code if there are NA values in the row	
  if(length(na_row) > 0){
    # Get the sample names with NA values
    na_samp <- pheno[na_row,1]

    # Identify the corresponding columns in the PCA data frame
    na_colm <- which(rownames(pca_df) %in% na_samp)

    # Remove rows with NA values from the phenotype data frame and corresponding columns from the PCA data frame
    pheno <- pheno[-na_row, ]
    pca_df <- pca_df[-na_colm, ]
  }

  # Assign class labels to PCA data:
  # Create a factor variable for class labels based on the "grade" column 
  pca_df$Class <- factor(pheno$grade)

  # Define a color palette for the 2D plot for multi-studies 
  my_pal <- c("#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#9A7D0A")

  # Execute the following code if the length of projname is greater than 1
  if(length(projname)>1){

    # Get the row names of the PCA data frame
    pca_rownames <- rownames(pca_df)

    # Create a batch variable for different studies
    study_list <- studies
    batch <-list()

    # Iterate over each element in study_list
    for(i in seq_along(study_list)){
      # Determine the number of samples from each study present in the PCA data frame
      l <- length(which(colnames(study_list[[i]]) %in% pca_rownames))
      # Create a list with repeated project names for each sample from the corresponding study
      batch[[i]] <- rep (projname [i], l)
    }
    # Combine the list into a single vector
    batch <- unlist(batch)
    # Assign the batch information as a factor variable in the PCA data frame
    pca_df$studies <- factor(batch)

    # Generate the PCA plot using ggplot2
    ggplot(data = pca_df, aes(x = PC1, y = PC2, color = Class, fill = Class)) +
      # Add points with shape based on the "Studies" variable
      geom_point(size = 3, aes( shape= studies)) +
      # Customize color scale for class labels
      scale_color_manual(values=c(my_pal)) +
      # Customize fill colors for class labels
      scale_fill_manual(values=c(paste(my_pal, "66", sep = ""))) +
      # Apply a classic theme to the plot
      theme_classic() +
      # Customize the appearance of the axes, text, and title
      theme(
	# Customize axis text font, size, and color
	axis.text = element_text(family = "Arial",size = 24 , colour = "black"),
	# Customize x-axis text font, color, and size
        axis.text.x = element_text(family = "Arial",colour = "black", size = 24),
	# Customize y-axis text font, color, and size
        axis.text.y = element_text(family = "Arial",colour = "black", size = 24),
	# Customize plot subtitle text font, size, color, and position
        plot.subtitle = element_text(family = "Arial",size = 24, colour = "black", hjust = 0.5),
	# Customize y-axis title font, size, and rotation angle
        axis.title.y = element_text(family = "Arial", size = 24, angle = 90),
	# Customize x-axis title font, size, and rotation angle
        axis.title.x = element_text(family = "Arial", size = 24, angle = 00),
	# Customize legend text size and font
        legend.text = element_text(size = 10, family = "Arial"), 
	# Customize legend title size and font
        legend.title = element_text(size = 20, family = "Arial")
      ) +
      # Set the subtitle of the plot
      labs(subtitle = "PCA plot- Merged studies")
    
  } else {
    
    # Plot a 3D plot for a single study:
    # Plotting using the plot3d function
    plot3d( 
	# Specify the coordinates for the 3D plot
        x=pca_df$PC1, y=pca_df$PC2, z=pca_df$PC3, 
        # Set color based on the class labels
        col = my_pal[pca_df$Class],
	# Set the point type as "s" for spheres 
        type = 's', 
	# Set the radius of the spheres
        radius = 4,
	# Set the size of the spheres
        size = 3,
	# Set the axis labels
        xlab="PC1", ylab="PC2", zlab="PC3"
    )
    
    # Add a legend to the plot
    legend3d(
	# Set the position of the legend
        "topright", 
	# Set the legend labels based on the class labels
        legend = levels(factor(pca_df$Class)), 
	# Set the legend title
        title = compare ,
	# Set the legend marker type
        pch = 16,
	# Set the legend marker colors
        col = my_pal, 
	# Set the legend text size
        cex=1.2, 
	# Set the inset position of the legend
        inset=c(0.02)
    )
    
    # Save the 3D plot as an HTML widget:
    # Combine project names into a single string
    projname  <- paste0(projname , collapse = "-") 
    # Save the 3D plot as an HTML widget using the saveWidget() function
    htmlwidgets::saveWidget(
        # Specify the dimensions of the widget
        rglwidget(width = 800, height = 800), 
        # Set the file name
        file = paste0("Plots/PCA_", projname , ".html"),
	# Specify the library directory
        libdir = "libs",
	# Set the title of the HTML file
	title = paste0("PCA - ", projname),
	# Specify whether the HTML file should be self-contained
        selfcontained = FALSE
    )
  }
}



###################################################################################################################################### 
##################################################### 6- VSNQuantilNorm Function######################################################
# VSNQuantilNorm: Performs VSN (Variance Stabilizing Normalization) followed by quantile normalization on the input gene expression matrix 

# Input (Arguments): 
#   - GEOmatrix: The gene expression matrix to be normalized 

# Output: 
#   - A data frame containing the normalized gene expression matrix 
 
# Function Description: 
# This function performs Variance Stabilizing Normalization (VSN) followed by quantile normalization on the  
# input gene expression matrix, `GEOmatrix`. VSN is used to stabilize the variance of gene expression values 
# across samples, and quantile normalization is applied to ensure that the distribution of gene expression values  
# is the same across samples. The function first sets the row names of the gene expression matrix and removes the  
# ID column. Then, it performs VSN normalization using the `normalizeVSN` function. After VSN normalization, the  
# function applies quantile normalization by ranking and sorting each column of the normalized matrix. It calculates  
# the mean across rows and applies an indexing function to retrieve the mean values using column indices. Finally,  
# the function returns a data frame containing the normalized gene expression matrix, where the first column is  
# named "ID" and row names are combined with the normalized matrix. 

VSNQuantilNorm <- function(GEOmatrix){
 # Set the seed for reproducibility of random processes 	
  set.seed(1234)
	
  # Reading the gene expression matrix and setting row names	
  rownames(GEOmatrix) <- GEOmatrix[,1]
  GEOmatrix <- GEOmatrix[,-1]
	
  # Normalize the gene expression data using the VSN (variance stabilization and normalization) method
  vsnnorm <- normalizeVSN(GEOmatrix)
	
  # Quantile normalization:
  # Rank each column of the normalized matrix
  df_rank <- apply(vsnnorm,2,rank,ties.method="min")
  # Sort each column of the normalized matrix
  df_sorted <- data.frame (apply(vsnnorm, 2, sort))
  # Calculate the mean across rows	
  df_mean <- apply(df_sorted, 1, mean)

  # Function to retrieve mean values using column indices
  index_to_mean <- function(my_index, my_mean){

    # Return the mean value at the specified index	
    return(my_mean[my_index])
  }

  # Apply the index_to_mean function to the ranked matrix to obtain the final normalized matrix	
  norm_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  # Convert the final normalized matrix into a data frame
  norm_final <- data.frame (norm_final)
  # Convert all columns of the data frame to numeric type	
  norm_finall <- data.frame (sapply(norm_final, function(x) as.numeric(as.character(x))), check.names= FALSE)
  # Combine the row names of the original matrix with the normalized matrix
  norm_finall <- cbind(rownames(GEOmatrix), norm_finall)
  # Rename the first column as "ID"	
  colnames(norm_finall)[1] <- "ID"

  # Return the normalized gene expression matrix
  return(norm_finall)
}


	
###################################################################################################################################### 
###################################################### 7- MergeStudies function ######################################################   
# MergeStudies: Merge multiple matrices 

# Input (Arguments): 
#   - studies: A list of matrices to be merged 

# Output: 
#   - The merged matrix

# Function Description: 
# This function merges multiple matrices provided as a list, 'studies'. It uses the 'merge' function to combine  
# the matrices, ensuring that all rows and columns from each matrix are included in the merged result by setting  
# the 'all=TRUE' argument. The function returns the merged matrix, which contains the combined data from all  
# input matrices. 
	
Mergestudies <- function(studies){
	
  # Merging two or more matrices using the 'merge' function with 'all=TRUE' to include all rows and columns
  study_list <- Reduce(function(x, y) merge(x, y, all= TRUE), studies)
  return(study_list)
}


	
###################################################################################################################################### 
#################################################### 8- StudyBatchEffect Function ####################################################  
# StudyBatchEffect: Perform batch effect correction on merged matrices 

# Input (Arguments): 
#   - studies: A list of matrices with batch effects

# Output: 
#   - The batch-corrected matrix 

# Function Description: 
# This function performs batch effect correction on a list of matrices provided in 'studies'.  
# It calculates the number of samples for each matrix and creates a vector of batch sizes based on the  
# number of samples in each matrix. The function then merges the matrices using the 'merge' function with  
# 'all=TRUE' to include all rows and columns. Next, it sets the row names of the merged matrix and performs  
# batch effect correction using the ComBat function. The corrected matrix is converted to a data frame, 
# ensuring columns are of numeric type. Finally, the row names are combined with the corrected matrix,  
# the first column is renamed as "ID", and the batch-corrected matrix is returned as the output of the function. 
	
StudyBatchEffect <- function(studies){
  # Calculating the number of samples for each matrix
  study_list <- studies
  # Create an empty list to store batches	
  batch <-list()

  # Iterate over the indices of the study_list
  for(i in seq_along(study_list)){
    # Subtract 1 to exclude the column containing row names
    batch[[i]] <- ncol(study_list[[i]]) - 1
  }

  # Flatten the list to create a vector of batch sizes
  batch <- unlist(batch)
  # Replicate the sequence of study indices based on batch sizes	
  batch <- rep (seq_along(study_list), batch)
	
  # Merging matrices using the 'merge' function with 'all=TRUE' to include all rows and columns
  study_list <- Reduce(function(x, y) merge(x, y, all= TRUE), study_list)

  # Setting the row names of the merged matrix to the first column
  rownames(study_list) <- study_list[,1]
  # Exclude the first column containing row names	
  study_list <- study_list[,-1]
	
  # Perform batch effect correction using the ComBat function
  corrected_batch <- ComBat(
	# Input data for batch effect correction
	dat=study_list,
	# Batch information for correction
        batch=batch,
	# No covariates specified for batch effect correction
        mod=NULL ,
	# Estimate empirical Bayes parameters
        par.prior= TRUE, 
	# Disable diagnostic plots
        prior.plots= FALSE)
	
  # Convert the corrected matrix to a data frame	
  corrected_batch <- data.frame (corrected_batch)

  # Convert columns of the data frame to numeric type
  corrected_batch <- data.frame (
	# Convert each column to numeric type using as.numeric(as.character(x))
	sapply(corrected_batch, function(x) as.numeric(as.character(x))),
	# Set check.names to FALSE to preserve column names
        check.names= FALSE, 
	# Set row.names to the row names of the original corrected_batch object
	row.names = rownames(corrected_batch))

  # Combine the row names with the corrected matrix	
  corrected_batch <- cbind(rownames(corrected_batch), corrected_batch)

  # Rename the first column as "ID"
  colnames(corrected_batch)[1] <- "ID"

  # Return the batch-corrected matrix as the output of the function	
  return(corrected_batch)
}


	
###################################################################################################################################### 
#################################################### 9- MakeDensityPlot Function ####################################################  
# MakeDensityPlot: Create a density plot of gene expression values 

# Input: 
#   - matrix: The gene expression matrix 
#   - studies: A list of matrices representing different studies 

# Output: 
#   - None (the plot is displayed) 

# Function Description: 
# This function creates a density plot to visualize the distribution of gene expression values across different studies.  
# It converts the gene expression matrix to a long format data frame, calculates the number of samples for each study,  
# and adds a study factor variable to the data frame. The function then plots the density plot using ggplot2,  
# with gene expression values on the x-axis and color-coded by studies. It customizes the appearance of the plot  
# and adds axis labels and a subtitle. The density plot provides insights into the distribution patterns and  
# potential differences in gene expression between studies. 

MakeDensityPlot <- function(matrix, studies){
  # Set the random seed for reproducibility
  set.seed(1234)
	
  # Convert the gene expression matrix to a long format data frame
  melted <- melt(matrix, id.vars = "ID")
	
  # Calculating the number of samples for each study:
  # Create an empty list named 'st'
  st <-list()


  # Loop over the sequence of indices of the 'studies' vector
  # The variable 'i' represents the current index in each iteration
  # The loop iterates over each element in 'studies'	
  for(i in seq_along(studies)){
    # Get the number of rows in the study matrix
    row_n <- nrow(studies[[i]])
    # Get the number of columns in the study matrix excluding the ID column
    col_n <- ncol(studies[[i]]) -1
    # Create a vector of study names repeated for each sample
    st[[i]] <- rep (names(studies)[i], row_n*col_n)
  }
	
  
  st <- unlist(st)
  # Add the study factor variable to the melted data frame	
  melted$studies <- factor(st)
	
  # Plotting a density plot:
  # Define a color palette for the plot
  my_pal <- c("#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#9A7D0A")
  # Create a ggplot object with the melted data as the input data
  # Specify x as gene expression values and color by studies
  ggplot(data = melted, aes(x = value, color = studies)) +
    # Add density plot with a specified line size
    geom_density(size = 1.5) +
    # Set axis labels	
    labs(y= 'Density', x= 'Value') +
    # Apply a classic theme to the plot
    theme_classic() +
    # Set manual color scale based on study names
    scale_color_manual(values=c(my_pal)) +
    # Set manual fill color scale
    scale_fill_manual(values=c(paste(my_pal, "66", sep = ""))) +
    # Customize the appearance of the axes, text, and title
    theme(
	# Customize axis text styles
	axis.text = element_text(family = "Arial",size = 24 , colour = "black", angle = 90),
	# Customize x-axis text style
        axis.text.x = element_text(family = "Arial",colour = "black", size = 12),
	# Customize y-axis text style
        axis.text.y = element_text(family = "Arial",colour = "black", size = 24),
	# Customize subtitle
        plot.subtitle = element_text(family = "Arial",size = 24, colour = "black", hjust = 0.5),
	# Customize y-axis title
        axis.title.y = element_text(family = "Arial", size = 24, angle = 90),
	# Customize x-axis title
        axis.title.x = element_text(family = "Arial", size = 24, angle = 00)
    ) +
    # Add a subtitle to the plot
    labs(subtitle = "Density plot")
} 


	
###################################################################################################################################### 
###################################################### 10- MakeMDSPlot Function ######################################################  
# MakeMDSPlot: Create a multidimensional scaling (MDS) plot 

# Input (Arguments): 
#   - matrix: The input matrix containing gene expression data 
#   - studies: A list of matrices representing different studies 
#   - projname: A vector of project names corresponding to each study 
#   - compare: A comparison variable 
#   - nstudies: The number of studies 

# Output: 
#   - None (the plot is displayed) 

# Function Description: 
# This function creates an MDS plot to visualize the relationship between samples based on their gene expression data.  
# The plot helps to identify clusters and patterns in the data, providing insights into similarities and differences among samples. 

# Define the MakeMDSPlot function with input parameters
MakeMDSPlot <- function(matrix, studies, projname, compare,  nstudies){
  # Set the random seed for reproducibility
  set.seed(1234)

  # Check if the matrix has an "ID" column, and if so, set row names and remove the ID column	
  if(colnames(matrix)[1] == "ID"){
    # Set row names to the values in the "ID" column
    rownames(matrix) <- matrix[,1]
    # Remove the "ID" column from the matrix
    matrix <- matrix[,-1]
  }

  # Transpose the matrix and convert it to a data frame	
  matrix <- data.frame (t(matrix))


	
  # Compute MDS using the dist() and cmdscale() functions:
  # Assign the input matrix to the variable mds
  mds <- matrix %>%
    # Compute the distance matrix using the dist() function
    dist() %>%          
    # Perform classical multidimensional scaling (MDS) on the distance matrix using the cmdscale() function.
    cmdscale() %>%
    # Convert the MDS output to a data frame and assign it to the variable mds.
    data.frame ()
  # Set column names as "Dim1" and "Dim2"
  colnames(mds) <- c("Dim1", "Dim2")

	
  # Perform clustering using k-means on the MDS coordinates
  # Perform k-means clustering on mds with nstudies clusters and assign the resulting cluster assignments to the variable clust.
  clust <- kmeans(mds, nstudies)$cluster %>%
    # Convert the cluster assignments to a factor data type.
    as.factor()

  # Add clustering results to the MDS data frame
  # Update the 'mds' data frame
  mds <- mds %>%
    # Add a new column called "Clusters" and assign cluster assignments from 'clust'
    mutate(clusters = clust)
	
  # Retrieve the row names of the MDS data frame
  mds_rownames <- rownames(mds)

  # Assign study names to each sample based on column names of study matrices:
  	
  study_list <- studies
  # Create an empty list called 'batch'	
  batch <-list()
  # Loop over the indices of 'study_list'	
  for(i in seq_along(study_list)){
    # Determine the number of samples in the current study
    l <- length(which(colnames(study_list[[i]]) %in% mds_rownames))
    # Create a vector of project names for the samples in the current study
    batch[[i]] <- rep (projname [i], l)
  }
	
  # Combine the project names into a single vector
  batch <- unlist(batch)
  # Add the project names as a factor to the MDS data frame	
  mds$studies <- factor(batch) 

  # MDS plot customization:
  # This section creates an MDS (Multidimensional Scaling) plot using the ggscatter function from the ggpubr package and customizes various plot elements:
  # Define a vector of color codes for plotting
  my_pal <- c("#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#9A7D0A")
  # Create the scatter plot using ggscatter from the ggpubr package
  g <- ggscatter(mds, x = "Dim1", y = "Dim2", 
                 # Specify the label, shape, color, palette, size, and other options for ggscatter:
                 # Set the labels for the points as row names of the MDS matrix
                 label = rownames(mds),
	         # Use the "Studies" variable for shaping the points
	         shape = "Studies", 
                 # Set the theme of the ggplot to "classic"
                 ggtheme = theme_classic(),
	         # Use the "Clusters" variable for coloring the points
                 color = "Clusters",
	         # Set the color palette for the points
                 palette = my_pal,
	         # Set the size of the points
                 size = 1.5, 
	         # Display ellipses around the points
                 ellipse = TRUE,
	         # Set the type of ellipse to "convex"
                 ellipse.type = "convex",
	         # Enable point label repulsion to avoid overlapping labels
                 repel = TRUE)
  # Customize the appearance of the axes, text, and title
  # Add the theme and labels to the ggplot object
  g + theme(
    # Customize axis text
    axis.text = element_text(family = "Arial",size = 24 , colour = "black"),
    # Customize x-axis text
    axis.text.x = element_text(family = "Arial",colour = "black", size = 24),
    # Customize y-axis text
    axis.text.y = element_text(family = "Arial",colour = "black", size = 24),
    # Customize subtitle
    plot.subtitle = element_text(family = "Arial",size = 24, colour = "black", hjust = 0.5),
    # Customize y-axis title
    axis.title.y = element_text(family = "Arial", size = 24, angle = 90),
    # Customize x-axis title
    axis.title.x = element_text(family = "Arial", size = 24, angle = 00)
    ) +
    # Set the subtitle of the plot as "MDS plot"
    labs(subtitle = "MDS plot")
}



###################################################################################################################################### 
##################################################### 11- CochranQTest Function ######################################################
# CochranQTest: Perform Cochran's Q Test on multiple studies 

# Input: 
#   - studies: A list of matrices representing different studies
#   - groups: A vector specifying the comparison groups (e.g., grade1 vs grade3) 
#   - compare: A string indicating the type of comparison (e.g., "grade") 

# Output: 
#   - A list containing the mean and variance for each comparison mode 

# Function Description: 
# Cochran's Q Test: Perform statistical analysis to compare proportions across multiple groups 

# This function performs Cochran's Q Test on multiple studies to compare proportions across different groups. 
# It takes a list of matrices representing different studies, a vector specifying the comparison groups,  
# and a string indicating the type of comparison. The function calculates the mean and variance for each comparison mode  
# and returns them as a list. Cochran's Q Test is a statistical test used to determine if there are significant differences  
# in proportions between multiple groups. It helps to assess the homogeneity of proportions across the groups and  
# identify any significant departures from the null hypothesis of equal proportions. 
	
CochranQTest <- function(studies, groups, compare ){
  # Set the random seed for reproducibility
  set.seed(1234)
	
  # Reading matrix and phenotype:
  # Initialize a list to store the mean values for each study
  d.adj.Split <-list()
  # List to store variance values for each study
  var.d.adj.Split <-list()

  # Iterate over each study in the 'studies' list
  for(i in seq_along(studies)){
    # Get the matrix for the current study
    matrix <- studies[[i]]

    # Check if the matrix has an "ID" column, and if so, handle it accordingly
    # (e.g., setting row names and removing the ID column)
    if(colnames(matrix)[1] == "ID"){
      # Set row names to the values in the first column
      rownames(matrix) <- matrix[,1]
      # Remove the first column from the matrix
      matrix <- matrix[,-1]
    }
	
    # Read the phenotype file for the corresponding comparison type:
    # Read the phenotype file
    pheno <- read_csv (paste0("phenotype/phenotype_", compare , ".csv"))
    # Convert the read data into a data frame
    pheno <- data.frame (pheno)

    # Select the samples from the phenotype file that match the column names in the matrix:
    # Identify the matching column indices
    paired_pheno <- which(pheno$Sample %in% colnames(matrix))
    # Select the rows that correspond to the matched samples
    pheno <- pheno[paired_pheno, ]

    # Identify the rows with missing values in the "grade" column
    na_row <- which(is.na(pheno$grade))


	
    # Check if there are any rows with missing values in the "grade" column
    # and handle the missing values accordingly (e.g., removing corresponding rows and columns)
    # Check if the length of 'na_row' is greater than zero (if there are any missing values in the 'na_row' variable)
    if(length(na_row)> 0){
      # Get the sample IDs with missing values
      na_samp <- pheno[na_row,1]
      # Identify the matching column indices in the matrix
      na_colm <- which(colnames(matrix) %in% na_samp)
      # Remove the rows with missing values from the phenotype data
      pheno <- pheno[-na_row, ]
      # Remove the corresponding columns from the matrix
      matrix <- matrix[,  -na_colm]
    }

    # Calculating the number of each comparison mode (grade1 vs grade3):
    # Empty list to store the factor levels (0 and 1) for each comparison group
    g_list <-list()
    # Empty list to store the sample IDs for each comparison group
    sampleid <-list()

    # Iterate over the comparison groups and perform operations for each group
    for(j in seq_along(groups)){
      # Identify the rows with the specified group 
      g_row <- which(pheno $grade %in% groups[j])
      # Get the sample IDs for the specified group 
      s_id <- pheno[g_row,1]
      # Store the sample IDs for the comparison group 
      sampleid[[j]] <- s_id
      # Create the factor levels for the comparison group 
      g_list[[j]] <- factor(rep(j - 1, length(g_row)))
    }
	
    # Applying those modes and calculating the mean and variance for Cochran' Q Test 
    # Unlist the factor levels for all comparison groups to create a single vector of factor levels 
    spl <- unlist(g_list)
    # Unlist the sample IDs for all comparison groups 
    spsam <- unlist(sampleid)
    # Select the columns in the matrix corresponding to the sample IDs 
    matrix <- matrix[,spsam]
    # Convert the matrix to a matrix object
    matrix <- as.matrix(matrix)
    # Calculate the mean difference for each comparison mode 
    d.Split <- getdF(matrix, spl)
    # Calculate the overall mean difference
    mean_d <- dstar(d.Split, length(spl))
    # Store the mean difference for the current study 
    d.adj.Split[[i]] <- mean_d
    # Calculate the variance 
    var_d <- sigmad(mean_d, sum(spl==0), sum(spl==1))
    # Store the variance for the current study 
    var.d.adj.Split[[i]] <- var_d
  }
	
  # Create a list containing the mean and variance values 
  output <-list(d.adj.Split, var.d.adj.Split)

  # Return the output as a list containing mean differences and variances 	
  return(output)

  # Display a message indicating the completion of the analysis 	
  message ("Analysis has completed.")
}



###################################################################################################################################### 
###################################################### 12- MakeQQPlot Function #######################################################  
# MakeQQPlot: Generate a QQ plot for Cochran's Q Test 

# Input: 
#   - cohranq: Output from CochranQTest function (list of means and variances) 
#   - nstudies: Number of studies 

# Output: 
#   - QQ plot 

# Function Description: 
# This function generates a QQ plot for Cochran's Q Test. It takes the output from the CochranQTest function,  
# which is a list of means and variances, and the number of studies as input. The function calculates the quantiles  
# of the chi-square distribution and the quantiles of the sample, and then plots them on a QQ plot.  
# A QQ plot (quantile-quantile plot) is used to assess whether a dataset follows a specific theoretical distribution  
# by comparing the observed quantiles against the expected quantiles. In this case, the QQ plot helps to visualize  
# the goodness-of-fit between the observed Cochran's Q statistics and the expected chi-square distribution.  
# Deviations from the diagonal reference line indicate departures from the expected distribution. 
	
MakeQQPlot <- function(cohranq, nstudies){
  # Unpack the previous list into data frames for means and variances 
  # Extract the means from the output list 
  output1 <- cohranq[[1]] 
  # Extract the variances from the output list
  output2 <- cohranq[[2]]
  # Combine the means into a data frame	
  mymns <- data.frame (do.call(cbind, output1))
  # Combine the variances into a data frame
  myvars <- data.frame (do.call(cbind, output2))

  # Calculate Cochran's Q Test:
  # Compute Cochran's Q Test using means and variances
  my.Q <- f.Q(mymns, myvars)
  # Store the number of studies for later use
  nstudies <- nstudies

	
  # Calculate quantiles of the chi-square distribution
  chisqq <- qchisq(seq(0, .9999, .001), df = nstudies-1)
  # Calculate quantiles of the sample
  tmp <- quantile(my.Q, seq(0, .9999, .001))
  # Create a data frame for QQ plot data
  qq <- data.frame (chisqq = chisqq, tmp = tmp)

  # Plotting the QQ plot:
  # Set data and aesthetics for the plot
  ggplot(data = qq, aes(x = chisqq, y = tmp, color = "#1B9E77")) +
    # Add points to the plot
    geom_point() +
    # Add a diagonal reference line to the plot
    geom_abline() +
    # Apply a classic theme to the plot
    theme_classic() +
    # Set labels for x and y axes
    labs(x= 'quantiles of Chi square', y= 'quantiles of Sample') +
    # Customize the appearance of the axes, text, and title
    theme(
      # Customize axis text appearance
      axis.text = element_text(family = "Arial",size = 24 , colour = "black"),
      # Customize x-axis text appearance
      axis.text.x = element_text(family = "Arial",colour = "black", size = 24),
      # Customize y-axis text appearance
      axis.text.y = element_text(family = "Arial",colour = "black", size = 24),
      # Customize plot subtitle appearance
      plot.subtitle = element_text(family = "Arial",size = 24, colour = "black", hjust = 0.5),
      # Customize y-axis title appearance
      axis.title.y = element_text(family = "Arial", size = 24, angle = 90),
      # Customize x-axis title appearance
      axis.title.x = element_text(family = "Arial", size = 24, angle = 00),
      # Remove the legend
      legend.position="none"
    ) +
    # Set the plot subtitle
    labs(subtitle = "QQ Plot")
} 



###################################################################################################################################### 
##################################################### 13- FEMREMAnalysis Function##################################################### 
# FEMREMAnalysis: Perform FEM-REM meta-analysis 

# Input: 
#   - studies: A list of expression matrices from different studies 
#   - projname: Names of the projects or studies 
#   - compare: Comparison mode (e.g., "grade1_vs_grade3")
#   - groups: List of comparison groups 
#   - model: Model selection ("FEM" for fixed effects model, "REM" for random effects model) 
#   - FDR: False discovery rate threshold for filtering genes 
#   - nperm: Number of permutations for z-score calculation 

# Output: 
#   - Meta-analysis results 

# Function Description: 
# This function performs FEM-REM meta-analysis on expression matrices from different studies.  
# It takes as input a list of expression matrices, names of the projects or studies, a comparison mode,  
# a list of comparison groups, a model selection (FEM for fixed effects model or REM for random effects model),  
# a false discovery rate (FDR) threshold for gene filtering, and the number of permutations for z-score calculation. 

# The function first creates folders to store the meta-analysis results. It reads the phenotype file for the comparison  
# mode and the annotation files for each project. It then performs data preprocessing, including handling missing values 
# and filtering based on comparison groups. The expression matrices, class labels, and annotation data are prepared for analysis.   

# Next, the function determines the type of analysis based on the model selection and performs z-score calculation using 
# the zScoreFDR function. The two-sided scores are extracted and combined with annotation information.  
# Missing gene symbols are removed, and the data frame is subsetted, counted, and checked for duplicated gene symbols.   
# For duplicated gene symbols, average expression values and associated information are calculated and stored.  
# The averaged genes are combined with unique genes into a single data frame, sorted by z-score in descending order, and written to files. 

# The function also performs FDR filtering on the data frame, writes the filtered results to files, 
# and displays summary information for each project. The completion message is displayed when the analysis is finished. 
	
FEMREMAnalysis <- function(studies, projname , compare , groups, model, FDR, nperm){
  # Set the random seed for reproducibility
  set.seed(1234)
	
  # Creating folders for meta-analysis results:
  # Create MetaAnalysis folder
  dir.create("MetaAnalysis", showWarnings = FALSE)
  # Create filtered subfolder	
  dir.create("MetaAnalysis/filtered", showWarnings = FALSE)
  # Create volcano subfolder
  dir.create("MetaAnalysis/volcano", showWarnings = FALSE)
	
  # Reading phenotype file
  phenotype <- read_csv (paste0("phenotype/phenotype_", compare , ".csv"))
  phenotype <- data.frame (phenotype)

  # Reading annotation file for each project:
  # Initialize an empty list to store annotations for each project
  annotations <-list()

  # Iterate over each project to process annotations
  for(i in seq_along(projname)){
    # Read annotation file
    annot<- read_tsv(paste0(projname [i], "/", "annot.tsv"))
    # Convert to data frame
    annot<- data.frame (annot)
    # Set row names
    rownames(annot) <- annot[,1]
    # Remove the first column
    annot<- annot[,-1]
    # Select specific columns
    annot<- annot[ , c("Gene.symbol", "Chromosome.location", "GO.Function")]
    # Combine row names and data
    annot<- cbind(rownames(annot), annot)
    # Rename the first column
    colnames(annot)[1] <- "ID" 
    # Order by ID column
    annot<- annot[order(annot$ID), ]
    # Store annotation data in a list
    annotations[[projname [i]]] <- annot
  }

  # Perform intersection of annotations	
  an <- Reduce(intersect, annotations)
  # Combine annotations into a data frame
  annotations <- data.frame (do.call(cbind, an))

  # Create a list to store expression matrices
  matrices_list <-list()
  # Create a list to store class labels	
  Class <-list()

  # Iterate over each study to process data	
  for(i in seq_along(studies)){
    # Get the expression matrix for the current study
    matrix <- studies[[i]]
    # Check if the matrix has an "ID" column
    if(colnames(matrix)[1] == "ID"){
      # Set row names from the first column
      rownames(matrix) <- matrix[,1]
      # Remove the first column
      matrix <- matrix[,-1]
    }
	
    # Data preprocessing and cleaning for phenotype and gene expression matrices:
    # Identify the indices of samples in the phenotype matching matrix that are present in the column names of the matrix
    paired_pheno <- which(phenotype$Sample %in% colnames(matrix))
    # Subset the phenotype data frame to include only the rows corresponding to the identified samples
    pheno <- phenotype[paired_pheno, ]
    # Identify the indices of rows in the phenotype data frame where the grade column contains missing values (NA)
    na_row <- which(is.na(pheno$grade))
	
   # Check if there are any missing values in the rows of the phenotype data
   if(length(na_row)> 0){
      # Get the sample IDs with missing grade values
      na_samp <- pheno[na_row,1]
      # Get the corresponding column indices
      na_colm <- which(colnames(matrix) %in% na_samp)
      # Remove rows with missing grade values from phenotype
      pheno <- pheno[-na_row, ]
      # Remove columns with missing grade values from matrix
      matrix <- matrix[,  -na_colm]
    }
    
    # Calculating the number of each comparison mode (grade1 vs grade3):
    # Create a list to store comparison groups
    g_list <-list()
    # Create a list to store sample IDs
    sampleid <-list()

    # Iterate over each group to process comparisons
    for(j in seq_along(groups)){
      # Identify rows matching the current comparison group
      g_row <- which(pheno $grade %in% groups[j])
      # Get the sample IDs for the current group
      s_id <- pheno[g_row,1]
      # Store sample IDs in the list
      sampleid[[j]] <- s_id
      # Create factor levels for the current group
      g_list[[j]] <- factor(rep(j - 1, length(g_row)))
    }
	
    # Prepare factor levels, sample IDs, and expression matrix for analysis:
    # Unlist the factor levels for all comparison groups
    spl <- unlist(g_list)
    # Convert factor levels to numeric
    spl <- as.numeric(as.character(spl))
    # Store class labels for the current study
    Class[[i]] <- spl
    # Unlist the sample IDs for all comparison groups
    spsam <- unlist(sampleid)
    # Subset the matrix using the selected sample IDs
    matrix <- matrix[,spsam]
    # Convert to matrix class
    matrix <- as.matrix(matrix)
    # Create an ExpressionSet object
    matrix <- ExpressionSet(assayData=matrix)
    # Store the expression matrix in the list
    matrices_list[[i]] <- matrix
  }
 
  # Determine the type of analysis based on the model type (FEM or REM):
  # Check if the model type is FEM for Fixed Effects Model analysis:
  # - If the selected model is 'FEM' (Fixed Effects Model), set the useREM flag to FALSE.
  # - This flag is used to indicate the use of the Fixed Effects Model in subsequent analysis.	
  if(model == 'FEM'){
    # Set useREM flag to FALSE for FEM model (indicating the use of Fixed Effects Model)
    useREM <- FALSE
  } 

  # - If the selected model is 'REM' (Random Effects Model), set the useREM flag to TRUE.
  # - This flag is used to indicate the use of the Random Effects Model in subsequent analysis.
  if(model == 'REM'){
    # Set useREM flag to TRUE for REM model (indicating the use of Random Effects Model)
    useREM <- TRUE
  }

  # Perform z-score calculation	
  ScoresFDR <- zScoreFDR(matrices_list, Class, useREM = useREM, nperm = nperm)
  # Get the two-sided scores
  twosided <- data.frame (ScoresFDR[["two.sided"]])
  # Sort by row names
  twosided <- twosided[order(rownames(twosided)),]
	
  # Identify matching annotations and add annotation information:
  # Identify matching annotations
  anno_row <- which(annotations$ID %in% rownames(twosided))
  # Add annotation information
  twosided <- cbind(twosided,  annotations[anno_row, c(2, 3, 4)])

  # Removing missing gene symbols:
  # Identify missing gene symbols
  genemissing <- which(is.na(twosided$Gene.symbol) == TRUE)
  # Remove rows with missing gene symbols
  twosided <- twosided[-genemissing, ]

  # Subset, Count, and Identify Duplicated Gene Symbols:
  # Subset columns
  combination <- twosided[ ,c("zSco", "FDR", "Gene.symbol", "Chromosome.location", "GO.Function")]
  # Count duplicated gene symbols
  dupgene_freq <- data.frame (table(combination$Gene.symbol))
  # Get duplicated gene symbols
  only_dupgene <- dupgene_freq[which(dupgene_freq$Freq> 1), 1]
  # Identify rows with duplicated gene symbols
  only_dupgene_row <- which(combination$Gene.symbol %in% only_dupgene)
  # Subset rows with duplicated gene symbols
  duplicat_deg <- combination[only_dupgene_row,]

  # Calculate and store the average expression values and associated information for each gene symbol with duplicated entries 
  # Create an empty list to store averaged genes
  averaged_genes <-list()

  # Iterate over each duplicated gene symbol
  for(b in 1:length(only_dupgene)){
    # Identify rows for the current duplicated gene
    d_row <- which(duplicat_deg$Gene.symbol %in% only_dupgene[b])
    # Calculate column means
    col_mean <- rbind(colMeans(duplicat_deg[d_row, 1:2]))
    # Combine column means with other columns
    col_all <- cbind(col_mean, duplicat_deg[d_row[1],3:5])
    # Store the averaged gene information
    averaged_genes[[b]] <- col_all
  }
 
  # Combine averaged genes into a single data frame
  averaged_genes <- do.call(rbind, averaged_genes)

  # Combine the rows representing unique gene symbols with the rows representing averaged genes:
  # Subset rows with unique gene symbols
  unique_deg <- combination[-only_dupgene_row, ]
  # Combine unique genes and averaged genes
  all_unique_genes_combi <- rbind(unique_deg, averaged_genes)
  # Reset row names
  rownames(all_unique_genes_combi) <- NULL 

  # Sort the data frame by z-score in descending order
  all_unique_genes_combi <- all_unique_genes_combi[order(all_unique_genes_combi[,1], decreasing = TRUE),]
  # Write the result to a file
  write_tsv(all_unique_genes_combi, "MetaAnalysis/volcano/combination.tsv")
	
  # Filter the data frame based on FDR threshold
  all_unique_genes_combi <- all_unique_genes_combi[which(all_unique_genes_combi[, "FDR"] < FDR),]
  # Write the filtered result to a file
  write_tsv(all_unique_genes_combi, "MetaAnalysis/filtered/combination.tsv")
  
  # Seperating each study, averaging duplicated genes
  # Get column indices for FDR values
  sfdr <- grep("FDR_Ex_{0-9}", colnames(twosided))

  # Iterate over the indices of sfdr
  for(i in seq_along(sfdr)){
    # Get column indices for gene symbols
    gs_col <- grep("Gene.symbol|Chromosome.location|GO.Function", colnames(twosided))
    # Subset columns for the current study
    common <- twosided[,c(sfdr[i]-1, sfdr[i], gs_col)]

    # Rename columns with project names
    colnames(common)[1:2] <- gsub("Ex_.", projname [i],colnames(common)[1:2])

    # Perform operations on duplicated gene symbols if there are multiple rows in the common data frame:
    # Check if there are multiple rows in the common data frame
    if(nrow(common)>1){ 
      # Count duplicated gene symbols
      dupgene_freq <- data.frame (table(common$Gene.symbol))
      # Get duplicated gene symbols
      only_dupgene <- dupgene_freq[which(dupgene_freq$Freq> 1), 1]
      # Identify rows with duplicated gene symbols
      only_dupgene_row <- which(common$Gene.symbol %in% only_dupgene)
      # Subset rows with duplicated gene symbols
      duplicat_deg <- common [only_dupgene_row,]

      # Iterate over each duplicated gene symbol, calculate column means, and store the averaged gene information in a list:
      # Create a list to store averaged genes
      averaged_genes <-list()

      # Iterate over each duplicated gene symbol
      for(b in 1:length(only_dupgene)){
	# Identify rows for the current duplicated gene
        d_row <- which(duplicat_deg$Gene.symbol %in% only_dupgene[b])
	# Calculate column means
        col_mean <- rbind(colMeans(duplicat_deg[d_row, 1:2]))
	# Combine column means with other columns
        col_all <- cbind(col_mean, duplicat_deg[d_row[1],3:5])
	# Store the averaged gene information
        averaged_genes[[b]] <- col_all
      }

      # Combine averaged genes and unique genes into a single data frame and sort by z-score in descending order:
      # Combine averaged genes into a single data frame
      averaged_genes <- do.call(rbind, averaged_genes)
      # Subset rows with unique gene symbols
      unique_deg <- common [-only_dupgene_row, ]
      # Combine unique genes and averaged genes
      all_unique_genes <- rbind(unique_deg, averaged_genes)
      # Reset row names
      rownames(all_unique_genes) <- NULL 
      # Sort the genes by z-score in descending order
      all_unique_genes <- all_unique_genes[order(all_unique_genes[,1], decreasing = TRUE),]

      # Write the result to a tsv file
      write_tsv(all_unique_genes, paste0("MetaAnalysis/volcano/", projname [i], ".tsv"))

      # Apply FDR threshold
      all_unique_genes <- all_unique_genes[which(all_unique_genes[,2] < FDR),]
      # Write the filtered result to a file
      write_tsv(all_unique_genes, paste0("MetaAnalysis/filtered/", projname [i], ".tsv"))
    }
	
    # Count genes with "///"
    eslash_genes <- length(grep("///", all_unique_genes$Gene.symbol))

    # Display summary message for the current project
    message ("###################### Summary for: ", projname [i], " ############################")

    # Display summary information about the genes analysis
    message (cat("Number of missing gene symbols:", length(genemissing),
                "\nNumber of duplicated genes:", length(only_dupgene),
                "\nTotal number of genes after averaging duplicated genes and FDR filter:", NROW(all_unique_genes),
                "\nNumber of genes having ///:", eslash_genes))
	
    message ("##########################################################################")
  }

    # Display completion message for the analysis
    message ("+++++++++++++++++++++++ Analysis has completed. +++++++++++++++++++++++")
}


###################################################################################################################################### 
###################################################### 14- Makevolcano Function####################################################### 
# Function to generate a volcano plot

# Input: 
#   - projname: Name of the project (default: "GSE11121") 
#   - padjlevel: Adjusted p-value threshold (default: 0.05) 
#   - pSE_thresh: Positive effect size threshold (default: 1) 
#   - nSE_thresh: Negative effect size threshold (default: -1) 
#   - ntop: Number of top genes to label (default: 10) 
#   - voltype: Type of volcano plot ("ggplot" or "plotly") (default: "ggplot") 

# Output: 
#   - A volcano plot visualization 
#   - DEG files (upregulated, downregulated, and top genes) 

# Function Description: 
# This function generates a volcano plot for differential gene expression analysis.  
# It takes several arguments to customize the plot, including the project name, 
# adjusted p-value threshold, effect size thresholds, number of top genes to label,  
# and the type of volcano plot (ggplot or plotly). 

# The function reads the volcano plot data, calculates adjusted p-values and assigns DEG labels based on the thresholds.  
# It then generates the volcano plot using either ggplot or plotly, visualizing the relationship between effect size and statistical significance.  
# The plot includes color-coded points representing upregulated, downregulated, and non-significant genes. 
# Additionally, the function saves DEG files containing upregulated, downregulated, and top genes,  
# and displays the number of upregulated and downregulated genes. 
  	
Makevolcano <- function(projname  = "GSE11121",
                        # padjlevel: Adjusted p-value threshold (default: 0.05)
                        padjlevel = 0.05,
                        # pSE_thresh: Positive effect size threshold (default: 1)
                        pSE_thresh = 1,
                        # nSE_thresh: Negative effect size threshold (default: -1)
                        nSE_thresh = -1,
                        # ntop: Number of top genes to label (default: 10)
                        ntop = 10,
                        # voltype: Type of volcano plot ("ggplot" or "plotly") (default: "ggplot")
                        voltype = "ggplot"){

  # Create a directory for storing DEG files
  dir.create("MetaAnalysis/DEG", showWarnings = FALSE)
	
  # Reading the DEGs table:
  # Read the DEGs table from a TSV file
  volcano <- read_tsv(paste0("MetaAnalysis/volcano/", projname , ".tsv"))
  # Convert the table to a data frame	
  volcano <- data.frame (volcano)

  # Remove project name from column names
  colnames(volcano)[1:2] <- gsub(paste0("_", projname), "", colnames(volcano)[1:2])

  # Calculate adjusted p-values and assign DEG labels based on the threshold
  volcano <- volcano %>%
    # Calculate adjusted p-values using the negative logarithm of FDR values
    mutate(AdjustedPvalue = -log10(FDR)) %>%
    # Assign initial label of "NotSignificant" to all genes
    # (initially, indicating that they are not significantly differentially expressed
    # before applying further criteria to determine if a gene should be labeled as "Upregulated" or "Downregulated")
    mutate(DEG = "NotSignificant") %>%
    # Assign "Upregulated" label to genes that meet adjusted p-value and positive effect size thresholds
    mutate(DEG = ifelse(AdjustedPvalue> -log10(padjlevel) & zSco> pSE_thresh, "Upregulated", DEG)) %>%
    # Assign "Downregulated" label to genes that meet adjusted p-value and negative effect size thresholds
    mutate(DEG = ifelse(AdjustedPvalue> -log10(padjlevel) & zSco < nSE_thresh, "Downregulated", DEG))
  
  # Handle zero FDR values to avoid infinite AdjustedPvalue: 
  # Identify the indices of rows with FDR equal to 0
  zero_fdr <- which(volcano$FDR == 0)
  # Calculate the minimum non-zero FDR and adjust it to avoid infinite AdjustedPvalue
  m <- -log10(min(volcano[-zero_fdr, "FDR"])) + 0.5
  # Replace the AdjustedPvalue of zero FDR rows with the adjusted value
  volcano[zero_fdr, "AdjustedPvalue"] <- m

  # Define a vector of color codes to be used for plotting
  my_pal <- c("#1B9E77","#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#9A7D0A")

   Calculating the number of each Upregulated, Downregulated, and NotSignificant probes#
  # Select rows from the volcano data frame where DEG is "Upregulated"
  volcano_up <- volcano[which(volcano$DEG == "Upregulated"),]
  # Select rows from the volcano data frame where DEG is "Downregulated"
  volcano_down <- volcano[which(volcano$DEG == "Downregulated"),]
  # Sort the volcano_down data frame based on the zSco column in ascending order
  volcano_down <- volcano_down[order(volcano_down$zSco, decreasing = FALSE), ]
  # Select rows from the volcano data frame where DEG is "NotSignificant"
  volcano_NA <- volcano[which(volcano$DEG == "NotSignificant"),]
	
  # Combine the data frames of upregulated, downregulated, and non-significant genes into a single data frame
  volcano_all <- rbind(volcano_up, volcano_down, volcano_NA)
  # Replace the labels of the top probes with the word "Top"
  ntop_rep <- rep ("Top", ntop)
  # Create the DEG column in the volcano_all data frame by combining labels from different data frames
  volcano_all$DEG <- c(ntop_rep, volcano_up$DEG[-c(1:ntop)], ntop_rep, volcano_down$DEG[-c(1:ntop)], volcano_NA$DEG)
  # Create a column 'Top' in the 'volcano_all' data frame
  # Assign gene symbols to the 'Top' column for genes labeled as 'Top'
  volcano_all$Top <- ifelse(volcano_all$DEG == "Top", volcano_all$Gene.symbol, "")
	
  # Plotting either with plotly or ggplot:
  # Check if the specified voltype is "ggplot"
  if(voltype == "ggplot"){
    # Create a ggplot object for generating the volcano plot using ggplot2
    p <- ggplot(data = volcano_all, aes(x = zSco, y= AdjustedPvalue, color= DEG, fill = DEG, label= Top)) + 
      # Set the axis labels
      labs(x= 'Effect size', y= "-log10(FDR)") +
      # Add points to the plot
      geom_point(size = 1, shape = 21) +
      # Add text labels to the plot with overlap checking  
      geom_text(check_overlap = TRUE,vjust = 0.1, nudge_y = 0.1) +
      # Set manual color scale for DEG categories
      scale_color_manual(values=c(my_pal)) +
      # Set manual fill color scale for DEG categories
      scale_fill_manual(values=c(paste(my_pal, "66", sep = ""))) +
      # Apply a classic theme to the plot
      theme_classic() +
      # Customize the appearance of the axes, text, and title
      theme(
	# Customize axis text appearance
	axis.text = element_text(family = "Arial",size = 24 , colour = "black"),
	# Customize x-axis text appearance 
        axis.text.x = element_text(family = "Arial",colour = "black", size = 24),
	# Customize y-axis text appearance
        axis.text.y = element_text(family = "Arial",colour = "black", size = 24),
	# Customize plot subtitle appearance 
        plot.subtitle = element_text(family = "Arial",size = 24, colour = "black", hjust = 0.5),
	# Customize y-axis title appearance
        axis.title.y = element_text(family = "Arial", size = 24, angle = 90),
	# Customize x-axis title appearance
        axis.title.x = element_text(family = "Arial", size = 24, angle = 00),
	# Customize legend text appearance
        legend.text = element_text(size = 10, family = "Arial"), 
	# Customize legend title appearance
        legend.title = element_text(size = 20, family = "Arial")
      ) +
      # Set the plot subtitle
      labs(subtitle = paste0("Volcano plot-", projname))
    # Return the ggplot object
    p
	
  } else {
    # Create a ggplot object, add layers and annotation to the plot
    p <- ggplot(
        data = volcano_all, 
	aes(
	  # Set the x-axis variable as 'zSco' (effect size)
	  x = zSco, 
	  # Set the y-axis variable as 'AdjustedPvalue' (-log10(FDR))
	  y= AdjustedPvalue, 
	  # Map the 'DEG' variable to the color aesthetic
	  color= DEG, 
	  # Map the 'DEG' variable to the fill aesthetic
	  fill = DEG, 
	  # Create text labels by concatenating gene information
	  text = paste(
	    # Include the gene symbol
            "<br>Gene: " , Gene.symbol,
	    # Include the chromosome location
      	    "<br>Chromosome: " , Chromosome.location,
	    # Include the GO function
            "<br>GO function: " , GO.function 
	  )
	)
      ) +
	# Set the x and y-axis labels
	labs(x= 'Effect size', y= , title = paste0("Volcano plot-", projname)) +
	# Add points to the plot
        geom_point(size = 1, shape = 21) +
	 # Set the color aesthetics using DEG
	scale_color_manual(values=c(my_pal)) +
	# Set the fill aesthetics using DEG
	scale_fill_manual(values=c(paste(my_pal, "66", sep = ""))) +
	# Apply a classic theme to the plot
	theme_classic() +
	# Customize the appearance of the axes, text, and title
	theme(
	  # Customize axis text font family, size, and color
	  axis.text = element_text(family = "Arial",size = 24 , colour = "black"),
	  # Customize x-axis text font family, color, and size
	  axis.text.x = element_text(family = "Arial",colour = "black", size = 24),
	  # Customize y-axis text font family, color, and size
	  axis.text.y = element_text(family = "Arial",colour = "black", size = 24),
	  # Customize plot subtitle font family, size, color, and horizontal justification
	  plot.subtitle = element_text(family = "Arial",size = 24, colour = "black", hjust = 0.5),
	  # Customize y-axis title font family, size, and angle
	  axis.title.y = element_text(family = "Arial", size = 24, angle = 90),
	  # Customize x-axis title font family, size, and angle
	  axis.title.x = element_text(family = "Arial", size = 24, angle = 00)
	)

    # Convert the ggplot object to plotly with text tooltip
    p <- ggplotly(p, tooltip = "text") %>%
      # Customize the layout settings for the plotly plot
      layout(
	# Set the title of the plot
	title = paste0("Volcano plot-", projname),
	# Customize font family, color, and size
        font =list(family = "Arial", color = "black", size = 24),
	# Customize legend settings
        legend =list(family = "Arial", color = "black", size = 20, itemclick = "toggle")
      )
  }

 # Calculate the number of up-regulated and down-regulated genes
  write_tsv(rbind(volcano_up, volcano_down), paste0("MetaAnalysis/DEG/", projname , "_up_down.tsv"))
	
  # Create the top table (up and down) based on the specified number of top genes (ntop)
  Top <- rbind(volcano_up[1:ntop, ], volcano_down[1:ntop, ])
  Top <- Top[,-6]
	
  # Save the top table into the project folder
  write_tsv(Top, paste0("MetaAnalysis/DEG/", projname , "_top.tsv"))

  # Print the number of up-regulated and down-regulated genes
  message ("################# Significance #################")
  message (cat("Up-regulated genes:", nrow(volcano_up),
              "\nDown-regulated genes:", nrow(volcano_down)))
  message ("################################################")

  # Return the generated plot	
  return(p)
}


	
###################################################################################################################################### 
####################################################### 15- MakeVenna Function######################################################## 
# Function to generate a volcano plot 

# Arguments: 
#    - projname: Name of the project (default: "GSE11121") 
#    - padjlevel: Adjusted p-value threshold (default: 0.05) 
#    - pSE_thresh: Positive effect size threshold (default: 1) 
#    - nSE_thresh: Negative effect size threshold (default: -1) 
#    - ntop: Number of top genes to label (default: 10) 
#    - voltype: Type of volcano plot ("ggplot" or "plotly") (default: "ggplot") 
 
# Input: Dataframe containing the differential expression analysis results 

# Output: Volcano plot object (ggplot or plotly) 

# Function Description: 
# This function generates a volcano plot based on differential expression analysis results.  
# It takes several arguments to customize the plot, including the project name, adjusted p-value threshold,  
# effect size thresholds, number of top genes to label, and the type of volcano plot (ggplot or plotly). 

# The function reads the differential expression analysis results from a dataframe.  
# It then identifies the common genes across all projects and creates separate files for common genes and top common genes.  
# The function calculates the mean values of z-score and FDR for the common genes and creates a new dataframe.  
# The common genes, as well as the Venn diagram showing the overlapping genes between projects, are saved in the "common" folder. 

# Finally, the function generates a volcano plot using either ggplot or plotly, visualizing the relationship between effect size and statistical significance of genes.  
# The plot includes color-coded points representing upregulated, downregulated, and non-significant genes. The function returns the volcano plot object. 

# Function to generate a volcano plot based on differential expression analysis results 
MakeVenna <- function(common,
                      ntop = 10){

  # Create a directory named "common" inside the "MetaAnalysis" folder to store common gene files
  dir.create("MetaAnalysis/common", showWarnings = FALSE)
	
  # Reading all top-down files for each project 
  # Create an empty list to store the differential expression gene tables for each project 
  com_table <-list()

  # Iterate over each project in the 'common' vector 
  for(i in seq_along(common)){
    # Read the TSV file for the current project and convert it to a data frame 
    DEG_table <- data.frame (read_tsv(paste0("MetaAnalysis/DEG/", common [i], "_up_down.tsv")))
    # Store the DEG table in the com_table list, using the project name as the list element name 
    com_table[[common [i]]] <- DEG_table[,-6]
  }
	
  # Calculating up and down-regulated genes for each project: 
  # Create an empty list to store the up and down-regulated genes for each project 
  DEGslist <-list()

  # Iterate over each project in the com_table list 
  for(i in seq_along(com_table)){
    # Store the up and down-regulated genes in the DEGslist 
    DEGslist[[common [i]]] <- com_table[[i]][[3]]
  }
	
  # Finding common genes across all projects 
  com_prob <- Reduce(intersect, DEGslist)

	  

  # Writing common genes into separate files and generating top files 
  # Create an empty list to store the common genes across all studies   
  all_com <-list()
  # Loop over each study in the com_table list 	
  for(i in seq_along(com_table)){
    # Find the indices of common genes in the current project's DEG table 
    comn <- which(com_table[[names(com_table)[i]]][["Gene.symbol"]] %in% com_prob)
    # Extract the common genes from the DEG table 
    com_final_DEG  <- com_table[[names(com_table)[i]]][comn, ]
    # Write the common genes into a separate file 
    write_tsv(com_final_DEG, paste0("MetaAnalysis/common/",names(com_table)[i],"_common_DEGs.tsv"))
    # Sort the common genes based on gene symbol 
    sorted_com_final_DEG <- com_final_DEG[order(com_final_DEG$Gene.symbol),]
    # Store the sorted common genes in the all_com list 
    all_com[[names(com_table)[i]]] <- sorted_com_final_DEG

    # Extract the top upregulated genes from the common genes 
    com_final_up <- com_final_DEG[which(com_final_DEG$DEG == "Upregulated"), ]
    # Selecting the top ntop upregulated genes 
    com_final_up <- com_final_up[1:ntop, ]

    # Extract the top downregulated genes from the common genes 
    com_final_down <- com_final_DEG[which(com_final_DEG$DEG == "Downregulated"), ]
    # Select the top ntop rows from the com_final_down dataframe 
    com_final_down <- com_final_down[1:ntop, ]

    # Write the top common genes (upregulated and downregulated) into a separate file 
    write_tsv(rbind(com_final_up, com_final_down), paste0("MetaAnalysis/common/", names(com_table)[i],"_top_common_DEGs.tsv"))
  }

  # Separating the first two columns and remaining columns for each project 
  # Create an empty list to store the first two columns of each dataframe in all_com 
  firstsecond_col <-list()

  # Create an empty list to store the remaining columns (3 to 6) of each dataframe in all_com 
  remaining <-list()

  # Iterate over the indices of all_com list 
  for(i in seq_along(all_com)){
    # Extract the first two columns of the i-th dataframe in all_com and store it in firstsecond_col list 
    firstsecond_col[[i]] <- data.frame (all_com[[i]][,1:2])
    # Extract the first two columns of the i-th dataframe in all_com and store it in firstsecond_col list 
    remaining[[i]] <- data.frame (all_com[[i]][,3:6])
  }

  # Combine the first two columns of the data frames in the firstsecond_col list and calculate the mean values of z-score and FDR, and 
  # create a new data frame (sz_fdr_df) containing the mean values and combine it with the remaining columns: 
  # Combine the first two columns of the data frames in firstsecond_col list into a single data frame 
  firstsecond_col_df <- do.call(cbind, firstsecond_col)
  # Find the column indices containing "zSco" in their column names 
  sz <- grep("zSco", colnames(firstsecond_col_df))
  # Find the column indices containing "FDR" in their column names 
  fdr <- grep("FDR", colnames(firstsecond_col_df))
  # Calculate the mean of z-score for each row 
  sz_mean <- rowMeans(firstsecond_col_df[,sz])
  # Calculate the mean of FDR for each row 
  fdr_mean <- rowMeans(firstsecond_col_df[,fdr])
  # Create a new data frame with columns for zSco and FDR means 
  sz_fdr_df <- data.frame (zSco = sz_mean, FDR= fdr_mean)
  # Combine the sz_fdr_df with the remaining columns from remaining list 
  sz_fdr_df <- cbind(sz_fdr_df, remaining[[1]])

  # Write the sz_fdr_df to a TSV file 
  write_tsv(sz_fdr_df, paste0("MetaAnalysis/common/all_mean_common_DEGs.tsv"))
  
  # Display a message indicating that common genes between all datasets have been saved into the common folder 
  message ("+++++ Common genes between all datasets have saved into common folder. +++++")
	
  # Creating venn diagram: 
  # Create an empty list to store the DEGs for each project 
  DEGslist <-list()

  # Iterate over the projects in the com_table 
  for(i in seq_along(com_table)){
    # Extract the DEGs for each project and store them in the DEGslist 
    DEGslist[[common [i]]] <- com_table[[i]][3]
  }

  # Create an empty list to store the Venn diagram input for each project 
  venlist <-list()

  # Extract the common DEGs for each project and store them in the venlist: 
  # Iterate over the elements in the DEGslist 
  for(i in seq_along(DEGslist)){
    # Assign the common elements from DEGslist to the corresponding element in venlist 
    venlist[[common [i]]] <- DEGslist[[i]][[1]]
  }

  # Set the seed for reproducibility in generating the Venn diagram 
  set.seed(1234)

  # Define the color palette for the Venn diagram 
  my_pal <- c("#1B9E77","#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#9A7D0A")

  # Generate the Venn diagram using the specified parameters 
  venaa <- ggvenn(venlist, 
	          # Set the fill colors of the Venn diagram circles using the my_pal vector 
                  fill_color = my_pal,
	          # Set the stroke size of the Venn diagram circles to 0.5 
                  stroke_size = 0.5,
	          # Set the size of the set names in the Venn diagram to 4 
                  set_name_size = 4)

  # Return the Venn diagram object 
  return(venaa)
}



###################################################################################################################################### 
################################################### 16- Makenetworkanalyst Function###################################################  
# Function to prepare data for NetworkAnalyst platform 

# Arguments: 
#   - matrix: Input matrix containing gene expression data 
#   - projname: Name of the project 
#   - compare: Comparison type (default: "grade") 
#   - groups: List of groups for comparison 

# Input: Gene expression matrix, phenotype data, annotation data 

# Output: Modified matrix file for NetworkAnalyst 

# Function Description: 
# This function prepares the gene expression data for the NetworkAnalyst platform.  
# It takes a gene expression matrix, project name, comparison type, and groups for comparison as input.  
# The function performs several steps to process the data and generate a modified matrix file compatible with NetworkAnalyst. 

# The function first creates a directory to store the network analysis files for NetworkAnalyst.  
# It then checks if the first column of the matrix is labeled as "ID" and adjusts the matrix accordingly.  
# The phenotype file is read and processed based on the specified comparison type.  
# The annotation file is also read and processed to ensure compatibility. 

# Next, the function finds matching samples based on the provided groups and subsets the matrix accordingly.  
# The values in the matrix are rounded to two decimal places. Rows with missing gene symbols are removed from the matrix.  
# Group information is extracted from the phenotype data and stored in the grouplis list. 

# The function renames the matrix file according to the group names specified by the user.  
# It modifies the column names and the second row to adhere to the file format required by NetworkAnalyst.  
# Finally, the modified matrix file is saved in the "MetaAnalysis/networkanalyst" directory with the name "projname_network_full.txt". 

# The completion message is displayed when the file has been successfully saved. 
	
Makenetworkanalyst <- function(matrix, projname , compare  = "grade", groups){
	
  # Creating directory to store network analysis files for NetworkAnalyst 
  dir.create("MetaAnalysis/networkanalyst", showWarnings = FALSE)

  # Check if the first column of the matrix is labeled as "ID" 
  if(colnames(matrix)[1] == "ID"){
    # Assign the values of the first column as row names 
    rownames(matrix) <- matrix[,1]
    # Remove the first column from the matrix 
    matrix <- matrix[,-1]
  }
	
  # Reading and processing phenotype file 
  # Read the phenotype file based on the comparison type 
  pheno <- read_csv (paste0("phenotype/phenotype_", compare , ".csv"))
  # Convert the read data into a data frame 
  pheno <- data.frame (pheno)
  # Find the indices of the samples that match the columns of the matrix 
  paired_pheno <- which(pheno$Sample %in% colnames(matrix))
  # Subset the phenotype data based on the matched samples 
  pheno <- pheno[paired_pheno, ]

  # Reading and processing annotation file: 
  # Read the annotation file and create a data frame 
  annot<- read_tsv(paste0(projname , "/", "annot.tsv"))
  annot<- data.frame (annot)

  # Set the row names of the annotation data frame 
  rownames(annot) <- annot[,1]

  # Remove the first column from the annotation data frame 
  annot<- annot[,-1]

	
  # Finding matching grades or subtypes: 
  # Create an empty list to store sample IDs for each group 
  sampleid <-list()

  # Iterate over the groups and find matching samples in the phenotype data 
  for(i in seq_along(groups)){
    # Find rows in the phenotype data matching the group 
    g_row <- which(pheno $grade %in% groups[[i]])
    # Store the matched samples in the list 
    sampleid[[i]] <- pheno[g_row, ]
  }

  # Combine the sample IDs into a single data frame 
  sampleid <- do.call(rbind, sampleid)

  # Sort the sample IDs by the "Sample" column
  sampleid <- sampleid[order(sampleid$Sample),]

  # Subset the matrix using the sample IDs 
  matrix <- matrix[, sampleid$Sample]

  # Round the values in the matrix to two decimal places 
  matrix <- round(matrix, 2)

  # Remove rows with missing gene symbols 
  identical(rownames(matrix), rownames(annot))
  matrix <- cbind(annot$Gene.symbol, matrix)

  # Rename the first column as "Gene.symbol" 
  colnames(matrix)[1] <- "Gene.symbol"
  # Find the indices of missing gene symbols in the matrix 
  genemissing <- which(is.na(matrix$Gene.symbol) == TRUE)
  # Remove rows with missing gene symbols from the matrix 
  matrix <- matrix[-genemissing, ]

  # Create an empty list to store group information 
  grouplis <-list()

  # Iterating over each group in the 'groups' list 
  for(i in seq_along(groups)){
    # Combine group names into a single string separated by "|" 
    name <- paste0(groups[[i]], collapse = "|")
    # Find rows in the phenotype data matching the group 
    rownum <- grep(name, pheno$grade)
    # Store the matched sample names in the grouplis list 
    grouplis[[names(groups)[i]]] <- pheno$Sample[rownum]
  }

	
  # Rename the matrix file according to the group names specified by the user 
  updown_matrix <- rbind(colnames(matrix), matrix)

  # Iterating over each element in the 'grouplis' list 
  for(i in seq_along(grouplis)){
    # Combine sample names for each group into a single string separated by "|" 
    matname <- paste0(grouplis[[i]], collapse = "|")
    # Find column indices in the matrix matching the group 
    matcoln <- grep(matname, colnames(updown_matrix))
    # Rename the corresponding columns with group names 
    updown_matrix[1,matcoln] <- names(grouplis)[i]
  }

	
  # Modify column names and the second row according to the file format for NetworkAnalyst 
  # Find column indices with group names 
  changedcol <- grep(paste0(names(grouplis), collapse = "|"), updown_matrix[1,])
  # Subset the matrix with the modified column names 
  updown_matrix <- updown_matrix[, c(1, changedcol)]
  # Rename the first cell in the matrix as "#CLASS" 
  updown_matrix[1,1] <- "#CLASS"
  # Rename the first column as "#NAME" 
  colnames(updown_matrix)[1] <- "#NAME"

  # Save the modified matrix file for NetworkAnalyst 
  write.table(updown_matrix, paste0("MetaAnalysis/networkanalyst/", projname ,"_network_full.txt"), quote = FALSE, row.names = FALSE, sep = "\t" )

  # Displaying a message indicating that the file has been saved 
  message ("+++++ Your file has saved into networkanalyst folder. +++++")
} 



###################################################################################################################################### 
###################################################### 17- MakeHeatmap Function ###################################################### 
# Function to create a heatmap based on gene expression data 

# Arguments: 
#   - matrix: The gene expression data matrix 
#   - projname: The project name or identifier 
#   - compare: The comparison for phenotype 
#   - groups: A list of groups for sample selection 
#   - ntop: The number of top genes to consider 
#   - sorted: A logical value indicating whether the heatmap should be sorted 

# Input: 
#   - Gene expression data matrix 
#   - Phenotype data 
#   - Annotation data 
#   - Differential expression gene (DEG) file 

# Output: 
#   - Heatmap plot visualizing the gene expression data 

# Function Description: 
# This function generates a heatmap based on gene expression data.  
# It takes a gene expression data matrix, project name, comparison for phenotype, groups for sample selection,  
# number of top genes to consider, and a logical value indicating whether the heatmap should be sorted as input. 

# The function starts by checking if the matrix's first column is named "ID" and adjusts the matrix if necessary.  
# It reads the phenotype CSV file and converts it to a data frame.  
# The annotation file is also read and processed based on the project name. 
# The matrix is adjusted to match the IDs in the annotation data frame. 

# Next, the function performs averaging for duplicated genes by calculating the column-wise mean for rows with the same gene symbol.  
# The averaged genes are combined with unique non-duplicated genes into a single data frame. 
# The DEG file is read, and the top upregulated and downregulated genes are extracted. 

# The function selects the genes from the combined data frame that match the selected DEGs.  
# The network data frame is created, and the sample IDs are matched with the network data frame columns based on the groups provided.  
# The grades are extracted from the network data frame columns. 

# The function generates the heatmap plot using the network data frame.  
# If the heatmap should be sorted, the columns are ordered based on column names.  
# The plot includes row (gene) annotations, a bottom annotation indicating the grades, and a legend for expression levels.  
# The row title is positioned on the left side, and the column names are hidden. 

# The completion message is displayed when the heatmap plot is generated. 
	
MakeHeatmap <- function(matrix, projname , compare , groups, ntop,  sorted){

  # Set the random seed for reproducibility 
  set.seed(1234)

 # Check if the first column of the matrix is named "ID" 
 if(colnames(matrix)[1] == "ID"){
    # If so, adjust the matrix to treat the IDs as row names: 
    # Set row names of the matrix to the values in the first column 
    rownames(matrix) <- matrix[,1]
    # Remove the first column from the matrix 
    matrix <- matrix[,-1]
    # Sort the matrix based on the row names
    matrix <- matrix[order(rownames(matrix)),]
  }

	
  # Read the phenotype CSV file and convert it to a data frame: 
  # Read the phenotype CSV file based on the "compare" value 
  pheno <- read_csv (paste0("phenotype/phenotype_", compare , ".csv"))
  # Convert the read data to a data frame
  pheno <- data.frame (pheno)

	
  # Read annotation file based on the "projname" value:
  # If the condition is true (if projname is equal to "combination"), perform the following operations 
  if(projname  == "combination"){
    # For "combination" project, read from "GSE25055/annot.tsv"
    annot<- data.frame (read_tsv(paste0("GSE25055", "/", "annot.tsv")))
  } else{
    # Read the annotation TSV file based on the "projname" value
    annot<- data.frame (read_tsv(paste0(projname , "/", "annot.tsv")))
  }
	
  # Select only the "ID" and "Gene.symbol" columns from the annotation data frame
  annot<- annot[ , c("ID","Gene.symbol")]
  # Sort the annotation data frame based on the "ID" column
  annot<- annot[order(annot$ID), ]

 # Check if the row names of the matrix match the IDs in the annotation
 if(identical(rownames(matrix), annot$ID)){
    # If so, add the "Gene.symbol" column to the matrix and adjust column names
    matrix <- cbind(annot$Gene.symbol, matrix)
    # Set the column name of the matrix at index 1 to "Gene.symbol"
    colnames(matrix)[1] <- "Gene.symbol"
  }

  # Perform averaging for duplicated genes:
  # Create a data frame with frequency counts of duplicated genes
  dupgene_freq <- data.frame (table(matrix$Gene.symbol))
  # Create a data frame with frequency counts of duplicated genes
  only_dupgene <- dupgene_freq[which(dupgene_freq$Freq > 1), 1]
  # Find the row indices of the matrix where the gene symbol matches the duplicated genes
  only_dupgene_row <- which(matrix$Gene.symbol %in% only_dupgene)
  # Extract the duplicated genes and their corresponding rows from the matrix
  duplicat_deg <- matrix[only_dupgene_row,]

  # Create an empty list to store the averaged genes
  averaged_genes <-list()

  # Iterate over each duplicated gene
  for(i in 1:length(only_dupgene)){
    # Find the row indices of the duplicated gene in the duplicat_deg data frame
    d_row <- which(duplicat_deg$Gene.symbol %in% only_dupgene[i])
    # Calculate the column-wise mean for the selected rows
    col_mean <- rbind(colMeans(duplicat_deg[d_row, -1]))
    # Combine the averaged gene data with its corresponding gene symbol
    col_all <- cbind(duplicat_deg[d_row[1],1], col_mean)
    # Add the averaged gene to the list
    averaged_genes[[i]] <- col_all
  }

  # Combine all the averaged genes into a single data frame
  averaged_genes <- data.frame (do.call(rbind, averaged_genes))
  # Set the column name of the averaged genes at index 1 to "Gene.symbol"
  colnames(averaged_genes)[1] <- "Gene.symbol"
  # Extract the unique non-duplicated genes from the matrix
  unique_deg <- matrix[-only_dupgene_row, ]
  # Combine the unique genes and averaged genes into a single data frame
  all_unique_genes <- rbind(unique_deg, averaged_genes)
  # Reset the row names of the combined genes data frame
  rownames(all_unique_genes) <- NULL 
	
  # Reading DEG file
  DEG_table <- data.frame (read_tsv(paste0("MetaAnalysis/DEG/", projname , "_up_down.tsv")))

	
  # Finding up and down regulated gene and make original matrix of them
  # Extract the upregulated genes from the DEG table based on the "DEG" column value
  up <- DEG_table[which(DEG_table$DEG == "Upregulated"),]
  # Select the top "ntop" rows and the second and third columns (gene symbol and expression values)
  up <- up[1:ntop,2:3]

  # Extract the downregulated genes from the DEG table based on the "DEG" column value
  down <- DEG_table[which(DEG_table$DEG == "Downregulated"),]
  # Select the top "ntop" rows and the second and third columns (gene symbol and expression values)
  down <- down[1:ntop, 2:3]

  # Combine the upregulated and downregulated genes into a single gene list column	
  genelist <- rbind(up, down)[,2]
  # Find the row indices of the genes in the all_unique_genes data frame
  num <- which(all_unique_genes$Gene.symbol %in% genelist)
  # Extract the network data frame with only the selected genes
  network_df <- data.frame (all_unique_genes[num,])
  # Extract the gene symbols column from the network data frame
  Gene.symbol <- network_df$Gene.symbol
  # Remove the gene symbols column from the network data frame
  network_df <- network_df[,-1]
  # Convert the network data frame to a numeric matrix
  network_df <- as.matrix(sapply(network_df, as.numeric))
  # Set the row names of the network matrix to the gene symbols
  rownames(network_df) <- Gene.symbol

	
  # Finding grades and matching with colnames:
  # Create an empty list to store sample IDs for each group
  sampleid <-list()

  # Iterate over each group and find the corresponding samples in the phenotype data
  for(i in seq_along(groups)){
    # Find the row indices in pheno where grade matches the current group
    g_row <- which(pheno $grade %in% groups[[i]])
    # Store the corresponding sample IDs for the current group in the sampleid list
    sampleid[[i]] <- pheno[g_row, ]
  }

  # Combine the sample IDs from different groups into a single data frame
  sampleid <- do.call(rbind, sampleid)
  # Find the column indices of the samples in the network_df matrix
  sampleid_n <- which(sampleid$Sample %in% colnames(network_df))
  # Filter the sample IDs based on the column indices
  sampleid <- sampleid[sampleid_n,]
  # Sort the sample IDs based on the "grade" column
  sampleid <- sampleid[order(sampleid$grade),]
  # Add a prefix "Grade_" to the values in the "grade" column
  sampleid $grade <- paste0("Grade_", sampleid$grade) 
  # Extract the corresponding columns from the network_df matrix based on the sample IDs
  network_df <- network_df[, sampleid$Sample]
  # Set the column names of the network_df matrix to the modified "grade" values
  colnames(network_df) <- sampleid $grade 
	
  # Grade info
  # Extract the column names (grades) of the network_df matrix
  type <- colnames(network_df)

  # Create a HeatmapAnnotation object with a data frame and annotation height
  ha = HeatmapAnnotation(
    df = data.frame (Grade = type),
    annotation_height = unit(4, "mm")
  )

	
  # Heatmap plot
  # If the heatmap should be sorted, provide additional parameters for ordering
  if(sorted == TRUE){
  Heatmap(network_df,
	  # Set the properties for row names (font size)
          row_names_gp = gpar(fontsize = 7),
	  # Add the HeatmapAnnotation object as the bottom annotation
          bottom_annotation= ha,
	  # Define the order of columns based on column names
          column_order = colnames(network_df),
	  # Set the color of the border around the heatmap
          border_gp = gpar(col = "black"),
	  # Set the title for the row (gene) annotations
          row_title = "DEGs",
	  # Position the row title on the left side
          row_title_side = "left",
	  # Hide the column names in the heatmap
          show_column_names = FALSE,
	
          heatmap_legend_param =list(
                                      # Set the orientation of the legend (vertical)
	                              legend_direction = "vertical",
                                      # Set the height of the legend
                                      legend_height = unit(50, "mm"),
                                      # Set the width of each grid in the legend
                                      grid_width = unit(6, "mm"),
                                      # Set the height of each grid in the legend
                                      grid_height = unit(50, "cm"),
                                      # Set the title of the legend
                                      title = "Expression"))
  } else {
    # Create the heatmap plot without sorting if "sorted" is FALSE
    Heatmap(network_df,
	    # Set the properties for row names (font size)
            row_names_gp = gpar(fontsize = 7),
	    # Add the HeatmapAnnotation object as the bottom annotation
            bottom_annotation= ha,
	    # Set the color of the border around the heatmap
            border_gp = gpar(col = "black"),
	    # Set the title for the row (gene) annotations
            row_title = "DEGs",
	    # Position the row title on the left side
            row_title_side = "left",
	    # Hide the column names in the heatmap
            show_column_names = FALSE,
	
            heatmap_legend_param =list(
	                                # Set the orientation of the legend (vertical)
	                                legend_direction = "vertical",
	                                # Set the height of the legend
                                        legend_height = unit(50, "mm"),
	                                # Set the width of each grid in the legend
                                        grid_width = unit(6, "mm"),
	                                # Set the height of each grid in the legend
                                        grid_height = unit(50, "cm"),
	                                # Set the title of the legend
                                        title = "Expression"))
  }
}



###################################################################################################################################### 
########################################################## 18- GSEAAnalysis ########################################################## 
# Function to perform Gene Set Enrichment Analysis (GSEA) 

# Arguments: 
#   - projname: The project name or identifier 
#   - FDR: The false discovery rate (FDR) threshold for significance 

# Input: 
#   - DEG_table: Data frame containing differential expression gene (DEG) information 
#   - H: Data frame containing hallmark pathway gene sets 
#   - FC: Data frame containing gene names and logFC values 

# Output: 
#   - GSEA plot visualizing the enriched gene sets 

# Function Description: 
# This function performs Gene Set Enrichment Analysis (GSEA) using the fgseaSimple function. 
# It takes the project name and FDR threshold as input and generates a GSEA plot of enriched gene sets. 

# The function reads the DEG table specific to the project name and prepares gene scores and gene set lookup table. 
# GSEA is performed using fgseaSimple, and the results are visualized in a GSEA plot. 
# The significance of gene sets is determined based on FDR and NES values. 
# The plot includes pathway names, NES values, and colored bars representing significance. 

# Function to perform Gene Set Enrichment Analysis (GSEA)
GSEAAnalysis <- function(projname , FDR){
	
  # Read the DEG table file specific to the given project name and create a data frame
  DEG_table <- data.frame (read_tsv(paste0("MetaAnalysis/DEG/", projname , "_up_down.tsv")))
	
  # Removing /// from some genes
  DEG_table$Gene.symbol <- gsub("///.*", "", DEG_table$Gene.symbol)
	
  # Retrieve hallmark pathway gene sets specific to the "Homo sapiens" species
  H <- msigdbr(species = "Homo sapiens", category = "H")
  # Store the result of the following operations in H.symbol.ls
  H.symbol.ls <- H %>% 
    # Select the gs_name and gene_symbol columns from H
    select(gs_name, gene_symbol) %>% 
    # Group the data by gs_name
    group_by(gs_name) %>% 
    # Summarize the unique gene_symbols for each gs_name
    summarise(all.genes =list(unique(gene_symbol))) %>% 
    # Convert the summarized data into a lookup table format
    deframe()
	
  # Create a new data frame FC with columns 3 and 1 from DEG_table (gene names and logFC)
  FC <- DEG_table[,c(3,1)]
  
  # Create a vector FC.vec containing the values from the zSco column of the FC data frame
  FC.vec <- FC$zSco
  # Assign the Gene.symbol column values from the FC data frame as the names of the FC.vec vector
  names(FC.vec) <- FC$Gene.symbol
	
  # Perform Gene Set Enrichment Analysis (GSEA) using the fgseaSimple function
  gsea.H <- fgseaSimple(
                        # - The pathways argument specifies the gene set lookup table H.symbol.ls            
	                pathways = H.symbol.ls,
                        # - The stats argument provides the FC.vec vector with gene scores            
                        stats = FC.vec,
                        # - The scoreType is set to "std" for standard scoring            
                        scoreType = "std",
                        # - The nperm parameter is set to 1000 for the number of permutations to perform 
                        # in order to estimate the significance of the enrichment scores            
                        nperm=1000)

  # Create a vector my_pal with hexadecimal color codes representing a color palette
  my_pal <- c("#1B9E77","#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#9A7D0A")
	
  # Clean up pathway names in the gsea.H object:
  # Remove the prefix "HALLMARK_" from the pathway column in the gsea.H object
  gsea.H$pathway <- gsub("HALLMARK_","", gsea.H$pathway)
  # Replace underscores "_" with spaces " " in the pathway column of the gsea.H object
  gsea.H$pathway <- gsub("_"," ", gsea.H$pathway)
	
  # Calculating the Significance based on FDR and NES
  gsea.H <- gsea.H %>%
    # Update the Significance column in the gsea.H object with initial values of "NotSignificance"
    mutate(Significance = "NotSignificance") %>%
    # Update the Significance values based on the conditions:
    #   - If padj is less than or equal to FDR and NES is greater than 0, assign "Upregulated"
    mutate(Significance = ifelse(padj <= FDR & NES > 0, "Upregulated", Significance)) %>%
    #   - If padj is less than or equal to FDR and NES is less than 0, assign "Downregulated"
    mutate(Significance = ifelse(padj <= FDR & NES < 0, "Downregulated", Significance))

  # Filter the gsea.H object to retain only the rows where Significance is not equal to "NotSignificance
  gsea.H <- gsea.H[which(gsea.H$Significance != "NotSignificance"), ]
	
  # Create a ggplot visualization using the gsea.H object
  gsea.H %>% 
    # - The x-axis represents the pathway names reordered by NES values
    # - The y-axis represents the NES values
    ggplot(aes(x=reorder(pathway, NES), y=NES)) +
    # - Display the bars with a width of 0.5, coloring them based on the Significance column
    geom_col(width = 0.5, aes(color= Significance, fill = Significance)) +
    # - Add labels for the x-axis ("Gene set") and y-axis ("Normalized enrichment score (NES)")
    labs(x= 'Gene set', y= 'Gene set') +
    # - Apply a classic theme to the plot
    theme_classic() +
    # - Customize the colors of the bars and legend using the "blue" and "red" values:
    # Set the colors for data points using a manual color scale
    scale_color_manual(values=c("blue", "red")) +
    # Set the colors for filled areas using a manual fill scale
    scale_fill_manual(values=c("blue", "red")) +
    # - Flip the coordinates to have horizontal bars instead of vertical bars
    coord_flip() +

    # Customize the appearance of various plot elements (the axes, text, and title) using the theme function
    theme(
	  # - Modify the font family, size, and color of the axis text
	  axis.text = element_text(family = "Arial",size = 24 , colour = "black"),
	  # - Customize the font family, color, and size of the x-axis text
          axis.text.x = element_text(family = "Arial",colour = "black", size = 24),
	  # - Customize the font family, color, and size of the y-axis text 
          axis.text.y = element_text(family = "Arial",colour = "black", size = 9),
	  # - Adjust the appearance of the plot subtitle, including font family, size, color, and horizontal justification
          plot.subtitle = element_text(family = "Arial",size = 24, colour = "black", hjust = 0.5),
	  # - Customize the appearance of the y-axis title, including font family, size, and angle (90-degree rotation)
          axis.title.y = element_text(family = "Arial", size = 24, angle = 90),
	  # - Customize the appearance of the x-axis title, including font family, size, and angle (0-degree rotation)
          axis.title.x = element_text(family = "Arial", size = 24, angle = 00),
	  # - Modify the size and font family of the legend text
          legend.text = element_text(size = 10, family = "Arial"), 
	  # - Modify the size and font family of the legend title
          legend.title = element_text(size = 20, family = "Arial")) +
    # - Set the subtitle of the plot to include the project name (projname)
    labs(subtitle = paste0("GSEA - ", projname))
}


###################################################################################################################################### 
################################################### Running all functions ############################################################ 

# Downloading GEO data for specific GEO accession numbers
DownloadGEO("GSE25065")
DownloadGEO("GSE11121")
DownloadGEO("GSE7390")
DownloadGEO("GSE25055")
	
# Read the GEO matrices into the R environment for the specified GEO accession numbers
GSE11121 <- ReadGEO ("GSE11121")
GSE25055 <- ReadGEO ("GSE25055")
GSE25065 <- ReadGEO ("GSE25065")
GSE7390 <- ReadGEO ("GSE7390")
	
# Call the MakePhenotype function to create phenotype data for the specified GEO datasets
Makephenotype (projname  = c("GSE11121", "GSE25055", "GSE25065", "GSE7390"),
              compare  = "grade")

# Call the MakeBoxPlot function to create box plots for the specified GEO matrices and project names
MakeBoxPlot(GEOmatrix = GSE11121, projname  = "GSE11121")
MakeBoxPlot(GEOmatrix = GSE25065, projname  = "GSE25065")
MakeBoxPlot(GEOmatrix = GSE25055, projname  = "GSE25055")
MakeBoxPlot(GEOmatrix = GSE7390, projname  = "GSE7390")

# Call the MakePCA function to perform Principal Component Analysis (PCA) for the specified GEO matrices
# and project names, using the "grade" comparison
MakePCA(GEOmatrix = GSE11121, projname  = "GSE11121", compare  = "grade")
MakePCA(GEOmatrix = GSE25065, projname  = "GSE25065", compare  = "grade")
MakePCA(GEOmatrix = GSE25055, projname  = "GSE25055", compare  = "grade")
MakePCA(GEOmatrix = GSE7390, projname  = "GSE7390", compare  = "grade")

# Apply VSN (Variance Stabilizing Normalization) quantile normalization to the specified GEO matrices
normalized_GSE11121 <- VSNQuantilNorm(GEOmatrix = GSE11121)
normalized_GSE25055 <- VSNQuantilNorm(GEOmatrix = GSE25055)
normalized_GSE25065 <- VSNQuantilNorm(GEOmatrix = GSE25065)
normalized_GSE7390 <- VSNQuantilNorm(GEOmatrix = GSE7390)

# Generate box plots and perform PCA analysis for the specified normalized GEO matrices and project names after normalization
MakeBoxPlot(GEOmatrix =  normalized_GSE11121, projname  = "GSE11121")
MakePCA(GEOmatrix =  normalized_GSE11121, projname  = "GSE11121", compare  = "grade")

MakeBoxPlot(GEOmatrix = normalized_GSE25055, projname  = "GSE25055")
MakePCA(GEOmatrix = normalized_GSE25055, projname  = "GSE25055", compare  = "grade")

MakeBoxPlot(GEOmatrix = normalized_GSE25065, projname  = "GSE25065")
MakePCA(GEOmatrix = normalized_GSE25065, projname  = "GSE25065", compare  = "grade")

MakeBoxPlot(GEOmatrix = normalized_GSE7390, projname  = "GSE7390")
MakePCA(GEOmatrix = normalized_GSE7390, projname  = "GSE7390", compare  = "grade")

# Merging studies without batch correction (just for comparison)
merge_studies <- Mergestudies(list(normalized_GSE25055,  normalized_GSE11121, normalized_GSE25065, normalized_GSE7390))

# Perform PCA (before batch effect) on merged gene expression matrix using multiple studies and compare based on "grade"
# Important note : studies should match projname
MakePCA(GEOmatrix = merge_studies, 
        studies = list(normalized_GSE25055,  normalized_GSE11121, normalized_GSE25065, normalized_GSE7390), 
        projname  = c("GSE25055", "GSE11121", "GSE25065", "GSE7390"), 
        compare  = "grade")

# Batch effect correction
corbatch <- StudyBatchEffect(studies = list(normalized_GSE25055,  normalized_GSE11121, normalized_GSE25065, normalized_GSE7390))

# Perform PCA on batch-corrected data and compare based on "grade"
# Important note : studies should match projname 
MakePCA(GEOmatrix = corbatch, 
        studies = list(normalized_GSE25055,  normalized_GSE11121, normalized_GSE25065, normalized_GSE7390), 
        projname  = c("GSE25055", "GSE11121", "GSE25065", "GSE7390"), compare  = "grade")

# Create density plots before batch effect correction for gene expression in each study 
MakeDensityPlot(matrix = merge_studies, 
                studies = list(GSE25055 = normalized_GSE25055, 
                               GSE11121 =  normalized_GSE11121, 
                               GSE25065 = normalized_GSE25065,
                                GSE7390 = normalized_GSE7390))

# Density plot after batch effect correction
MakeDensityPlot(matrix = corbatch, 
                studies = list(GSE25055 = normalized_GSE25055, 
                               GSE11121 =  normalized_GSE11121, 
                               GSE25065 = normalized_GSE25065,
                                GSE7390 = normalized_GSE7390))

# MDS plot before batch effect correction
# Important note : studies should match projname 
MakeMDSPlot(matrix = merge_studies, 
            studies = list(normalized_GSE25055,  normalized_GSE11121, normalized_GSE25065, normalized_GSE7390), 
            projname  = c("GSE25055", "GSE11121", "GSE25065", "GSE7390"), 
            compare  = "grade",
            nstudies =  4)


# Create MDS (Multidimensional Scaling) plot, after batch effect correction, using merged gene expression data from 4 studies, comparing based on "grade"
# Important note : studies should match projname 
MakeMDSPlot(matrix = corbatch, 
            studies = list(normalized_GSE25055,  normalized_GSE11121, normalized_GSE25065, normalized_GSE7390), 
            projname  = c("GSE25055", "GSE11121", "GSE25065", "GSE7390"), 
            compare  = "grade",
            nstudies =  4)

# Perform Cochran's Q test on the studies, comparing groups "3" and "1" based on "grade"
cohranqtest <- CochranQTest(studies = list(GSE11121, GSE25065,GSE25055, GSE7390),
                            groups = c("3", "1"),
                            compare  = "grade")
	
# Create QQ plot using Cochran's Q test results from 4 studies
MakeQQPlot(cohranq = cohranqtest, nstudies =  4)

# Perform REM (Random Effects Model) Analysis on the studies, comparing groups "3" and "1" based on "grade" with FDR of 0.05 and 50 permutations.
FEMREMAnalysis(studies = list(normalized_GSE25055, normalized_GSE11121, normalized_GSE25065, normalized_GSE7390),
               projname  = c("GSE25055", "GSE11121","GSE25065", "GSE7390"),
               compare  = "grade",
               model = "REM",
               groups = c("3", "1"),
               FDR = 0.05,
               nperm = 50)

# Create a volcano plot for project "GSE11121" with an adjusted p-value threshold of 0.05, positive SE threshold of 1, negative SE threshold of -1, displaying the top 10 results, using the "ggplot" style
volcano_GSE11121 <- Makevolcano(projname  = "GSE11121",
                                padjlevel = 0.05,
                                pSE_thresh = 1,
                                nSE_thresh = -1,
                                ntop = 10,
                                voltype = "ggplot")

# Create a volcano plot for project "GSE7390" with an adjusted p-value threshold of 0.05, positive SE threshold of 1, negative SE threshold of -1, displaying the top 10 results, using the "ggplot" style
volcano_GSE7390 <- Makevolcano(projname  = "GSE7390",
                               padjlevel = 0.05,
                               pSE_thresh = 1,
                               nSE_thresh = -1,
                               ntop = 10,
                               voltype = "ggplot")

# Create a volcano plot for project "GSE25065" with an adjusted p-value threshold of 0.05, positive SE threshold of 1, negative SE threshold of -1, displaying the top 10 results, using the "ggplot" style
volcano_GSE25065 <- Makevolcano(projname  = "GSE25065",
                                padjlevel = 0.05,
                                pSE_thresh = 1,
                                nSE_thresh = -1,
                                ntop = 10,
                                voltype = "ggplot")

# Create a volcano plot for project "GSE25055" with an adjusted p-value threshold of 0.05, positive SE threshold of 1, negative SE threshold of -1, displaying the top 10 results, using the "ggplot" style
volcano_GSE25055 <- Makevolcano(projname  = "GSE25055",
                                padjlevel = 0.05,
                                pSE_thresh = 1,
                                nSE_thresh = -1,
                                ntop = 10,
                                voltype = "ggplot")

# Combination
# Create a volcano plot for the combined data "combination" with an adjusted p-value threshold of 0.05, positive SE threshold of 1, negative SE threshold of -1, displaying the top 10 results, using the "ggplot" style
volcano_combination <- Makevolcano(projname  = "combination",
                                padjlevel = 0.05,
                                pSE_thresh = 1,
                                nSE_thresh = -1,
                                ntop = 10,
                                voltype = "ggplot")

# Create a Venn diagram for the common elements among "GSE25055", "GSE11121", "GSE25065", and "GSE7390" studies, considering the top 10 elements
MakeVenna(common = c("GSE25055", "GSE11121","GSE25065", "GSE7390"), ntop = 10)

# Create a heatmap using the gene expression matrix from "normalized_GSE25055" for project "GSE25055", comparing based on "grade" and grouping samples into "grade1" and "grade3". Display the top 25 elements and do not sort the heatmap.
MakeHeatmap(matrix = normalized_GSE25055,
            projname  = "GSE25055",
            compare  = "grade", 
            groups =list(grade1 = "1", grade3 = "3"),
            ntop = 25,
            sorted = FALSE)

# Create a heatmap using the gene expression matrix from "normalized_GSE11121" for project "GSE11121", comparing based on "grade" and grouping samples into "grade1" and "grade3". Display the top 25 elements and do not sort the heatmap. 
MakeHeatmap(matrix =  normalized_GSE11121,
            projname  = "GSE11121",
            compare  = "grade", 
            groups =list(grade1 = "1", grade3 = "3"),
            ntop = 25,
            sorted = FALSE)

# Create a heatmap using the gene expression matrix from "normalized_GSE25065" for project "GSE25065", comparing based on "grade" and grouping samples into "grade1" and "grade3". Display the top 25 elements and do not sort the heatmap.
MakeHeatmap(matrix = normalized_GSE25065,
            projname  = "GSE25065",
            compare  = "grade", 
            groups =list(grade1 = "1", grade3 = "3"),
            ntop = 25,
            sorted = FALSE)

# Create a heatmap using the gene expression matrix from "normalized_GSE7390" for project "GSE7390", comparing based on "grade" and grouping samples into "grade1" and "grade3". Display the top 25 elements and do not sort the heatmap. 
MakeHeatmap(matrix = normalized_GSE7390,
            projname  = "GSE7390",
            compare  = "grade", 
            groups =list(grade1 = "1", grade3 = "3"),
            ntop = 25,
            sorted = FALSE)

# Combination
# Create a heatmap using the merged gene expression matrix for the "combination" project, comparing based on "grade" and grouping samples into "grade1" and "grade3". Display the top 25 elements and do not sort the heatmap. 
MakeHeatmap(matrix = merge_studies,
            projname  = "combination",
            compare  = "grade", 
            groups =list(grade1 = "1", grade3 = "3"),
            ntop = 25,
            sorted = FALSE)
	
# (Sorted grade) Create a heatmap using the merged gene expression matrix for the "combination" project, comparing based on "grade" and grouping samples into "grade1" and "grade3". Display the top 25 elements and sort the heatmap based on the values of "grade".
MakeHeatmap(matrix = merge_studies,
            projname  = "combination",
            compare  = "grade", 
            groups =list(grade1 = "1", grade3 = "3"),
            ntop = 25,
            sorted = TRUE)

# Perform Gene Set Enrichment Analysis (GSEA) for project "GSE25055" with a false discovery rate (FDR) of 0.05.
GSEAAnalysis(projname  = "GSE25055", FDR = 0.05)
# Perform Gene Set Enrichment Analysis (GSEA) for project "GSE11121" with a false discovery rate (FDR) of 0.05.
GSEAAnalysis(projname  = "GSE11121", FDR = 0.05)
# Perform Gene Set Enrichment Analysis (GSEA) for project "GSE25065" with a false discovery rate (FDR) of 0.05.
GSEAAnalysis(projname  = "GSE25065", FDR = 0.05)
# Perform Gene Set Enrichment Analysis (GSEA) for project "GSE7390" with a false discovery rate (FDR) of 0.05.
GSEAAnalysis(projname  = "GSE7390", FDR = 0.05)

# (Combination) Perform Gene Set Enrichment Analysis (GSEA) for the combined project "combination" with a false discovery rate (FDR) of 0.05.
GSEAAnalysis(projname  = "combination", FDR = 0.05)

# Perform network analysis using the gene expression matrix from "GSE11121" for project "GSE11121", comparing based on "grade" and grouping samples into "grade1" and "grade3".
Makenetworkanalyst(matrix = GSE11121,
                   projname  = "GSE11121", 
                   compare  = "grade", 
                   groups =list(grade1 = "1", grade3 = "3"))

# Perform network analysis using the gene expression matrix from "GSE25055" for project "GSE25055", comparing based on "grade" and grouping samples into "grade1" and "grade3".
Makenetworkanalyst(matrix =GSE25055,
                   projname  = "GSE25055", 
                   compare  = "grade", 
                   groups =list(grade1 = "1", grade3 = "3"))

# Perform network analysis using the gene expression matrix from "GSE25065" for project "GSE25065", comparing based on "grade" and grouping samples into "grade1" and "grade3".
Makenetworkanalyst(matrix = GSE25065,
                   projname  = "GSE25065", 
                   compare  = "grade", 
                   groups =list(grade1 = "1", grade3 = "3"))

# Perform network analysis using the gene expression matrix from "GSE7390" for project "GSE7390", comparing based on "grade" and grouping samples into "grade1" and "grade3". 
Makenetworkanalyst(matrix = GSE7390,
                   projname  = "GSE7390", 
                   compare  = "grade", 
                   groups =list(grade1 = "1", grade3 = "3"))
