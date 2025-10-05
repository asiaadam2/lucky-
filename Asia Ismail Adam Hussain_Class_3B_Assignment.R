setwd()
dir.create("Raw_Data")
dir.create("Scripts")
dir.create("Plots")
dir.create("Results")


if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

BiocManager::install(c("GEOquery","affy","arrayQualityMetrics"))

install.packages("dplyr")

library(GEOquery) 
library(affy) 
library(arrayQualityMetrics)
library(dplyr) 

gse_data <- getGEO("GSE3292", GSEMatrix = TRUE)

length(gse_data)
eset <- gse_data[[1]]

options(download.file.method = "libcurl")
expression_data <- exprs(gse_data$GSE3292_series_matrix.txt.gz)


feature_data <-  fData(gse_data$GSE3292_series_matrix.txt.gz)

phenotype_data <-  pData(gse_data$GSE3292_series_matrix.txt.gz)

sum(is.na(phenotype_data$source_name_ch1)) 


getGEOSuppFiles("GSE3292", baseDir = "Raw_Data", makeDirectory = TRUE)


untar("Raw_Data/GSE3292_RAW.tar", exdir = "Raw_Data/CEL_Files")


raw_data <- ReadAffy(celfile.path = "Raw_Data/CEL_Files")

raw_data 


arrayQualityMetrics(expressionset = raw_data,
                    outdir = "Results/QC_Raw_Data",
                    force = TRUE,
                    do.logtransform = TRUE)


list.file("Results/QC_Raw_Data")
#based on the arrayQualityMetrics reports, 2 arrays were flagged as outliner.


normalized_data <- rma(raw_data)
  
arrayQualityMetrics(expressionset = normalized_data[, c(17,26)],
                    outdir = "Results/QC_Normalized_Data",
                    force = TRUE)

arrayQualityMetrics(expressionset = normalized_data,
                    outdir = "Results/QC_Normalized_Data",
                    force = TRUE)

processed_data <- as.data.frame(exprs(normalized_data))
dim(processed_data)

install.packages("matrixStats")
library(matrixStats)
row_median <- rowMedians(as.matrix(processed_data))

hist(row_median,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")
threshold <- 3.5 
abline(v = threshold, col = "black", lwd = 2)


indx <- row_median > threshold 
filtered_data <- processed_data[indx, ] 
dim(filtered_data)

colnames(filtered_data) <- rownames(phenotype_data)
processed_data <- filtered_data 


class(phenotype_data$source_name_ch1)

groups <- factor(phenotype_data$source_name_ch1,
                 levels = c(" normal mucosa", "tumor"),
                 label = c("normal", "cancer"))

class(groups)
levels(groups)




boxplot(filtered_data, main = "boxplot of normalized and filtered data",
        xlab = "samples",
        ylab = "expression values",
        col = "lightblue",
        outline = FALSE,
        las = 2)
dev.off()
setwd()
dir.create("Plots", showWarnings = FALSE)

png("Plots/boxplot_expression.png", width = 1200, height = 800)
png("Results/")

pdf("Plots/boxplot_expression.pdf", width = 12, height = 8)

boxplot(processed_data, main = "expression boxplot",
        col = "lightblue",
        las = 2)
dev.off()

pdf("Plots/pca_plot.pdf", width = 12, height = 8)
pca <- prcomp(t(filtered_data), scale. = TRUE)
plot(pca$x[,1], pca$x[,2],
     xlab = "pc1",
     ylab = "pc2",
     main = "pca of normalized & filtered data",
     pch = 19,
     col = "blue")
grid()
dev.off()
