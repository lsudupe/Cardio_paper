### SCRIPT: Spatial object creation, add metadata, plots

## 12.08.24 Laura Sudupe , git @lsudupe

library(stats)
library(Seurat)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(tidyverse)


# Define dir
data.dir <- "./GSE265828_RAW/"

dpi3 <- Load10X_Spatial(data.dir = paste0(data.dir, "./dpi3/"),
                                filename = "filtered_feature_bc_matrix.h5",
                                assay = "Spatial",
                                slice = "dpi3",
                                filter.matrix = TRUE)

dpi5 <- Load10X_Spatial(data.dir = paste0(data.dir,"./dpi5m/"),
                                  filename = "filtered_feature_bc_matrix.h5",
                                  assay = "Spatial",
                                  slice = "dpi5_male",
                                  filter.matrix = TRUE)

dpi3 <- SCTransform(dpi3, assay = "Spatial", verbose = FALSE)
dpi5 <- SCTransform(dpi5, assay = "Spatial", verbose = FALSE)


# Filter 
dpi3 <- subset(x = dpi3, 
               subset= (nCount_Spatial >= 500) & 
                 (nFeature_Spatial >= 250))

dpi5 <- subset(x = dpi5, 
                 subset= (nCount_Spatial >= 500) & 
                        (nFeature_Spatial >= 250))

# Read meta
meta_dpi3 <- read.csv("./meta_dpi3_process.csv", row.names = 1, stringsAsFactors = FALSE)
meta_dpi5 <- read.csv("./meta_dpi5_process.csv", row.names = 1, stringsAsFactors = FALSE)


# Add meta
dpi3 <- AddMetaData(dpi3, meta_dpi3)
dpi5 <- AddMetaData(dpi5, meta_dpi5)

# Colors
area_colors <- c("IZ" = "darkblue", "BZ1" = "skyblue", "BZ2" = "darkgreen", "RZ" = "yellow")
sample_colors <- c("dpi3" = "gray", "dpi5_female" = "pink", "dpi5_male" = "lightblue")

# Spatial areas plot
SpatialDimPlot(dpi3, group.by ="new_area", cols = area_colors, alpha = 0.6)
SpatialDimPlot(dpi5, group.by ="new_area", cols = area_colors, alpha = 0.6)


## Pairwise comparation
dpi3_meta <- dpi3@meta.data
colnames(dpi3_meta)[colnames(dpi3_meta) == "Ratio.D2.D1"] <- "Ratio D2/D1"
dpi5_meta <- dpi5@meta.data
colnames(dpi5_meta)[colnames(dpi5_meta) == "Ratio.D2.D1"] <- "Ratio D2/D1"

# Function to process metadata aand generate dataframes with results and annotations
process_metadata <- function(metadata, sample_name) {
  # Filter and select the columns
  new_df <- metadata %>%
    dplyr::select(
      new_area,
      CC9,
      CC5,
      CC8,
      CC11,
      `Ratio D2/D1`,
      Dynamics.1,
      Dynamics.2
    )
  
  # List of column to analyze
  columns_to_analyze <- c("CC9", "CC5", "CC8", "CC11", "Ratio D2/D1", "Dynamics.1", "Dynamics.2")
  
  # Inizialaize dataframe to store data
  results <- data.frame(variable = character(), group1 = character(), group2 = character(), p_value = numeric(), stringsAsFactors = FALSE)
  
  # Wilcoxon test and p-valores
  for (col in columns_to_analyze) {
    areas <- unique(new_df$new_area)
    for (i in 1:(length(areas) - 1)) {
      for (j in (i + 1):length(areas)) {
        group1_data <- new_df[new_df$new_area == areas[i], col, drop = TRUE]
        group2_data <- new_df[new_df$new_area == areas[j], col, drop = TRUE]
        test_result <- wilcox.test(group1_data, group2_data, exact = FALSE, correct = FALSE)
        results <- rbind(results, data.frame(variable = col, group1 = areas[i], group2 = areas[j], p_value = test_result$p.value))
      }
    }
  }
  
  results$p_adjusted <- p.adjust(results$p_value, method = "BH")
  results$Significant <- ifelse(results$p_adjusted < 0.05, 1, 0)
  
  results <- results %>%
    mutate(comparison = paste(group1, "vs", group2, sep = "_"))
  
  annotation_df <- results %>%
    select(variable, comparison, Significant) %>%
    pivot_wider(names_from = comparison, values_from = Significant, values_fill = list(Significant = 0))
  
  df_parte1 <- annotation_df %>% slice(1:4)
  df_parte2 <- annotation_df %>% slice(5:nrow(annotation_df))
  
  # Create table with significance 
  star_table <- results %>%
    mutate(stars = cut(
      p_adjusted, 
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
      labels = c("***", "**", "*", ""), 
      right = FALSE
    )) %>%
    select(variable, comparison, stars) %>%
    pivot_wider(names_from = comparison, values_from = stars, values_fill = list(stars = ""))
  
  star_table_part1 <- star_table %>% slice(1:4)
  star_table_part2 <- star_table %>% slice(5:nrow(star_table))
  
  # Generate dataframes of adjusted p-values
  results_df <- results %>%
    select(variable, comparison, p_adjusted) %>%
    pivot_wider(names_from = comparison, values_from = p_adjusted)
  
  pvalue_df1 <- results_df %>% slice(1:4)
  pvalue_df2 <- results_df %>% slice(5:nrow(results_df))
  
  list(df1 = df_parte1, df2 = df_parte2, star_table1 = star_table_part1, star_table2 = star_table_part2, pvalue_df1 = pvalue_df1, pvalue_df2 = pvalue_df2)
}

# Lists of metadata 
metadatas <- list(
  dpi3 = dpi3_meta,
  dpi5 = dpi5_meta
)

# List of results
results_df1_list <- list()
results_df2_list <- list()
star_df1 <- list()
star_df2 <- list()
pvalue_df1 <- list()
pvalue_df2 <- list()

for (sample_name in names(metadatas)) {
  processed <- process_metadata(metadatas[[sample_name]], sample_name)
  results_df1_list[[sample_name]] <- processed$df1
  results_df2_list[[sample_name]] <- processed$df2
  star_df1[[sample_name]] <- processed$star_table1
  star_df2[[sample_name]] <- processed$star_table2
  pvalue_df1 <- processed$pvalue_df1
  pvalue_df2 <- processed$pvalue_df2
}

## Median
# Function to prepare metadata and calculate median values
calculate_mean_values <- function(metadata, sample_name) {
  # Select relevant columns
  new_df <- metadata %>%
    dplyr::select(
      new_area,
      CC9,
      CC5,
      CC8,
      CC11,
      `Ratio D2/D1`,
      Dynamics.1,
      Dynamics.2
    )
  
  # Caluclate the median per area and signature
  mean_values_df <- new_df %>%
    group_by(new_area) %>%
    summarise(across(everything(), list(mean = ~mean(.x, na.rm = TRUE)))) %>%
    pivot_longer(-new_area, names_to = "signature", values_to = "mean_value") %>%
    mutate(signature = str_replace(signature, "_mean", ""),
           sample = sample_name)
  
  return(mean_values_df)
}

# Metadata list to process
metadatas <- list(
  dpi3 = dpi3_meta,
  dpi5 = dpi5_meta
)

# LIst to store dataframes
mean_values_list <- list()

# Iterate metadata
for (sample_name in names(metadatas)) {
  mean_values_list[[sample_name]] <- calculate_mean_values(metadatas[[sample_name]], sample_name)
}

# Combine all dataframes inone
mean_values_df <- bind_rows(mean_values_list)

# Transform dataframe to wide
df_wide <- pivot_wider(mean_values_df, 
                       names_from = c("sample", "new_area"), 
                       names_sep = "_", 
                       values_from = "mean_value", 
                       values_fill = list(mean_value = NA))


# Prepara rows and columns for the matrix
signatures <- c("CC9", "CC5", "CC8", "CC11", "Ratio D2/D1", "Dynamics.1", "Dynamics.2")
rownames(df_wide) <- signatures

# Separate df_wide in two dataframes
dpi3_meta_df <- select(df_wide, contains("dpi3"))
dpi5_meta_df <- select(df_wide, contains("dpi5"))


colnames(dpi3_meta_df) <- gsub("dpi3_", "", colnames(dpi3_meta_df))
colnames(dpi5_meta_df) <- gsub("dpi5_", "", colnames(dpi5_meta_df))

rownames(dpi3_meta_df) <- signatures
rownames(dpi5_meta_df) <- signatures

# Create daratframes
median_list <- list(
  dpi3_meta_df = dpi3_meta_df,
  dpi5_meta_df = dpi5_meta_df
)


## Heatmaps
area_colors <- c("IZ" = "red", "BZ1" = "skyblue", "BZ2" = "darkgreen", "RZ" = "yellow")
sample_colors <- c("dpi3" = "gray", "dpi5" = "lightblue")

# dpi3
dpi3_median <- median_list[["dpi3_meta_df"]]
dpi3_median <- t(dpi3_median)

row.names(dpi3_median) <- c("BZ1", "BZ2", "IZ", "RZ")
dpi3_median <- dpi3_median[c("RZ", "BZ1", "BZ2", "IZ"), ]

dpi3_median1 <- dpi3_median[, 1:4]
dpi3_median2 <- dpi3_median[, 5:ncol(dpi3_median)]

# dpi3 df1
dpi3_df1_s <- results_df1_list[["dpi3"]]

top_annotation <- HeatmapAnnotation(df=data.frame(Area = c("RZ","BZ1","BZ2","IZ"), 
                                                  Sample = c("dpi3")),
                                    col=list(Area = area_colors,
                                             Sample = sample_colors))
right_annotation <- rowAnnotation(df=data.frame(BZ1vsRZ=c("significant","significant","significant","significant"), 
                                                BZ1vsIZ=c("significant","significant","significant","significant"),
                                                BZ2vsRZ=c("significant","significant","significant","significant"), 
                                                BZ2vsIZ=c("significant","significant","significant","significant"),
                                                RZvsIZ=c("significant","significant","significant","significant"), 
                                                BZ1vsBZ2=c("no-significant","no-significant","no-significant","no-significant")),
                                  col = list(BZ1vsRZ=c("significant"="black", "no-significant"="white"),
                                             BZ1vsIZ=c("significant"="black", "no-significant"="white"),
                                             BZ2vsRZ=c("significant"="black", "no-significant"="white"),
                                             BZ2vsIZ=c("significant"="black", "no-significant"="white"),
                                             RZvsIZ=c("significant"="black", "no-significant"="white"),
                                             BZ1vsBZ2=c("significant"="black", "no-significant"="white")))
dpi3_median1 <- t(dpi3_median1)

#pdf("./dp3_CC9_CC5_CC8_CC11.pdf", width = 8)
ComplexHeatmap::Heatmap(dpi3_median1,
                        top_annotation = top_annotation, 
                        right_annotation = right_annotation,
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE)
#dev.off()

# dpi3 df2
dpi3_df2_s <- results_df2_list[["dpi3"]]
top_annotation <- HeatmapAnnotation(df=data.frame(Area = c("RZ","BZ1","BZ2","IZ"), 
                                                  Sample = c("dpi3")),
                                    col=list(Area = area_colors,
                                             Sample = sample_colors))
right_annotation <- rowAnnotation(df=data.frame(BZ1vsRZ=c("no-significant","significant","significant"), 
                                                BZ1vsIZ=c("significant","no-significant","significant"),
                                                BZ2vsRZ=c("no-significant","significant","significant"), 
                                                BZ2vsIZ=c("significant","significant","significant"),
                                                RZvsIZ=c("significant","significant","significant"), 
                                                BZ1vsBZ2=c("no-significant","no-significant","no-significant")),
                                  col = list(BZ1vsRZ=c("significant"="black", "no-significant"="white"),
                                             BZ1vsIZ=c("significant"="black", "no-significant"="white"),
                                             BZ2vsRZ=c("significant"="black", "no-significant"="white"),
                                             BZ2vsIZ=c("significant"="black", "no-significant"="white"),
                                             RZvsIZ=c("significant"="black", "no-significant"="white"),
                                             BZ1vsBZ2=c("significant"="black", "no-significant"="white")))
dpi3_median2 <- t(dpi3_median2)

#pdf("./dp3_D1_D2_ratio.pdf", width = 8)
ComplexHeatmap::Heatmap(dpi3_median2,
                        top_annotation = top_annotation, 
                        right_annotation = right_annotation,
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE)
#dev.off()

# dpi5
dpi5_median <- median_list[["dpi5_meta_df"]]
dpi5_median <- t(dpi5_median)

row.names(dpi5_median) <- c("BZ1", "BZ2", "IZ", "RZ")
dpi5_median <- dpi5_median[c("RZ", "BZ1", "BZ2", "IZ"), ]

dpi5_median1 <- dpi5_median[, 1:4]
dpi5_median2 <- dpi5_median[, 5:ncol(dpi5_median)]

# dpi5 df1
dpi5_df1_s <- results_df1_list[["dpi5"]]
top_annotation <- HeatmapAnnotation(df=data.frame(Area = c("RZ","BZ1","BZ2","IZ"), 
                                                  Sample = c("dpi5")),
                                    col=list(Area = area_colors,
                                             Sample = sample_colors))
right_annotation <- rowAnnotation(df=data.frame(BZ1vsRZ=c("significant","significant","significant","significant"), 
                                                BZ1vsIZ=c("significant","significant","significant","significant"),
                                                BZ2vsRZ=c("significant","significant","significant","significant"), 
                                                BZ2vsIZ=c("significant","significant","significant","significant"),
                                                RZvsIZ=c("significant","significant","significant","significant"), 
                                                BZ1vsBZ2=c("no-significant","no-significant","no-significant","no-significant")),
                                  col = list(BZ1vsRZ=c("significant"="black", "no-significant"="white"),
                                             BZ1vsIZ=c("significant"="black", "no-significant"="white"),
                                             BZ2vsRZ=c("significant"="black", "no-significant"="white"),
                                             BZ2vsIZ=c("significant"="black", "no-significant"="white"),
                                             RZvsIZ=c("significant"="black", "no-significant"="white"),
                                             BZ1vsBZ2=c("significant"="black", "no-significant"="white")))
dpi5_median1 <- t(dpi5_median1)

#pdf("./dpi5_CC9_CC5_CC8_CC11.pdf", width = 8)
ComplexHeatmap::Heatmap(dpi5_median1,
                        top_annotation = top_annotation, 
                        right_annotation = right_annotation,
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE)
#dev.off()

# dpi5 df2
dpi5_df2_s <- results_df2_list[["dpi5"]]
top_annotation <- HeatmapAnnotation(df=data.frame(Area = c("RZ","BZ1","BZ2","IZ"), 
                                                  Sample = c("dpi5")),
                                    col=list(Area = area_colors,
                                             Sample = sample_colors))
right_annotation <- rowAnnotation(df=data.frame(BZ1vsRZ=c("no-significant","significant","significant"), 
                                                BZ1vsIZ=c("significant","no-significant","significant"),
                                                BZ2vsRZ=c("no-significant","significant","significant"), 
                                                BZ2vsIZ=c("significant","significant","significant"),
                                                RZvsIZ=c("significant","significant","significant"), 
                                                BZ1vsBZ2=c("no-significant","no-significant","no-significant")),
                                  col = list(BZ1vsRZ=c("significant"="black", "no-significant"="white"),
                                             BZ1vsIZ=c("significant"="black", "no-significant"="white"),
                                             BZ2vsRZ=c("significant"="black", "no-significant"="white"),
                                             BZ2vsIZ=c("significant"="black", "no-significant"="white"),
                                             RZvsIZ=c("significant"="black", "no-significant"="white"),
                                             BZ1vsBZ2=c("significant"="black", "no-significant"="white")))
dpi5_median2 <- t(dpi5_median2)

#pdf("./dpi5_GC1_GC2_ratio.pdf", width = 8)
ComplexHeatmap::Heatmap(dpi5_median2,
                        top_annotation = top_annotation, 
                        right_annotation = right_annotation,
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE)
#dev.off()

# Spatial ratio plot
bl <- colorRampPalette(c("navy","royalblue","lightskyblue"))(200)                      
re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)

l <- c(min(dpi3@meta.data[["Ratio.D2.D1"]]), max(dpi3@meta.data[["Ratio.D2.D1"]]))
p1 <- SpatialFeaturePlot(dpi3, features = c("Ratio.D2.D1"), combine = FALSE, ncol = 2)
fix.p1 <- scale_fill_gradientn(colours=c(bl,"white", re),
                               breaks=l,
                               labels=c("Min","Max"),
                               na.value = "grey98",
                               limits = l)
p2 <- lapply(p1, function (x) x + fix.p1)
print(CombinePlots(p2))

l <- c(min(dpi5@meta.data[["Ratio.D2.D1"]]), max(dpi5@meta.data[["Ratio.D2.D1"]]))
p1 <- SpatialFeaturePlot(dpi5, features = c("Ratio.D2.D1"), combine = FALSE, ncol = 2)
fix.p1 <- scale_fill_gradientn(colours=c(bl,"white", re),
                               breaks=l,
                               labels=c("Min","Max"),
                               na.value = "grey98",
                               limits = l)
p2 <- lapply(p1, function (x) x + fix.p1)
print(CombinePlots(p2))


