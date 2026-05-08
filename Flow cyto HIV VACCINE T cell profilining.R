###############################################################
# HIGH-DIMENSIONAL CyTOF ANALYSIS OF MVA-HIV VACCINE RESPONSES
# CD4+ T-cell longitudinal profiling in vaccinated macaques
#
# Workflow:
# 1. Import FCS files
# 2. Arcsinh transformation
# 3. Data preprocessing + downsampling
# 4. t-SNE visualization
# 5. CellVizR UMAP + clustering
# 6. Quality control
# 7. Marker visualization
# 8. Statistical comparisons
# 9. Volcano plots + heatmaps
# 10. Longitudinal cluster dynamics
# 11. Metacluster generation
###############################################################

############################
# LOAD REQUIRED LIBRARIES
############################

library(flowCore)
library(CellVizR)
library(Rtsne)
library(uwot)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(devtools)

############################
# SET WORKING DIRECTORY
############################

setwd("/Users/Margaux/Desktop/Flow Cytometry Analysis/data/T-cells_CD4p/")

############################
# IMPORT FCS FILES
############################

# Retrieve all FCS files from directory
fcs_files <- list.files(
  pattern = "\\.fcs$",
  full.names = TRUE
)

# Read all FCS files into a flowSet object
flow_set <- read.flowSet(fcs_files)

############################
# ARCSINH TRANSFORMATION
############################

# Retrieve marker names
marker_names <- colnames(flow_set)

# Create arcsinh transformation
transformation <- arcsinhTransform()

# Apply transformation to all markers
transform_list <- transformList(
  from = marker_names,
  tfun = transformation
)

transformed_flowset <- transform_list %on% flow_set

############################
# CONVERT FLOWSET TO DATAFRAME
############################

data <- data.frame()

for (i in seq_along(transformed_flowset)) {
  
  # Extract sample name
  sample_name <- basename(fcs_files[i])
  sample_name <- gsub("\\.fcs$", "", sample_name)
  
  # Extract expression matrix
  sample_data <- data.frame(
    exprs(transformed_flowset[[i]])
  )
  
  # Add sample identifier
  sample_data$sample <- sample_name
  
  # Merge into master dataframe
  data <- rbind(data, sample_data)
}

# Rename columns with marker names
marker_labels <- markernames(transformed_flowset)

colnames(data) <- c(
  "Time",
  marker_labels,
  "sample"
)

############################
# DOWNSAMPLING
############################

# Set seed for reproducibility
set.seed(123)

# Downsample to 25,000 cells
downsampled_data <- data[
  sample(nrow(data), 25000),
]

############################
# REMOVE TECHNICAL CHANNELS
############################

exclude_markers <- c(
  "Time",
  "Cell_length",
  "Dead",
  "Beads",
  "Cells",
  "Cells.1",
  "CD197",
  "FoxP3",
  "CD62L",
  "CD185",
  "CD45RA",
  "CD279",
  "MIP1b",
  "CD195",
  "CD107a",
  "CD95",
  "CD127"
)

# Remove technical markers
downsampled_data <- downsampled_data[
  ,
  !colnames(downsampled_data) %in% exclude_markers
]

# Create matrix for t-SNE
data_no_tech_channels <- downsampled_data[
  ,
  !colnames(downsampled_data) %in% c("sample")
]

############################
# t-SNE DIMENSIONALITY REDUCTION
############################

set.seed(123)

tsne_result <- Rtsne(
  data_no_tech_channels,
  dims = 2,
  num_threads = 4
)

# Store t-SNE coordinates
downsampled_data$tSNE_1 <- tsne_result$Y[,1]
downsampled_data$tSNE_2 <- tsne_result$Y[,2]

############################
# SAMPLE METADATA EXTRACTION
############################

# Extract condition/timepoint from sample names
downsampled_data$condition <- gsub(
  "_.*",
  "",
  downsampled_data$sample
)

############################
# GLOBAL t-SNE VISUALIZATION
############################

ggplot(
  downsampled_data,
  aes(x = tSNE_1, y = tSNE_2)
) +
  geom_point(
    aes(color = sample),
    alpha = 0.5
  ) +
  theme_classic()

############################
# FACETED t-SNE BY SAMPLE
############################

ggplot(
  downsampled_data,
  aes(x = tSNE_1, y = tSNE_2)
) +
  geom_point(
    aes(color = sample),
    alpha = 0.5
  ) +
  facet_wrap(~sample) +
  theme_classic()

############################
# FACETED t-SNE BY CONDITION
############################

ggplot(
  downsampled_data,
  aes(x = tSNE_1, y = tSNE_2)
) +
  geom_point(
    aes(color = sample),
    alpha = 0.5
  ) +
  facet_wrap(~condition) +
  theme_classic()

############################
# MARKER EXPRESSION t-SNE
############################

# Identify markers to visualize
all_columns <- colnames(downsampled_data)

exclude_columns <- c(
  "tSNE_1",
  "tSNE_2",
  "sample",
  "condition"
)

markers_plot <- all_columns[
  !(all_columns %in% exclude_columns)
]

# Loop through all markers
for (marker in markers_plot) {
  
  p <- ggplot(
    downsampled_data,
    aes(x = tSNE_1, y = tSNE_2)
  ) +
    geom_point(
      aes_string(color = marker),
      alpha = 0.5
    ) +
    facet_wrap(~condition) +
    scale_color_viridis_c() +
    ggtitle(paste("Marker:", marker)) +
    theme_classic()
  
  print(p)
}

############################
# IMPORT DATA INTO CELLVIZR
############################

DataCell <- import(
  fcs_files,
  filetype = "fcs",
  transform = "arcsinh",
  exclude.markers = exclude_markers,
  d.method = "uniform",
  parameters.method = list(
    "target.percent" = 0.05
  )
)

############################
# CREATE METADATA
############################

metadata <- data.frame(
  individual = gsub(
    "\\.fcs",
    "",
    basename(fcs_files)
  )
)

# Extract condition
metadata$condition <- sapply(
  strsplit(metadata$individual, "_"),
  function(x) x[1]
)

# Extract timepoint
metadata$timepoint <- sapply(
  strsplit(metadata$individual, "_"),
  function(x) x[1]
)

rownames(metadata) <- metadata$individual

# Assign metadata to DataCell
DataCell <- assignMetadata(
  DataCell,
  metadata = metadata
)

############################
# QUALITY CONTROL
############################

# Marker range QC
QCR <- QCMarkerRanges(fcs_files)
plot(QCR)

# Cell count QC
plotCellCounts(
  DataCell,
  stats = c(
    "min",
    "median",
    "mean",
    "q75",
    "max"
  ),
  samples = NULL,
  sort = TRUE
)

# Small cluster QC
QCS <- QCSmallClusters(
  DataCell,
  th.size = 50,
  plot.device = TRUE
)

# Uniformity QC
QCU <- QCUniformClusters(
  DataCell,
  uniform.test = "both",
  th.pvalue = 0.05,
  th.IQR = 2,
  plot.device = TRUE
)

############################
# UMAP MANIFOLD GENERATION
############################

DataCell <- generateManifold(
  DataCell,
  method = "UMAP"
)

############################
# K-MEANS CLUSTERING
############################

# 200-cluster model
DataCell <- identifyClusters(
  DataCell,
  space = "manifold",
  method = "kmeans",
  centers = 200
)

# 100-cluster validation model
DataCell_100 <- identifyClusters(
  DataCell,
  space = "manifold",
  method = "kmeans",
  centers = 100
)

############################
# CLUSTER VISUALIZATION
############################

plotClustersCounts(
  DataCell,
  clusters = NULL,
  sort = TRUE
)

# Density overlay
plotManifold(
  DataCell,
  markers = "density"
)

# Cluster overlay
plotManifold(
  DataCell,
  markers = "clusters"
) +
  guides(color = "none")

############################
# SAMPLE GROUPS
############################

# Post-prime samples
pp_samples <- c(
  "PPD03_BB078_CD4",
  "PPD03_BB231_CD4",
  "PPD03_BC621_CD4",
  "PPD03_BD619_CD4",
  "PPD03_BD620_CD4",
  
  "PPD08_BB078_CD4",
  "PPD08_BB231_CD4",
  "PPD08_BC621_CD4",
  "PPD08_BD619_CD4",
  "PPD08_BD620_CD4",
  
  "PPD14_BB078_CD4",
  "PPD14_BB231_CD4",
  "PPD14_BC621_CD4",
  "PPD14_BD619_CD4",
  "PPD14_BD620_CD4",
  
  "PPD57_BB078_CD4",
  "PPD57_BB231_CD4",
  "PPD57_BC621_CD4",
  "PPD57_BD619_CD4",
  "PPD57_BD620_CD4"
)

# Post-boost samples
pb_samples <- c(
  "PBD03_BB078_CD4",
  "PBD03_BB231_CD4",
  "PBD03_BC621_CD4",
  "PBD03_BD619_CD4",
  "PBD03_BD620_CD4",
  
  "PBD08_BB078_CD4",
  "PBD08_BB231_CD4",
  "PBD08_BC621_CD4",
  "PBD08_BD619_CD4",
  "PBD08_BD620_CD4",
  
  "PBD14_BB078_CD4",
  "PBD14_BB231_CD4",
  "PBD14_BC621_CD4",
  "PBD14_BD619_CD4",
  "PBD14_BD620_CD4",
  
  "PBD28_BB078_CD4",
  "PBD28_BB231_CD4",
  "PBD28_BC621_CD4",
  "PBD28_BD619_CD4",
  "PBD28_BD620_CD4"
)

############################
# MARKER UMAP VISUALIZATION
############################

for (m in markers_plot) {
  
  p_prime <- plotManifold(
    DataCell,
    markers = m,
    samples = pp_samples
  ) +
    ggtitle(paste(m, "- Post Prime"))
  
  p_boost <- plotManifold(
    DataCell,
    markers = m,
    samples = pb_samples
  ) +
    ggtitle(paste(m, "- Post Boost"))
  
  grid.arrange(
    p_prime,
    p_boost,
    ncol = 2
  )
}

############################
# HEATMAP VISUALIZATION
############################

# Expression heatmap
hm.exp <- plotHmExpressions(DataCell)

grid.arrange(hm.exp)

# Metacluster heatmap
hm.exp_10 <- plotHmExpressions(
  DataCell,
  metaclusters = 10
)

grid.arrange(hm.exp_10)

############################
# MARKER DENSITY FOR SPECIFIC CLUSTER
############################

plotMarkerDensity(
  DataCell,
  clusters = "1"
)

############################
# STATISTICAL COMPARISONS
############################

# Retrieve all unique timepoints
unique_tps <- unique(
  DataCell@metadata$timepoint
)

# Store samples by timepoint
samples_by_timepoint <- list()

for (tp in unique_tps) {
  
  samples_by_timepoint[[tp]] <- selectSamples(
    DataCell,
    timepoint = tp
  )
}

############################
# DEFINE COMPARISONS
############################

comparisons_to_run <- list(
  
  list(
    "test" = "PPD03",
    "ref"  = "BPD19",
    "name" = "PPD03 vs Baseline"
  ),
  
  list(
    "test" = "PBD03",
    "ref"  = "PPD03",
    "name" = "PBD03 vs PPD03"
  ),
  
  list(
    "test" = "PBD28",
    "ref"  = "PPD57",
    "name" = "PBD28 vs PPD57"
  ),
  
  list(
    "test" = "PPD57",
    "ref"  = "PPD03",
    "name" = "PPD57 vs PPD03"
  ),
  
  list(
    "test" = "PBD28",
    "ref"  = "PBD03",
    "name" = "PBD28 vs PBD03"
  )
)

############################
# COMPUTE STATISTICS
############################

DataCell@statistic <- data.frame()

for (comp in comparisons_to_run) {
  
  test_samples <- samples_by_timepoint[[comp$test]]
  ref_samples  <- samples_by_timepoint[[comp$ref]]
  
  DataCell <- computeStatistics(
    DataCell,
    condition = test_samples,
    ref.condition = ref_samples,
    comparison = comp$name,
    test.statistics = "t.test",
    paired = FALSE
  )
}

############################
# VOLCANO PLOTS
############################

for (comp in comparisons_to_run) {
  
  p <- plotVolcano(
    DataCell,
    comparison = comp$name,
    th.pv = 1.3,
    th.fc = 1.5,
    plot.text = TRUE
  ) +
    labs(
      title = paste(
        "Comparison:",
        comp$name
      )
    )
  
  print(p)
}

############################
# STATISTICAL HEATMAP
############################

hm.stats <- plotHmStatistics(
  DataCell,
  clusters = NULL,
  statistics = "pvalue"
)

grid.arrange(hm.stats)

############################
# ORDER TIMEPOINTS
############################

DataCell@metadata$timepoint <- factor(
  DataCell@metadata$timepoint,
  levels = c(
    "BPD19",
    "PPD03",
    "PPD08",
    "PPD14",
    "PPD57",
    "PBD03",
    "PBD08",
    "PBD14",
    "PBD28"
  )
)

############################
# LONGITUDINAL CLUSTER DYNAMICS
############################

clusters_of_interest <- c(
  "63",
  "2",
  "16",
  "47",
  "96",
  "15",
  "77",
  "88",
  "111",
  "21",
  "55",
  "28",
  "176",
  "49"
)

for (clust in clusters_of_interest) {
  
  p <- plotBoxplot(
    DataCell,
    clusters = clust,
    samples = NULL,
    value.y = "percentage",
    observation = "timepoint",
    test.statistics = "t.test"
  )
  
  print(p)
}

############################
# METACLUSTER GENERATION
############################

DataCell <- createMetaclusters(
  DataCell,
  clusters = c(
    "63",
    "2",
    "16",
    "47",
    "96",
    "15",
    "77",
    "88",
    "111",
    "164"
  ),
  metacluster.name = "Metacluster1"
)

############################
# METACLUSTER BOXPLOT
############################

plotBoxplot(
  DataCell,
  clusters = "Metacluster1",
  samples = NULL,
  value.y = "percentage",
  observation = "timepoint",
  test.statistics = "t.test"
)

###############################################################
# END OF ANALYSIS PIPELINE
###############################################################