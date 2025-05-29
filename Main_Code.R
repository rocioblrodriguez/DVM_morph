# packages:
# "tidyr", "dplyr", "readr", "stringr", "data.table", "ggplot2", "readxl", 
# "ggpubr", "gridExtra", "RColorBrewer", "colorspace", "FactoMineR", "factoextra",
# "gginnards", "cellWise", "corrplot", "vegan", "morphr", "purrr", "imager",
# "ggrepel", "cowplot", "Nmisc", "bestNormalize"

#### 1. Load Data ####

# metadata
download.file(
  url      = "https://www.dropbox.com/scl/fi/crlvzq8w6gmy0g0okflq5/Final_Samples_Gradients.csv?rlkey=cjfl9n3xgxwm5x2tyk94gpm40&st=1k69598w&dl=1&raw=1",
  destfile = "/tmp/metadata.csv",
  method   = "auto"
)

metadata <- readr::read_csv("/tmp/metadata.csv") |>
  janitor::clean_names() |>
  dplyr::mutate(
    cruise = dplyr::case_when(
      cruise == "New Horizon" ~ "nh1208",
      cruise == "Oceanus" ~ "oc473",
      TRUE ~ cruise
    ),
    cast = paste0("m", cast)
  ) |>
pointblank::expect_col_vals_not_null(vars(cruise, cast))

# env_data
download.file(
  url      = "https://www.dropbox.com/scl/fi/eysmmohpo9fkrya9uleuo/Environmental_data.csv?rlkey=en2jo0uozdm73ta20ko7mee7m&st=946l6s0y&dl=1&raw=1",
  destfile = "/tmp/env_data.csv",
  method   = "auto"
)

env_data <- readr::read_csv("/tmp/env_data.csv") |>
  janitor::clean_names() |>
  dplyr::mutate(cruise_moc_net = paste(cruise, tow, paste0("n", net), sep = "_")) |>
  dplyr::select(-cruise, -net)

# note that `nh1208_m19_n9` and `oc473_m14_n4` are present in the environmental but not ecotaxa data!

# data (corrected dates)
download.file(
  url      = "https://www.dropbox.com/scl/fi/anvu2t0zvn9f6skqm6tgy/ecotaxa_export_5421_20250307_2215_rm.tsv?rlkey=qloai0j5lyd79ky05a17tabjx&st=275m9wk2&dl=1&raw=1",
  destfile = "/tmp/data.tsv",
  method   = "auto"
)

data <- ecotaxaLoadR::load_eco_taxa(file_path = "/tmp/data.tsv") 

data_full <- dplyr::left_join(
  x  = data,
  y  = env_data,
  by = c("cruise_moc_net")
) |>
  dplyr::left_join(
    y  = metadata[, c("cruise", "cast", "d_n")],
    by = c(
      "cruise" = "cruise",
      "moc"    = "cast"
    )
  ) |>
  pointblank::expect_row_count_match(count = nrow(data))


# standardize selected environmental variables

env_vars <- c(
  "avg_t",
  "avg_sal",
  "avg_o2",
  "avg_fluor"
)

data_full[, env_vars] <- vegan::decostand(
  x      = data_full[, env_vars],
  method = "standardize"
)

# apply Yeo-Johnson transformation with standardization

morphological_vars <- c(
  "object_area", "object_mean", "object_stddev", "object_mode", "object_min", "object_max", "object_x", "object_y", "object_xm", "object_ym", "object_perim", "object_bx", "object_by", "object_width", "object_height", "object_major", "object_minor", "object_angle", "object_circ", "object_feret", "object_intden", "object_median", "object_skew", "object_kurt", "object_area", "object_xstart", "object_ystart", "object_area_exc", "object_fractal", "object_skelarea", "object_slope", "object_histcum1", "object_histcum2", "object_histcum3", "object_xmg5", "object_ymg5", "object_nb1", "object_nb2", "object_nb3", "object_compentropy", "object_compmean", "object_compslope", "object_compm1", "object_compm2", "object_compm3", "object_symetrieh", "object_symetriev", "object_symetriehc", "object_symetrievc", "object_convperim", "object_convarea", "object_fcons", "object_thickr", "object_tag", "object_esd", "object_elongation", "object_range", "object_meanpos", "object_centroids", "object_cv", "object_sr", "object_perimareaexc", "object_feretareaexc", "object_perimferet", "object_perimmajor", "object_circex", "object_cdexc"
  )

YJ_transformed <- apply(
  data_full[, c(morphological_vars)],
  2,
  function(x) bestNormalize::yeojohnson(x, standardize = TRUE)$x.t
) |>
  as.data.frame()

data_full <- data_full |>
  dplyr::select(-tidyselect::all_of(morphological_vars)) |>
  dplyr::bind_cols(YJ_transformed)

# dataset is ready for analysis

data_full <- data_full |>
  dplyr::mutate(
    region = case_when(
      cruise == "nh1208" & station %in% c(7, 11) ~ "PSAE",
      cruise == "nh1208" & station %in% c(15, 21, 23) ~ "NPPF",
      cruise == "nh1208" & station %in% c(27, 34) ~ "NPTG",
      cruise == "oc473"  & station %in% c(5, 8) ~ "NASW",
      cruise == "oc473"  & station %in% c(13, 21) ~ "GFST",
      cruise == "oc473"  & station %in% c(26, 31) ~ "NADR",
      TRUE ~ NA_character_
    )
  )

# taxa:
# "Actinopterygii", "Annelida", "Bryozoa", "Cephalochordata", "Chaetognatha",
# "Cnidaria<Metazoa", "Hydrozoa", "Siphonophorae", "Amphipoda", "Calanoida",
# "Cyclopoida", "Harpacticoida", "Mormonilla", "Decapoda", "Euphausiacea",
# "calyptopsis<Euphausiacea", "Ostracoda", "Echinodermata", "Harosa",
# "Foraminifera", "Heteropoda", "Mollusca", "Cephalopoda",
# "Gastropoda<Mollusca", "Cavoliniidae", "Creseidae", "Gymnosomata",
# "Limacinidae", "Oikopleura", "Pseudothecosomata", "Doliolida", "Salpida"

filtered <- data_full |>
  dplyr::filter(
    grepl(
      pattern     = "nh1208",
      x           = cruise_moc_net,
      ignore.case = TRUE
    ),
    grepl(
      pattern     = "Calanoida",
      x           = object_annotation_category,
      ignore.case = TRUE
    ),
    object_annotation_status == "validated"
  )

## image processing ####

# images_directory <- "imgs/"
# Cop_data <- Cop_data %>% mutate(img_path = str_c(images_directory, id, ".jpg"))

morphological_drops <- c(
  "object_ystart",
  "object_bx",
  "object_by",
  "object_xstart",
  "object_angle"
)

prune_morpho_columns <- function(
  df,
  morph_cols,
  drop_cols
) {

  morpho <- df |>
    dplyr::select(tidyselect::all_of(morph_cols))

  na_cols <- names(morpho)[colSums(is.na(morpho)) > 0]

  if (length(na_cols) > 0) {
    message(
      "The following columns were dropped due to NA values: ",
      paste(na_cols, collapse = ", ")
    )
  } else {
    message("No columns were dropped due to NA values.")
  }

  morpho <- morpho |>
    dplyr::select(dplyr::where(~ !any(is.na(.)))) |>
    dplyr::select(-tidyselect::all_of(drop_cols))

  return(morpho)

}

pruned <- prune_morpho_columns(
  filtered,
  morphological_vars,
  morphological_drops
  )

# compute correlation matrix
cor_matrix <- cor(
  x   = pruned,
  use = "pairwise.complete.obs"
  )  

# find highly correlated variables
highly_correlated_morph_cols <- caret::findCorrelation(
  x      = cor_matrix,
  cutoff = 0.99,
  names  = TRUE
  )  

# weighting
highly_correlated_morph_data <- filtered[, highly_correlated_morph_cols]

highly_correlated_cor_matrix <- cor(
  x = highly_correlated_morph_data,
  use = "pairwise.complete.obs"
)

ggcorrplot::ggcorrplot(
  corr = highly_correlated_cor_matrix,
  lab  = TRUE,
  type = "upper"
)

abundance_sqrt    <- sqrt(filtered$abundance)
weighted_data     <- sweep(highly_correlated_morph_data, 1, abundance_sqrt, `*`)
weighted_data     <- weighted_data[, colSums(is.na(weighted_data)) < nrow(weighted_data)]
weighted_env_cols <- c(env_vars, "object_depth_max", "dic_umol_kg", "d_n") # assume dic_umol_kg ~ CO2
weighted_data     <- cbind(filtered[, weighted_env_cols], weighted_data)

#### 8. PCA and K-Means Clustering Analysis ####
# Perform PCA on the weighted morphological data
#res_pca <- PCA(weighted_data, ncp = 10, graph = FALSE)  # ncp = 10 limits to the first 10 PCs

res_pca <- FactoMineR::PCA(
  X          = weighted_data,
  scale.unit = TRUE,
  quanti.sup = 1:6,
  quali.sup  = 7,
  graph      = FALSE,
  ncp        = 20
)

pamk_best <- fpc::pamk(
  data      = highly_correlated_morph_data,
  criterion = "multiasw",
  usepam    = FALSE
)

source("find_optimal_k.R")

result <- find_optimal_k(
  highly_correlated_morph_data,
  k_range = 2:15
)

print(result$plot)


# Define calanoida clusters
#set.seed(123)
#write_csv(weighted_data, "weighted_data.csv")
#sil<-fviz_nbclust(weighted_data, kmeans,  method="silhouette")
#sil + geom_vline(xintercept = 6, linetype="dashed", 
 #               color = "red")+
#theme(text = element_text(size=15))+
#labs( tag = "A")
# Load necessary library
library(ggdendro)

# Compute distance matrix using K-means cluster centers
#dist_centroids <- dist(kmeans_result$centers, method = "euclidean")

# Perform hierarchical clustering on centroids
#hclust_centroids <- hclust(dist_centroids, method = "ward.D2")

# Plot dendrogram for cluster centroids
#plot(hclust_centroids, main = "Dendrogram of K-means Cluster Centroids", 
 #    xlab = "Cluster", sub = "", cex = 0.8)
# Sample 10,000 rows correctly (rows only, keeping all columns)
#set.seed(123)  # Ensure reproducibility
sample_data <- morphological_data[sample(nrow(morphological_data), 20000), ]

# Compute the optimal number of clusters using the Elbow Method
fviz_nbclust(sample_data, kmeans, method = "wss", k.max = 15)

# Perform hierarchical clustering
# distance_matrix <- dist(sample_data)  # Compute Euclidean distance matrix
distance_matrix <- parallelDist::parDist(as.matrix(sample_data), threads = parallel::detectCores() - 1)

hclust_result <- hclust(distance_matrix, method = "ward.D2")  # Ward's method for clustering
# Assign cluster labels to k-means results
Cop_data$Cluster <- factor(kmeans_result$cluster, labels = labels)

# Convert hierarchical clustering to a dendrogram
dend <- as.dendrogram(hclust_result)

# Assign k-means cluster labels to dendrogram leaves
cluster_assignment <- kmeans_result$cluster
names(cluster_assignment) <- rownames(morphological_data)  # Ensure labels match

# Reorder the dendrogram based on k-means clusters
dend <- reorder(dend, wts = cluster_assignment)

# Assign colors to dendrogram branches based on k-means clusters
dend <- color_branches(dend, k = 10, col = cols[as.character(labels[kmeans_result$cluster])])

# Plot the dendrogram with colored branches
plot(dend, main = "Hierarchical Dendrogram Colored by K-means Clusters",
     ylab = "Height", cex = 1.2)

# Add a legend to indicate cluster colors
legend("topright", legend = names(cols), fill = cols, border = "black", cex = 1.2, title = "Clusters")

# Perform K-means clustering on the weighted morphological data
set.seed(123)
kmeans_result <- kmeans(morphological_data, centers = 6, nstart = 25, iter.max = 100)
# Create a vector with the labels in the desired order
labels <- c("ST", "SD","MD","MT", "LT","LD")

# Match the order of labels to the clusters in your result
# Assuming your clusters originally are ordered as 1: ST, 2: MT, 3: LT, 4: SD, 5: MD, 6: LD
# But you want them as 1: LD, 2: ST, 3: MT, 4: SD, 5: MD, 6: LT
# We rearrange labels to match this new order
ordered_labels <- labels[c(1,2,3,4,5,6)]

# Use factor to assign these labels to the cluster ids
Cop_data$Cluster <- kmeans_result$cluster

Cop_data$Cluster <- factor(kmeans_result$cluster, labels = ordered_labels)

# **1. Calculate the mean transparency value for each cluster**
cluster_transparency <- Cop_data %>%
  group_by(Cluster) %>%
  summarise(mean_transparency = mean(mean, na.rm = TRUE)) %>%
  arrange(desc(mean_transparency))  # Sort from most transparent to darkest

# **2. Extract the ordered cluster names**
ordered_clusters <- cluster_transparency$Cluster

# **3. Reorder `Cluster` factor levels based on transparency ranking**
Cop_data$Cluster <- factor(Cop_data$Cluster, levels = ordered_clusters)

# **4. Define the color mapping based on transparency ranking**
cols <- c("ST" = "#C6DBEF",  # Most Transparent
          "MT" = "#9ECAE1",  # Mid Transparent
          "LT" = "#6BAED6",  # Less Transparent
          "SD" = "#4292C6",  # A Bit Dark
          "MD" = "#08519C",  # Darker
          "LD" = "#08306B")  # Darkest

# **Ensure colors are mapped in the new order**
cols <- cols[ordered_clusters]

# **5. PCA Variable Extraction & Processing**
pca.vars <- rbind(res_pca$var$coord, res_pca$quanti.sup$coord, res_pca$quali.sup$coord) %>% as.data.frame()

# Filter important variables
pca.vars1 <- pca.vars[c("mean", "major", "stddev", "Temp", "Sal", "Fluor", "O2", "depth_max","CO2", "D", "N"), ]

# Rename the row names for PCA variables
rownames(pca.vars1) <- c("TRANSPARENT", "LARGE", "DARK", "Temp", "Sal", "Fluor", "O2", "Depth","CO2","D", "N" )

# Create the "SMALL" row as the inverse of "LARGE"
small_row <- -pca.vars1["LARGE", ]
pca.vars1 <- rbind(pca.vars1, SMALL = small_row)

# **6. Define Region Labels for Plot Titles**
Region_labs <- c(
  NPPF = "North Pacific Polar Front (NPPF)",
  NPTG = "North Pacific Tropical Gyre (NPTG)",
  PSAE = "Eastern Pacific Subarctic Gyres (PSAE)",
  GFST = "Gulf Stream (GFST)",
  NASW = "Northwest Atlantic Subtropical Gyre (NASW)",
  NADR = "North Atlantic Drift (NADR)"
)

cluster_order <- c("ST","MT","LT","SD", "MD", "LD"  )

# Convert Cluster column to a factor and specify the levels order
#Cop_data$Cluster <- factor(Cop_data$Cluster, levels = cluster_order)
#Cop_data$Cluster <- as.factor(Cop_data$Cluster)

# **7. Create an Empty List to Store PCA Plots**
pca_plots_list <- list()

# Get a vector of unique regions
unique_regions <- unique(Cop_data$Region)

# **8. Loop Through Each Region to Generate PCA Plots**
for (region in unique_regions) {
  
  # Filter data for the current region
  region_data <- Cop_data[Cop_data$Region == region, ]
  
  # Merge PCA scores with region_data
  region_data <- region_data %>% 
    mutate(Dim.1 = res_pca$ind$coord[rownames(region_data), 1],
           Dim.2 = res_pca$ind$coord[rownames(region_data), 2])
  
  # Remove missing PCA coordinates
  region_data <- region_data %>% drop_na(Dim.1, Dim.2)
  
  # Select environmental variables
  env_vars <- pca.vars1[rownames(pca.vars1) %in% c("Temp", "Sal", "O2", "Fluor", "Depth", "D", "N"), , drop = FALSE]
  other_vars <- pca.vars1[!rownames(pca.vars1) %in% c("Temp", "Sal", "O2", "Fluor", "Depth", "D", "N"), , drop = FALSE]
  
  # Convert rownames to a column for labeling
  env_vars <- env_vars %>% mutate(label = rownames(env_vars))
  other_vars <- other_vars %>% mutate(label = rownames(other_vars))
  
  # **9. Create PCA Plot for Each Region**
  pca_plot <- ggplot(data = region_data, aes(x = Dim.1, y = Dim.2)) +
    geom_point(aes(fill = Cluster), size=6, stroke=0.3, shape=21) +
    theme_classic() +
    xlab(paste0("PC1 (", round(res_pca$eig[1, 2], 1), "%)")) +
    ylab(paste0("PC2 (", round(res_pca$eig[2, 2], 1), "%)")) +
    scale_fill_manual(values = cols, name = "Clusters") +
    ggtitle(Region_labs[region]) +
    scale_x_continuous(limits = c(-20, 30)) +
    scale_y_continuous(limits = c(-20, 20)) +
    
    # Environmental Variables Labels (Red)
    geom_label(data = env_vars, aes(x = Dim.1 * 15, y = Dim.2 * 15, label = label), 
               fill = "white",       
               color = "red", 
               size = 4, 
               fontface = "bold") +
    
    # Other Variable Labels (Black)
    geom_label_repel(
      data = other_vars, 
      aes(x = Dim.1 * 15, y = Dim.2 * 15, label = label), 
      fill = "white",       
      color = "black",      
      size = 4,             
      fontface = "bold"
    ) +
    
    theme(panel.background = element_rect(fill='transparent'), 
          plot.background = element_rect(fill='transparent', color=NA),
          text=element_text(size=25),
          legend.position = "bottom",
          panel.grid.major = element_blank(), #remove major gridlines
          panel.grid.minor = element_blank(), #remove minor gridlines
          legend.background = element_rect(fill='transparent'), #transparent legend bg
          legend.box.background = element_rect(fill='transparent'))
  
  # Save the plot in the list
  pca_plots_list[[region]] <- pca_plot
  ggsave(paste0("PCA_Cluster_", region, ".png"), plot = pca_plot, width = 7, height = 6, bg='transparent')
}

# **10. Arrange PCA Plots for Selected Regions**
PCA_Regions <- ggarrange(
  pca_plots_list$PSAE, pca_plots_list$NPPF, pca_plots_list$NPTG, 
  nrow = 1, common.legend = TRUE, legend = "right"
)


PCA_Regions
ggsave("PCA_Regions.png", PCA_Regions, width = 20, height = 6, bg="transparent")

#### All regions ####
# Combine data for all regions
all_regions_data <- filtered |>
  dplyr::mutate(
    Dim.1 = res_pca$ind$coord[rownames(filtered), 1],
    Dim.2 = res_pca$ind$coord[rownames(filtered), 2]
  ) |>
  tidyr::drop_na(Dim.1, Dim.2)

# Environmental and other variables for plotting
env_vars <- pca.vars1[rownames(pca.vars1) %in% c("Temp", "Sal", "O2", "Fluor", "Depth","CO2", "D", "N"), , drop = FALSE] %>%
  mutate(label = rownames(.))

other_vars <- pca.vars1[!rownames(pca.vars1) %in% c("Temp", "Sal", "O2", "Fluor", "Depth","CO2", "D", "N"), , drop = FALSE] %>%
  mutate(label = rownames(.))

# PCA Plot for all regions combined
combined_pca_plot <- ggplot(data = all_regions_data, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(fill = Cluster), size = 6, stroke = 0.3, shape = 21) +
  scale_fill_manual(values = cols, name = "Clusters") +
  geom_label(data = env_vars, aes(x = Dim.1 * 15, y = Dim.2 * 15, label = label), 
             fill = "white", color = "red", size = 4, fontface = "bold") +
  geom_label_repel(data = other_vars, aes(x = Dim.1 * 15, y = Dim.2 * 15, label = label), 
                   fill = "white", color = "black", size = 4, fontface = "bold") +
  theme_classic() +
  theme(panel.background = element_rect(fill='transparent'), 
        plot.background = element_rect(fill='transparent', color=NA),
        text=element_text(size=30),
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent'))+
  xlab(paste0("PC1 (", round(res_pca$eig[1, 2], 1), "%)")) +
  ylab(paste0("PC2 (", round(res_pca$eig[2, 2], 1), "%)")) +
  ggtitle("PCA Combined Regions") +
  scale_x_continuous(limits = c(-20, 30)) +
  scale_y_continuous(limits = c(-20, 20))

  theme(
    axis.text = element_text(color = "white"),
    axis.title = element_text(color = "white"),
    plot.title = element_text(color = "white", hjust = 0.5),
    legend.title = element_text(color = "white"),
    legend.text = element_text(color = "white"),
    axis.line = element_line(color = "white"),
    axis.ticks = element_line(color = "white"),
    panel.background = element_rect(fill = 'transparent'), 
    plot.background = element_rect(fill = 'transparent', color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent'),
    text=element_text(size=20),
  ) 
# Save plot
combined_pca_plot

ggsave("PCA_Cluster_All_Regions.png", combined_pca_plot, width = 7, height = 6, bg = "transparent")
#### 2. Create an empty list to store combined plots for each region ####
combined_plots_list <- list()

#### 3. Loop through each unique region ####
for (region in unique(Cop_data$Region)) {
  
  # Filter the data by region
  region_data <- Cop_data %>% filter(Region == region)
  
  # Summarize the data to get the total density per 'moc' and then calculate the average
  summarized_data <- region_data %>%
    group_by(D.N) %>%
    summarize(total_density = sum(density, na.rm = TRUE)) %>%
    ungroup() %>%
    summarize(avg_density = mean(total_density, na.rm = TRUE))
  
  # Print or use the average density as needed
  print(paste("Average density for region", region, ":", summarized_data$avg_density))
  
  # Create the day plot for the region
  day_plot <- ggplot(region_data %>% filter(D.N == "D"), aes(x = net, y = Abundance, fill = Cluster)) +
    geom_bar(stat="identity", position = "stack") +
    scale_fill_manual(values = cols) +
    coord_flip() +
    labs(x = "Depth (m)", y = NULL, fill = "Cluster") + # Depth (m) label stays
    theme_classic() +
    theme(legend.position = "none") +
    scale_y_reverse(limits = c(60000,0)) +
    scale_x_discrete(
      breaks = unique(region_data$net), 
      labels = c("1000", "800", "600", "400", "200", "100", "50", "25")) +
    ggtitle(paste0("DAY"))+
    
    theme(#legend.position = "none",
          panel.background = element_rect(fill='transparent'), 
          plot.background = element_rect(fill='transparent', color=NA),
          panel.grid.major = element_blank(), #remove major gridlines
          panel.grid.minor = element_blank(), #remove minor gridlines
          legend.background = element_rect(fill='transparent'), #transparent legend bg
          legend.box.background = element_rect(fill='transparent'), 
          plot.title = element_text(hjust = 0.5))
  
  # Create the night plot for the region
  night_plot <- ggplot(region_data %>% filter(D.N == "N"), aes(x = net, y = Abundance, fill = Cluster )) +
    geom_bar(stat="identity", position = "stack") +
    scale_fill_manual(values = cols) +
    coord_flip() +
    labs(x = "Depth (m)", y = NULL, fill = "Cluster") + # Depth (m) label stays
    theme_classic() +
    theme(legend.position = "none",
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.title.y = element_blank(), 
          axis.line.y = element_blank(), 
          plot.title = element_text(hjust = 1)) +
    scale_x_discrete(
      breaks = unique(region_data$net), 
      labels = c("1000", "800", "600", "400", "200", "100", "50", "25")) +
    scale_y_continuous(limits = c(0,60000))+
    ggtitle(paste0("NIGHT"))+
    
    theme(panel.background = element_rect(fill='transparent'), 
          plot.background = element_rect(fill='transparent', color=NA),
          panel.grid.major = element_blank(), #remove major gridlines
          panel.grid.minor = element_blank(), #remove minor gridlines
          legend.background = element_rect(fill='transparent'), #transparent legend bg
          legend.box.background = element_rect(fill='transparent'), 
          plot.title = element_text(hjust = 0.5))
  
  # Combine the day and night plots side by side
  combined_plot <- plot_grid(day_plot, night_plot, align = "h", ncol = 2, rel_widths = c(1, 1))
  
  # Add the shared x-axis label
  combined_plot_with_label <- ggdraw() + 
    draw_plot(combined_plot, x = 0, y = 0.05, width = 1, height = 0.95) + 
    draw_label(expression(Abundance~(ind~m^{-2})), 
               y = 0.04, # Position the x-axis label at the bottom
               fontface = "bold", 
               size = 12)
  
  # Save the combined plot for the region as a PNG file
  ggsave(filename = paste0("combined_day_night_density_plot_", region, ".png"), 
         plot = combined_plot_with_label, 
         width = 7, height = 5, bg="transparent")
  
  # Save the combined plot in the list, with the region name as the key
  combined_plots_list[[region]] <- combined_plot_with_label
}

#### 4. View the list of combined plots ####
combined_plots_list

#### 4. View the list of combined plots ####
abundance_plot<-ggarrange(combined_plots_list$NPTG, combined_plots_list$PSAE, combined_plots_list$NPPF, nrow = 1, common.legend = TRUE, legend = "right")
ggsave("abundance_plot_Regions.png", abundance_plot, width=17, height=6, bg="white")
# Calculate density difference between day and night
density_difference <- Cop_data %>%
  filter(D.N %in% c("D", "N")) %>%
  group_by(Cluster, net) %>%
  summarise(day_density = sum(density[D.N == "D"]),
            night_density = sum(density[D.N == "N"]),
            density_diff = day_density - night_density) %>%
  ungroup()

#### 2. Create an empty list to store Abundance Difference plots for each region ####
Abundance_diff_plots_list <- list()

#### 3. Loop through each unique region ####
for (region in unique(Cop_data$Region)) {
  
  # Filter the data by region
  region_data <- Cop_data %>% filter(Region == region)
  
  # Calculate Abundance difference (Day - Night)
  Abundance_difference <- region_data %>%
    group_by(net, Cluster) %>%
    summarize(
      day_Abundance = sum(Abundance[D.N == "D"], na.rm = TRUE),
      night_Abundance = sum(Abundance[D.N == "N"], na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(Abundance_diff = day_Abundance - night_Abundance)
  
  # Plot Abundance difference
  Abundance_diff_plot <- ggplot(Abundance_difference, aes(x = net, y = Abundance_diff, fill = as.factor(Cluster))) +
    geom_col(position = "stack") +
    scale_fill_manual(values = cols) +
    coord_flip() +
    labs(
      x = "Net Depth (m)", 
      y = expression("Abundance Difference (Day - Night) (" * ind ~ m^{-2} * ")"), 
      fill = "Cluster"
    ) +
    theme_classic() +
    scale_y_continuous(limits = c(-40000, 40000)) + # Use normal y-axis as we have positive and negative differences
    scale_x_discrete(
      breaks = unique(Abundance_difference$net), 
      labels = c("1000", "800", "600", "400", "200", "100", "50", "25")) +
    #ggtitle(paste0("Difference in Abundance (Day - Night) - ", region))+
    theme(legend.position = "none",
          panel.background = element_rect(fill='transparent'), 
          plot.background = element_rect(fill='transparent', color=NA),
          panel.grid.major = element_blank(), #remove major gridlines
          panel.grid.minor = element_blank(), #remove minor gridlines
          legend.background = element_rect(fill='transparent'), #transparent legend bg
          legend.box.background = element_rect(fill='transparent'))
  
  # Save Abundance difference plot for the region as a PNG file
  ggsave(filename = paste0("Abundance_difference_plot_", region, ".png"), 
         plot = Abundance_diff_plot, 
         width = 8, height = 6)
  
  # Save the plot in the list, with the region name as the key
  Abundance_diff_plots_list[[region]] <- Abundance_diff_plot
}

#### 4. View the list of Abundance difference plots ####
Abundance_diff_plots_list
Differenece_plot<-ggarrange(Abundance_diff_plots_list$Subarctic, Abundance_diff_plots_list$Central, Abundance_diff_plots_list$Subtropical, nrow = 1, common.legend = TRUE, legend = "right")
ggsave("Differenece_plot_Regions.png", Differenece_plot, width=17, height=5, bg="white")

# Loop through each region to create Abundance per Cluster plots
for (region in unique(Cop_data$Region)) {
  
  # Filter the data by region
  region_data <- Cop_data %>% filter(Region == region)
  
  # Plot Abundance per cluster
  dens_plot_cl <- ggplot(region_data, aes(x = net, y = Abundance, group = interaction(as.factor(Cluster), D.N), 
                                          colour = as.factor(Cluster), linetype = D.N, fill = as.factor(Cluster), alpha = D.N)) +
    stat_summary(geom = "line", fun = sum, size = 1.5) +
    stat_summary(geom = "area", fun = sum) +
    facet_wrap(~Cluster, scales = "free_y", nrow = 1) +
    scale_x_discrete(
      limits = rev,
      breaks = unique(region_data$net), 
      labels = c("1000", "800", "600", "400", "200", "100", "50", "25"), 
      name = "Depth (m)") +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    scale_linetype_manual(values = c("solid", "dashed"),
                          breaks = c("D", "N"), 
                          labels = c("Day", "Night"), 
                          name = "Day/night") +
    scale_alpha_manual(values = c(0.1, 0.3), 
                       name = "Day/night", 
                       breaks = c("D", "N"), 
                       labels = c("Day", "Night")) +
    theme_classic() +
    theme(text = element_text(size = 25, face = "bold"), 
          panel.border = element_rect(colour = "black", fill = NA, size = 3),
          #legend.position = "none", 
          axis.text.y = element_text(angle = 90, hjust = 1),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_text(angle = 180, vjust = 1),
          plot.background = element_rect(fill = "transparent"),
          legend.key.width = unit(3, "line")) +
    ylab(expression(bold("Abundance ind" ~ m^-3))) +
    guides(alpha = guide_legend(nrow = 1), 
           linetype = guide_legend(nrow = 1))
  
  # Save Abundance per cluster plot for the region
  ggsave(paste0("abund_plot_cal_", region, ".png"), plot = dens_plot_cl, width = 40, height = 10, bg = "white")
}

####Plot transparency ~ size with marginal distribution#####
# Step 1: Filter `data_full` the same way as `Cop_data`
data_filtered <- data_full %>%
  filter(Cruise == "nh1208", 
         annotation_category %in% c("Calanoida"),
         annotation_status == "validated")

# Assign esd_mm and Mean1 from filtered data correctly
Cop_data <- Cop_data %>%
  mutate(
    esd_mm = data_filtered$esd * 0.0053,
    Mean1 = data_filtered$mean,
    Cluster = factor(Cluster, levels = cluster_order),
    Type = if_else(Cluster %in% c("ST", "MT", "SD"), "Transparent", "Dark")
  )



# Initialize a list to store the plots
plots_list <- list()


for (region in unique(Cop_data$Region)) {
  
  region_data <- Cop_data %>% filter(Region == region)
  
  # Convert to quantile-based binning
  region_data <- region_data %>%
    group_by(Cluster)%>%
    # Step 1: Define quantile-based breaks from 5% to 95%
    mutate(
      mean1_breaks = list(
        quantile(Mean1, probs = seq(0.05, 0.95, length.out = 50), na.rm = TRUE)
      ),
      esd_breaks = list(
        quantile(esd_mm, probs = seq(0.05, 0.95, length.out = 50), na.rm = TRUE)
      )
    ) %>%
    # Step 2: Convert each row to the correct bin
    rowwise() %>%
    mutate(
      Mean1_bin = cut(
        Mean1,
        breaks = mean1_breaks,
        include.lowest = TRUE,
        labels = mean1_breaks[-1]  # upper-bound labels for each bin
      ),
      esd_mm_bin = cut(
        esd_mm,
        breaks = esd_breaks,
        include.lowest = TRUE,
        labels = esd_breaks[-1]
      )
    ) %>%
    ungroup() %>%
    # Step 3: Convert bin factors → numeric (optional)
    mutate(
      Mean1_bin = as.numeric(as.character(Mean1_bin)),
      esd_mm_bin = as.numeric(as.character(esd_mm_bin))
    ) %>%
    # (Optional) drop the temporary columns with break vectors
    select(-mean1_breaks, -esd_breaks)
  
  # Total abundance by cluster
  cluster_abundance <- region_data %>%
    group_by(Cluster) %>%
    dplyr::summarise(total_abundance = sum(Abundance, na.rm = TRUE), .groups = 'drop')
  
  # Abundance-based summary per bin for Mean1
  binned_abundance_mean <- region_data %>%
    group_by(Cluster, Mean1_bin, Type) %>%
    summarise(
      total_abundance = sum(Abundance, na.rm = TRUE), 
      mean_Mean1 = mean(Mean1, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Abundance-based summary per bin for esd_mm
  binned_abundance_esd <- region_data %>%
    group_by(Cluster, esd_mm_bin, Type) %>%
    summarise(
      total_abundance = sum(Abundance, na.rm = TRUE),
      mean_esd_mm = mean(esd_mm, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Cluster summary statistics for plots
  df_summary <- region_data %>%
    group_by(Cluster) %>%
    summarise(
      mean_Mean1 = mean(Mean1, na.rm = TRUE),
      mean_esd_mm = mean(esd_mm, na.rm = TRUE),
      q05_Mean1 = quantile(Mean1, 0.05, na.rm = TRUE),
      q95_Mean1 = quantile(Mean1, 0.95, na.rm = TRUE),
      q05_esd_mm = quantile(esd_mm, 0.05, na.rm = TRUE),
      q95_esd_mm = quantile(esd_mm, 0.95, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Base summary plot
  plot_summary <- ggplot(df_summary, aes(x = mean_esd_mm, y = mean_Mean1, fill = Cluster, color=Cluster)) +
    geom_errorbar(aes(ymin = q05_Mean1, ymax = q95_Mean1), width = 0.08, size = 1.5) +
    geom_errorbarh(aes(xmin = q05_esd_mm, xmax = q95_esd_mm), height = 2.2, size = 1.5) +
    geom_point(size = 10,  shape=21) +
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    theme_classic(base_size = 16) +
    scale_y_continuous(expand = c(0,0), limits = c(100,230)) +
    scale_x_continuous(expand = c(0,0), limits = c(0,5)) +
    labs(x = "ESD (mm)", y = "Transparency") +
    theme(panel.background = element_rect(fill='transparent'), 
          plot.background = element_rect(fill='transparent', color=NA),
          legend.position = "none",
          text=element_text(size=40))
  
  # Transparent Mean1 abundance plot
  plot_mean_transparent <- ggplot(binned_abundance_mean, 
                                  aes(x = as.numeric(Mean1_bin), y = total_abundance, fill = Cluster)) +
    stat_smooth(
      geom = "area",
      alpha = 0.3,            # transparency for the area fill
      aes(fill = Cluster,          # fill color
      color = Cluster),        # outline color
      size = 1,               # thickness of the outline
      position = "identity"
    )+
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    scale_x_continuous(expand = c(0, 0), limits = c(100, 230)) +
    theme_classic() +
    coord_flip() +
    theme(panel.background = element_rect(fill='transparent'), 
          plot.background = element_rect(fill='transparent', color=NA),
          legend.position = "none")
  
  # Transparent esd_mm abundance plot
  plot_esd_transparent <- ggplot(binned_abundance_esd, 
    aes(x = as.numeric(mean_esd_mm), y = total_abundance, fill = Cluster)) +
    stat_smooth(
      geom = "area",
      alpha = 0.3,            # transparency for the area fill
      aes(fill = Cluster,          # fill color
          color = Cluster),        # outline color
      size = 1,               # thickness of the outline
      position = "identity"
    )+
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 5)) +
    theme_classic() +
    theme(panel.background = element_rect(fill='transparent'), 
          plot.background = element_rect(fill='transparent', color=NA),
          legend.position = "none")
  
  # Combine plots
  combined_plot <- plot_summary %>%
    insert_xaxis_grob(plot_esd_transparent, grid::unit(.3, "null"), position = "top") %>%
    insert_yaxis_grob(plot_mean_transparent, grid::unit(.2, "null"), position = "right")
  
  plots_list[[region]] <- ggdraw(combined_plot)
  
  ggsave(filename = paste0("marginal_plot_cop", region, ".png"), plot = plots_list[[region]], width = 10, height = 10, bg = "transparent")
}

plots_list



#marginal_plot<-ggarrange(plots_list$PSAE, plots_list$NPPF, plots_list$NPTG, nrow = 1)
#ggsave("marginal_plot_regions.png", marginal_plot, width=24, height = 8, bg="white")


#### WMD & DVM ####
calc_plot_wmd_dvm <- function(data, colors) {
  
  regions <- unique(data$Region)
  wmd_plots_list <- list()
  dvm_plots_list <- list()
  
  for (region in regions) {
    
    # Filter data by region
    region_data <- data %>% filter(Region == region)
    
    # Calculate WMD
    wmd_stats <- region_data %>%
      mutate(Zm = (depth_max + depth_min) / 2) %>%
      group_by(Cluster, D.N) %>%
      summarise(
        total_abundance=sum(Abundance),
        # Weighted Mean Depth
        WMD = sum(Abundance * Zm, na.rm = TRUE) / total_abundance,
        
        # Weighted Variance & SD
        var_w = sum(Abundance * ((Zm - WMD)^2), na.rm = TRUE) / total_abundance,
        SD = sqrt(var_w),
        
        # Weighted SE using total organisms (sum(Abundance))
        SE = SD / sqrt(sum(Abundance, na.rm = TRUE)),
        
        # 95% CI
        #lower_ci = WMD - 1.96 * SE,
        #upper_ci = WMD + 1.96 * SE,
        .groups = "drop"
      )%>%
      mutate(lower_ci = WMD - 1.96 * SE,
             upper_ci = WMD + 1.96 * SE,
             )
    
    # Add a new column for fill color: White for Day, Cluster Color for Night
    wmd_stats <- wmd_stats %>%
      mutate(
        Cluster = as.character(Cluster),  # Ensure Cluster is character
        fill_color = ifelse(D.N == "D", "white", cols[Cluster]) # Map cluster colors for Night
      )
    wmd_stats$Cluster <- factor(wmd_stats$Cluster, levels = c("ST", "MT", "LT", "SD", "MD", "LD"))
    
    # Plot WMD
    WMD_plot <- ggplot(wmd_stats, aes(x = Cluster, y = WMD, fill = fill_color, color = Cluster, shape = D.N)) +
      geom_point(size = 8, stroke=1.5) +
      
      # Shape styles: Day (D) = 25 (triangle down), Night (N) = 24 (triangle up)
      scale_shape_manual(values = c("D" = 25, "N" = 24), breaks = c("D","N"),
                         labels=c("Day", "Night")) +
      
      # Fill: Uses the predefined fill_color column
      scale_fill_identity() +
      
      # Outline color: Cluster colors
      scale_color_manual(values = cols, guide = FALSE) +
      
      labs(x = "Cluster", y = "WMD (m)") +
      #scale_y_continuous(position = "right") +
      
      theme_classic() +
      theme(
        panel.background = element_rect(fill = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color = NA),
        text = element_text(size = 30),
        legend.key.height = unit(1, "cm"),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.position = "none"
      )
    
    wmd_plots_list[[region]] <- WMD_plot
    ggsave(paste0("WMD_plot_", region, ".png"), WMD_plot, width = 7, height = 6)
    
    # Calculate DVM
    dvm_stats <- wmd_stats %>%
      group_by(Cluster) %>%
      summarise(
        DVM_Day = mean(WMD[D.N == "D"], na.rm = TRUE),
        DVM_Night = mean(WMD[D.N == "N"], na.rm = TRUE),
        DVM = DVM_Day - DVM_Night,
        CI_Lower = DVM - 1.96 * sqrt(mean(SE[D.N == "D"], na.rm = TRUE)^2 + mean(SE[D.N == "N"], na.rm = TRUE)^2),
        CI_Upper = DVM + 1.96 * sqrt(mean(SE[D.N == "D"], na.rm = TRUE)^2 + mean(SE[D.N == "N"], na.rm = TRUE)^2),
        .groups = 'drop'
      )
    
    # Plot DVM
    DVM_plot <- ggplot(dvm_stats, aes(x = Cluster, y = DVM, fill = Cluster)) +
      geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.1, color="black") +
      geom_point(stat = "identity", color = "black", size=5, shape=21) +
      scale_fill_manual(values = colors, guide = FALSE) +
      labs( x = "Cluster", y = "DVM amplitude (m)") +
      theme_classic()+
      theme(
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.position = "none",
        text = element_text(size=20)
        
      )
    
    dvm_plots_list[[region]] <- DVM_plot
    ggsave(paste0("DVM_plot_", region, ".png"), DVM_plot, width = 7, height = 6)
  }
  
  list(WMD = wmd_plots_list, DVM = dvm_plots_list)
}

plots<-calc_plot_wmd_dvm(Cop_data, cols)
plots$WMD$PSAE
plots$DVM$PSAE

calc_plot_dvm_st <- function(
    data,
    colors,
    error_width_esd_list,
    error_height_esd_list,
    error_width_transparency_list,
    error_height_transparency_list
) {
  
  regions <- unique(data$Region)
  dvm_plots_list <- list()
  
  for (region in regions) {
    
    # Retrieve specific error-bar width/height for ESD plots (region-specific)
    current_width_esd <- error_width_esd_list[[region]]
    current_height_esd <- error_height_esd_list[[region]]
    
    # Retrieve specific error-bar width/height for Transparency plots (region-specific)
    current_width_transparency <- error_width_transparency_list[[region]]
    current_height_transparency <- error_height_transparency_list[[region]]
    
    # Filter data by region
    region_data <- data %>% filter(Region == region)
    
    # Calculate WMD
    wmd_stats <- region_data %>%
      mutate(Zm = (depth_max + depth_min) / 2) %>%
      group_by(Cluster, D.N) %>%
      summarise(
        WMD = sum(Abundance * Zm, na.rm = TRUE) / sum(Abundance, na.rm = TRUE),
        SD = sqrt(
          sum(Abundance * (Zm^2), na.rm = TRUE) / sum(Abundance, na.rm = TRUE) - WMD^2
        ),
        SE = SD / sqrt(n()),
        lower_ci = WMD - 1.96 * SE,
        upper_ci = WMD + 1.96 * SE,
        .groups = 'drop'
      )
    
    # Calculate DVM
    dvm_stats <- wmd_stats %>%
      group_by(Cluster) %>%
      summarise(
        DVM_Day = mean(WMD[D.N == "D"], na.rm = TRUE),
        DVM_Night = mean(WMD[D.N == "N"], na.rm = TRUE),
        DVM = DVM_Day - DVM_Night,
        CI_Lower = DVM - 1.96 * sqrt(
          mean(SE[D.N == "D"], na.rm = TRUE)^2 + mean(SE[D.N == "N"], na.rm = TRUE)^2
        ),
        CI_Upper = DVM + 1.96 * sqrt(
          mean(SE[D.N == "D"], na.rm = TRUE)^2 + mean(SE[D.N == "N"], na.rm = TRUE)^2
        ),
        .groups = 'drop'
      )
    
    # Morphology stats
    df_summary <- region_data %>%
      group_by(Cluster) %>%
      summarise(
        mean_Mean1 = mean(Mean1, na.rm = TRUE),
        mean_esd_mm = mean(esd_mm, na.rm = TRUE),
        q05_Mean1 = quantile(Mean1, 0.05, na.rm = TRUE),
        q95_Mean1 = quantile(Mean1, 0.95, na.rm = TRUE),
        q05_esd_mm = quantile(esd_mm, 0.05, na.rm = TRUE),
        q95_esd_mm = quantile(esd_mm, 0.95, na.rm = TRUE),
        .groups = 'drop'
      )
    
    # Merge DVM + morphology
    merged_data <- dvm_stats %>%
      left_join(df_summary, by = "Cluster")
    
    # ESD vs. DVM
    plot_dvm_esd <- ggplot(merged_data, aes(x = mean_esd_mm, y = DVM, color = Cluster)) +
      geom_point(size = 8) +
      # Vertical errorbars
      geom_errorbar(
        aes(ymin = CI_Lower, ymax = CI_Upper),
        width = current_width_esd,
        size = 1.5
      ) +
      # Horizontal errorbars
      geom_errorbarh(
        aes(xmin = q05_esd_mm, xmax = q95_esd_mm),
        height = current_height_esd,
        size = 1.5
      ) +
      scale_color_manual(values = colors) +
      theme_classic(base_size = 16) +
      labs(
        x = "ESD (mm)",
        y = "DVM amplitude (m)"
      ) +
      theme(
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.position = "none", 
        text = element_text(size=30)
      )
    
    dvm_plots_list[[paste0(region, "_esd")]] <- plot_dvm_esd
    ggsave(
      filename = paste0("DVM_ESD_plot_", region, ".png"),
      plot = plot_dvm_esd,
      width = 7,
      height = 6,
      bg = "transparent"
    )
    
    # Transparency vs. DVM
    plot_dvm_transparency <- ggplot(merged_data, aes(x = mean_Mean1, y = DVM, color = Cluster)) +
      geom_point(size = 8) +
      geom_errorbar(
        aes(ymin = CI_Lower, ymax = CI_Upper),
        width = current_width_transparency,
        size = 1
      ) +
      geom_errorbarh(
        aes(xmin = q05_Mean1, xmax = q95_Mean1),
        height = current_height_transparency,
        size = 1
      ) +
      scale_color_manual(values = colors) +
      theme_classic(base_size = 16) +
      labs(
        x = "Transparency",
        y = "DVM amplitude (m)"
        
      ) +
      theme(
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.position = "none",
        text = element_text(size=30)
        
      )
    
    dvm_plots_list[[paste0(region, "_transparency")]] <- plot_dvm_transparency
    ggsave(
      filename = paste0("DVM_Transparency_plot_", region, ".png"),
      plot = plot_dvm_transparency,
      width = 7,
      height = 6,
      bg = "transparent"
    )
  }
  
  list(DVM = dvm_plots_list)
}

# Example Usage:
esd_widths <- list("PSAE" = 0.08, "NPPF" = 0.03, "NPTG" = 0.04)
esd_heights <- list("PSAE" = 2, "NPPF" = 4, "NPTG" = 6)
trans_widths <- list("PSAE" = 2, "NPPF" = 2, "NPTG" = 2)
trans_heights <- list("PSAE" = 3, "NPPF" = 4, "NPTG" = 5)
result <- calc_plot_dvm_st(
  Cop_data,
  cols,
  error_width_esd_list = esd_widths,
  error_height_esd_list = esd_heights,
  error_width_transparency_list = trans_widths,
  error_height_transparency_list = trans_heights
)


#########################PCs~WMD&DVM FOR Each Region####################
# Extract the loadings (contributions) of morphological variables on each PC
# Assuming 'res_pca' is the result from your PCA analysis
loadings <- as.data.frame(res_pca$var$coord[, 1:5])  # Extract loadings for PC1 to PC10

# Optionally, rename the columns for clarity
colnames(loadings) <- paste0("PC", 1:5)

# Add variable names as a column for easier visualization
loadings$Variable <- rownames(loadings)

# Reshape the data for plotting (long format)
library(reshape2)
loadings_long <- melt(loadings, id.vars = "Variable")


# Create the heatmap to show influence of variables on PCs
ggplot(loadings_long, aes(x = variable, y = Variable, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
                       limit = c(-1, 1), space = "Lab", 
                       name="Contribution") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        text = element_text(size = 14)) +
  labs(x = "Principal Component", y = "Morphological Variable", 
       title = "Contribution of Morphological Variables to PCs (PC1 to PC10)")

# Create bins for PC1 and PC2
breaks_pc1 <- 50  # Specify the number of bins for PC1
breaks_pc2 <- 50  # Specify the number of bins for PC2

Cop_data$PC1<-res_pca$ind$coord[, 1]
Cop_data$PC2<-res_pca$ind$coord[, 2]
Cop_data$logPC1<-log10(Cop_data$PC1)
Cop_data$logPC2<-log10(Cop_data$PC2)
Cop_data$PC3<-res_pca$ind$coord[, 3]
Cop_data$PC4<-res_pca$ind$coord[, 4]
Cop_data$PC5<-res_pca$ind$coord[, 5]
Cop_data$logPC3<-log10(Cop_data$PC3)
Cop_data$logPC4<-log10(Cop_data$PC4)
Cop_data$logPC5<-log10(Cop_data$PC5)
# Generate pretty breaks from the esd_mm data
pc1_bin <- pretty(Cop_data$PC1, n = breaks_pc1)

Cop_data$bin_pc1 <- cut(Cop_data$PC1, breaks = 51, pc1_bin[-1])

# Generate pretty breaks from the esd_mm data
PC2_bin <- pretty(Cop_data$PC2, n = breaks_pc2)

Cop_data$bin_PC2 <- cut(Cop_data$PC2, breaks = 70, PC2_bin[-1])


# Calculate WMD and DVM by PC1 bins
result_pc1 <- Cop_data %>%
  mutate(Zm = (depth_max + depth_min) / 2) %>%
  group_by(bin_pc1, D.N) %>%
  summarise(
    total_abundance=sum(Abundance),
    # Weighted Mean Depth
    WMD = sum(Abundance * Zm, na.rm = TRUE) / total_abundance,
    
    # Weighted Variance & SD
    var_w = sum(Abundance * ((Zm - WMD)^2), na.rm = TRUE) / total_abundance,
    SD = sqrt(var_w),
    
    # Weighted SE using total organisms (sum(Abundance))
    SE = SD / sqrt(sum(Abundance, na.rm = TRUE)),
    
    # 95% CI
    #lower_ci = WMD - 1.96 * SE,
    #upper_ci = WMD + 1.96 * SE,
    .groups = "drop"
  )%>%
  mutate(lower_ci = WMD - 1.96 * SE,
         upper_ci = WMD + 1.96 * SE,
  )
# Calculate DVM by PC1 bins
DVM_Cop_datac1 <- result_pc1 %>%
  group_by(bin_pc1) %>%
  mutate(D.N_count = n_distinct(D.N)) %>%
  filter(all(c("D", "N") %in% D.N)) %>%
  summarise(
    DVM_Day = mean(WMD[D.N == "D"], na.rm = TRUE),
    DVM_Night = mean(WMD[D.N == "N"], na.rm = TRUE),
    DVM = DVM_Day - DVM_Night,
    CI_Lower = DVM - 1.96 * sqrt(mean(SE[D.N == "D"], na.rm = TRUE)^2 + mean(SE[D.N == "N"], na.rm = TRUE)^2),
    CI_Upper = DVM + 1.96 * sqrt(mean(SE[D.N == "D"], na.rm = TRUE)^2 + mean(SE[D.N == "N"], na.rm = TRUE)^2),
    .groups = 'drop'
  )

# Repeat the same process for PC2 bins
result_pc2 <- Cop_data %>%
  mutate(Zm = (depth_max + depth_min) / 2) %>%
  group_by(bin_pc2, D.N) %>%
  summarise(
    total_abundance=sum(Abundance),
    # Weighted Mean Depth
    WMD = sum(Abundance * Zm, na.rm = TRUE) / total_abundance,
    
    # Weighted Variance & SD
    var_w = sum(Abundance * ((Zm - WMD)^2), na.rm = TRUE) / total_abundance,
    SD = sqrt(var_w),
    
    # Weighted SE using total organisms (sum(Abundance))
    SE = SD / sqrt(sum(Abundance, na.rm = TRUE)),
    
    # 95% CI
    #lower_ci = WMD - 1.96 * SE,
    #upper_ci = WMD + 1.96 * SE,
    .groups = "drop"
  )%>%
  mutate(lower_ci = WMD - 1.96 * SE,
         upper_ci = WMD + 1.96 * SE,
  )
# Calculate DVM by PC2 bins
DVM_Cop_datac2 <- result_pc2 %>%
  group_by(bin_pc2) %>%
  mutate(D.N_count = n_distinct(D.N)) %>%
  filter(all(c("D", "N") %in% D.N)) %>%
  summarise(
    DVM_Day = mean(WMD[D.N == "D"], na.rm = TRUE),
    DVM_Night = mean(WMD[D.N == "N"], na.rm = TRUE),
    DVM = DVM_Day - DVM_Night,
    CI_Lower = DVM - 1.96 * sqrt(mean(SE[D.N == "D"], na.rm = TRUE)^2 + mean(SE[D.N == "N"], na.rm = TRUE)^2),
    CI_Upper = DVM + 1.96 * sqrt(mean(SE[D.N == "D"], na.rm = TRUE)^2 + mean(SE[D.N == "N"], na.rm = TRUE)^2),
    .groups = 'drop'
  )

result_pc1$Region <- factor(result_pc1$Region, levels = c("PSAE", "NPPF", "NPTG"))
DVM_Cop_datac1$Region <- factor(DVM_Cop_datac1$Region, levels = c("PSAE", "NPPF", "NPTG"))

result_pc1_PSAE<-subset(result_pc1, Region=="PSAE")
result_pc1_NPPF<-subset(result_pc1, Region=="NPPF")
result_pc1_NPTG<-subset(result_pc1, Region=="NPTG")
result_pc1$bin_pc1<-as.numeric(as.character(result_pc1$bin_pc1))
# WMD Plot by PC1
WMD_plot_pc1 <- ggplot(result_pc1, aes(x = bin_pc1, y = WMD)) +
  # Map both shape and fill to D.N
  geom_point(aes(shape = D.N, fill=D.N, size=total_abundance), color = "black") +
  
  # Error bars for Day (D)
  geom_errorbar(
  data = filter(result_pc1, D.N == "D"),
   aes(ymin = WMD, ymax = upper_ci),
  width = 0, color = "black"
  ) +
  
  # Error bars for Night (N)
  geom_errorbar(
    data = filter(result_pc1, D.N == "N"),
    aes(ymin = lower_ci, ymax = WMD),
    width = 0, color = "black"
  ) +
  
  # Classic theme
  theme_classic() +
  
  # Assign shapes 25 (Day) and 24 (Night), or whichever you prefer
  scale_shape_manual(values = c("D" = 25, "N" = 24)) +
  
  # Day fill = "white", Night fill = "black" (can adjust)
  scale_fill_manual(values = c("D" = "white", "N" = "black")) +
  
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    #legend.position = "none",
    text = element_text(size = 30)
  ) +
  labs(x = "PC1", y = "WMD (m)")+
  scale_size_continuous(breaks = c(5000, 10000, 15000, 20000, 25000))
  #scale_x_continuous(limits = c(0,55))+
  #scale_y_continuous(limits = c(0,700))

DVM_Cop_datac1_PSAE<-subset(DVM_Cop_datac1, Region=="PSAE")
DVM_Cop_datac1_NPPF<-subset(DVM_Cop_datac1, Region=="NPPF")
DVM_Cop_datac1_NPTG<-subset(DVM_Cop_datac1, Region=="NPTG")
# DVM Amplitude Plot by PC1
DVM_plot_pc1 <- ggplot(DVM_Cop_datac1, aes(x = bin_pc1, y = DVM)) +
  geom_point(size = 6, shape = 21, fill="black") +
  #geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 4, raw=TRUE), se = FALSE, color = "black", linetype = 2) +
  theme_classic() +
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    legend.position = "none",
    text = element_text(size=30)
    
  )+
  labs(x = "PC1", y = "DVM Amplitude (m)") +
  #scale_x_continuous(limits = c(0,45))+
  #scale_y_continuous(limits = c(-300,300))+
  facet_wrap(~Region, scales = "free")

# Print the plots
print(WMD_plot_pc1)
print(DVM_plot_pc1)

ggsave("WMD_plot_by_PC1_byregion.png", WMD_plot_pc1, width = 16, height = 6)
ggsave("DVM_plot_by_PC1_byregion.png", DVM_plot_pc1, width = 12, height = 6)

WMD_DVM_by_Region_PC1<-ggarrange(WMD_plot_pc1, DVM_plot_pc1, nrow=2, common.legend = TRUE)
ggsave("WMD_DVM_by_Region_PC1.png", WMD_DVM_by_Region_PC1, width=17, height = 13, bg="white")


result_pc2$Region <- factor(result_pc2$Region, levels = c("PSAE", "NPPF", "NPTG"))
DVM_Cop_datac2$Region <- factor(DVM_Cop_datac2$Region, levels = c("PSAE", "NPPF", "NPTG"))

# WMD Plot by PC2
WMD_plot_pc2 <- ggplot(result_pc2, aes(x = bin_pc2, y = WMD)) +
  # Map both shape and fill to D.N
  geom_point(aes(shape = D.N, fill = D.N, size=total_abundance),color = "black") +
  
  # Error bars for Day (D)
  geom_errorbar(
    data = filter(result_pc2, D.N == "D"),
    aes(ymin = WMD, ymax = upper_ci),
    width = 0, color = "black"
  ) +
  
  # Error bars for Night (N)
  geom_errorbar(
    data = filter(result_pc2, D.N == "N"),
    aes(ymin = lower_ci, ymax = WMD),
    width = 0, color = "black"
  ) +
  
  # Classic theme
  theme_classic() +
  
  # Assign shapes 25 (Day) and 24 (Night), or whichever you prefer
  scale_shape_manual(values = c("D" = 25, "N" = 24)) +
  
  # Day fill = "white", Night fill = "black" (can adjust)
  scale_fill_manual(values = c("D" = "white", "N" = "black")) +
  
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    #legend.position = "none",
    text = element_text(size = 30)
  ) +
  labs(x = "PC2", y = "WMD (m)")
#scale_x_continuous(limits = c(0,55))+
#scale_y_continuous(limits = c(0,700))


# DVM Amplitude Plot by PC2
DVM_plot_pc2 <- ggplot(DVM_Cop_datac2, aes(x = bin_pc2, y = DVM)) +
  geom_point(size = 6, shape = 21, fill="black") +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 4, raw = TRUE), se = FALSE, color = "black", linetype = 2) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    legend.position = "none",
    text = element_text(size=30)
    
  )+
  labs(x = "PC2", y = "DVM Amplitude (m)") +
  #ggtitle("DVM Amplitude by PC2") +
  facet_wrap(~Region, scales = "free")

# Print the plots
print(WMD_plot_pc2)
print(DVM_plot_pc2)


ggsave("WMD_plot_by_PC2.png", WMD_plot_pc2, width = 10, height = 8)
ggsave("DVM_plot_by_PC2_byregion.png", DVM_plot_pc2, width = 12, height = 6)

WMD_DVM_by_Region_PC2<-ggarrange(WMD_plot_pc2, DVM_plot_pc2, nrow=2, common.legend = TRUE)
ggsave("WMD_DVM_by_Region_PC2.png", WMD_DVM_by_Region_PC2, width=17, height = 13, bg="white")


# Function to save WMD and DVM plots for each region
save_region_plots <- function(data, region_name, pc_column, plot_title_suffix) {
  # Filter data for the specific region
  region_data <- data %>% filter(Region == region_name)
  
  # WMD Plot
  WMD_plot <- ggplot(region_data, aes(x = !!sym(pc_column), y = WMD, shape = D.N)) +
    geom_point(size = 3, color = "black") +
    geom_errorbar(data = filter(region_data, D.N == "D"), aes(ymin = WMD, ymax = upper_ci), width = 0, color = "black") +
    geom_errorbar(data = filter(region_data, D.N == "N"), aes(ymin = lower_ci, ymax = WMD), width = 0, color = "black") +
    theme_minimal() +
    scale_shape_manual(labels = c("Day", "Night"), values = c(1, 16)) +
    theme(text = element_text(size = 20, face = "bold"),
          panel.border = element_rect(colour = "black", fill = NA, size = 2),
          legend.position = "bottom") +
    labs(x = pc_column, y = "Weighted Mean Depth (m)") 
  
  # DVM Amplitude Plot
  DVM_plot <- ggplot(region_data, aes(x = !!sym(pc_column), y = DVM)) +
    geom_point(size = 3, shape = 17) +
    geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0) +
    geom_smooth(method = "lm", formula = y ~ poly(x, 4, raw = TRUE), se = FALSE, color = "black", linetype = 2) +
    theme_minimal() +
    theme(text = element_text(size = 20, face = "bold"),
          panel.border = element_rect(colour = "black", fill = NA, size = 2),
          legend.position = "bottom") +
    labs(x = pc_column, y = "DVM Amplitude") +
    ggtitle(paste("DVM Amplitude by", pc_column, "for", region_name))
  
  # Save the plots
  ggsave(paste0("WMD_plot_by_", pc_column, "_", region_name, ".png"), WMD_plot, width = 10, height = 8)
  ggsave(paste0("DVM_plot_by_", pc_column, "_", region_name, ".png"), DVM_plot, width = 10, height = 8)
  
  # Combine WMD and DVM plots
  combined_plot <- ggarrange(WMD_plot, DVM_plot, nrow = 2, common.legend = TRUE)
  ggsave(paste0("WMD_DVM_by_", pc_column, "_", region_name, ".png"), combined_plot, width = 17, height = 13, bg = "white")
}

# List of regions
regions <- unique(result_pc1$Region)

# Save plots for each region for PC1
for (region in regions) {
  save_region_plots(result_pc1, region, "bin_pc1", "PC1")
}

# Save plots for each region for PC2
for (region in regions) {
  save_region_plots(result_pc2, region, "bin_pc2", "PC2")
}







#########################Size######################################################

library(dplyr)
library(ggplot2)
library(ggpubr)

#---------------------------------------------------------------------
# 1) CREATE BIN COLUMN FOR ESD (size)
#---------------------------------------------------------------------
# Suppose 'Cop_data' has columns:
#   - Region, D.N (day/night), Abundance (or # of organisms), 
#   - depth_min, depth_max, ESD_mm (size)
#   - We want 50 bins for ESD
ESD_breaks <- 50

Cop_data <- Cop_data %>%
  mutate(
    # create bin_ESD from the ESD_mm column
    bin_ESD = cut(esd_mm, breaks = ESD_breaks, labels = FALSE, include.lowest = TRUE)
  )

#---------------------------------------------------------------------
# 2) CALCULATE WMD FOR (bin_ESD, D.N, Region)
#---------------------------------------------------------------------
result_ESD <- Cop_data %>%
  # Calculate midpoint depth for each row
  mutate(Zm = (depth_min + depth_max) / 2) %>%
  group_by(bin_ESD, D.N, Region) %>%
  summarise(
    total_abundance = sum(Abundance, na.rm = TRUE),
    
    # Weighted Mean Depth
    WMD = sum(Abundance * Zm, na.rm = TRUE) / total_abundance,
    
    # Weighted variance, standard deviation
    var_w = sum(Abundance * (Zm - WMD)^2, na.rm = TRUE) / total_abundance,
    SD = sqrt(var_w),
    
    # Weighted standard error (if each organism is an independent measurement)
    SE = SD / sqrt(total_abundance),
    .groups = 'drop'
  ) %>%
  mutate(
    # 95% Confidence Interval for WMD
    lower_ci = WMD - 1.96 * SE,
    upper_ci = WMD + 1.96 * SE
  )

#---------------------------------------------------------------------
# 3) CALCULATE DVM BY ESD BIN
#---------------------------------------------------------------------
# We want day (D) vs night (N) Weighted Mean Depth difference
DVM_ESD <- result_ESD %>%
  group_by(bin_ESD, Region) %>%
  # Only keep bins that have both Day and Night
  filter(all(c("D", "N") %in% D.N)) %>%
  summarise(
    DVM_Day = mean(WMD[D.N == "D"], na.rm = TRUE),
    DVM_Night = mean(WMD[D.N == "N"], na.rm = TRUE),
    DVM = DVM_Day - DVM_Night,
    # Combine day & night SE by sqrt( SE_day^2 + SE_night^2 )
    SE_Day = mean(SE[D.N == "D"], na.rm = TRUE),
    SE_Night = mean(SE[D.N == "N"], na.rm = TRUE),
    SE_DVM = sqrt(SE_Day^2 + SE_Night^2),
    CI_Lower = DVM - 1.96 * SE_DVM,
    CI_Upper = DVM + 1.96 * SE_DVM,
    .groups = 'drop'
  )

#---------------------------------------------------------------------
# 4) FUNCTIONS TO MAKE & SAVE PLOTS
#---------------------------------------------------------------------
# We define 2 functions that produce ggplot objects: one for WMD, one for DVM.

make_wmd_plot <- function(data, region_name, bin_column, x_label) {
  # Filter data by region
  region_data <- data %>% filter(Region == region_name)
  
  p <- ggplot(region_data, aes(x = !!sym(bin_column), y = WMD, shape = D.N, fill = D.N)) +
    geom_point(size = 6, color = "black") +
    geom_errorbar(
      data = filter(region_data, D.N == "D"),
      aes(ymin = WMD, ymax = upper_ci),
      width = 0, color = "black"
    ) +
    geom_errorbar(
      data = filter(region_data, D.N == "N"),
      aes(ymin = lower_ci, ymax = WMD),
      width = 0, color = "black"
    ) +
    # Day = shape 25 with white fill, Night = shape 24 with black fill
    scale_shape_manual(values = c("D" = 25, "N" = 24)) +
    scale_fill_manual(values = c("D" = "white", "N" = "black")) +
    theme_classic(base_size = 16) +
    labs(
      x = x_label,
      y = "WMD (m)"
      
    ) +
    theme(
      panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      legend.position = "none",
      text = element_text(size=30)
    )
  return(p)
}

make_dvm_plot <- function(data, region_name, bin_column, x_label) {
  region_data <- data %>% filter(Region == region_name)
  
  p <- ggplot(region_data, aes(x = !!sym(bin_column), y = DVM)) +
    geom_point(size = 6, shape = 21, fill = "black", color = "black") +
    geom_smooth(method = "loess", se=FALSE, linetype=2, color="black")+
    #geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0) +
    #geom_smooth(
     # method = "lm", formula = y ~ poly(x, 4, raw = TRUE),
      #se = FALSE, color = "black", linetype = 2
    #) +
    theme_classic(base_size = 16) +
    labs(
      x = x_label,
      y = "DVM (m)"
    ) +
    theme(
      panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      legend.position = "none",
      text = element_text(size=30)
    )
  return(p)
}

#---------------------------------------------------------------------
# 5) LOOP OVER REGIONS, SAVE PLOTS TO LISTS & TO DISK
#---------------------------------------------------------------------
regions <- c("PSAE", "NPPF", "NPTG")  # or unique(result_ESD$Region)

# Lists to store the final plots
wmd_plots_list <- list()
dvm_plots_list <- list()

for(r in regions) {
  
  # (A) Make WMD plot
  wmd_plot <- make_wmd_plot(result_ESD, r, "bin_ESD", "ESD (mm)")
  wmd_plots_list[[r]] <- wmd_plot
  
  # Save WMD plot
  ggsave(
    filename = paste0("WMD_ESD_loess_", r, ".png"),
    plot = wmd_plot,
    width = 7,
    height = 6
  )
  
  # (B) Make DVM plot
  dvm_plot <- make_dvm_plot(DVM_ESD, r, "bin_ESD", "ESD (mm)")
  dvm_plots_list[[r]] <- dvm_plot
  
  # Save DVM plot
  ggsave(
    filename = paste0("DVM_ESD_", r, ".png"),
    plot = dvm_plot,
    width = 7,
    height = 6
  )
  
  # (C) Combine WMD + DVM in one PNG
  # combined_plot <- ggarrange(wmd_plot, dvm_plot, nrow = 2, common.legend = TRUE)
  #ggsave(
  # filename = paste0("WMD_DVM_ESD_", r, ".png"),
  #plot = combined_plot,
  #width = 7,
  #height = 8
  #)
}




#########################TRANSPARENCY######################################################

library(dplyr)
library(ggplot2)
library(ggpubr)

#---------------------------------------------------------------------
# 1) CREATE BIN COLUMN FOR Mean1 (size)
#---------------------------------------------------------------------
# Suppose 'Cop_data' has columns:
#   - Region, D.N (day/night), Abundance (or # of organisms), 
#   - depth_min, depth_max, Mean1_mm (size)
#   - We want 50 bins for Mean1
Mean1_breaks <- 50

Cop_data <- Cop_data %>%
  mutate(
    # create bin_Mean1 from the Mean1_mm column
    bin_Mean1 = cut(Mean1, breaks = Mean1_breaks, labels = FALSE, include.lowest = TRUE)
  )

#---------------------------------------------------------------------
# 2) CALCULATE WMD FOR (bin_Mean1, D.N, Region)
#---------------------------------------------------------------------
result_Mean1 <- Cop_data %>%
  # Calculate midpoint depth for each row
  mutate(Zm = (depth_min + depth_max) / 2) %>%
  group_by(bin_Mean1, D.N, Region) %>%
  summarise(
    total_abundance = sum(Abundance, na.rm = TRUE),
    
    # Weighted Mean Depth
    WMD = sum(Abundance * Zm, na.rm = TRUE) / total_abundance,
    
    # Weighted variance, standard deviation
    var_w = sum(Abundance * (Zm - WMD)^2, na.rm = TRUE) / total_abundance,
    SD = sqrt(var_w),
    
    # Weighted standard error (if each organism is an independent measurement)
    SE = SD / sqrt(total_abundance),
    .groups = 'drop'
  ) %>%
  mutate(
    # 95% Confidence Interval for WMD
    lower_ci = WMD - 1.96 * SE,
    upper_ci = WMD + 1.96 * SE
  )

#---------------------------------------------------------------------
# 3) CALCULATE DVM BY Mean1 BIN
#---------------------------------------------------------------------
# We want day (D) vs night (N) Weighted Mean Depth difference
DVM_Mean1 <- result_Mean1 %>%
  group_by(bin_Mean1, Region) %>%
  # Only keep bins that have both Day and Night
  filter(all(c("D", "N") %in% D.N)) %>%
  summarise(
    DVM_Day = mean(WMD[D.N == "D"], na.rm = TRUE),
    DVM_Night = mean(WMD[D.N == "N"], na.rm = TRUE),
    DVM = DVM_Day - DVM_Night,
    # Combine day & night SE by sqrt( SE_day^2 + SE_night^2 )
    SE_Day = mean(SE[D.N == "D"], na.rm = TRUE),
    SE_Night = mean(SE[D.N == "N"], na.rm = TRUE),
    SE_DVM = sqrt(SE_Day^2 + SE_Night^2),
    CI_Lower = DVM - 1.96 * SE_DVM,
    CI_Upper = DVM + 1.96 * SE_DVM,
    .groups = 'drop'
  )

#---------------------------------------------------------------------
# 4) FUNCTIONS TO MAKE & SAVE PLOTS
#---------------------------------------------------------------------
# We define 2 functions that produce ggplot objects: one for WMD, one for DVM.

make_wmd_plot <- function(data, region_name, bin_column, x_label) {
  # Filter data by region
  region_data <- data %>% filter(Region == region_name)
  
  p <- ggplot(region_data, aes(x = !!sym(bin_column), y = WMD, shape = D.N, fill = D.N)) +
    geom_point(size = 6, color = "black") +
    geom_errorbar(
      data = filter(region_data, D.N == "D"),
      aes(ymin = WMD, ymax = upper_ci),
      width = 0, color = "black"
    ) +
    geom_errorbar(
      data = filter(region_data, D.N == "N"),
      aes(ymin = lower_ci, ymax = WMD),
      width = 0, color = "black"
    ) +
    # Day = shape 25 with white fill, Night = shape 24 with black fill
    scale_shape_manual(values = c("D" = 25, "N" = 24)) +
    scale_fill_manual(values = c("D" = "white", "N" = "black")) +
    theme_classic(base_size = 16) +
    labs(
      x = x_label,
      y = "WMD (m)"
      
    ) +
    theme(
      panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      legend.position = "none",
      text = element_text(size=30)
    )
  return(p)
}

make_dvm_plot <- function(data, region_name, bin_column, x_label) {
  region_data <- data %>% filter(Region == region_name)
  
  p <- ggplot(region_data, aes(x = !!sym(bin_column), y = DVM)) +
    geom_point(size = 6, shape = 21, fill = "black", color = "black") +
    #geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0) +
    geom_smooth(
      method = "lm", formula = y ~ poly(x, 4, raw = TRUE),
      se = FALSE, color = "black", linetype = 2
    ) +
    theme_classic(base_size = 16) +
    labs(
      x = x_label,
      y = "DVM (m)"
    ) +
    theme(
      panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      legend.position = "none",
      text = element_text(size=30)
    )
  return(p)
}

#---------------------------------------------------------------------
# 5) LOOP OVER REGIONS, SAVE PLOTS TO LISTS & TO DISK
#---------------------------------------------------------------------
regions <- c("PSAE", "NPPF", "NPTG")  # or unique(result_Mean1$Region)

# Lists to store the final plots
wmd_plots_list <- list()
dvm_plots_list <- list()

for(r in regions) {
  
  # (A) Make WMD plot
  wmd_plot <- make_wmd_plot(result_Mean1, r, "bin_Mean1", "Transparency")
  wmd_plots_list[[r]] <- wmd_plot
  
  # Save WMD plot
  ggsave(
    filename = paste0("WMD_Mean1_", r, ".png"),
    plot = wmd_plot,
    width = 7,
    height = 6
  )
  
  # (B) Make DVM plot
  dvm_plot <- make_dvm_plot(DVM_Mean1, r, "bin_Mean1", "Transparency")
  dvm_plots_list[[r]] <- dvm_plot
  
  # Save DVM plot
  ggsave(
    filename = paste0("DVM_Mean1_", r, ".png"),
    plot = dvm_plot,
    width = 7,
    height = 6
  )
  
  # (C) Combine WMD + DVM in one PNG
  # combined_plot <- ggarrange(wmd_plot, dvm_plot, nrow = 2, common.legend = TRUE)
  #ggsave(
  # filename = paste0("WMD_DVM_Mean1_", r, ".png"),
  #plot = combined_plot,
  #width = 7,
  #height = 8
  #)
}