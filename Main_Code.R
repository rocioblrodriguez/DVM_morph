#### 1. Load and Install Required Packages ####
packages <- c("tidyr", "dplyr", "readr", "stringr", "data.table", "ggplot2", "readxl", 
              "ggpubr", "gridExtra", "RColorBrewer", "colorspace", "FactoMineR", 
              "factoextra", "gginnards", "cellWise", "corrplot", "vegan", "morphr", 
              "purrr", "imager", "ggrepel", "cowplot", "Nmisc", "bestNormalize")

package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})


#### 1. Load Data ####
#data <-read.table("/Users/rociorodriguez/Desktop/Calanoida data/ecotaxa_export_5421_20241212_2208.tsv", header = TRUE, sep = "\t")
# Read the data
data <- read.table("/Users/rociorodriguez/Desktop/export_5421_20250421_1659/ecotaxa_export_5421_20250421_1659.tsv", 
                   header = TRUE, sep = "\t")
env_data <- read.csv("/Users/rociorodriguez/Desktop/export_5421_20240919_0518/Environmental_data.csv")
metadata <- read.csv("/Users/rociorodriguez/Desktop/export_5421_20240919_0518/Final_Samples_Gradients.csv")

#### 2. Data Preprocessing ####
# Split 'object_id' and calculate Abundance and abundance.m2
data <- data %>%
  separate(object_id, c("cruise", "moc", "net", "fraction"), extra = "drop", remove = FALSE) %>%
  mutate(density = (acq_sub_part / sample_tot_vol),
         Abundance = (acq_sub_part / sample_tot_vol)* (object_depth_max - object_depth_min))

# Rename columns in environmental data and standardize 'Net' for merging
env_data <- env_data %>%
  rename(
    Temp = `avg..t`,
    Sal = `avg.sal`,
    Fluor = `avg.fluor`,
    O2 = `avg.o2`,
    Depth_max = `Max.Depth`,
    Depth_min = `Min.Depth`
  ) %>%
  mutate(Net = paste0("n", Net))

# Rename and modify columns in metadata for merging
metadata <- metadata %>%
  rename(moc = Cast, lat = Latitude, lon = Longitude) %>%
  mutate(moc = paste0("m", moc),
         Cruise = ifelse(Cruise == "New Horizon", "nh1208", Cruise))

#### 3. Rename 'cruise' to 'Cruise' in the main dataset ####
colnames(data)[colnames(data) == "cruise"] <- "Cruise"

#### 4. Merge All Dataframes ####
# Merge environmental data with the main copepod dataset
data_merged <- data %>%
  left_join(env_data, by = c("Cruise" = "Cruise", "moc" = "Tow", "net" = "Net")) %>%
  rename_with(~ gsub("^object_", "", .x))

# Merge metadata with the transformed copepod + environmental data
data_full <- merge(data_merged, metadata[, c("moc", "Cruise", "D.N")], 
                   by = c("moc", "Cruise"), all.x = TRUE)

#write_csv(data_full, "data_full.csv")



#### 5. Standardize Environmental Variables ####
# Standardize selected environmental variables
data_full[, c("Temp", "Sal", "O2", "Fluor")] <- decostand(data_full[, c("Temp", "Sal", "O2", "Fluor")], method = "standardize")


# Select morphological variables for transformation
morphological_vars <- data_full[, 23:89]  # Adjust indices as per your dataset

# Apply Yeo-Johnson transformation with standardization
YJ_transformed <- apply(morphological_vars, 2, function(x) yeojohnson(x, standardize = TRUE)$x.t)

# Combine transformed morphological variables back into the dataset
data_t2 <- data_full %>%
  select(-one_of(colnames(morphological_vars))) %>%
  bind_cols(as.data.frame(YJ_transformed))


# Final dataset `data_transformed` is ready for analysis

# Create 'Region' column based on station values
data_t2 <- data_t2 %>%
  mutate(Region = case_when(
    # Stations from Cruise "nh1208"
    Cruise == "nh1208" & Station %in% c(7, 11) ~ "PSAE",
    Cruise == "nh1208" & Station %in% c(15, 21, 23) ~ "NPPF",
    Cruise == "nh1208" & Station %in% c(27, 34) ~ "NPTG",
    
    # Stations from Cruise "oc473"
    Cruise == "oc473" & Station %in% c(5, 8) ~ "NASW",
    Cruise == "oc473" & Station %in% c(13, 21) ~ "GFST",
    Cruise == "oc473" & Station %in% c(26, 31) ~ "NADR",
    
    # Default case
    TRUE ~ NA_character_
  ))
# Subset the data to only include rows where the "annotation_category" column is equal to "Calanoida"
#Cop_data <- subset(data_t2, annotation_category == "Calanoida" & Cruise =="nh1208")
# Variables
Cruise <- "nh1208"
Taxa <- c("Actinopterygii", "Annelida", "Bryozoa", "Cephalochordata", "Chaetognatha",
          "Cnidaria<Metazoa", "Hydrozoa", "Siphonophorae", "Amphipoda", "Calanoida", 
          "Cyclopoida", "Harpacticoida", "Mormonilla", "Decapoda", "Euphausiacea", 
          "calyptopsis<Euphausiacea", "Ostracoda", "Echinodermata", "Harosa", "Foraminifera",  
          "Heteropoda", "Mollusca", "Cephalopoda", "Gastropoda<Mollusca", 
          "Cavoliniidae", "Creseidae", "Gymnosomata", "Limacinidae", "Oikopleura", 
          "Pseudothecosomata", "Doliolida", "Salpida")
Status<-"validated"

# Filtering the data
 
Cop_data <- data_t2 %>%
  filter(Cruise == "nh1208", 
         annotation_category %in% c("Calanoida"),
         annotation_status == Status)

## Image Processing ####

# Add image paths for each data entry
#images_directory <- "imgs/"
#Cop_data <- Cop_data %>% mutate(img_path = str_c(images_directory, id, ".jpg"))
library(corrr)   # For correlation analysis
library(caret)   # For finding highly correlated variables
# Subset columns 116 to 182
df_subset <- Cop_data[, 116:182]
# Remove columns with NA values
Cop_data_cleaned <- df_subset[, colSums(is.na(df_subset)) == 0]
# Convert any factor/character columns to numeric if needed
# (If they're actually categories, consider whether correlation is meaningful.)
cols_to_drop <- c("ystart", "bx", "by", "xstart", "angle")
Cop_data_pruned <- Cop_data_cleaned[
  , !(names(Cop_data_cleaned) %in% cols_to_drop)
]
cor_matrix <- cor(Cop_data_pruned, use = "pairwise.complete.obs")  # Compute correlation matrix

highly_correlated <- findCorrelation(cor_matrix, cutoff = 0.99, names = TRUE)  # Find highly correlated variables
selected_morpho_data <- Cop_data[, colnames(morphological_vars) %in% highly_correlated]
#### 7. Compute Abundance and Prepare for PCA and Clustering ####
morphological_variables <- c('mean', 'stddev', 'mode', 'min', 'skew', 'sr', 'kurt', 'histcum3', "nb1", 'fcons',
                             'meanpos', 'height', 'major', 'max', 'perim.', 'width', 'minor', 'xstart', 'ystart',
                             'cdexc', 'fractal', 'angle', 'circ.', 'symetrieh', 'thickr', 'elongation', 
                             'perimferet', 'perimmajor')

morphological_variables <- c("mean", "mode", "min", "max", "angle", "circ.", 
  "kurt", "fractal", "nb1", "symetrievc", "fcons", 
  "thickr",  "elongation", "meanpos", "sr", 
  "feretareaexc", "perimferet", "cdexc",  "major", "stddev")


# Multiply morphological data by the square root of Abundance
morphological_data <-Cop_data[,highly_correlated]
abundance_sqrt <- sqrt(Cop_data$Abundance)
weighted_data <- sweep(morphological_data, 1, abundance_sqrt, `*`)
weighted_data <- weighted_data[, colSums(is.na(weighted_data)) < nrow(weighted_data)]
weighted_data1<-cbind(Cop_data[, c("depth_max", "Temp", "Sal", 
                                           "O2", "Fluor",  "CO2","D.N")], weighted_data)
#### 8. PCA and K-Means Clustering Analysis ####
# Perform PCA on the weighted morphological data
#res.pca <- PCA(weighted_data, ncp = 10, graph = FALSE)  # ncp = 10 limits to the first 10 PCs

res.pca <- PCA(weighted_data1, 
               scale.unit = TRUE, 
               quanti.sup = 1:6, 
               quali.sup = 7,
               graph = FALSE, ncp = 20)

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
distance_matrix <- dist(sample_data)  # Compute Euclidean distance matrix
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
dend <- color_branches(dend, k = 6, col = cols[as.character(labels[kmeans_result$cluster])])

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
pca.vars <- rbind(res.pca$var$coord, res.pca$quanti.sup$coord, res.pca$quali.sup$coord) %>% as.data.frame()

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
    mutate(Dim.1 = res.pca$ind$coord[rownames(region_data), 1],
           Dim.2 = res.pca$ind$coord[rownames(region_data), 2])
  
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
    xlab(paste0("PC1 (", round(res.pca$eig[1, 2], 1), "%)")) +
    ylab(paste0("PC2 (", round(res.pca$eig[2, 2], 1), "%)")) +
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
all_regions_data <- Cop_data %>%
  mutate(Dim.1 = res.pca$ind$coord[rownames(Cop_data), 1],
         Dim.2 = res.pca$ind$coord[rownames(Cop_data), 2]) %>%
  drop_na(Dim.1, Dim.2)

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
  xlab(paste0("PC1 (", round(res.pca$eig[1, 2], 1), "%)")) +
  ylab(paste0("PC2 (", round(res.pca$eig[2, 2], 1), "%)")) +
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
# Assuming 'res.pca' is the result from your PCA analysis
loadings <- as.data.frame(res.pca$var$coord[, 1:5])  # Extract loadings for PC1 to PC10

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

Cop_data$PC1<-res.pca$ind$coord[, 1]
Cop_data$PC2<-res.pca$ind$coord[, 2]
Cop_data$logPC1<-log10(Cop_data$PC1)
Cop_data$logPC2<-log10(Cop_data$PC2)
Cop_data$PC3<-res.pca$ind$coord[, 3]
Cop_data$PC4<-res.pca$ind$coord[, 4]
Cop_data$PC5<-res.pca$ind$coord[, 5]
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



