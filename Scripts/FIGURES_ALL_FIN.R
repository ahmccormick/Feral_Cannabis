
########################################################################################################################################################################
# Figure 1
########################################################################################################################################################################

#######################
# Figure 1A
#######################
#MAP
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(readr)

#/Users/annamccormick/R/Feral_Cannabis_EGS/VCF_overlays/RESULTS2_noLD
# your custom palette
# 1) Read your points
df <- read_csv("/Users/annamccormick/R/Feral_Cannabis_EGS/VCF_overlays/RESULTS2_noLD/Ren_Soorni_Aina_PCA_with_FusionPop_coords_noLD.csv",
               col_types = cols(
                 sample.id          = col_character(),
                 EV1                = col_double(),
                 EV2                = col_double(),
                 Fusion_Populations = col_character(),
                 Latitude           = col_double(),
                 Longitude          = col_double()
               ))

my_cols <- c(
  "Aina_cluster_1"                = "darkred",      
  "Aina_cluster_2"                = "darkblue",     
  "Aina_cluster_3"                = "darkgreen",    
  "Aina_cluster_4"                = "purple4",      
  "Aina_cluster_5"                = "darkorange4",  
  "Ren_basal"            = "orange",       
  "Ren_drug-type"        = "red",          
  "Ren_drug-type feral"  = "blue",         
  "Ren_hemp-type"        = "green",        
  "Soorni_Population_1"     = "purple",       
  "Soorni_Population_2"     = "#8DA0CB"       
)

pop_labels <- c(
  "Aina_cluster_1"       = "Aina Cluster 1",
  "Aina_cluster_2"       = "Aina Cluster 2",
  "Aina_cluster_3"       = "Aina Cluster 3",
  "Aina_cluster_4"       = "Aina Cluster 4",
  "Aina_cluster_5"       = "Aina Cluster 5",
  "Ren_basal"            = "Ren Basal",
  "Ren_drug-type"        = "Ren Drug-type",
  "Ren_drug-type feral"  = "Ren Drug-type feral",
  "Ren_hemp-type"        = "Ren Hemp-type",
  "Soorni_Population_1"  = "Soorni Population 1",
  "Soorni_Population_2"  = "Soorni Population 2"
)

# 2) Load world map
world <- ne_countries(scale = "medium", returnclass = "sf")

p_map <- ggplot(data = world) +
  geom_sf(fill = "gray95", colour = "gray70", size = 0.2) +
  geom_point(
    data = df,
    aes(x = Longitude, y = Latitude, color = Population),
    size = 2, alpha = 0.8
  ) +
  scale_color_manual(
    name   = "Population",
    values = my_cols,
    labels = pop_labels,
    na.value = "black"
  ) +
  coord_sf(xlim = c(-180, 180), ylim = c(10, 70), expand = FALSE) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(
    panel.background = element_rect(fill = "aliceblue"),
    panel.grid.minor = element_blank()
  )
p_map
ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/Figure1A.pdf",
  plot     = p_map,
  width    = 12,
  height   = 10,
  units    = "in"
)

#######################
# Figure 1B
#######################
#ONLY AINA
library(dplyr)
library(ggplot2)
library(maps)

# 0) your data
df_clean <- read.csv(
  "/Users/annamccormick/R/Feral_Cannabis_EGS/VCF/filtering/10chr_PCA_with_coord_and_clusters.csv",
  stringsAsFactors = FALSE
) %>%
  filter(
    !is.na(HCPC_cluster),
    !is.na(Longitude),
    !is.na(Latitude)
  )

# 1) state polygon
states <- map_data("state")

# 2) state label positions
#    maps::state.center gives a data.frame of x/y for each state in state.name
state_labs <- data.frame(
  state = tolower(state.name),
  long  = state.center$x,
  lat   = state.center$y,
  stringsAsFactors = FALSE
) %>% 
  # keep only those in our lon/lat window
  filter(long >= -100, long <= -75, lat >= 40, lat <= 50)

# 3) hulls + counts (as before)
counts <- df_clean %>%
  count(HCPC_cluster) %>%
  arrange(HCPC_cluster) %>%
  mutate(label = paste0("Cluster ", HCPC_cluster, " (n=", n, ")"))

cluster_labels <- setNames(counts$label, counts$HCPC_cluster)

hulls <- df_clean %>%
  group_by(HCPC_cluster) %>%
  slice(chull(Longitude, Latitude)) %>%
  ungroup()

# 4) plot
ggplot() +
  # state fill
  geom_polygon(
    data   = states,
    aes(x = long, y = lat, group = group),
    fill   = "grey95", colour = "grey70", size = 0.3
  ) +
  # hulls
  geom_polygon(
    data   = hulls,
    aes(
      x     = Longitude,
      y     = Latitude,
      group = HCPC_cluster,
      fill  = factor(HCPC_cluster)
    ),
    alpha  = 0.15, colour = NA
  ) +
  # points
  geom_point(
    data = df_clean,
    aes(x = Longitude, y = Latitude, colour = factor(HCPC_cluster)),
    size  = 2, alpha = 0.8
  ) +
  # state labels
  geom_text(
    data    = state_labs,
    aes(x = long, y = lat, label = tools::toTitleCase(state)),
    size    = 3,
    colour  = "grey40",
    fontface = "plain"
  ) +
  coord_quickmap(xlim = c(-100, -75), ylim = c(40, 50)) +
  scale_fill_brewer(
    "HCPC cluster",
    palette = "Set1",
    labels  = cluster_labels
  ) +
  scale_colour_brewer(
    "HCPC cluster",
    palette = "Set1",
    labels  = cluster_labels
  ) +
  theme_minimal() +
  labs(
    title = "Feral Cannabis Sample Locations by HCPC Cluster",
    x     = "Longitude",
    y     = "Latitude"
  ) +
  theme(
    panel.background = element_rect(fill = "aliceblue"),
    panel.grid       = element_line(colour = "white"),
    legend.position  = "right"
  )

#labels fix
my_cols <- c(
  "Aina_cluster_1" = "darkred",
  "Aina_cluster_2" = "darkblue",
  "Aina_cluster_3" = "darkgreen",
  "Aina_cluster_4" = "purple4",
  "Aina_cluster_5" = "darkorange4"
)

df_clean <- df_clean %>%
  mutate(
    Aina_cluster = paste0("Aina_cluster_", HCPC_cluster)
  )

hulls <- df_clean %>%
  group_by(Aina_cluster) %>%
  slice(chull(Longitude, Latitude)) %>%
  ungroup()


counts <- df_clean %>%
  count(Aina_cluster) %>%
  arrange(Aina_cluster) %>%
  mutate(label = paste0(Aina_cluster, " (n=", n, ")"))

cluster_labels <- setNames(counts$label, counts$Aina_cluster)

my_shapes <- c(
  "Aina_cluster_1" = 16,  # filled circle
  "Aina_cluster_2" = 17,  # filled triangle
  "Aina_cluster_3" = 15,  # filled square
  "Aina_cluster_4" = 18,  # filled diamond
  "Aina_cluster_5" = 8    # star
)

counts <- df_clean %>%
  count(Aina_cluster) %>%
  arrange(Aina_cluster) %>%
  mutate(
    cluster_num = gsub("Aina_cluster_", "", Aina_cluster),
    label = paste0("Aina Cluster ", cluster_num, " (n=", n, ")")
  )

cluster_labels <- setNames(counts$label, counts$Aina_cluster)

p <- ggplot() +
  geom_polygon(
    data   = states,
    aes(x = long, y = lat, group = group),
    fill   = "grey95",
    colour = "grey70",
    size   = 0.3
  ) +
  
  geom_point(
    data = df_clean,
    aes(
      x      = Longitude,
      y      = Latitude,
      colour = Aina_cluster,
      shape  = Aina_cluster
    ),
    size  = 2.4,
    alpha = 0.85
  ) +
  
  geom_text(
    data     = state_labs,
    aes(x = long, y = lat, label = tools::toTitleCase(state)),
    size     = 3,
    colour   = "grey40"
  ) +
  
  coord_quickmap(xlim = c(-100, -75), ylim = c(40, 50)) +
  
  scale_colour_manual(
    name   = "Aina cluster",
    values = my_cols,
    labels = cluster_labels,
    drop   = FALSE
  ) +
  scale_shape_manual(
    name   = "Aina cluster",
    values = my_shapes,
    labels = cluster_labels,
    drop   = FALSE
  ) +
  
  theme_minimal() +
  labs(
    title = "",
    x     = "Longitude",
    y     = "Latitude"
  ) +
  theme(
    panel.background = element_rect(fill = "aliceblue"),
    panel.grid       = element_line(colour = "white"),
    legend.position  = "right"
  )

p

ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/Figure1B.pdf",
  plot     = p,
  width    = 12,
  height   = 4.5,
  units    = "in"
)

#######################
# Figure 1C
#######################
library(dplyr)
library(ggplot2)
library(maps)
pca_df <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_VCF_remake/VCF_overlays/PCA_merged.csv")

# define your palette
my_cols <- c(
  "1"                = "darkred",  # red
  "2"                = "darkblue",  # blue
  "3"                = "darkgreen",  # green
  "4"                = "purple4",  # purple
  "5"                = "darkorange4",  # orange
  "basal"            = "orange",  # brown
  "drug-type"        = "red",  # pink
  "drug-type feral"  = "blue",  # grey
  "hemp-type"        = "green",  # teal
  "Population_1"     = "purple",  # yellow
  "Population_2"     = "#8DA0CB"   # lavender
)

# ensure clusters match the palette keys (as character/factor)
pca_df$HCPC_cluster <- factor(pca_df$HCPC_cluster, levels = names(cluster_cols))

pca_df$HCPC_cluster <- factor(
  pca_df$HCPC_cluster,
  levels = names(my_cols)
)


ggplot(pca_df, aes(x = EV1, y = EV2, color = HCPC_cluster)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(
    name   = "HCPC Group",
    values = my_cols
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x     = "PC1 (4.17%)",
    y     = "PC2 (1.22%)"
  ) +
  theme(
    panel.grid.minor = element_blank()
  )

legend_labels <- c(
  "1"               = "Aina Cluster 1",
  "2"               = "Aina Cluster 2",
  "3"               = "Aina Cluster 3",
  "4"               = "Aina Cluster 4",
  "5"               = "Aina Cluster 5",
  "basal"           = "Ren Basal",
  "drug-type"       = "Ren Drug-type",
  "drug-type feral" = "Ren Drug-type (feral)",
  "hemp-type"       = "Ren Hemp-type",
  "Population_1"    = "Soorni Population 1",
  "Population_2"    = "Soorni Population 2"
)

p <- ggplot(pca_df, aes(EV1, EV2)) +
  stat_ellipse(
    aes(color = HCPC_cluster),
    type = "norm",
    level = 0.95,
    linewidth = 0.6,
    na.rm = TRUE
  ) +
  geom_point(
    aes(color = HCPC_cluster),
    size = 2,
    alpha = 0.8
  ) +
  scale_color_manual(
    values = my_cols,
    labels = legend_labels,
    name   = "Population"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x = "PC1 (4.17%)",
    y = "PC2 (1.22%)"
  ) +
  theme(panel.grid.minor = element_blank())

p

ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/Figure1A.pdf",
  plot     = p,
  width    = 10,
  height   = 6,
  units    = "in"
)

########################################################################################################################################################################
# Figure 2
########################################################################################################################################################################

####################################################################
#Importing the contree file and adding colours and aligning tip labels
####################################################################
library(ape)
library(ggtree)
library(dplyr)
library(ggplot2)

# Load the tree
tree <- read.tree("/Users/annamccormick/R/Feral_Cannabis_EGS/newTree/scratch_round1/merged_Aina_Ren_Soorni_n909_13389.min1.phy.varsites.phy.contree")
class(tree)

# Load metadata
metadata <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/newTree/scratch_round1/Ren_Soorni_Aina_metadata_info.csv")
colnames(metadata)

# Rename key columns for clarity
# Assuming columns are: sample.id, FIN_name, [other metadata]
colnames(metadata)[which(colnames(metadata) == "sample.id")] <- "Current"
colnames(metadata)[which(colnames(metadata) == "FIN_name")]  <- "Full"


# Replace tree tip labels with FIN_name, but preserve unmatched tips ("hops")
tree$tip.label <- ifelse(
  tree$tip.label %in% metadata$Current,
  metadata$Full[match(tree$tip.label, metadata$Current)],
  tree$tip.label  # keep labels like "hops"
)

# Assign groups for coloring
metadata <- metadata %>%
  mutate(
    Group = case_when(
      grepl("Aina_cluster_1", Full)        ~ "Aina_cluster_1",
      grepl("Aina_cluster_2", Full)        ~ "Aina_cluster_2",
      grepl("Aina_cluster_3", Full)        ~ "Aina_cluster_3",
      grepl("Aina_cluster_4", Full)        ~ "Aina_cluster_4",
      grepl("Aina_cluster_5", Full)        ~ "Aina_cluster_5",
      grepl("Ren_basal", Full)             ~ "Ren_basal",
      grepl("Ren_drug-type feral", Full)   ~ "Ren_drug-type feral",
      grepl("Ren_drug-type", Full)         ~ "Ren_drug-type",
      grepl("Ren_hemp-type", Full)         ~ "Ren_hemp-type",
      grepl("Soorni_Population_1", Full)   ~ "Soorni_Population_1",
      grepl("Soorni_Population_2", Full)   ~ "Soorni_Population_2",
      TRUE ~ "Other"
    )
  )

# Join tree labels with metadata in FIN_name
tree_data <- data.frame(label = tree$tip.label) %>%
  left_join(metadata %>% distinct(Full, .keep_all = TRUE), by = c("label" = "Full"))

#Ladderize and root the tree 
tree <- ladderize(tree, right = FALSE)
root_sample <- "hops"
tree <- root(tree, outgroup = root_sample, resolve.root = TRUE)

# Save processed tree
#write.tree(tree, file = "/Users/annamccormick/R/Feral_Cannabis_EGS/newTree/scratch_round1/Aina_Ren_Soorni_n909_13389_processed.contree")


####################################################################
# Define your color palette
my_cols <- c(
  "Aina_cluster_1"        = "darkred",      
  "Aina_cluster_2"        = "darkblue",     
  "Aina_cluster_3"        = "darkgreen",    
  "Aina_cluster_4"        = "purple4",      
  "Aina_cluster_5"        = "darkorange4",  
  "Ren_basal"             = "orange",       
  "Ren_drug-type"         = "red",          
  "Ren_drug-type feral"   = "blue",         
  "Ren_hemp-type"         = "green",        
  "Soorni_Population_1"   = "purple",       
  "Soorni_Population_2"   = "#8DA0CB",
  "Other"                 = "gray50",        # fallback group
  "hops"                  = "black"          # manual assignment if needed
)

# OPTIONAL: manually assign 'hops' group if it’s not in the metadata
tree_data$Group[is.na(tree_data$Group) & tree_data$label == "hops"] <- "hops"

# Plot with ggtree
ggtree(tree) %<+% tree_data +
  geom_tiplab(size = 2) +
  geom_tippoint(aes(color = Group), size = 2) +
  scale_color_manual(values = my_cols) +
  theme_tree2() +
  scale_x_continuous(
    trans = "log1p",
    breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 10, 15),
    labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "1", "10", "15"),
    expand = c(0.2, 0.5)
  ) +
  theme(legend.position = "right")

#move labels to far right
ggtree(tree) %<+% tree_data +
  geom_tiplab(size = 0.7, align = TRUE, linesize = 0.2) +  # align tips + smaller font
  geom_tippoint(aes(color = Group), size = 1.8) +         # adjust dot size to match
  scale_color_manual(values = my_cols) +
  theme_tree2() +
  scale_x_continuous(
    trans = "log1p",
    breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 10, 15),
    labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "1", "10", "15"),
    expand = c(0.2, 0.5)
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
#ggsave("/Users/annamccormick/R/Feral_Cannabis_EGS/newTree/scratch_round1/final_tree.pdf", width = 12, height = 14)

#
ggtree(tree) %<+% tree_data +
  geom_tiplab(size = 0.5, align = TRUE, linesize = 0.2, offset = 0.01) +
  geom_tippoint(aes(color = Group), size = 1.8) +
  geom_text2(aes(label = ifelse(!isTip, label, "")), hjust = -0.2, size = 1.5) +
  scale_color_manual(values = my_cols) +
  theme_tree2() +
  scale_x_continuous(
    trans = "log1p",
    breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 1),
    labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "1"),
    expand = c(0.01, 0.8)
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

# Save the final figure
#ggsave("/Users/annamccormick/R/Feral_Cannabis_EGS/newTree/scratch_round1/final_tree.pdf", 
       #width = 49, height = 14)

################################################
#adding in bootstrap good 95% support branche markers on tree
library(ape)
library(ggtree)
library(dplyr)
library(ggplot2)

# Load the tree
tree <- read.tree("/Users/annamccormick/R/Feral_Cannabis_EGS/newTree/scratch_round1/merged_Aina_Ren_Soorni_n909_13389.min1.phy.varsites.phy.contree")
class(tree)

# Load metadata
metadata <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/newTree/scratch_round1/Ren_Soorni_Aina_metadata_info.csv")
colnames(metadata)

# Rename key columns for clarity
# Assuming columns are: sample.id, FIN_name, [other metadata]
colnames(metadata)[which(colnames(metadata) == "sample.id")] <- "Current"
colnames(metadata)[which(colnames(metadata) == "FIN_name")]  <- "Full"


# Replace tree tip labels with FIN_name, but preserve unmatched tips ("hops")
tree$tip.label <- ifelse(
  tree$tip.label %in% metadata$Current,
  metadata$Full[match(tree$tip.label, metadata$Current)],
  tree$tip.label  # keep labels like "hops"
)

# Assign groups for coloring
metadata <- metadata %>%
  mutate(
    Group = case_when(
      grepl("Aina_cluster_1", Full)        ~ "Aina_cluster_1",
      grepl("Aina_cluster_2", Full)        ~ "Aina_cluster_2",
      grepl("Aina_cluster_3", Full)        ~ "Aina_cluster_3",
      grepl("Aina_cluster_4", Full)        ~ "Aina_cluster_4",
      grepl("Aina_cluster_5", Full)        ~ "Aina_cluster_5",
      grepl("Ren_basal", Full)             ~ "Ren_basal",
      grepl("Ren_drug-type feral", Full)   ~ "Ren_drug-type feral",
      grepl("Ren_drug-type", Full)         ~ "Ren_drug-type",
      grepl("Ren_hemp-type", Full)         ~ "Ren_hemp-type",
      grepl("Soorni_Population_1", Full)   ~ "Soorni_Population_1",
      grepl("Soorni_Population_2", Full)   ~ "Soorni_Population_2",
      TRUE ~ "Other"
    )
  )

# Join tree labels with metadata in FIN_name
tree_data <- data.frame(label = tree$tip.label) %>%
  left_join(metadata %>% distinct(Full, .keep_all = TRUE), by = c("label" = "Full"))

#Ladderize and root the tree 
tree <- ladderize(tree, right = FALSE)
root_sample <- "hops"
tree <- root(tree, outgroup = root_sample, resolve.root = TRUE)

####################################################################
# Define your color palette
my_cols <- c(
  "Aina_cluster_1"        = "darkred",      
  "Aina_cluster_2"        = "darkblue",     
  "Aina_cluster_3"        = "darkgreen",    
  "Aina_cluster_4"        = "purple4",      
  "Aina_cluster_5"        = "darkorange4",  
  "Ren_basal"             = "orange",       
  "Ren_drug-type"         = "red",          
  "Ren_drug-type feral"   = "blue",         
  "Ren_hemp-type"         = "green",        
  "Soorni_Population_1"   = "purple",       
  "Soorni_Population_2"   = "#8DA0CB",
  "Other"                 = "gray50",        # fallback group
  "hops"                  = "black"          # manual assignment if needed
)

# OPTIONAL: manually assign 'hops' group if it’s not in the metadata
tree_data$Group[is.na(tree_data$Group) & tree_data$label == "hops"] <- "hops"


# Extract bootstrap supports (these are stored as `node.label` for internal nodes only)
node_supports <- as.numeric(tree$node.label)

# Create a data frame with node numbers and support values
# Internal nodes are numbered from (Ntip + 1) to (Ntip + Nnode)
support_df <- data.frame(
  node = (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode),
  support = node_supports
)

# Label which nodes are well supported (≥95%)
support_df <- support_df %>%
  mutate(well_supported = support >= 95)



library(ggtree)
# Build base tree
p <- ggtree(tree) %<+% tree_data

# Get ggtree's internal data and merge with support_df
tree_df <- p$data
tree_df <- left_join(tree_df, support_df, by = "node")

# Now plot with internal support overlay
p <- p +
  geom_tiplab(size = 0.01, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = Group), size = 1.8) + 
  geom_point2(
    data = tree_df %>% filter(!isTip & support >= 95),
    aes(subset = !isTip & support >= 95),
    size = 1.5, shape = 21, fill = "pink", color = "pink"
  ) +
  scale_color_manual(values = my_cols) +
  theme_tree2() +
  scale_x_continuous(
    trans = "log1p",
    breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 10, 15),
    labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "1", "10", "15"),
    expand = c(0.2, 0.5)
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

# Print
print(p)

# Save the final figure
#ggsave("/Users/annamccormick/R/Feral_Cannabis_EGS/newTree/scratch_round1/tree_n909_supports_95.pdf", width = 20, height = 49)
ggsave("/Users/annamccormick/R/Feral_Cannabis_EGS/newTree/scratch_round1/tree_n909_supports_95_dimensions.pdf", width = 25, height = 15)

########################
# Bootstrap % checks
########################
library(ape)

# Load the consensus tree
tree <- read.tree("/Users/annamccormick/R/Feral_Cannabis_EGS/newTree/scratch_round1/merged_Aina_Ren_Soorni_n909_13389.min1.phy.varsites.phy.contree")

# Extract node support values
supports <- as.numeric(tree$node.label)

# Summary
summary(supports)

# Count how many nodes have support ≥ 95%
high_support <- sum(supports >= 95, na.rm = TRUE)
total_nodes  <- length(supports)
percent_high_support <- 100 * high_support / total_nodes

cat(sprintf("Nodes with ≥95%% support: %d out of %d (%.1f%%)\n",
            high_support, total_nodes, percent_high_support))

##############
######circular TREE
library(ggtree)
library(dplyr)
library(ggplot2)
library(ggrepel)

# ---- 1) Make plotting tree and shorten hops terminal branch ----
tree_plot <- tree
hop_tip <- which(tree_plot$tip.label == "hops")
hop_edge_idx <- which(tree_plot$edge[, 2] == hop_tip)
tree_plot$edge.length[hop_edge_idx] <- tree_plot$edge.length[hop_edge_idx] * 0.2

# ---- 2) Build base plot + attach metadata ----
p <- ggtree(tree_plot, layout = "circular") %<+% tree_data +
  geom_tree(linewidth = 0.2)

# ---- 3) Parse node supports robustly ----
node_supports <- suppressWarnings(
  as.numeric(gsub("^\\s*([0-9.]+).*$", "\\1", tree_plot$node.label))
)

support_df <- data.frame(
  node = (Ntip(tree_plot) + 1):(Ntip(tree_plot) + tree_plot$Nnode),
  support = node_supports
)

p$data <- left_join(p$data, support_df, by = "node")

# ---- 4) Decide threshold depending on scale ----
thr <- if (max(p$data$support, na.rm = TRUE) <= 1.01) 0.95 else 95

# Optional: choose how to format labels (probabilities vs percent)
p$data <- p$data %>%
  mutate(
    support_label = case_when(
      is.na(support) ~ NA_character_,
      thr < 1.01     ~ sprintf("%.2f", support),         # e.g., 0.97
      TRUE          ~ sprintf("%d", round(support))      # e.g., 97
    )
  )

# ---- 5) Plot: tips + pink dots + numeric labels for high supports ----
p_circ <- p +
  geom_tippoint(aes(color = Group), size = 0.9, alpha = 0.9) +
  
  # Pink markers for high support internal nodes
  geom_point2(
    aes(subset = !isTip & !is.na(support) & support >= thr),
    size = 0.9, shape = 21, fill = "pink", color = "pink", stroke = 0.2
  ) +
  
  # Numeric labels for the same high support nodes
  geom_label_repel(
    data = p$data %>% filter(!isTip, !is.na(support), support >= thr),
    aes(x = x, y = y, label = support_label),
    size = 2.2,
    label.size = 0.15,
    fill = "white",
    min.segment.length = 0,
    segment.size = 0.2,
    max.overlaps = 200
  ) +
  
  scale_color_manual(values = my_cols, drop = FALSE) +
  theme_void() +
  theme(legend.position = "right")

p_circ

library(ggtree)
library(dplyr)
library(ggplot2)
library(ggrepel)

# after you've built p and joined support into p$data (as you already do)
# and after thr is defined...

# choose what "main divisions" means: inner 12% of root-to-tip distance
x_cut <- quantile(p$data$x[!p$data$isTip], 0.005, na.rm = TRUE)

label_df <- p$data %>%
  filter(!isTip, !is.na(support), support >= thr, x <= x_cut) %>%
  mutate(lbl = if (thr < 1.01) sprintf("%.2f", support) else sprintf("%d", round(support)))

p_circ_main <- p +
  geom_tippoint(aes(color = Group), size = 1.5, alpha = 0.9) + #leaf tip circle size
  
  # keep pink dots for all >=95 if you want, OR restrict them similarly
  geom_point2(
    aes(subset = !isTip & !is.na(support) & support >= thr),
    size = 2.0, shape = 21, fill = "pink", color = "pink", stroke = 0.2
  ) +
  
  # print only the main-division supports
  geom_label_repel(
    data = label_df,
    aes(x = x, y = y, label = lbl),
    size = 3,
    label.size = 0.15,
    fill = "white",
    nudge_x = -0.6,
    segment.size = 0.2,
    min.segment.length = 0,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = my_cols, drop = FALSE) +
  theme_void() +
  theme(legend.position = "right")

p_circ_main


ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/tree1a.pdf",
  plot     = p_circ_main,
  width    = 15,
  height   = 15,
  units    = "in",
  device   = "pdf"
)

########################################################################################################################################################################
# Figure 3
########################################################################################################################################################################

####################
# Figure 3A
####################

#######################################################################################################################################
#only 9 above 0.5
#######################################################################################################################################
# Load data
gebv_df <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/newAina_Ren_Soorni_EGS/rrblup_GEAV_values_Ren_Soorni_Aina_new.csv")
str(gebv_df)

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# Optional: pivot data if needed (only if not already in long format)
gebv_long <- gebv_df %>%
  pivot_longer(cols = starts_with("bio_"), names_to = "bio_variable", values_to = "GEBV")

# Convert bio_names to a named vector
bio_names <- c(
  bio_01 = "Annual Mean Temperature (bio 01)",
  bio_02 = "Mean Diurnal Range (bio 02)",
  bio_03 = "Isothermality (bio 03)",
  bio_04 = "Temperature Seasonality (bio 04)",
  bio_05 = "Max Temperature of Warmest Month (bio 05)",
  bio_06 = "Min Temperature of Coldest Month (bio 06)",
  bio_07 = "Temperature Annual Range (bio 07)",
  bio_08 = "Mean Temperature of Wettest Quarter (bio 08)",
  bio_09 = "Mean Temperature of Driest Quarter (bio 09)",
  bio_10 = "Mean Temperature of Warmest Quarter (bio 10)",
  bio_11 = "Mean Temperature of Coldest Quarter (bio 11)",
  bio_12 = "Annual Precipitation (bio 12)",
  bio_13 = "Precipitation of Wettest Month (bio 13)",
  bio_14 = "Precipitation of Driest Month (bio 14)",
  bio_15 = "Precipitation Seasonality (bio 15)",
  bio_16 = "Precipitation of Wettest Quarter (bio 16)",
  bio_17 = "Precipitation of Driest Quarter (bio 17)",
  bio_18 = "Precipitation of Warmest Quarter (bio 18)",
  bio_19 = "Precipitation of Coldest Quarter (bio 19)"
)

# Custom colors
my_cols2 <- c(
  "Aina_cluster_1"       = "darkred",
  "Aina_cluster_2"       = "darkblue",
  "Aina_cluster_3"       = "darkgreen",
  "Aina_cluster_4"       = "purple4",
  "Aina_cluster_5"       = "darkorange4",
  "Ren_basal"            = "orange",
  "Ren_drug-type"        = "red",
  "Ren_drug-type feral"  = "blue",
  "Ren_hemp-type"        = "green",
  "Soorni_Population_1"  = "purple",
  "Soorni_Population_2"  = "#8DA0CB"
)

# Ensure ordering of BIO variables
bio_levels <- paste0("bio_", str_pad(1:19, 2, pad = "0"))

# Define variables you want to keep
selected_bios <- c("bio_02", "bio_05", "bio_08", "bio_09", "bio_12", "bio_14","bio_17", "bio_18", "bio_19")

# Reorder and label
gebv_long <- gebv_long %>%
  filter(bio_variable %in% selected_bios) %>%
  mutate(
    bio_variable = factor(bio_variable, levels = bio_levels),
    bio_label = factor(bio_names[bio_variable], levels = bio_names[bio_levels])
  )

# Plot
p<- ggplot(gebv_long, aes(x = Population, y = GEBV, fill = Population)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.4) +
  scale_fill_manual(values = my_cols2) +
  facet_wrap(~ bio_label, scales = "free_y", ncol = 3) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "none"
  ) +
  labs(
    title = "",
    x = "Population",
    y = "GEAV"
  )

ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/Fig2B.pdf",
  plot = p,
  width = 15,
  height = 10,
  units = "in"
)




####################
# Figure 3B
####################
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(forcats)
library(scales)

# -----------------------------
# 0) INPUT FILE
# -----------------------------
# Use your local path, OR (if you want) the uploaded path:
# file_path <- "/mnt/data/FigureS6_koppen_info.xlsx"
file_path <- "/Users/annamccormick/R/Feral_Cannabis_EGS/matching/FigureS6_koppen_info.xlsx"

# -----------------------------
# 1) Köppen–Geiger colours (your mapping)
# IMPORTANT: names must be characters "1","2",... not numeric
# -----------------------------
kg_colors <- c(
  "1"  = "#0000FF",  # Af
  "2"  = "#0078FF",  # Am
  "3"  = "#46A6FF",  # Aw
  
  "4"  = "#FF0000",  # BWh
  "5"  = "#FF9696",  # BWk
  "6"  = "#F5A300",  # BSh
  "7"  = "#FFDC64",  # BSk
  
  "8"  = "#FFFF00",  # Csa
  "9"  = "#C8FF00",  # Csb
  "10" = "#96FF00",  # Csc
  "11" = "#00FF00",  # Cwa
  "12" = "#00C800",  # Cwb
  "13" = "#009600",  # Cwc
  "14" = "#96FF96",  # Cfa
  "15" = "#64C864",  # Cfb
  "16" = "#329632",  # Cfc
  
  "17" = "#00FFFF",  # Dsa
  "18" = "#37C8FF",  # Dsb
  "19" = "#007DFF",  # Dsc
  "20" = "#0050C8",  # Dsd
  "21" = "#B4E6FF",  # Dwa
  "22" = "#64B4FF",  # Dwb
  "23" = "#3296FF",  # Dwc
  "24" = "#0064C8",  # Dwd
  "25" = "#C8C8FF",  # Dfa
  "26" = "#9696FF",  # Dfb
  "27" = "#6464FF",  # Dfc
  "28" = "#3232C8",  # Dfd
  
  "29" = "#BEBEBE",  # ET
  "30" = "#808080"   # EF
)

# -----------------------------
# 2) Full legend labels (your mapping)
# -----------------------------
kg_labels <- c(
  "1"  = "Af - Tropical, rainforest",
  "2"  = "Am - Tropical, monsoon",
  "3"  = "Aw - Tropical, savanna",
  "4"  = "BWh - Arid, desert, hot",
  "5"  = "BWk - Arid, desert, cold",
  "6"  = "BSh - Arid, steppe, hot",
  "7"  = "BSk - Arid, steppe, cold",
  "8"  = "Csa - Temperate, dry summer, hot summer",
  "9"  = "Csb - Temperate, dry summer, warm summer",
  "10" = "Csc - Temperate, dry summer, cold summer",
  "11" = "Cwa - Temperate, dry winter, hot summer",
  "12" = "Cwb - Temperate, dry winter, warm summer",
  "13" = "Cwc - Temperate, dry winter, cold summer",
  "14" = "Cfa - Temperate, no dry season, hot summer",
  "15" = "Cfb - Temperate, no dry season, warm summer",
  "16" = "Cfc - Temperate, no dry season, cold summer",
  "17" = "Dsa - Cold, dry summer, hot summer",
  "18" = "Dsb - Cold, dry summer, warm summer",
  "19" = "Dsc - Cold, dry summer, cold summer",
  "20" = "Dsd - Cold, dry summer, very cold winter",
  "21" = "Dwa - Cold, dry winter, hot summer",
  "22" = "Dwb - Cold, dry winter, warm summer",
  "23" = "Dwc - Cold, dry winter, cold summer",
  "24" = "Dwd - Cold, dry winter, very cold winter",
  "25" = "Dfa - Cold, no dry season, hot summer",
  "26" = "Dfb - Cold, no dry season, warm summer",
  "27" = "Dfc - Cold, no dry season, cold summer",
  "28" = "Dfd - Cold, no dry season, very cold winter",
  "29" = "ET - Polar, tundra",
  "30" = "EF - Polar, frost"
)

# -----------------------------
# 3) Read + standardise columns
# -----------------------------
df <- read_excel(file_path) %>%
  rename_with(str_trim)
df

df2 <- df %>%
  dplyr::rename(
    Present = `Present Climate Class`,
    Future  = `Future Climate Class`
  ) %>%
  dplyr::mutate(
    Population = stringr::str_trim(as.character(Population)),
    Present    = stringr::str_trim(as.character(Present)),
    Future     = stringr::str_trim(as.character(Future))
  )

# -----------------------------
# 4) Long format + keep only valid KG codes
# -----------------------------
long <- df2 %>%
  select(Population, Present, Future) %>%
  pivot_longer(
    cols = c(Present, Future),
    names_to = "Scenario",
    values_to = "KG_class"
  ) %>%
  filter(!is.na(Population), !is.na(KG_class), KG_class != "") %>%
  mutate(
    KG_class = as.character(KG_class),
    Scenario = factor(Scenario, levels = c("Present", "Future"))
  )

# -----------------------------
# 5) Drop legend entries that do not occur in the data
#    (THIS removes unused KG classes from the legend)
# -----------------------------
used_classes <- sort(unique(long$KG_class))

kg_colors_used <- kg_colors[used_classes]
kg_labels_used <- kg_labels[used_classes]

# -----------------------------
# 6) Force Population order (Cluster_1, Cluster_2, ...)
#    Option A: explicit manual order (recommended for papers)
# -----------------------------
pop_order <- c(
  "Cluster_1", "Cluster_2", "Cluster_3", "Cluster_4", "Cluster_5",
  "Ren_basal", "Ren_drug-type feral", "Ren_hemp-type",
  "Soorni_Population_1", "Soorni_Population_2"
)

# Keep anything not listed (if present) and append it to the end:
pop_order_final <- c(pop_order, setdiff(unique(long$Population), pop_order))

long <- long %>%
  mutate(Population = factor(Population, levels = pop_order_final))

# -----------------------------
# 7) Plot: proportions
# -----------------------------
p_prop <- ggplot(long, aes(x = Population, fill = KG_class)) +
  geom_bar(position = "fill", width = 0.9, color = NA) +
  facet_wrap(~Scenario, ncol = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(
    name   = "Köppen–Geiger climate class",
    values = kg_colors_used,
    breaks = used_classes,      # legend order = only used classes
    labels = kg_labels_used,
    drop   = TRUE
  ) +
  labs(x = "Population", y = "Percent of samples") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

p_prop


# -----------------------------
# 6b) Pretty display labels for populations
# -----------------------------
pop_labels <- c(
  "Cluster_1" = "Aina Cluster 1",
  "Cluster_2" = "Aina Cluster 2",
  "Cluster_3" = "Aina Cluster 3",
  "Cluster_4" = "Aina Cluster 4",
  "Cluster_5" = "Aina Cluster 5",
  
  "Ren_basal"           = "Ren basal",
  "Ren_drug-type feral" = "Ren drug-type feral",
  "Ren_hemp-type"       = "Ren hemp-type",
  
  "Soorni_Population_1" = "Soorni Population 1",
  "Soorni_Population_2" = "Soorni Population 2"
)

# Apply labels (any populations not listed will have _ replaced with spaces)
long <- long %>%
  mutate(
    Population_label = ifelse(
      Population %in% names(pop_labels),
      pop_labels[as.character(Population)],
      str_replace_all(as.character(Population), "_", " ")
    )
  )

p_prop <- ggplot(long, aes(x = Population_label, fill = KG_class)) +
  geom_bar(position = "fill", width = 0.9, color = NA) +
  facet_wrap(~Scenario, ncol = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(
    name   = "Köppen–Geiger climate class",
    values = kg_colors_used,
    breaks = used_classes,
    labels = kg_labels_used,
    drop   = TRUE
  ) +
  labs(x = "Population", y = "Percent of samples") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )
p_prop

ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/FigS6_KoppenG_proportions.pdf",
  plot     = p_prop,
  width    = 8,
  height   = 6,
  units    = "in"
)


########################################################################################################################################################################
# Figure 4 (Supplemental Figure 7)
########################################################################################################################################################################
library(vcfR)
library(dplyr)

# read VCF
vcf <- read.vcfR("/Users/annamccormick/R/Feral_Cannabis_EGS/GEA/merged_on_shared.vcf.gz")

# get genotypes
gt <- extract.gt(vcf)  # matrix samples x SNPs

# load popmap
popmap <- read.table("/Users/annamccormick/R/Feral_Cannabis_EGS/GEA/popmap.txt", header=FALSE, col.names=c("sample","pop"))
pop_order <- unique(popmap$pop)


# function to convert genotypes to allele counts
parse_ref_alt_counts <- function(g) {
  if (is.null(g) || length(g) == 0) return(c(REF = 0L, ALT = 0L))
  # drop missing
  g <- g[!is.na(g)]
  g <- g[g != "./." & g != ".|."]
  if (length(g) == 0) return(c(REF = 0L, ALT = 0L))
  
  # split "a/b" or "a|b" -> c("a","b")
  alle <- strsplit(g, "[/|]", perl = TRUE)
  alle <- unlist(alle, use.names = FALSE)
  
  # keep only 0/1 alleles (ignore "2","3" etc. if present)
  alle <- alle[alle %in% c("0","1")]
  if (length(alle) == 0) return(c(REF = 0L, ALT = 0L))
  
  c(REF = as.integer(sum(alle == "0")),
    ALT = as.integer(sum(alle == "1")))
}

# 2) For ONE SNP, aggregate by population and return a data.frame with pop, REF, ALT
geno2counts_safe <- function(geno_col, samples, popmap) {
  df <- data.frame(sample = samples, gt = geno_col, stringsAsFactors = FALSE)
  df <- merge(df, popmap, by = "sample")   # attach pop to each sample genotype
  
  # compute counts per pop in the fixed order; fill zeros if a pop has no samples at this SNP
  out_list <- lapply(pop_order, function(p) {
    g <- df$gt[df$pop == p]
    c(pop = p, parse_ref_alt_counts(g))
  })
  out <- do.call(rbind, out_list)
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  out$REF <- as.integer(out$REF); out$ALT <- as.integer(out$ALT)
  out
}


# correct the sample vector
samples_vec <- colnames(gt)

# apply over ROWS (margin = 1), one SNP at a time
allele_counts <- apply(gt, 1, geno2counts_safe, samples = samples_vec, popmap = popmap)


library(purrr)
# allele_counts is a list, one element per SNP
# Convert to a big matrix
pop_order <- unique(popmap$pop)

combine_counts <- function(df) {
  # reorder by pop order
  df <- df[match(pop_order, df$pop), ]
  as.vector(t(df[, c("REF","ALT")]))
}

geno_mat <- map_dfr(allele_counts, combine_counts)

# Write to file
# if rows == 2 * pops and cols == #SNPs, transpose to SNPs x (2*pops)
if (nrow(geno_mat) == 2L*length(pop_order)) {
  geno_mat <- t(as.matrix(geno_mat))
}

# sanity checks
stopifnot(ncol(geno_mat) == 2L * length(pop_order))   ######### 22 for 11 pops for REF and ALT for each population #########
cat("SNPs (rows):", nrow(geno_mat), "  Columns:", ncol(geno_mat), "\n")
stopifnot(!anyNA(geno_mat))

# write BayPass allele-count file
write.table(
  geno_mat,
  file = "/Users/annamccormick/R/Feral_Cannabis_EGS/GEA/geno_counts.txt",
  quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE
)


######################
#Environmental data
######################

#popmap <- read.table("/Users/annamccormick/R/Feral_Cannabis_EGS/GEA/popmap.txt", header=FALSE, col.names=c("sample","pop"))
pop_order <- unique(popmap$pop)

env_df <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/GEA/n909_climate_data_all_copy.csv",
                   header=TRUE, stringsAsFactors=FALSE)

# Ensure exact same pops and order:
env_df <- env_df[match(pop_order, env_df$pop), ]
stopifnot(all(env_df$pop == pop_order))

# Drop the 'pop' column and write numeric matrix:
env_mat <- as.matrix(env_df[ , setdiff(names(env_df), "pop"), drop=FALSE])

# Optional: scale covariates yourself; BayPass by default scales unless you pass -nocovscaling
env_mat <- scale(env_mat)

write.table(
  env_mat,
  file = "/Users/annamccormick/R/Feral_Cannabis_EGS/GEA/env.txt",
  quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE
)

######################
#Remove the Ren_Drug_Type population
######################

# reload your geno matrix
geno_mat <- read.table("/Users/annamccormick/R/Feral_Cannabis_EGS/GEA/geno_counts.txt")

# which population is missing?
drop_pop <- "Ren_drug_type"

# reconstruct the order you used originally
pop_order <- unique(read.table("/Users/annamccormick/R/Feral_Cannabis_EGS/GEA/popmap.txt",
                               header=FALSE, col.names=c("sample","pop"))$pop)

# figure out which columns belong to that pop
drop_idx <- which(pop_order == drop_pop)
cols_to_drop <- (2*drop_idx - 1):(2*drop_idx)

# drop them
geno_mat_reduced <- geno_mat[ , -cols_to_drop ]

# save reduced counts
write.table(geno_mat_reduced,
            file="/Users/annamccormick/R/Feral_Cannabis_EGS/GEA/geno_counts_reduced.txt",
            quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE)

# also update pop_order
pop_order_reduced <- pop_order[-drop_idx]


######################
#Remove from ENV file also
######################
# load your environment file
env_mat <- as.matrix(read.table("/Users/annamccormick/R/Feral_Cannabis_EGS/GEA/env.txt"))

# drop the matching row
env_mat_reduced <- env_mat[-drop_idx, ]

# save new env file
write.table(env_mat_reduced,
            file="/Users/annamccormick/R/Feral_Cannabis_EGS/GEA/env_reduced.txt",
            quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE)

######################
#need to transpose env matrix and rewrite
######################

# pop order (reduced)
popmap <- read.table("/Users/annamccormick/R/Feral_Cannabis_EGS/GEA/popmap.txt",
                     header = FALSE, col.names = c("sample","pop"), sep = "")
pop_order <- unique(popmap$pop)
drop_pop <- "Ren_drug_type"
pop_order_reduced <- setdiff(pop_order, drop_pop)

# read your per-pop CSV (has a 'pop' column + 19 covariate columns)
env_df <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/GEA/n909_climate_data_all_copy.csv",
                   stringsAsFactors = FALSE)

# keep only reduced pops and enforce the exact order
env_df <- env_df[match(pop_order_reduced, env_df$pop), ]
stopifnot(all(env_df$pop == pop_order_reduced))

# take only covariate columns
env_mat <- as.matrix(env_df[ , setdiff(names(env_df), "pop"), drop = FALSE])

# handle any remaining NAs (simple mean-impute; replace if you prefer)
for (j in seq_len(ncol(env_mat))) {
  if (anyNA(env_mat[, j])) env_mat[is.na(env_mat[, j]), j] <- mean(env_mat[, j], na.rm = TRUE)
}


# TRANSPOSE so rows = covariates, COLUMNS = populations (this is the key change)
env_mat_baypass <- t(env_mat)

# sanity checks: columns MUST equal #pops (10), no NAs
stopifnot(ncol(env_mat_baypass) == length(pop_order_reduced))
stopifnot(!anyNA(env_mat_baypass))

# write the BayPass covariate file (no headers)
write.table(env_mat_baypass,
            file = "/Users/annamccormick/R/Feral_Cannabis_EGS/GEA/env_reduced.txt",
            quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)


#############################
#extra needed for WZA - Build per-SNP metadata (CSV)
#############################
library(vcfR) 
library(dplyr) 
library(readr)

vcf <- read.vcfR("/Users/annamccormick/R/Feral_Cannabis_EGS/GEA/merged_on_shared.vcf.gz")
fix <- as.data.frame(getFIX(vcf)) |>
  transmute(chr = as.character(CHROM), pos = as.integer(as.character(POS)))

# 50 kb windows per chromosome
w <- 50000L
meta <- fix |>
  group_by(chr) |>
  mutate(window = paste(chr, floor((pos-1L)/w), sep=":")) |>
  ungroup()

write_csv(meta, "/Users/annamccormick/R/Feral_Cannabis_EGS/GEA/alleleFrequencyMetadata.csv")


###########################
# Make one env file per covariate (single line, 10 numbers)
###########################
#correlationsFromBaypass.py only reads the first line of the env file (one covariate). 
#Create 19 one-line files, one per WorldClim variable, with pops in the same order you used for the counts:

# env_per_pop_df: columns = pop, BIO1, BIO2, ... BIO19 in the exact pop order you used
# After dropping the one pop with no env:
# pop order (reduced)
popmap <- read.table("/Users/annamccormick/R/Feral_Cannabis_EGS/GEA/popmap.txt",
                     header = FALSE, col.names = c("sample","pop"), sep = "")
pop_order <- unique(popmap$pop)
drop_pop <- "Ren_drug_type"
pop_order_reduced <- setdiff(pop_order, drop_pop)
writeLines(pop_order_reduced, "/Users/annamccormick/R/Feral_Cannabis_EGS/GEA/pop_order_reduced.txt")

pops <- scan("/Users/annamccormick/R/Feral_Cannabis_EGS/GEA/pop_order_reduced.txt", what="character")  # or build from your R object

####################
# make one-line env files from env_reduced.txt (rows=covariates, cols=pops) ---
gea_dir <- "/Users/annamccormick/R/Feral_Cannabis_EGS/GEA"
E <- as.matrix(read.table(file.path(gea_dir, "env_reduced.txt")))
stopifnot(ncol(E) == 10)    # 10 pops after dropping one

# Name your 19 traits however you like:
trait_names <- sprintf("BIO%02d", seq_len(nrow(E)))  # BIO01..BIO19
for (i in seq_len(nrow(E))) {
  line <- paste(E[i, ], collapse = " ")
  writeLines(line, file.path(gea_dir, paste0("env_", trait_names[i], ".txt")))
}
####################

# 1) metadata rows == SNP rows in counts
meta <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/WZA/alleleFrequencyMetadata.csv")
n_snps_counts <- length(readLines("/Users/annamccormick/R/Feral_Cannabis_EGS/WZA/geno_counts_reduced.txt"))
stopifnot(nrow(meta) == n_snps_counts)

# 2) window column exists
stopifnot("window" %in% names(meta))


#############
#run in terminal 
# run Kendall's tau per SNP for, say, BIO01
python3 correlationsFromBaypass.py \
--geno alleleFrequencyData.txt.gz \
--env env_BIO01.txt \
--meta alleleFrequencyMetadata.csv \
--output BIO01_singleSNP_GEA.tsv

# aggregate with WZA over your 'window' column
python3 general_WZA_script.py \
--correlations BIO01_singleSNP_GEA.tsv \
--summary_stat K_tau_p \
--window window \
--output BIO01_wza.csv \
--retain chr pos

###########
#PLOTTING RESULTS
###########
library(readr)
library(dplyr)
library(ggplot2)

wza <- read_csv("/Users/annamccormick/R/Feral_Cannabis_EGS/WZA/BIO01_wza.csv")

# Require Z_pVal column from WZA output
stopifnot("Z_pVal" %in% names(wza))

# Clean + compute -log10 p strictly from Z_pVal
wza <- wza |>
  mutate(
    Z_pVal = as.numeric(Z_pVal),
    Z_pVal = pmin(pmax(Z_pVal, .Machine$double.xmin), 1),  # clamp to avoid Inf
    ml10p  = -log10(Z_pVal)
  ) |>
  filter(is.finite(ml10p))

# Plot using only Z_pVal
library(ggplot2)
ggplot(wza, aes(pos/1e6, ml10p)) +
  geom_point(size = 1.2, alpha = 0.8) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "red", linewidth = 0.6) +
  facet_wrap(~ chr, scales = "free_x") +
  labs(x = "Position (Mbp)", y = expression(-log[10](Z[pVal]))) +
  theme_classic()

###########
# All results
###########
library(readr)
library(dplyr)
library(ggplot2)
library(purrr)

wza_dir <- "/Users/annamccormick/R/Feral_Cannabis_EGS/WZA_all/WZA_out"  # <-- adjust if needed
files <- list.files(wza_dir, pattern = "_wza\\.csv$", full.names = TRUE)

# read & stack; tag each with its trait (BIO01..BIO19)
df <- map_dfr(files, ~{
  tr <- sub("_wza\\.csv$", "", basename(.x))     # e.g., "BIO01"
  read_csv(.x, show_col_types = FALSE) |> mutate(trait = tr)
})

# clean: use only Z_pVal
stopifnot("Z_pVal" %in% names(df))
df <- df |>
  mutate(
    Z_pVal = pmin(pmax(as.numeric(Z_pVal), .Machine$double.xmin), 1),
    ml10p  = -log10(Z_pVal),
    chr    = as.character(chr),
    pos    = as.numeric(pos)
  ) |>
  filter(is.finite(ml10p))

# nice ordering for traits & chromosomes
df$trait <- factor(df$trait, levels = sprintf("BIO%02d", 1:19), ordered = TRUE)
chr_levels <- unique(df$chr)[order(as.numeric(gsub("\\D+", "", unique(df$chr))), unique(df$chr))]
df$chr <- factor(df$chr, levels = chr_levels)

# write a multi-page PDF: 1 page per trait, each faceted by chromosome
pdf(file.path(wza_dir, "WZA_byChr_ALL_traits.pdf"), width = 10, height = 6)
for (tr in levels(df$trait)) {
  p <- ggplot(filter(df, trait == tr), aes(pos/1e6, ml10p)) +
    geom_point(size = 0.9, alpha = 0.8) +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "red", linewidth = 0.5) +
    facet_wrap(~ chr, scales = "free_x") +
    labs(title = paste(tr, "WZA (−log10 Z_pVal)"),
         x = "Position (Mbp)", y = expression(-log[10](Z[pVal]))) +
    theme_classic(base_size = 11)
  print(p)
}
dev.off()

##############
#  PLOT
##############
library(readr)
library(dplyr)
library(purrr)
library(ggplot2)
library(grid)   # for unit()

# 1) Load & stack all WZA outputs
wza_dir <- "/Users/annamccormick/R/Feral_Cannabis_EGS/WZA_all/WZA_out"  # adjust if needed
files   <- list.files(wza_dir, pattern = "_wza\\.csv$", full.names = TRUE)

df <- map_dfr(files, ~{
  tr <- sub("_wza\\.csv$", "", basename(.x))     # e.g., BIO01
  read_csv(.x, show_col_types = FALSE) |> mutate(trait = tr)
})

# 2) Use Z_pVal only → make -log10
stopifnot("Z_pVal" %in% names(df))
df <- df |>
  mutate(
    Z_pVal = pmin(pmax(as.numeric(Z_pVal), .Machine$double.xmin), 1),
    ml10p  = -log10(Z_pVal),
    chr    = as.character(chr),
    pos    = as.numeric(pos)
  ) |>
  filter(is.finite(ml10p))

# 3) Order traits BIO01..BIO19 and chromosomes 1..10 (natural order)
trait_levels <- sprintf("BIO%02d", 1:19)
df$trait <- factor(df$trait,
                   levels = intersect(trait_levels, unique(df$trait)),
                   ordered = TRUE)

chr_levels <- df |>
  distinct(chr) |>
  mutate(n = suppressWarnings(as.numeric(gsub("[^0-9]", "", chr)))) |>
  arrange(n, chr) |>
  pull(chr)
df$chr <- factor(df$chr, levels = chr_levels, ordered = TRUE)


p<- ggplot(df, aes(pos/1e6, ml10p)) +
  geom_point(size=0.35, alpha=0.75) +
  geom_hline(yintercept=-log10(0.01), linetype="dashed", color="red", linewidth=0.4) +
  facet_grid(trait ~ chr, scales = "free_y") +   # <-- free y per panel
  labs(x="Position (Mbp)", y=expression(-log[10](Z[pVal]))) +
  theme_classic(base_size = 9)

p
ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/FigS7.pdf",
  plot = p,
  width = 15,
  height = 10,
  units = "in"
)

############################################
#counting number of snps past threshold

snps_bio1_11 <- df %>%
  filter(
    trait %in% sprintf("BIO%02d", 1:11),
    ml10p >= thr
  ) %>%
  distinct(trait, chr, pos) %>%   # ensures unique SNPs
  nrow()

snps_bio1_11


snps_bio1_11_unique_loci <- df %>%
  filter(trait %in% sprintf("BIO%02d", 1:11), ml10p >= thr) %>%
  distinct(chr, pos) %>%
  nrow()

snps_bio1_11_unique_loci



snps_bio12_19 <- df %>%
  filter(
    trait %in% sprintf("BIO%02d", 12:19),
    ml10p >= thr
  ) %>%
  distinct(trait, chr, pos) %>%
  nrow()

snps_bio12_19


snps_bio12_19_unique_loci <- df %>%
  filter(trait %in% sprintf("BIO%02d", 12:19),
         ml10p >= thr) %>%
  distinct(chr, pos) %>%
  nrow()

snps_bio12_19_unique_loci


##########################################################################################
#red dots for bios appearing 3 times or more and above line for TEMPERATURE
##########################################################################################
thr <- -log10(0.01)  # p ≤ 0.01

# Keep only temperature traits
df_temp <- df %>%
  filter(trait %in% sprintf("BIO%02d", 1:11))

# Count in how many temp traits each locus passes threshold
temp_replication <- df_temp %>%
  filter(ml10p >= thr) %>%
  distinct(chr, pos, trait) %>%
  group_by(chr, pos) %>%
  summarise(
    n_temp_traits = n(),
    traits = paste(trait, collapse = ", "),
    .groups = "drop"
  )

# Loci that replicate in ≥3 temperature variables
rep3 <- temp_replication %>%
  filter(n_temp_traits >= 3) %>%
  select(chr, pos, n_temp_traits)

df_temp2 <- df_temp %>%
  left_join(rep3, by = c("chr", "pos")) %>%
  mutate(
    is_rep3 = !is.na(n_temp_traits),
    highlight = is_rep3 & (ml10p >= thr)
  )

p_temp_rep3 <- ggplot(df_temp2, aes(pos/1e6, ml10p)) +
  geom_point(size = 0.35, alpha = 0.6) +
  geom_hline(yintercept = thr, linetype = "dashed", color = "red", linewidth = 0.4) +
  geom_point(
    data = df_temp2 %>% filter(highlight),
    size = 0.7, alpha = 0.9, color = "red"
  ) +
  facet_grid(trait ~ chr, scales = "free_y") +
  labs(
    x = "Position (Mbp)",
    y = expression(-log[10](Z[pVal])),
    title = "BIO01–BIO11: loci replicated across ≥3 BIO variables (highlighted when significant)"
  ) +
  theme_classic(base_size = 9)

p_temp_rep3



ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/BIO01_BIO11_WZA_manhattan_red.pdf",
  plot     = p_temp_rep3,
  width    = 10,
  height   = 8,
  units    = "in"
)


######################################################
# only plot bios wanted for TEMPERATURE
######################################################
# Keep only selected temperature variables
df_temp_sel <- df_temp2 %>%
  filter(trait %in% c("BIO01", "BIO06", "BIO09", "BIO11"))

p_temp_rep3_sel <- ggplot(df_temp_sel, aes(pos/1e6, ml10p)) +
  geom_point(size = 0.35, alpha = 0.6) +
  geom_hline(
    yintercept = thr,
    linetype = "dashed",
    color = "red",
    linewidth = 0.4
  ) +
  geom_point(
    data = df_temp_sel %>% filter(highlight),
    size = 0.7,
    alpha = 0.9,
    color = "red"
  ) +
  facet_grid(trait ~ chr, scales = "free_y") +
  labs(
    x = "Position (Mbp)",
    y = expression(-log[10](Z[pVal])),
    title = "BIO01, BIO06, BIO09, BIO11: replicated temperature-associated loci (≥3 traits)"
  ) +
  theme_classic(base_size = 9)

p_temp_rep3_sel

ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/BIO01_BIO11_WZA_manhattan_red_subset.pdf",
  plot     = p_temp_rep3_sel,
  width    = 15,
  height   = 4,
  units    = "in"
)

##########################
# representative 
##########################
# Figure 4A
##########################
# Keep only selected temperature variables

df_temp_sel <- df_temp2 %>%
  filter(trait %in% c("BIO01"))

p_temp_rep3_sel <- ggplot(df_temp_sel, aes(pos/1e6, ml10p)) +
  # background points (exclude highlighted)
  geom_point(
    data = df_temp_sel %>% filter(!highlight),
    size = 0.6, alpha = 1.0, colour = "grey52"
  ) +
  geom_hline(
    yintercept = thr,
    linetype = "dashed",
    color = "red2",
    linewidth = 0.4
  ) +
  
  # halo for highlighted (draw first)
  geom_point(
    data = df_temp_sel %>% filter(highlight),
    size = 3.2, alpha = 1, colour = "white", shape = 16
  ) +
  # red highlighted points (draw on top)
  geom_point(
    data = df_temp_sel %>% filter(highlight),
    size = 2.2, alpha = 1, colour = "red", shape = 16
  ) +
  
  facet_grid(trait ~ chr, scales = "free_y") +
  labs(
    x = "Position (Mbp)",
    y = expression(-log[10](Z[pVal])),
    title = "BIO01, BIO06, BIO09, BIO11: replicated temperature-associated loci (≥3 traits)"
  ) +
  theme_classic(base_size = 9)

p_temp_rep3_sel

ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/BIO01_BIO11_WZA_manhattan_red_subset_bio1.pdf",
  plot     = p_temp_rep3_sel,
  width    = 16,
  height   = 4,
  units    = "in"
)

########################################################################
# Precipitation WZA: BIO12–BIO19
# Highlight loci that are significant (ml10p >= thr) AND replicated across >= 3 precip traits
##############

library(readr)
library(dplyr)
library(purrr)
library(ggplot2)

# ---- settings ----
wza_dir <- "/Users/annamccormick/R/Feral_Cannabis_EGS/WZA_all/WZA_out"
files   <- list.files(wza_dir, pattern = "_wza\\.csv$", full.names = TRUE)

thr <- -log10(0.01)  # WZA significance threshold: p <= 0.01  (ml10p >= 2)

# ---- 1) load & stack all traits ----
df <- map_dfr(files, ~{
  tr <- sub("_wza\\.csv$", "", basename(.x))  # e.g., "BIO01"
  read_csv(.x, show_col_types = FALSE) |> mutate(trait = tr)
})

stopifnot("Z_pVal" %in% names(df))

df <- df |>
  mutate(
    Z_pVal = pmin(pmax(as.numeric(Z_pVal), .Machine$double.xmin), 1),
    ml10p  = -log10(Z_pVal),
    chr    = as.character(chr),
    pos    = as.numeric(pos)
  ) |>
  filter(is.finite(ml10p), !is.na(chr), !is.na(pos))

# ---- 2) order traits/chromosomes (optional but nice) ----
trait_levels <- sprintf("BIO%02d", 1:19)
df$trait <- factor(df$trait,
                   levels = intersect(trait_levels, unique(df$trait)),
                   ordered = TRUE)

chr_levels <- df |>
  distinct(chr) |>
  mutate(n = suppressWarnings(as.numeric(gsub("[^0-9]", "", chr)))) |>
  arrange(n, chr) |>
  pull(chr)
df$chr <- factor(df$chr, levels = chr_levels, ordered = TRUE)

# ---- 3) subset precipitation traits ----
df_precip <- df %>%
  filter(trait %in% sprintf("BIO%02d", 12:19))

# ---- 4) replication count: in how many precip traits is each locus significant? ----
precip_replication <- df_precip %>%
  filter(ml10p >= thr) %>%
  distinct(chr, pos, trait) %>%       # one hit per trait
  group_by(chr, pos) %>%
  summarise(
    n_precip_traits = n(),
    traits = paste(trait, collapse = ", "),
    .groups = "drop"
  )

# loci replicated across >= 3 precipitation variables
rep3_precip <- precip_replication %>%
  filter(n_precip_traits >= 3) %>%
  select(chr, pos, n_precip_traits)

# ---- 5) flag the points to highlight (ONLY if above threshold) ----
df_precip2 <- df_precip %>%
  left_join(rep3_precip, by = c("chr", "pos")) %>%
  mutate(
    is_rep3 = !is.na(n_precip_traits),
    highlight = is_rep3 & (ml10p >= thr)   # ensures nothing below red line is highlighted
  )

# ---- 6) plot: base points black, replicated significant loci blue ----
p_precip_rep3 <- ggplot(df_precip2, aes(pos/1e6, ml10p)) +
  geom_point(size = 0.35, alpha = 0.60) +
  geom_hline(yintercept = thr, linetype = "dashed", color = "red", linewidth = 0.4) +
  geom_point(
    data = df_precip2 %>% filter(highlight),
    size = 0.75, alpha = 0.95, color = "blue"
  ) +
  facet_grid(trait ~ chr, scales = "free_y") +
  labs(
    x = "Position (Mbp)",
    y = expression(-log[10](Z[pVal])),
    title = "BIO12–BIO19: precipitation-associated loci replicated across ≥3 BIO variables (blue = significant & replicated)"
  ) +
  theme_classic(base_size = 9)

p_precip_rep3


ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/BIO12_BIO19_WZA_manhattan_rep3_blue.pdf",
  plot     = p_precip_rep3,
  width    = 10,
  height   = 8,
  units    = "in"
)


library(readr)

write_csv(
  precip_replication %>% filter(n_precip_traits >= 3),
  "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/BIO12_BIO19_rep3_precipitation_GEA_SNPs.csv"
)


############
#Sepcific Bios for ppt 
df_precip_sel <- df_precip2 %>%
  filter(trait %in% c("BIO14", "BIO15", "BIO17"))

p_precip_rep3_sel <- ggplot(df_precip_sel, aes(pos / 1e6, ml10p)) +
  geom_point(size = 0.35, alpha = 0.6) +
  geom_hline(
    yintercept = thr,
    linetype = "dashed",
    color = "red",
    linewidth = 0.4
  ) +
  geom_point(
    data = df_precip_sel %>% filter(highlight),
    size = 0.8,
    alpha = 0.95,
    color = "blue"
  ) +
  facet_grid(trait ~ chr, scales = "free_y") +
  labs(
    x = "Position (Mbp)",
    y = expression(-log[10](Z[pVal])),
    title = "BIO14, BIO15, BIO17: replicated precipitation-associated loci (≥3 traits)"
  ) +
  theme_classic(base_size = 9)

p_precip_rep3_sel

ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/BIO14_BIO15_BIO17_WZA_precipitation_rep3_blue.pdf",
  plot     = p_precip_rep3_sel,
  width    = 15,
  height   = 4,
  units    = "in"
)


##########################
#one representative 
##########################
# Figure 4B
##########################
df_precip_sel <- df_precip2 %>%
  filter(trait %in% c("BIO14"))

p_precip_rep3_sel <- ggplot(df_precip_sel, aes(pos/1e6, ml10p)) +
  geom_point(
    data = df_precip_sel %>% filter(!highlight),
    size = 0.6, alpha = 1.0, colour = "grey52"
  ) +
  geom_hline(yintercept = thr, linetype = "dashed", color = "red", linewidth = 0.4) +
  
  # halo (draw first)
  geom_point(
    data = df_precip_sel %>% filter(highlight),
    size = 3.2, alpha = 1, colour = "white", shape = 16
  ) +
  # colored highlight (draw on top)
  geom_point(
    data = df_precip_sel %>% filter(highlight),
    size = 2.2, alpha = 1, colour = "blue", shape = 16
  ) +
  
  facet_grid(trait ~ chr, scales = "free_y") +
  theme_classic(base_size = 9)
p_precip_rep3_sel

ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/BIO14_BIO15_BIO17_WZA_precipitation_rep3_blue_bio14.pdf",
  plot     = p_precip_rep3_sel,
  width    = 15,
  height   = 4,
  units    = "in"
)
######



########################################################################################################################################################################
# SUPPLEMENTAL FIGURES
########################################################################################################################################################################


########################################################################################################################################################################
# Figure S1
########################################################################################################################################################################
library("FactoMineR")
library("factoextra")

#aina alone
pca <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_VCF_remake/VCF_overlays/PCA_merged.csv")

#aina+ren+soorni
pca <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/newAina_Ren_Soorni_EGS/new_Aina_n760_with_core_flag.csv")


#Select 3rd and 4th principal components
#pca2 <- pca[3:4]
pca2 <- pca[3:6]

head(pca2)
# Standardizing the PCA data before clustering
my_data <- scale(pca2)

#elbow 
a <- fviz_nbclust(my_data, kmeans, method = "wss") + ggtitle("the Elbow Method")
a
ggsave("/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/Aina_Soorni_Ren_elbow.pdf", plot = a, width = 6, height = 4, device = "pdf")

#silhouette
b<- fviz_nbclust(my_data, kmeans, method = "silhouette") + ggtitle("The Silhouette Plot")
b
ggsave("/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/Aina_Soorni_Ren_silhouette.pdf", plot = b, width = 6, height = 4, device = "pdf")

###################
#FIGURE S4B
###################
#fastSTRUCTURE
###################

library(MetBrewer)
library(ggplot2)
library(tidyr)
library(dplyr)

install.packages("ggnewscale")
library(ggnewscale)

##############
#K2
##############
# Load the Hokusai palette
hiroshige_palette <- met.brewer("Hiroshige", 8)
k_colors <- hiroshige_palette[c(1, 6)]  # Colors for K1 and K2

# Load your data
test <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/fastSTRUCTURE_new/merged_cannaK2.csv")

# Sort by K1 and K2 percentages
test <- test %>%
  arrange(desc(K1), desc(K2))

# Convert to long format after sorting
long_data <- pivot_longer(test, cols = c("K1", "K2"), names_to = "category", values_to = "value")

# Reorder the 'Name' factor based on the sorted data
long_data$Name <- factor(long_data$Name, levels = unique(test$Name))

# Assign specific colors to each Type
type_colors <- c(
  "Ren_basal"           = "orange",
  "Ren_hemp-type"       = "green",
  "Ren_drug-type feral" = "blue",
  "Ren_drug-type"       = "red",
  "Soorni_Population_1" = "#9370DB",
  "Soorni_Population_2" = "#008080",
  "Aina_cluster_1"      = "darkred",
  "Aina_cluster_2"      = "darkblue",
  "Aina_cluster_3"      = "darkgreen",
  "Aina_cluster_4"      = "purple4",
  "Aina_cluster_5"      = "darkorange4"
)


# Plot with the colored bar for Type underneath
final_plot <- ggplot() +
  # Plot the K1 and K2 bars with their own fill colors
  geom_bar(data = long_data, aes(x = Name, y = value, fill = category), stat = "identity") +
  scale_fill_manual(name = "K Categories", values = k_colors) +  # Colors for K1 and K2
  
  # Add the colored bars underneath for Type using a separate geom_tile
  new_scale_fill() +  # Necessary to separate the two scales
  geom_tile(data = test, aes(x = Name, y = -0.05, fill = Type), height = 0.1) +
  scale_fill_manual(name = "Type", values = type_colors) +  # Colors for Type
  
  labs(title = "", x = "", y = "Percent Identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 1),  # Make the text smaller
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 6),       # Smaller legend text
    legend.title = element_text(size = 8),      # Smaller legend title
    legend.key.size = unit(0.4, 'cm'),          # Smaller legend keys
    legend.spacing.y = unit(0.2, 'cm')          # Reduce vertical spacing between legend items
  )

print(final_plot)

# Save the plot as a PDF with reduced height
ggsave(filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/fastSTRUCTURE_new/fastSTRUCTURE_plot_k2.pdf", 
       plot = final_plot, 
       width = 14,    # Width of the plot (adjust as needed)
       height = 3)   # Reduced height (adjust as needed)

##############
#K3
##############
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggnewscale)
library(MetBrewer)

hiroshige_palette <- met.brewer("Hiroshige", 8)
k_colors <- hiroshige_palette[c(1, 4, 6)]

# Load your data
test <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/fastSTRUCTURE_new/merged_cannaK3.csv")

# Create a new column for the dominant K category
test <- test %>%
  mutate(dominant_K = case_when(
    K1 >= K2 & K1 >= K3 ~ "K1",
    K2 >= K1 & K2 >= K3 ~ "K2",
    K3 >= K1 & K3 >= K2 ~ "K3"
  ))

# Sort by the dominant K category and then by the actual values in descending order
test <- test %>%
  arrange(dominant_K, desc(K1), desc(K2), desc(K3))

# Convert to long format after sorting
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3"), names_to = "category", values_to = "value")

# Reorder the 'Name' factor based on the sorted data
long_data$Name <- factor(long_data$Name, levels = unique(test$Name))

# Assign specific colors to each Type and K category
type_colors <- c(
  "Ren_basal"           = "orange",
  "Ren_hemp-type"       = "green",
  "Ren_drug-type feral" = "blue",
  "Ren_drug-type"       = "red",
  "Soorni_Population_1" = "#9370DB",
  "Soorni_Population_2" = "#008080",
  "Aina_cluster_1"      = "darkred",
  "Aina_cluster_2"      = "darkblue",
  "Aina_cluster_3"      = "darkgreen",
  "Aina_cluster_4"      = "purple4",
  "Aina_cluster_5"      = "darkorange4"
)

# Hiroshige palette for K categories
k_colors <- hiroshige_palette[c(1, 4, 6)]

# Plot with the colored bar for Type underneath
final_plot <- ggplot() +
  # Plot the K1 and K2 bars with their own fill colors
  geom_bar(data = long_data, aes(x = Name, y = value, fill = category), stat = "identity") +
  scale_fill_manual(name = "K Categories", values = k_colors) +  # Colors for K1 and K2
  
  # Add the colored bars underneath for Type using a separate geom_tile
  new_scale_fill() +  # Necessary to separate the two scales
  geom_tile(data = test, aes(x = Name, y = -0.05, fill = Type), height = 0.1) +
  scale_fill_manual(name = "Type", values = type_colors) +  # Colors for Type
  
  labs(title = "", x = "", y = "Percent Identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 1),  # Make the text smaller
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 6),       # Smaller legend text
    legend.title = element_text(size = 8),      # Smaller legend title
    legend.key.size = unit(0.4, 'cm'),          # Smaller legend keys
    legend.spacing.y = unit(0.2, 'cm')          # Reduce vertical spacing between legend items
  )

print(final_plot)

# Save the plot as a PDF with reduced height
ggsave(filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/fastSTRUCTURE_new/fastSTRUCTURE_plot_k3.pdf", 
       plot = final_plot, 
       width = 14,    # Width of the plot (adjust as needed)
       height = 3)   # Reduced height (adjust as needed)

##############
#K4
##############
#load fastSTRUCTURE data

hiroshige_palette <- met.brewer("Hiroshige", 8)
#selected_colors <- hiroshige_palette[c(1, 4,6, 8)]
k_colors <- hiroshige_palette[c(1,4,6,8)]  # Colors for K

test <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/fastSTRUCTURE_new/merged_cannaK4.csv")

# Create a new column for the dominant K category
test <- test %>%
  mutate(dominant_K = case_when(
    K1 >= K2 & K1 >= K3 & K1 >= K4 ~ "K1",
    K2 >= K1 & K2 >= K3 & K2 >= K4 ~ "K2",
    K3 >= K1 & K3 >= K2 & K3 >= K4 ~ "K3",
    K4 >= K1 & K4 >= K2 & K4 >= K3 ~ "K4"
  ))

# Sort by the dominant K category and then by the actual values in descending order
test <- test %>%
  arrange(dominant_K, desc(K1), desc(K2), desc(K3), desc(K4))

# Convert to long format after sorting
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3", "K4"), names_to = "category", values_to = "value")

# Reorder the 'Name' factor based on the sorted data
long_data$Name <- factor(long_data$Name, levels = unique(test$Name))

# Assign specific colors to each Type and K category
type_colors <- c(
  "Ren_basal"           = "orange",
  "Ren_hemp-type"       = "green",
  "Ren_drug-type feral" = "blue",
  "Ren_drug-type"       = "red",
  "Soorni_Population_1" = "#9370DB",
  "Soorni_Population_2" = "#008080",
  "Aina_cluster_1"      = "darkred",
  "Aina_cluster_2"      = "darkblue",
  "Aina_cluster_3"      = "darkgreen",
  "Aina_cluster_4"      = "purple4",
  "Aina_cluster_5"      = "darkorange4"
)

# Hiroshige palette for K categories (adjust the palette if needed)
k_colors <- hiroshige_palette[c(1, 4, 6, 8)]  # Choose appropriate colors for K1 to K4

# Plot with the colored bar for Type underneath
final_plot <- ggplot() +
  # Plot the K1 and K2 bars with their own fill colors
  geom_bar(data = long_data, aes(x = Name, y = value, fill = category), stat = "identity") +
  scale_fill_manual(name = "K Categories", values = k_colors) +  # Colors for K1 and K2
  
  # Add the colored bars underneath for Type using a separate geom_tile
  new_scale_fill() +  # Necessary to separate the two scales
  geom_tile(data = test, aes(x = Name, y = -0.05, fill = Type), height = 0.1) +
  scale_fill_manual(name = "Type", values = type_colors) +  # Colors for Type
  
  labs(title = "", x = "", y = "Percent Identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 1),  # Make the text smaller
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 6),       # Smaller legend text
    legend.title = element_text(size = 8),      # Smaller legend title
    legend.key.size = unit(0.4, 'cm'),          # Smaller legend keys
    legend.spacing.y = unit(0.2, 'cm')          # Reduce vertical spacing between legend items
  )

print(final_plot)

# Save the plot as a PDF with reduced height
ggsave(filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/fastSTRUCTURE_new//fastSTRUCTURE_plot_k4.pdf", 
       plot = final_plot, 
       width = 14,    # Width of the plot (adjust as needed)
       height = 3)   # Reduced height (adjust as needed)

##############
#K5
##############
hiroshige_palette <- met.brewer("Hiroshige", 8)

test <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/fastSTRUCTURE_new/merged_cannaK5.csv")

# Create a new column for the dominant K category
test <- test %>%
  mutate(dominant_K = case_when(
    K1 >= K2 & K1 >= K3 & K1 >= K4 & K1 >= K5 ~ "K1",
    K2 >= K1 & K2 >= K3 & K2 >= K4 & K2 >= K5 ~ "K2",
    K3 >= K1 & K3 >= K2 & K3 >= K4 & K3 >= K5 ~ "K3",
    K4 >= K1 & K4 >= K2 & K4 >= K3 & K4 >= K5 ~ "K4",
    K5 >= K1 & K5 >= K2 & K5 >= K3 & K5 >= K4 ~ "K5"
  ))

# Sort by the dominant K category and then by the actual values in descending order
test <- test %>%
  arrange(dominant_K, desc(K1), desc(K2), desc(K3), desc(K4), desc(K5))

# Convert to long format after sorting
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3", "K4", "K5"), names_to = "category", values_to = "value")

# Reorder the 'Name' factor based on the sorted data
long_data$Name <- factor(long_data$Name, levels = unique(test$Name))

# Assign specific colors to each Type and K category
type_colors <- c(
  "Ren_basal"           = "orange",
  "Ren_hemp-type"       = "green",
  "Ren_drug-type feral" = "blue",
  "Ren_drug-type"       = "red",
  "Soorni_Population_1" = "#9370DB",
  "Soorni_Population_2" = "#008080",
  "Aina_cluster_1"      = "darkred",
  "Aina_cluster_2"      = "darkblue",
  "Aina_cluster_3"      = "darkgreen",
  "Aina_cluster_4"      = "purple4",
  "Aina_cluster_5"      = "darkorange4"
)
# Hiroshige palette for K categories (adjust the palette if needed)
k_colors <- hiroshige_palette[c(1, 3, 5, 6, 8)]  # Choose appropriate colors for K1 to K5

# Plot with the colored bar for Type underneath
final_plot <- ggplot() +
  # Plot the K1 and K2 bars with their own fill colors
  geom_bar(data = long_data, aes(x = Name, y = value, fill = category), stat = "identity") +
  scale_fill_manual(name = "K Categories", values = k_colors) +  # Colors for K1 and K2
  
  # Add the colored bars underneath for Type using a separate geom_tile
  new_scale_fill() +  # Necessary to separate the two scales
  geom_tile(data = test, aes(x = Name, y = -0.05, fill = Type), height = 0.1) +
  scale_fill_manual(name = "Type", values = type_colors) +  # Colors for Type
  
  labs(title = "", x = "", y = "Percent Identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 1),  # Make the text smaller
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 6),       # Smaller legend text
    legend.title = element_text(size = 8),      # Smaller legend title
    legend.key.size = unit(0.4, 'cm'),          # Smaller legend keys
    legend.spacing.y = unit(0.2, 'cm')          # Reduce vertical spacing between legend items
  )

print(final_plot)

# Save the plot as a PDF with reduced height
ggsave(filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/fastSTRUCTURE_new//fastSTRUCTURE_plot_k5.pdf", 
       plot = final_plot, 
       width = 14,    # Width of the plot (adjust as needed)
       height = 3)   # Reduced height (adjust as needed)

##############
#K6
##############
#load fastSTRUCTURE data

# Load your data
test <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/fastSTRUCTURE_new/merged_cannaK6.csv")

# Create a new column for the dominant K category
test <- test %>%
  mutate(dominant_K = case_when(
    K1 >= K2 & K1 >= K3 & K1 >= K4 & K1 >= K5 & K1 >= K6 ~ "K1",
    K2 >= K1 & K2 >= K3 & K2 >= K4 & K2 >= K5 & K2 >= K6 ~ "K2",
    K3 >= K1 & K3 >= K2 & K3 >= K4 & K3 >= K5 & K3 >= K6 ~ "K3",
    K4 >= K1 & K4 >= K2 & K4 >= K3 & K4 >= K5 & K4 >= K6 ~ "K4",
    K5 >= K1 & K5 >= K2 & K5 >= K3 & K5 >= K4 & K5 >= K6 ~ "K5",
    K6 >= K1 & K6 >= K2 & K6 >= K3 & K6 >= K4 & K6 >= K5 ~ "K6"
  ))

# Sort by the dominant K category and then by the actual values in descending order
test <- test %>%
  arrange(dominant_K, desc(K1), desc(K2), desc(K3), desc(K4), desc(K5), desc(K6))

# Convert to long format after sorting
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3", "K4", "K5", "K6"), names_to = "category", values_to = "value")

# Reorder the 'Name' factor based on the sorted data
long_data$Name <- factor(long_data$Name, levels = unique(test$Name))

# Assign specific colors to each Type and K category
type_colors <- c(
  "Ren_basal"           = "orange",
  "Ren_hemp-type"       = "green",
  "Ren_drug-type feral" = "blue",
  "Ren_drug-type"       = "red",
  "Soorni_Population_1" = "#9370DB",
  "Soorni_Population_2" = "#008080",
  "Aina_cluster_1"      = "darkred",
  "Aina_cluster_2"      = "darkblue",
  "Aina_cluster_3"      = "darkgreen",
  "Aina_cluster_4"      = "purple4",
  "Aina_cluster_5"      = "darkorange4"
)

# Hiroshige palette for K categories (adjust the palette if needed)
k_colors <- hiroshige_palette[c(1, 3, 5, 6, 7, 8)]  # Choose appropriate colors for K1 to K6

# Plot with the colored bar for Type underneath
final_plot <- ggplot() +
  # Plot the K1 and K2 bars with their own fill colors
  geom_bar(data = long_data, aes(x = Name, y = value, fill = category), stat = "identity") +
  scale_fill_manual(name = "K Categories", values = k_colors) +  # Colors for K1 and K2
  
  # Add the colored bars underneath for Type using a separate geom_tile
  new_scale_fill() +  # Necessary to separate the two scales
  geom_tile(data = test, aes(x = Name, y = -0.05, fill = Type), height = 0.1) +
  scale_fill_manual(name = "Type", values = type_colors) +  # Colors for Type
  
  labs(title = "", x = "", y = "Percent Identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 1),  # Make the text smaller
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 6),       # Smaller legend text
    legend.title = element_text(size = 8),      # Smaller legend title
    legend.key.size = unit(0.4, 'cm'),          # Smaller legend keys
    legend.spacing.y = unit(0.2, 'cm')          # Reduce vertical spacing between legend items
  )

print(final_plot)

# Save the plot as a PDF with reduced height
ggsave(filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/fastSTRUCTURE_new/fastSTRUCTURE_plot_k6.pdf", 
       plot = final_plot, 
       width = 14,    # Width of the plot (adjust as needed)
       height = 3)   # Reduced height (adjust as needed)





########################################################################################################################################################################
# Figure S2
########################################################################################################################################################################
library(tess3r)
#packageVersion("tess3r")
#citation("tess3r")
library(fields)
library(rworldmap)

## 1) Read genotypes (rrBLUP/HMP-like: markers in rows, samples in columns)
gd1 <- read.delim("/Users/annamccormick/R/Feral_Cannabis_EGS/newAina_Ren_Soorni_EGS/merged_on_shared_rrBLUP_format_hmp.txt",
                  sep = "\t", header = TRUE, check.names = FALSE,
                  stringsAsFactors = FALSE)

meta_idx <- 1:4                      # rs., allele, chrom, pos
snp_ids  <- gd1[[1]]
geno     <- as.matrix(gd1[ , -meta_idx, drop = FALSE])  # samples are columns
storage.mode(geno) <- "numeric"
X <- t(geno)                          # individuals x SNPs
if (all(na.omit(unique(as.vector(X))) %in% c(-1,0,1))) X <- X + 1
rownames(X) <- colnames(gd1)[-meta_idx]
colnames(X) <- snp_ids

## 2) Read coordinates and match to X (keep SRR exactly; normalize others)
envdat <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/newAina_Ren_Soorni_EGS/n909_climate_data_all.csv", stringsAsFactors = FALSE)

canon <- function(x){
  x <- trimws(x); x <- tolower(x)
  x <- gsub("-", ".", x); x <- gsub("[ _]+", ".", x); gsub("\\.+",".", x)
}
geno_id_real <- rownames(X)
env_id_real  <- envdat$sample.id

geno_key <- ifelse(grepl("^SRR", geno_id_real),
                   canon(geno_id_real),
                   canon(sub("\\.\\d+$","", geno_id_real)))
env_key  <- canon(env_id_real)

coord_tbl <- data.frame(key = env_key,
                        lon = envdat$Longitude, lat = envdat$Latitude)
idx      <- match(geno_key, coord_tbl$key)
coord    <- coord_tbl[idx, c("lon","lat")]
have_xy  <- complete.cases(coord)

X_sub     <- X[have_xy, , drop = FALSE]
coord_sub <- as.matrix(coord[have_xy, ])
rownames(coord_sub) <- rownames(X_sub)

## 3) Light QC on SNPs (≤20% missing; polymorphic)
miss_rate <- colMeans(is.na(X_sub))
keep_miss <- miss_rate <= 0.20
is_mono   <- function(v){ u <- unique(v[!is.na(v)]); length(u) <= 1 }
keep_poly <- !apply(X_sub, 2, is_mono)
X_filt    <- X_sub[, keep_miss & keep_poly, drop = FALSE]

## 4) Run tess3 over K = 1..14 (3 reps each)
K_range <- 1:14
fit <- tess3(X = X_filt, coord = coord_sub,
             K = K_range, rep = 3, lambda = 1,
             ploidy = 2, max.iteration = 5000, openMP.core = 4)

## 5) Cross-entropy & pick K (lower is better)
ce_mat <- sapply(fit, function(x) x$crossentropy)   # rows = reps, cols = K
best_ce        <- apply(ce_mat, 2, min, na.rm = TRUE)
best_rep_index <- apply(ce_mat, 2, which.min)

plot(K_range, best_ce, type = "b", xlab = "K", ylab = "Cross-entropy")

# choose K (e.g., elbow); change this if you prefer K=11
Kstar   <- 6
repstar <- best_rep_index[which(K_range == Kstar)]

## 6) Get Q and plot like the vignette (kriged, NA+Eurasia window)
Q   <- tess3r::qmatrix(fit, K = Kstar, rep = repstar)
pal <- grDevices::colorRampPalette(c(
  "#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00",
  "#ffff33","#a65628","#f781bf","#999999"
))(ncol(Q))

wm      <- getMap("low")
win_all <- c(-130, 150, 10, 60)  # xmin, xmax, ymin, ymax

pdf(
  file = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/FigS5_k6.pdf",
  width = 12,
  height = 5
)

plot(Q,
     coord       = coord_sub,
     method      = "map.max",
     interpol    = tess3r::FieldsKrigModel(10),
     window      = win_all,
     map.polygon = wm,
     background  = TRUE,
     resolution  = c(500, 300),
     cex         = 0.4,
     col.palette = pal,
     main        = sprintf("Ancestry coefficients (K = %d)", Kstar),
     xlab        = "Longitude", ylab = "Latitude")

dev.off()
#################
# choose K (e.g., elbow); change this if you prefer K=11
Kstar   <- 4
repstar <- best_rep_index[which(K_range == Kstar)]

## 6) Get Q and plot like the vignette (kriged, NA+Eurasia window)
Q   <- tess3r::qmatrix(fit, K = Kstar, rep = repstar)
pal <- grDevices::colorRampPalette(c(
  "#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00",
  "#ffff33","#a65628","#f781bf","#999999"
))(ncol(Q))

wm      <- getMap("low")
win_all <- c(-130, 150, 10, 60)  # xmin, xmax, ymin, ymax


pdf(
  file = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/FigS5_k4.pdf",
  width = 12,
  height = 5
)
plot(Q,
     coord       = coord_sub,
     method      = "map.max",
     interpol    = tess3r::FieldsKrigModel(10),
     window      = win_all,
     map.polygon = wm,
     background  = TRUE,
     resolution  = c(500, 300),
     cex         = 0.4,
     col.palette = pal,
     main        = sprintf("Ancestry coefficients (K = %d)", Kstar),
     xlab        = "Longitude", ylab = "Latitude")

dev.off()
#################
# choose K (e.g., elbow); change this if you prefer K=11
Kstar   <- 5
repstar <- best_rep_index[which(K_range == Kstar)]

## 6) Get Q and plot like the vignette (kriged, NA+Eurasia window)
Q   <- tess3r::qmatrix(fit, K = Kstar, rep = repstar)
pal <- grDevices::colorRampPalette(c(
  "#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00",
  "#ffff33","#a65628","#f781bf","#999999"
))(ncol(Q))

wm      <- getMap("low")
win_all <- c(-130, 150, 10, 60)  # xmin, xmax, ymin, ymax

pdf(
  file = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/FigS5_k5.pdf",
  width = 12,
  height = 5
)

plot(Q,
     coord       = coord_sub,
     method      = "map.max",
     interpol    = tess3r::FieldsKrigModel(10),
     window      = win_all,
     map.polygon = wm,
     background  = TRUE,
     resolution  = c(500, 300),
     cex         = 0.4,
     col.palette = pal,
     main        = sprintf("Ancestry coefficients (K = %d)", Kstar),
     xlab        = "Longitude", ylab = "Latitude")

dev.off()


#################
# choose K (e.g., elbow); change this if you prefer K=11
Kstar   <- 2
repstar <- best_rep_index[which(K_range == Kstar)]

## 6) Get Q and plot like the vignette (kriged, NA+Eurasia window)
Q   <- tess3r::qmatrix(fit, K = Kstar, rep = repstar)
pal <- grDevices::colorRampPalette(c(
  "#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00",
  "#ffff33","#a65628","#f781bf","#999999"
))(ncol(Q))

wm      <- getMap("low")
win_all <- c(-130, 150, 10, 60)  # xmin, xmax, ymin, ymax


pdf(
  file = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/FigS5_k2.pdf",
  width = 12,
  height = 5
)
plot(Q,
     coord       = coord_sub,
     method      = "map.max",
     interpol    = tess3r::FieldsKrigModel(10),
     window      = win_all,
     map.polygon = wm,
     background  = TRUE,
     resolution  = c(500, 300),
     cex         = 0.4,
     col.palette = pal,
     main        = sprintf("Ancestry coefficients (K = %d)", Kstar),
     xlab        = "Longitude", ylab = "Latitude")

dev.off()

#################
# choose K (e.g., elbow); change this if you prefer K=11
Kstar   <- 3
repstar <- best_rep_index[which(K_range == Kstar)]

## 6) Get Q and plot like the vignette (kriged, NA+Eurasia window)
Q   <- tess3r::qmatrix(fit, K = Kstar, rep = repstar)
pal <- grDevices::colorRampPalette(c(
  "#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00",
  "#ffff33","#a65628","#f781bf","#999999"
))(ncol(Q))

wm      <- getMap("low")
win_all <- c(-130, 150, 10, 60)  # xmin, xmax, ymin, ymax


pdf(
  file = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/FigS5_k3.pdf",
  width = 12,
  height = 5
)
plot(Q,
     coord       = coord_sub,
     method      = "map.max",
     interpol    = tess3r::FieldsKrigModel(10),
     window      = win_all,
     map.polygon = wm,
     background  = TRUE,
     resolution  = c(500, 300),
     cex         = 0.4,
     col.palette = pal,
     main        = sprintf("Ancestry coefficients (K = %d)", Kstar),
     xlab        = "Longitude", ylab = "Latitude")

dev.off()




########################################################################################################################################################################
# Figure S3
########################################################################################################################################################################

#############
#PI
#############
#use vcftools on merged vcf files using a population .txt file
#past eg fpr pi for one sample
#bcftools view -S BSh_population_file.txt -Oz -o BSh_population.vcf.gz merged_Ren_Soorni_snps.vcf.gz

#getting separate pi for each new vcf file
#vcftools --gzvcf BSh_population.vcf.gz --site-pi --out BSh_diversity


#checking sample IDs
library(vcfR)
library(dplyr)

# read VCF
vcf <- read.vcfR("/Users/annamccormick/R/Feral_Cannabis_EGS/pi/merged_on_shared.vcf.gz")

# Option 1: quick sample names
vcf@gt %>% colnames()   # includes "FORMAT" + all sample IDs

# Option 2: cleaner — sample names only (excluding FORMAT)
samples <- colnames(vcf@gt)[-1]
samples


#######
#bcftools view -S Aina_cluster1.txt -Oz -o Aina_cluster1.vcf.gz merged_on_shared.vcf.gz
#bcftools view -S Aina_cluster2.txt -Oz -o Aina_cluster2.vcf.gz merged_on_shared.vcf.gz
#bcftools view -S Aina_cluster3.txt -Oz -o Aina_cluster3.vcf.gz merged_on_shared.vcf.gz
#bcftools view -S Aina_cluster4.txt -Oz -o Aina_cluster4.vcf.gz merged_on_shared.vcf.gz
#bcftools view -S Aina_cluster5.txt -Oz -o Aina_cluster5.vcf.gz merged_on_shared.vcf.gz
#bcftools view -S Ren_basal.txt -Oz -o Ren_basal.vcf.gz merged_on_shared.vcf.gz
#bcftools view -S Ren_hemp.txt -Oz -o Ren_hemp.vcf.gz merged_on_shared.vcf.gz
#bcftools view -S Ren_dt_feral.txt -Oz -o Ren_dt_feral.vcf.gz merged_on_shared.vcf.gz
#bcftools view -S Ren_drug_type.txt -Oz -o Ren_drug_type.vcf.gz merged_on_shared.vcf.gz
#bcftools view -S Soorni_pop1.txt -Oz -o Soorni_pop1.vcf.gz merged_on_shared.vcf.gz
#bcftools view -S Soorni_pop2.txt -Oz -o Soorni_pop2.vcf.gz merged_on_shared.vcf.gz

#Step 2 with vcf tools for pi
#vcftools --gzvcf Aina_cluster1.vcf.gz --site-pi --out Aina_cluster1_diversity
#vcftools --gzvcf Aina_cluster2.vcf.gz --site-pi --out Aina_cluster2_diversity
#vcftools --gzvcf Aina_cluster3.vcf.gz --site-pi --out Aina_cluster3_diversity
#vcftools --gzvcf Aina_cluster4.vcf.gz --site-pi --out Aina_cluster4_diversity
#vcftools --gzvcf Aina_cluster5.vcf.gz --site-pi --out Aina_cluster5_diversity
#vcftools --gzvcf Ren_basal.vcf.gz --site-pi --out Ren_basal_diversity
#vcftools --gzvcf Ren_hemp.vcf.gz --site-pi --out Ren_hemp_diversity
#vcftools --gzvcf Ren_dt_feral.vcf.gz --site-pi --out Ren_dt_feral_diversity
#vcftools --gzvcf Ren_drug_type.vcf.gz --site-pi --out Ren_drug_type_diversity
#vcftools --gzvcf Soorni_pop1.vcf.gz --site-pi --out Soorni_pop1_diversity
#vcftools --gzvcf Soorni_pop2.vcf.gz --site-pi --out Soorni_pop2_diversity

library(ggplot2)
library(dplyr)

# Load the populations
# Load the new π results
aina1_pi <- read.table("~/R/Feral_Cannabis_EGS/pi/Aina_cluster1_diversity.sites.pi", header = TRUE) %>% 
  mutate(Population = "Aina Cluster 1, n=20")

aina2_pi <- read.table("~/R/Feral_Cannabis_EGS/pi/Aina_cluster2_diversity.sites.pi", header = TRUE) %>% 
  mutate(Population = "Aina Cluster 2, n=342")

aina3_pi <- read.table("~/R/Feral_Cannabis_EGS/pi/Aina_cluster3_diversity.sites.pi", header = TRUE) %>% 
  mutate(Population = "Aina Cluster 3, n=326")

aina4_pi <- read.table("~/R/Feral_Cannabis_EGS/pi/Aina_cluster4_diversity.sites.pi", header = TRUE) %>% 
  mutate(Population = "Aina Cluster 4, n=40")

aina5_pi <- read.table("~/R/Feral_Cannabis_EGS/pi/Aina_cluster5_diversity.sites.pi", header = TRUE) %>% 
  mutate(Population = "Aina Cluster 5, n=14")

ren_basal_pi <- read.table("~/R/Feral_Cannabis_EGS/pi/Ren_basal_diversity.sites.pi", header = TRUE) %>% 
  mutate(Population = "Ren Basal, n=14")

ren_hemp_pi <- read.table("~/R/Feral_Cannabis_EGS/pi/Ren_hemp_diversity.sites.pi", header = TRUE) %>% 
  mutate(Population = "Ren Hemp, n=35")

ren_dt_feral_pi <- read.table("~/R/Feral_Cannabis_EGS/pi/Ren_dt_feral_diversity.sites.pi", header = TRUE) %>% 
  mutate(Population = "Ren Drug-type Feral, n=17")

ren_drug_type_pi <- read.table("~/R/Feral_Cannabis_EGS/pi/Ren_drug_type_diversity.sites.pi", header = TRUE) %>% 
  mutate(Population = "Ren Drug-type, n=16")

soorni1_pi <- read.table("~/R/Feral_Cannabis_EGS/pi/Soorni_pop1_diversity.sites.pi", header = TRUE) %>% 
  mutate(Population = "Soorni Pop 1, n=49")

soorni2_pi <- read.table("~/R/Feral_Cannabis_EGS/pi/Soorni_pop2_diversity.sites.pi", header = TRUE) %>% 
  mutate(Population = "Soorni Pop 2, n=18")

# Combine all population data into one dataframe
pi_data <- bind_rows(
  aina1_pi, aina2_pi, aina3_pi, aina4_pi, aina5_pi,
  ren_basal_pi, ren_hemp_pi, ren_dt_feral_pi, ren_drug_type_pi,
  soorni1_pi, soorni2_pi
)


# Define chromosome mapping with both number and NC identifier
chromosome_mapping <- c(
  "NC_044371.1" = "1 - NC_044371.1",
  "NC_044375.1" = "2 - NC_044375.1",
  "NC_044372.1" = "3 - NC_044372.1",
  "NC_044373.1" = "4 - NC_044373.1",
  "NC_044374.1" = "5 - NC_044374.1",
  "NC_044377.1" = "6 - NC_044377.1",
  "NC_044378.1" = "7 - NC_044378.1",
  "NC_044379.1" = "8 - NC_044379.1",
  "NC_044376.1" = "9 - NC_044376.1",
  "NC_044370.1" = "10 - NC_044370.1"
)

# Add Chromosome column with labels and ordered levels
pi_data <- pi_data %>%
  mutate(Chromosome = factor(chromosome_mapping[CHROM], 
                             levels = c("1 - NC_044371.1", "2 - NC_044375.1", "3 - NC_044372.1", 
                                        "4 - NC_044373.1", "5 - NC_044374.1", "6 - NC_044377.1", 
                                        "7 - NC_044378.1", "8 - NC_044379.1", "9 - NC_044376.1", 
                                        "10 - NC_044370.1")))

# Plot by chromosome with combined labels
# Define custom colors for each population to improve visibility
custom_colors <- c(
  "Aina Cluster 1, n=20" = "red3",   # red
  "Aina Cluster 2, n=342" = "darkblue",   # blue
  "Aina Cluster 3, n=326" = "darkgreen",   # green
  "Aina Cluster 4, n=40" = "purple",   # purple
  "Aina Cluster 5, n=14" = "darkorange4",   # orange
  "Ren Basal, n=14" = "orange",        # brown
  "Ren Hemp, n=35" = "green",         # pink
  "Ren Drug-type Feral, n=17" = "blue", # gray
  "Ren Drug-type, n=16" = "red",    # yellow
  "Soorni Pop 1, n=49" = "#66c2a5",     # teal
  "Soorni Pop 2, n=18" = "#8da0cb"      # light blue-violet
)

p <- ggplot(pi_data, aes(x = POS, y = PI, color = Population)) +
  geom_smooth(alpha = 0.5, size = 0.8, se = FALSE) +
  labs(title = "",
       x = "Genomic Position (Mb)",
       y = "Nucleotide Diversity (pi)") +
  facet_wrap(~ Chromosome, ncol = 2, nrow = 5, scales = "free_x") +  # Set two columns and five rows
  theme_minimal() +
  scale_color_manual(values = custom_colors) +
  scale_x_continuous(breaks = seq(0, max(pi_data$POS), by = 1e7),  # Set tick marks every 10 Mb
                     labels = function(x) paste0(x / 1e6, " Mb")) +  # Convert to Mb for readability
  theme(
    axis.text.x = element_text(size = 18, angle = 90, hjust = 1),  # Increase font size for x-axis text
    axis.text.y = element_text(size = 20),  # Increase font size for y-axis text
    axis.title.x = element_text(size = 20, face = "plain"),  # Increase font size for x-axis title
    axis.title.y = element_text(size = 20, face = "plain"),  # Increase font size for y-axis title
    plot.title = element_text(size = 20, face = "plain", hjust = 0.5),  # Increase font size for plot title
    strip.text = element_text(size = 20, face = "plain"),  # Increase font size for facet labels
    legend.text = element_text(size = 20),  # Increase font size for legend text
    legend.title = element_text(size = 20, face = "plain")  # Increase font size for legend title
  ) +
  theme(panel.grid.major.x = element_line(color = "grey80"),
        panel.grid.minor.x = element_blank())  # Remove minor grid lines

p

ggsave("/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/FigS16.pdf", plot = p, width = 20, height = 22)



########################################################################################################################################################################
# Figure S4
########################################################################################################################################################################

######################################################
# Prediction Accuracy (Aina+Ren+Soorni)
######################################################
library(rrBLUP)  # RR-BLUP package for Ridge Regression Best Linear Unbiased Prediction
library(hibayes) # hibayes package for Bayesian models
library(dplyr)

# Load custom functions for k-fold cross-validation
source("/Users/annamccormick/R/Feral_Cannabis_EGS/EGS/xval_kfold_functions.R")

# Read genotype data in rrBLUP format
data <- read.delim("/Users/annamccormick/R/Feral_Cannabis_EGS/newAina_Ren_Soorni_EGS/merged_on_shared_rrBLUP_format_hmp.txt", sep = "\t", header = TRUE)

# Read environmental data from CSV
#envdat <- read.csv('/Users/annamccormick/R/Feral_Cannabis_EGS/EGS_fusion_Ren_Soorni_Aina/n909_climate_data_all.csv', header = TRUE, stringsAsFactors = FALSE)
envdat <- read.csv('/Users/annamccormick/R/Feral_Cannabis_EGS/newAina_Ren_Soorni_EGS/n909_climate_data_all.csv', header = TRUE, stringsAsFactors = FALSE)


# Read training dataset 
trainingset <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/newAina_Ren_Soorni_EGS/new_n191_worldclim_data.csv", header = TRUE, stringsAsFactors = FALSE)

###
head(colnames(data), 20)

#first few sample IDs in your envdat
head(envdat$sample.id, 20)

# overlap
common_ids <- intersect(colnames(data), envdat$sample.id)
length(common_ids)


# Set the row names for the environmental and training datasets
row.names(envdat) <- envdat$sample.id
row.names(trainingset) <- trainingset$sample.id

# Filter genotype data to include only columns matching Sample_ID in envdat
gd2 <- data[, c("rs.", "allele", "chrom", "pos", colnames(data)[colnames(data) %in% envdat$sample.id])]

# Set the SNP IDs as the row names for the genotype data
row.names(gd2) <- gd2$rs.

# Remove the first four columns (rs., allele, chrom, pos) to create the genotype matrix
g.in <- gd2[, -c(1:4)]  # Isolate genotype data
g.in <- as.matrix(g.in)
g.in <- t(g.in)  # Transpose to align individuals as rows
cat("Dimensions of g.in", dim(g.in), "\n")

########################
#Remove snps with greater then 20% missingness
#g.in <- g.in[, colSums(is.na(g.in)) == 0]  # Remove SNPs with missing data
g.in <- g.in[, colMeans(is.na(g.in)) <= 0.20]
# Check dimensions of the cleaned genotype matrix
cat("Dimensions of g.in after removing SNPs with missing data:", dim(g.in), "\n")


## where to insert imputation with mean appropriatly - matching the lables used this time
# missingness BEFORE imputation

# overall missing‐data rate, before train/pred split #All data: 871622 / 8090880 missing (10.77%)
total_missing <- sum(is.na(g.in))
total_cells   <- prod(dim(g.in))
cat(sprintf(
  "All data: %d / %d missing (%.2f%%)\n",
  total_missing, total_cells, 100 * total_missing / total_cells
))

# 2) Impute remaining missing values with the SNP (column) mean
marker_means <- colMeans(g.in, na.rm = TRUE)
for (j in seq_len(ncol(g.in))) {
  nas <- is.na(g.in[, j])
  if (any(nas)) {
    g.in[nas, j] <- marker_means[j]
  }
}
cat("Total missing after imputation:", sum(is.na(g.in)), "\n")



### Load Environmental Data ###
# Select unique identifier and environmental data columns for training set
# Select column 1 plus bio columns
colnames(envdat)[4:22]

y.trainset <- trainingset[, c(1, 5:23)]
# Move sample.id into row names
rownames(y.trainset) <- y.trainset$sample.id


# Select unique identifier and environmental data columns for the full dataset
y.in <- envdat[, c(1, 4:ncol(envdat))]

# Convert environmental data to matrix format for RR-BLUP
y.trainset.mat <- as.matrix(y.trainset[, -1])  # Exclude the first column (Sample_ID)
y.in.mat <- as.matrix(y.in[, -1])  # Exclude the first column (Sample_ID)
#y.in.mat <- y.in.mat[, -ncol(y.in.mat)]
dim(y.in.mat)

# Ensure row alignment between genotype and environmental data
# Ensure the same for the training set

# Check dimensions before cross-validation
cat("Dimensions of g.in:", dim(g.in), "\n")
cat("Dimensions of y.in.mat:", dim(y.in.mat), "\n")
cat("Dimensions of y.trainset.mat:", dim(y.trainset.mat), "\n")

length(intersect(rownames(g.in), rownames(y.in.mat)))  # Should be 909
length(intersect(rownames(g.in), rownames(y.trainset.mat)))  # Should be 191


########################
# RR-BLUP ###
########################
# k.xval function from external script, running 10-fold cross-validation with 50 repetitions
xval_k10_rrblup <- k.xval(g.in = g.in, y.in = y.in.mat, y.trainset = y.trainset.mat, k.fold = 10, reps = 50)

saveRDS(xval_k10_rrblup, "/Users/annamccormick/R/Feral_Cannabis_EGS/newAina_Ren_Soorni_EGS/xval_rrblup_kfold_10.RData")
########################################################################

########################
# Gaussian Kernel
########################
K <- A.mat(g.in) #Compute relationship matrix using Gaussian Kernel
k_dist <- dist(K) #Calculate distance matrix from relationship matrix

#Run k-fold cross-validation for Gaussian Kernel
xval_k10_GAUSS <- k.xval.GAUSS(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k_dist = k_dist, k.fold = 10, reps = 50)
saveRDS(xval_k10_GAUSS, "/Users/annamccormick/R/Feral_Cannabis_EGS/newAina_Ren_Soorni_EGS/xval_GAUSS_kfold_10.RData")

########################
#Exponential Kernel
########################
#Compute relationship matrix using Exponential Kernel
K.Exp=Kernel_computation(X=g.in, name="exponential", degree=NULL, nL=NULL)

#Set row and column names for the Exponential Kernel matrix
row.names(K.Exp) <- rownames(g.in)
colnames(K.Exp) <- rownames(g.in)

#Calculate distance matrix from Exponential Kernel relationship matrix
exp_dist <- dist(K.Exp) # Calculate Relationship Matrix\

#Run k-fold cross-validation for Exponential Kernel
xval_k10_EXP <- k.xval.EXP(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k_dist = exp_dist, k.fold = 10, reps = 50)
saveRDS(xval_k10_EXP, "/Users/annamccormick/R/Feral_Cannabis_EGS/newAina_Ren_Soorni_EGS/xval_EXP_kfold_10.RData")

########################
#BayesCPi
########################
#Run k-fold cross-validation for BayesCPi model
xval_k10_BayesCpi <- k.xval.BayesCpi(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k.fold = 10, reps = 50, niter=3000,nburn=1200)
saveRDS(xval_k10_BayesCpi, "/Users/annamccormick/R/Feral_Cannabis_EGS/newAina_Ren_Soorni_EGS/xval_BayesCpi_kfold_10.RData")


#############################################################################################################################
#PLOTS
#####
library(ggplot2)
setwd("/Users/annamccormick/R/Feral_Cannabis_EGS/newAina_Ren_Soorni_EGS//")

# RR-BLUP - K-Fold Cross-Validation
#Load cross-validation results for rrBLUP
rrblup_kfold10 <- readRDS("xval_rrblup_kfold_10.RData")
rrblup_kfold10 <- rrblup_kfold10$xval.result
rrblup_kfold10$r.mean <- as.numeric(rrblup_kfold10$r.mean)

# Gaussian Kernel - K-fold
#Load cross-validation results for Gaussian Kernel
gauss_kfold_10 <- readRDS('xval_GAUSS_kfold_10.RData')
gauss_kfold_10 <- gauss_kfold_10$xval.result
gauss_kfold_10$r.mean <- as.numeric(gauss_kfold_10$r.mean)

# Exponential Kernel - K-fold
#Load cross-validation results for Exponential Kernel
EXP_kfold_10 <- readRDS('xval_EXP_kfold_10.RData')
EXP_kfold_10 <- EXP_kfold_10$xval.result
EXP_kfold_10$r.mean <- as.numeric(EXP_kfold_10$r.mean)

# BayesCpi - K-fold
#Load cross-validation results for BayesCpi model
bayescpi_kfold_10 <- readRDS('xval_BayesCpi_kfold_10.RData')
bayescpi_kfold_10 <- bayescpi_kfold_10$xval.result
bayescpi_kfold_10$r.mean <- as.numeric(bayescpi_kfold_10$r.mean)


# Organize all dataframes for merging
## Rename model names
rrblup_kfold10$model <- "rrBLUP"
gauss_kfold_10$model <- "Gaussian Kernel"
EXP_kfold_10$model <- "Exponential Kernel"
bayescpi_kfold_10$model <- "BayesCpi"

## Input xval type
rrblup_kfold10$xval <- "Ten-Fold"
gauss_kfold_10$xval <- "Ten-Fold"
EXP_kfold_10$xval <- "Ten-Fold"
bayescpi_kfold_10$xval <- "Ten-Fold"


#Combine all model results into a single list
model_list <- list(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10,bayescpi_kfold_10)
#model_list <- list(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10)

#Remove any NA values from the model results
model_list1 <- lapply(model_list, na.omit)

#Combine all models into a single dataframe
all_models <- do.call("rbind", model_list1)

#Convert standard deviation values to numeric
all_models$r.sd <- as.numeric(all_models$r.sd)

#Ensure cross-validation type is a factor with specified levels
all_models$xval <- factor(all_models$xval, levels = c("Ten-Fold"))

#Filter for specific traits of interest
all_bio <- all_models[all_models$trait %in% c('bio_01', 'bio_02', 'bio_03', 'bio_04', 'bio_05','bio_06',
                                              'bio_07','bio_08','bio_09','bio_010','bio_11',
                                              'bio_12','bio_13','bio_14','bio_15','bio_16','bio_17','bio_18','bio_19'),]

#Plot 
all_bio$trait <- factor(all_bio$trait, levels = paste0("bio", 1:19))

p<- ggplot(all_models, aes(y = r.mean, x = model, color = model)) +
  theme_bw() +
  geom_errorbar(aes(x = model, ymin = r.mean-r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 1) + 
  geom_hline(yintercept = 0.5, color="red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1)

p

ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/FigS1.pdf",
  plot = p,
  width = 15,
  height = 10,
  units = "in"
)





########################################################################################################################################################################
# Figure S5
########################################################################################################################################################################

####################
#EGS results all 19 for Aina+Ren+Soorni
####################
#Selecting representatives from each population for the Aina dataset - at random
#core covers 75 pops
#have 5 pop with no lat long
#have 5 pop with lat long but not in core - so selected one genotype randomly from these
#####################################
library(readxl)
library(dplyr)

#df <- read_excel("/Users/annamccormick/R/Feral_Cannabis_EGS/EGS_fusion_Ren_Soorni_Aina/85_pop_selections.xlsx", sheet = 1)
df <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/newAina_Ren_Soorni_EGS/new_Aina_n760_with_core_flag.csv")
head(df)
str(df)


# Step 1: How many populations are covered by core samples?
core_covered_pops <- df %>%
  filter(is_core_310 == TRUE) %>%
  distinct(pop_base)

num_covered <- nrow(core_covered_pops)
cat("Number of populations covered by the core:", num_covered, "\n")

# Step 2: For those covered populations, choose one representative at random
core_representatives <- df %>%
  filter(is_core_310 == TRUE) %>%
  group_by(pop_base) %>%
  slice_sample(n = 1) %>%
  ungroup()

# Step 3: Identify populations with no core sample
all_pops <- df %>% distinct(pop_base)

uncovered_pops <- all_pops %>%
  filter(!pop_base %in% core_covered_pops$pop_base)

num_uncovered <- nrow(uncovered_pops)
cat("Number of populations NOT covered by the core:", num_uncovered, "\n")

# Step 4: Randomly pick one sample from each uncovered population
uncovered_representatives <- df %>%
  filter(pop_base %in% uncovered_pops$pop_base) %>%
  group_by(pop_base) %>%
  slice_sample(n = 1) %>%
  ungroup()

# Step 5: Combine into one full representative list
all_representatives <- bind_rows(core_representatives, uncovered_representatives)

cat("Total representative individuals:", nrow(all_representatives), "\n")
write.csv(all_representatives, file = "/Users/annamccormick/R/Feral_Cannabis_EGS/newAina_Ren_Soorni_EGS/core_representatives.csv", row.names = FALSE)


##################################################################################################################################################################
library(dplyr)

# Exclude selected representative samples from the original data
non_representatives <- df %>%
  anti_join(all_representatives, by = "sample.id")

write.csv(non_representatives %>% select(sample.id),
          file = "/Users/annamccormick/R/Feral_Cannabis_EGS/EGS_fusion_Ren_Soorni_Aina/non_representative_sample_ids.csv",
          row.names = FALSE)
#5 sample IDs lose - pulled from the all_representatives where the lat/long was NA


######################################################
# EGS
######################################################
# Step 5: Genomic Selection with RR-BLUP
######################################################
library(rrBLUP)
library(dplyr)

# Load genotype and environmental data
gd1 <- read.delim("/Users/annamccormick/R/Feral_Cannabis_EGS/newAina_Ren_Soorni_EGS/merged_on_shared_rrBLUP_format_hmp.txt", sep = "\t", header = TRUE)


#envdat <- read.csv('/Users/annamccormick/R/Feral_Cannabis_EGS/EGS_fusion_Ren_Soorni_Aina/n191_worldclim_data.csv', head = TRUE)
envdat <- read.csv('/Users/annamccormick/R/Feral_Cannabis_EGS/newAina_Ren_Soorni_EGS/new_n191_worldclim_data.csv', head = TRUE)

# Set row names for genotype data (use 'rs.' column which contains SNP IDs)
row.names(gd1) <- gd1$rs.
head(rownames(gd1))
head(gd1$rs.)

# Remove the non-genotype columns ('rs.', 'allele', 'chrom', 'pos') for analysis
gd3 <- gd1[, -c(1:4)]  # SNP columns only
g.in <- as.matrix(gd3)  # Convert to matrix
g.in <- t(g.in)  # Transpose g.in to align with g.train structure (SNPs as columns, individuals as rows)
print(g.in [1:10, 1:10])  # Visualize 

# Set row names for environmental data using 'Sample_ID'
row.names(envdat) <- envdat$sample.id

# Subset environmental data for analysis (start from the 5th column)
y.in.rr <- envdat[, 5:23]  #14-32


minicore_entries <- which(envdat$is_core_310 == TRUE)  # Subset to core lines where 'Core' is TRUE
y.in.rr <- y.in.rr[minicore_entries, ]  
y.in.mat <- as.matrix(y.in.rr)

print(y.in.mat [1:10, 1:10])  # Visualize 
# Now fix the rownames to match the dots in g.in
#rownames(y.in.mat) <- gsub("-", "\\.", rownames(y.in.mat))
#print(y.in.mat [1:10, 1:10])  # Visualize 

##########################
# Training set
train <- row.names(y.in.mat)  # Names of training lines

# fix labelling issue # replace every "-" with "."  
train_dots <- gsub("-", "\\.", train)
train <- train_dots

# 2. now you can subset
#g.train <- g.in[ train_dots, ] # Subset rows (individuals) directly
g.train <- g.in[ train, ] # Subset rows (individuals) directly

print(g.train[1:10, 1:10])  # Visualize part of g.train for debugging
dim(g.train)


# Prediction set 
#pred    <- setdiff(row.names(g.in), train_dots) 
pred    <- setdiff(row.names(g.in), train)
g.pred  <- g.in[pred, ]   
dim(g.pred)

########################################################################################################################
# g.train is 310 × M (M = number of SNPs before filtering)
miss_rate <- colMeans(is.na(g.train))
summary(miss_rate)
hist(miss_rate, breaks=50,
     main="Missingness per SNP", xlab="Proportion missing")

# Compute per‑SNP missing rates
miss_rate <- colMeans(is.na(g.train))

# Identify SNPs to keep less then ≤ 20% missing. If greater then 20% missingness discarded
keep_snps <- names(miss_rate)[miss_rate <= 0.20]
cat("After missingness filter, SNP count =", length(keep_snps), "\n")

# Subset both training and prediction matrices to those SNPs
g.train <- g.train[, keep_snps]
g.pred  <- g.pred[,  keep_snps]

dim(g.train)
dim(g.pred)



############################################################
# 2) Impute the remaining missing genotypes with the marker mean
############################################################
# missingness BEFORE imputation

# overall missing‐data rate, before train/pred split
total_missing <- sum(is.na(g.in))
total_cells   <- prod(dim(g.in))
cat(sprintf(
  "All data: %d / %d missing (%.2f%%)\n",
  total_missing, total_cells, 100 * total_missing / total_cells
))


# Function to mean‑impute each column
mean_impute <- function(mat) {
  for (j in seq_len(ncol(mat))) {
    na_idx <- is.na(mat[,j])
    if (any(na_idx)) {
      mu <- mean(mat[,j], na.rm = TRUE)
      # if your algorithm wants integer codes (–1,0,1), you could round(mu)
      mat[na_idx, j] <- mu
    }
  }
  mat
}

g.train <- mean_impute(g.train)
g.pred  <- mean_impute(g.pred)

# sanity check
stopifnot(!any(is.na(g.train)))
stopifnot(!any(is.na(g.pred)))
cat("No missing values remain in g.train or g.pred\n")

############################################################
# RR‑BLUP loop
# List of traits to analyze
traits <- colnames(y.in.mat)

# Initialize object for storing results
gebv_df <- data.frame(matrix(nrow = nrow(g.in), ncol = length(traits)))
colnames(gebv_df) <- traits
row.names(gebv_df) <- row.names(g.in)

# RR-BLUP loop for each trait
for (t in 1:length(traits)) {
  trait <- traits[t]
  
  # Set up training set for the trait
  y.train <- as.matrix(y.in.mat[train, trait])  
  
  # Run RR-BLUP model
  solve.out <- mixed.solve(y = y.train, Z = g.train, SE = FALSE, return.Hinv = FALSE)
  u.hat <- solve.out$u
  
  # Calculate GEBVs for both prediction and training sets
  GEBV <- g.pred %*% u.hat
  GEBV_train <- g.train %*% u.hat
  
  # Store results in the combined dataframe
  pred_rows <- match(row.names(g.pred), row.names(g.in))  # Indices for pred rows
  train_rows <- match(train, row.names(g.in))  # Indices for train rows
  
  gebv_df[pred_rows, t] <- GEBV  # Predictions for test lines
  gebv_df[train_rows, t] <- GEBV_train  # Predictions for training lines
}

# Final output: GEBVs for all individuals
write.csv(gebv_df, '/Users/annamccormick/R/Feral_Cannabis_EGS/newAina_Ren_Soorni_EGS/rrblup_GEAV_values_Ren_Soorni_Aina_new.csv')


#######################################################################################################################################
#adding back in metadata
#######################################################################################################################################
gebv_df <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/newAina_Ren_Soorni_EGS/rrblup_GEAV_values_Ren_Soorni_Aina_new.csv")

str(gebv_df)

library(ggplot2)
library(dplyr)
library(tidyr)

# Convert bio_names to a named vector if not already
bio_names <- c(
  bio_01 = "Annual Mean Temperature (bio 01)",
  bio_02 = "Mean Diurnal Range (bio 02)",
  bio_03 = "Isothermality (bio 03)",
  bio_04 = "Temperature Seasonality (bio 04)",
  bio_05 = "Max Temperature of Warmest Month (bio 05)",
  bio_06 = "Min Temperature of Coldest Month (bio 06)",
  bio_07 = "Temperature Annual Range (bio 07)",
  bio_08 = "Mean Temperature of Wettest Quarter (bio 08)",
  bio_09 = "Mean Temperature of Driest Quarter (bio 09)",
  bio_10 = "Mean Temperature of Warmest Quarter (bio 10)",
  bio_11 = "Mean Temperature of Coldest Quarter (bio 11)",
  bio_12 = "Annual Precipitation (bio 12)",
  bio_13 = "Precipitation of Wettest Month (bio 13)",
  bio_14 = "Precipitation of Driest Month (bio 14)",
  bio_15 = "Precipitation Seasonality (bio 15)",
  bio_16 = "Precipitation of Wettest Quarter (bio 16)",
  bio_17 = "Precipitation of Driest Quarter (bio 17)",
  bio_18 = "Precipitation of Warmest Quarter (bio 18)",
  bio_19 = "Precipitation of Coldest Quarter (bio 19)"
)

my_cols2 <- c(
  "Aina_cluster_1"               = "darkred",
  "Aina_cluster_2"               = "darkblue",
  "Aina_cluster_3"               = "darkgreen",
  "Aina_cluster_4"               = "purple4",
  "Aina_cluster_5"               = "darkorange4",
  "Ren_basal"                    = "orange",
  "Ren_drug-type"                = "red",
  "Ren_drug-type feral"          = "blue",
  "Ren_hemp-type"                = "green",
  "Soorni_Population_1"          = "purple",
  "Soorni_Population_2"          = "#8DA0CB"
)

# Ensure bio_variable is ordered numerically
library(dplyr)
library(tidyr)
library(stringr)


gebv_long <- gebv_df %>%
  pivot_longer(
    cols = starts_with("bio_"),       # bio_01 ... bio_19
    names_to = "bio_variable",
    values_to = "GEBV"
  )

# Define proper order of bio variables
bio_levels <- paste0("bio_", str_pad(1:19, 2, pad = "0"))

# Add factor ordering + descriptive labels
gebv_long <- gebv_long %>%
  mutate(
    bio_variable = factor(bio_variable, levels = bio_levels),
    bio_label    = factor(bio_names[bio_variable], levels = bio_names[bio_levels])
  )


ggplot(gebv_long, aes(x = Population, y = GEBV, fill = Population)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.4) +
  scale_fill_manual(values = my_cols2) +
  facet_wrap(~ bio_label, scales = "free_y", ncol = 3) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "none"
  ) +
  labs(
    title = "",
    x = "Population",
    y = "GEAV"
  )


########################################################################################################################################################################
# Figure S6
########################################################################################################################################################################

####################
# Figure S6A
####################
#Koppen-Geiger by dataset
library(raster)
library(dplyr)
library(ggplot2)

# Load your sample data
data <- read.csv('/Users/annamccormick/R/Feral_Cannabis_EGS/EGS_fusion_Ren_Soorni_Aina/n191_lat_long_koppen_geiger.csv')

# Load the Köppen-Geiger raster
climate_raster <- raster("/Users/annamccormick/R/Feral_Cannabis_EGS/Beck_KG_V1/Beck_KG_V1_present_0p083.tif")

# Extract climate class for each lat-long point
coordinates <- data.frame(lon = data$Longitude, lat = data$Latitude)
data$Climate_Class <- raster::extract(climate_raster, coordinates)


# Define climate class labels
climate_labels <- c(
  `1` = "Af - Tropical, rainforest",
  `2` = "Am - Tropical, monsoon",
  `3` = "Aw - Tropical, Savanna",
  `4` = "BWh - Arid, desert, hot",
  `5` = "BWk - Arid, desert, cold",
  `6` = "BSh - Arid, steppe, hot",
  `7` = "BSk - Arid, steppe, cold",
  `8` = "Csa - Temperate, dry summer, hot summer",
  `9` = "Csb - Temperate, dry summer, warm summer",
  `10` = "Csc -Temperate, dry summer, cold summer",
  `11` = "Cwa - Temperate, dry winter, hot summer",
  `12` = "Cwb - Temperate, dry winter, warm summer",
  `13` = "Cwc - Temperate, dry winter, cold summer",
  `14` = "Cfa - Temperate, no dry season, hot summer",
  `15` = "Cfb - Temperate, no dry season, warm summer",
  `16` = "Cfc - Temperate, no dry season, cold summer",
  `17` = "Dsa - Cold, dry summer, hot summer",
  `18` = "Dsb - Cold, dry summer, warm summer",
  `19` = "Dsc - Cold, dry summer, cold summer",
  `20` = "Dsd - Cold, dry summer, very cold winter ",
  `21` = "Dwa - Cold, dry winter, hot summer",
  `22` = "Dwb - Cold, dry winter, warm summer",
  `23` = "Dwc - Cold, dry winter, cold summer",
  `24` = "Dwd - Cold, dry winter, very cold winter",
  `25` = "Dfa - Cold, no dry season, hot summer",
  `26` = "Dfb - Cold, no dry season, warm summer",
  `27` = "Dfc - Cold, no dry season, cold summer",
  `28` = "Dfd - Cold, no dry season, very cold winter",
  `29` = "ET - Polar, tundra",
  `30` = "EF - Polar, frost"
)

# Assign readable labels
data$Climate_Label <- factor(climate_labels[as.character(data$Climate_Class)], levels = climate_labels)


# Count and normalize by climate class
climate_counts <- data %>%
  count(Climate_Label, Dataset) %>%
  group_by(Climate_Label) %>%
  mutate(Proportion = n / sum(n)) %>%
  ungroup()


p <- ggplot(climate_counts, aes(x = Climate_Label, y = Proportion, fill = Dataset)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c(
    "Aina"   = "darkred",
    "Ren"    = "darkgreen",
    "Soorni" = "darkblue"
  )) +
  labs(
    x = "Köppen-Geiger Climate Class",
    y = "Proportion of Samples",
    fill = "Dataset"
  ) +
  coord_flip() +  # Rotate the plot 90 degrees
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9),     # Adjust y-axis (climate labels) text size
    axis.title = element_text(size = 12),
    legend.position = "right"
  )

p


ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/Fig2A.pdf",
  plot = p,
  width = 8,
  height = 6,
  units = "in"
)

####################
# Figure S6B
####################
# simple present to future Köppen–Geiger
################

library(terra)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# Paths
kg_dir  <- "/Users/annamccormick/R/Feral_Cannabis_EGS/koppen_geiger_tif"
pts_csv <- "/Users/annamccormick/R/Feral_Cannabis_EGS/EGS_fusion_Ren_Soorni_Aina/n191_lat_long_koppen_geiger.csv"

names(read_csv(pts_csv, show_col_types = FALSE))
# 1) Read points
pts_df <- read_csv(pts_csv, show_col_types = FALSE) %>%
  dplyr::rename(
    sample = Sample_ID,
    lon    = Longitude,
    lat    = Latitude
  )

pts_v <- vect(pts_df, geom = c("lon","lat"), crs = "EPSG:4326")

# 2) Pick just two rasters: present (1991–2020) and 2041–2070 (SSP585)
r_present <- rast(file.path(kg_dir, "1991_2020/koppen_geiger_0p5.tif"))          # adjust if needed
r_future  <- rast(file.path(kg_dir, "2041_2070/ssp585/koppen_geiger_0p5.tif"))   # adjust if needed

# 3) Extract climate class at points
vals_present <- terra::extract(r_present, pts_v, method = "near")[, 2]
vals_future  <- terra::extract(r_future,  pts_v, method = "near")[, 2]

# 4) Combine into one table
out_tbl <- pts_df %>%
  mutate(
    kg_present = vals_present,
    kg_future  = vals_future,
    changed    = kg_present != kg_future
  )

# 5) Quick summaries
table(out_tbl$changed)
prop.table(table(out_tbl$changed)) * 100

out_changes <- out_tbl %>%
  filter(changed) %>%
  count(kg_present, kg_future, sort = TRUE)

# Save if needed
#dir.create(file.path(kg_dir, "KG_outputs"), showWarnings = FALSE, recursive = TRUE)
#write_csv(out_tbl, file.path(kg_dir, "KG_outputs/kg_present_vs_future.csv"))

################
# Plot setup
################

# Climate class labels (your coding 1–30)
climate_labels <- c(
  `1` = "Af - Tropical, rainforest",
  `2` = "Am - Tropical, monsoon",
  `3` = "Aw - Tropical, savanna",
  `4` = "BWh - Arid, desert, hot",
  `5` = "BWk - Arid, desert, cold",
  `6` = "BSh - Arid, steppe, hot",
  `7` = "BSk - Arid, steppe, cold",
  `8` = "Csa - Temperate, dry summer, hot summer",
  `9` = "Csb - Temperate, dry summer, warm summer",
  `10` = "Csc - Temperate, dry summer, cold summer",
  `11` = "Cwa - Temperate, dry winter, hot summer",
  `12` = "Cwb - Temperate, dry winter, warm summer",
  `13` = "Cwc - Temperate, dry winter, cold summer",
  `14` = "Cfa - Temperate, no dry season, hot summer",
  `15` = "Cfb - Temperate, no dry season, warm summer",
  `16` = "Cfc - Temperate, no dry season, cold summer",
  `17` = "Dsa - Cold, dry summer, hot summer",
  `18` = "Dsb - Cold, dry summer, warm summer",
  `19` = "Dsc - Cold, dry summer, cold summer",
  `20` = "Dsd - Cold, dry summer, very cold winter",
  `21` = "Dwa - Cold, dry winter, hot summer",
  `22` = "Dwb - Cold, dry winter, warm summer",
  `23` = "Dwc - Cold, dry winter, cold summer",
  `24` = "Dwd - Cold, dry winter, very cold winter",
  `25` = "Dfa - Cold, no dry season, hot summer",
  `26` = "Dfb - Cold, no dry season, warm summer",
  `27` = "Dfc - Cold, no dry season, cold summer",
  `28` = "Dfd - Cold, no dry season, very cold winter",
  `29` = "ET - Polar, tundra",
  `30` = "EF - Polar, frost"
)

# Köppen–Geiger map-style colour scheme (named by your class codes)
# (This matches the categorical look of the Beck/Peel-style KG maps: A=blue, B=red/orange,
#  C=greens, D=cyan/blue-violet, E=greys)
kg_colors <- c(
  `1`  = "#0000FF",  # Af
  `2`  = "#0078FF",  # Am
  `3`  = "#46A6FF",  # Aw
  
  `4`  = "#FF0000",  # BWh
  `5`  = "#FF9696",  # BWk
  `6`  = "#F5A300",  # BSh
  `7`  = "#FFDC64",  # BSk
  
  `8`  = "#FFFF00",  # Csa
  `9`  = "#C8FF00",  # Csb
  `10` = "#96FF00",  # Csc
  `11` = "#00FF00",  # Cwa
  `12` = "#00C800",  # Cwb
  `13` = "#009600",  # Cwc
  `14` = "#96FF96",  # Cfa
  `15` = "#64C864",  # Cfb
  `16` = "#329632",  # Cfc
  
  `17` = "#00FFFF",  # Dsa
  `18` = "#37C8FF",  # Dsb
  `19` = "#007DFF",  # Dsc
  `20` = "#0050C8",  # Dsd
  `21` = "#B4E6FF",  # Dwa
  `22` = "#64B4FF",  # Dwb
  `23` = "#3296FF",  # Dwc
  `24` = "#0064C8",  # Dwd
  `25` = "#C8C8FF",  # Dfa
  `26` = "#9696FF",  # Dfb
  `27` = "#6464FF",  # Dfc
  `28` = "#3232C8",  # Dfd
  
  `29` = "#BEBEBE",  # ET
  `30` = "#808080"   # EF
)

# 1) Long format
df <- out_tbl %>%
  dplyr::select(
    sample,
    `1991–2020`          = kg_present,
    `2041–2070 (SSP585)` = kg_future
  ) %>%
  tidyr::pivot_longer(
    -sample,
    names_to = "group",
    values_to = "code"
  ) %>%
  dplyr::filter(!is.na(code))


# 2) Factor order by numeric code
codes_present <- sort(unique(df$code))
df <- df %>% mutate(code = factor(code, levels = codes_present))

# 3) Legend labels: "code – description"
legend_labs <- setNames(
  sapply(as.character(codes_present), function(z) {
    nm <- climate_labels[[z]]
    if (is.null(nm) || is.na(nm)) paste0(z, " – (unmapped)") else paste0(z, " – ", nm)
  }),
  as.character(codes_present)
)

# 4) Count per group × class
comp <- df %>% count(group, code, name = "n")

# 5) Subset KG colours to only those classes that appear (prevents warnings)
kg_colors_use <- kg_colors[levels(comp$code)]

# 6) Plot
p<- ggplot(comp, aes(x = code, y = n, fill = code)) +
  geom_col() +
  geom_text(aes(label = paste0("n=", n)),
            vjust = -0.25, size = 3.3, show.legend = FALSE) +
  facet_wrap(~ group, ncol = 1, scales = "free_y") +
  scale_fill_manual(
    values = kg_colors_use,
    breaks = levels(comp$code),
    labels = legend_labs,
    name   = "Köppen–Geiger climate class"
  ) +
  labs(
    x = "Climate class code",
    y = "Number of samples",
    title = "Köppen–Geiger classes: 1991–2020 vs 2041–2070 (SSP585)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

p
ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/FigS8.pdf",
  plot = p,
  width = 15,
  height = 10,
  units = "in"
)

########################################################################################################################################################################
# Figure S7
########################################################################################################################################################################
#SNPRelate
library(SNPRelate)
setwd("/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_VCF_remake/")

vcf.fn <- "/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_VCF_remake/Cannabis_sativa_Aina_n760_filtered.vcf.gz"

#Imports variants & convert to Genlight object
snpgdsVCF2GDS(vcf.fn, "SNPs_mergedAll_2.gds", method="copy.num.of.ref")
snpgdsSummary("SNPs_mergedAll_2.gds")
genofile <- snpgdsOpen("SNPs_mergedAll_2.gds")
set.seed(1234)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,autosome.only = F) 
names(snpset)
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2, autosome.only = F)

# === NEW SECTION: Keep only SNPs from NC_ chromosomes ===
# Extract SNP metadata
snp.id <- read.gdsn(index.gdsn(genofile, "snp.id"))
snp.chrom <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
snp.annot <- data.frame(id = snp.id, chrom = snp.chrom)
rownames(snp.annot) <- as.character(snp.annot$id)

# Flatten snpset and match to annotation
snpset.id <- unlist(snpset)
snp.annot.sub <- snp.annot[as.character(snpset.id), ]

# Keep only SNPs on chromosomes starting with "NC_"
nc_snps <- rownames(snp.annot.sub)[grepl("^NC_", snp.annot.sub$chrom)]
snpset.id <- snpset.id[snpset.id %in% nc_snps]
# === END OF NEW SECTION ===

# PCA
pca <- snpgdsPCA(genofile, snp.id = snpset.id, num.thread = 2, autosome.only = FALSE)

# PCA results to dataframe
pc.percent <- pca$varprop * 100
head(round(pc.percent, 2))
# 1.95 1.38 1.22 1.17 0.73 0.68
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],
                  EV2 = pca$eigenvect[,2],
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)

# Save PCA table
write.csv(tab, "/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_VCF_remake/PCA_filtered.csv")

#########################
# Figure S7A - plot
#########################
cluster_cols <- c(
  "1" = "#E41A1C",  # red
  "2" = "#377EB8",  # blue
  "3" = "#4DAF4A",  # green
  "4" = "#984EA3",  # purple
  "5" = "#FF7F00"   # orange
)

pca_df <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_VCF_remake/PCA_filtered.csv")

# ensure clusters match the palette keys (as character/factor)
pca_df$HCPC_cluster <- factor(pca_df$HCPC_cluster, levels = names(cluster_cols))
scale_color_manual(
  values = cluster_cols,
  labels = c(
    "1" = "Aina Cluster 1",
    "2" = "Aina Cluster 2",
    "3" = "Aina Cluster 3",
    "4" = "Aina Cluster 4",
    "5" = "Aina Cluster 5"
  ),
  na.value = "grey70",
  drop = FALSE
)

p<- ggplot(pca_df, aes(x = EV1, y = EV2, color = HCPC_cluster)) +
  stat_ellipse(aes(color = HCPC_cluster),
               type = "norm", level = 0.95, linewidth = 0.6, na.rm = TRUE) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_manual(values = cluster_cols, na.value = "grey70", drop = FALSE) +
  theme_minimal() +
  labs(
    title = "3,762 SNP",
    x = "Eigenvector 1 (1.95)",
    y = "Eigenvector 2 (1.38)",
    color = "HCPC cluster"
  )

p

ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/FigS9A.pdf",
  plot = p,
  width = 15,
  height = 10,
  units = "in"
)

##fix labels
#########################
# Figure S7A – PCA plot
#########################

library(ggplot2)
library(dplyr)

#-----------------------------
# 1) Define colour palette
#-----------------------------
cluster_cols <- c(
  "1" = "darkred",
  "2" = "darkblue",
  "3" = "darkgreen",
  "4" = "purple4",
  "5" = "darkorange4"
)

#-----------------------------
# 2) Load PCA data
#-----------------------------
pca_df <- read.csv(
  "/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_VCF_remake/PCA_filtered.csv"
)

# Ensure cluster is a factor with fixed order
pca_df$HCPC_cluster <- factor(
  pca_df$HCPC_cluster,
  levels = names(cluster_cols)
)

#-----------------------------
# 3) PCA plot
#-----------------------------
p <- ggplot(pca_df, aes(x = EV1, y = EV2, color = HCPC_cluster)) +
  
  stat_ellipse(
    type = "norm",
    level = 0.95,
    linewidth = 0.6,
    na.rm = TRUE
  ) +
  
  geom_point(size = 2, alpha = 0.9) +
  
  scale_color_manual(
    values = cluster_cols,
    labels = c(
      "1" = "Aina Cluster 1",
      "2" = "Aina Cluster 2",
      "3" = "Aina Cluster 3",
      "4" = "Aina Cluster 4",
      "5" = "Aina Cluster 5"
    ),
    na.value = "grey70",
    drop = FALSE
  ) +
  
  theme_minimal() +
  
  labs(
    title = "3,762 SNP",
    x = "Eigenvector 1 (1.95)",
    y = "Eigenvector 2 (1.38)",
    color = "Aina cluster"
  )

print(p)

ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/FigS9A.pdf",
  plot = p,
  width = 8,
  height = 6,
  units = "in"
)


#######################
# Figure S7B
#######################
#ONLY AINA
library(dplyr)
library(ggplot2)
library(maps)

# 0) your data
df_clean <- read.csv(
  "/Users/annamccormick/R/Feral_Cannabis_EGS/VCF/filtering/10chr_PCA_with_coord_and_clusters.csv",
  stringsAsFactors = FALSE
) %>%
  filter(
    !is.na(HCPC_cluster),
    !is.na(Longitude),
    !is.na(Latitude)
  )

# 1) state polygon
states <- map_data("state")

# 2) state label positions
#    maps::state.center gives a data.frame of x/y for each state in state.name
state_labs <- data.frame(
  state = tolower(state.name),
  long  = state.center$x,
  lat   = state.center$y,
  stringsAsFactors = FALSE
) %>% 
  # keep only those in our lon/lat window
  filter(long >= -100, long <= -75, lat >= 40, lat <= 50)

# 3) hulls + counts (as before)
counts <- df_clean %>%
  count(HCPC_cluster) %>%
  arrange(HCPC_cluster) %>%
  mutate(label = paste0("Cluster ", HCPC_cluster, " (n=", n, ")"))

cluster_labels <- setNames(counts$label, counts$HCPC_cluster)

hulls <- df_clean %>%
  group_by(HCPC_cluster) %>%
  slice(chull(Longitude, Latitude)) %>%
  ungroup()

aina_cols <- c(
  "1" = "darkred",
  "2" = "darkblue",
  "3" = "darkgreen",
  "4" = "purple4",
  "5" = "darkorange4"
)


# 4) plot
p <- ggplot() +
  geom_polygon(
    data   = states,
    aes(x = long, y = lat, group = group),
    fill   = "grey95", colour = "grey70", size = 0.3
  ) +
  geom_polygon(
    data   = hulls,
    aes(
      x     = Longitude,
      y     = Latitude,
      group = HCPC_cluster,
      fill  = factor(HCPC_cluster)
    ),
    alpha  = 0.15, colour = NA
  ) +
  geom_point(
    data = df_clean,
    aes(x = Longitude, y = Latitude, colour = factor(HCPC_cluster)),
    size  = 2, alpha = 0.8
  ) +
  geom_text(
    data    = state_labs,
    aes(x = long, y = lat, label = tools::toTitleCase(state)),
    size    = 3,
    colour  = "grey40"
  ) +
  coord_quickmap(xlim = c(-100, -75), ylim = c(40, 50)) +
  scale_fill_manual(
    values = aina_cols,
    labels = c(
      "1" = "Aina Cluster 1",
      "2" = "Aina Cluster 2",
      "3" = "Aina Cluster 3",
      "4" = "Aina Cluster 4",
      "5" = "Aina Cluster 5"
    ),
    drop = FALSE
  ) +
  scale_colour_manual(
    values = aina_cols,
    labels = c(
      "1" = "Aina Cluster 1",
      "2" = "Aina Cluster 2",
      "3" = "Aina Cluster 3",
      "4" = "Aina Cluster 4",
      "5" = "Aina Cluster 5"
    ),
    drop = FALSE
  ) +
  theme_minimal() +
  labs(
    title  = "Feral Cannabis Sample Locations by Aina Cluster",
    x      = "Longitude",
    y      = "Latitude",
    fill   = "Aina cluster",
    colour = "Aina cluster"
  ) +
  theme(
    panel.background = element_rect(fill = "aliceblue"),
    panel.grid       = element_line(colour = "white"),
    legend.position  = "right"
  )

p

ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/FigS9B.pdf",
  plot = p,
  width = 10,
  height = 8,
  units = "in"
)

########################################################################################################################################################################
# Figure S8
########################################################################################################################################################################
# FIGURE S8A/B
##################
library("FactoMineR")
library("factoextra")

#aina alone
pca <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_VCF_remake/VCF_overlays/PCA_merged.csv")

#aina+ren+soorni
#pca <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/newAina_Ren_Soorni_EGS/new_Aina_n760_with_core_flag.csv")


#Select 3rd and 4th principal components
#pca2 <- pca[3:4]
pca2 <- pca[3:6]

head(pca2)
# Standardizing the PCA data before clustering
my_data <- scale(pca2)

#elbow 
a <- fviz_nbclust(my_data, kmeans, method = "wss") + ggtitle("the Elbow Method")
a
ggsave("/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/Aina_elbow.pdf", plot = a, width = 6, height = 4, device = "pdf")

#silhouette
b<- fviz_nbclust(my_data, kmeans, method = "silhouette") + ggtitle("The Silhouette Plot")
b
ggsave("/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/Aina_silhouette.pdf", plot = b, width = 6, height = 4, device = "pdf")

##################
# FIGURE S8C
##################
#fastSTRUCTURE 
#AINA ALONE

library(MetBrewer)
library(ggplot2)
library(tidyr)
library(dplyr)

#install.packages("ggnewscale")
library(ggnewscale)

##############
#K2
##############
# Load the Hokusai palette
hiroshige_palette <- met.brewer("Hiroshige", 8)
k_colors <- hiroshige_palette[c(1, 6)]  # Colors for K1 and K2

# Load your data
test <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/fastSTRUCTURE/canna_k2.csv")

# Sort by K1 and K2 percentages
test <- test %>%
  arrange(desc(K1), desc(K2))

# Convert to long format after sorting
long_data <- pivot_longer(test, cols = c("K1", "K2"), names_to = "category", values_to = "value")

# Reorder the 'Name' factor based on the sorted data
long_data$Name <- factor(long_data$Name, levels = unique(test$Name))

# Assign specific colors to each Type
type_colors <- c(
  "Aina_cluster_1"      = "darkred",
  "Aina_cluster_2"      = "darkblue",
  "Aina_cluster_3"      = "darkgreen",
  "Aina_cluster_4"      = "purple4",
  "Aina_cluster_5"      = "darkorange4"
)


# Plot with the colored bar for Type underneath
final_plot <- ggplot() +
  # Plot the K1 and K2 bars with their own fill colors
  geom_bar(data = long_data, aes(x = Name, y = value, fill = category), stat = "identity") +
  scale_fill_manual(name = "K Categories", values = k_colors) +  # Colors for K1 and K2
  
  # Add the colored bars underneath for Type using a separate geom_tile
  new_scale_fill() +  # Necessary to separate the two scales
  geom_tile(data = test, aes(x = Name, y = -0.05, fill = Type), height = 0.1) +
  scale_fill_manual(name = "Type", values = type_colors) +  # Colors for Type
  
  labs(title = "", x = "", y = "Percent Identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 1),  # Make the text smaller
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 6),       # Smaller legend text
    legend.title = element_text(size = 8),      # Smaller legend title
    legend.key.size = unit(0.4, 'cm'),          # Smaller legend keys
    legend.spacing.y = unit(0.2, 'cm')          # Reduce vertical spacing between legend items
  )

print(final_plot)

# Save the plot as a PDF with reduced height
ggsave(filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/fastSTRUCTURE/fastSTRUCTURE_plot_k2.pdf", 
       plot = final_plot, 
       width = 14,    # Width of the plot (adjust as needed)
       height = 3)   # Reduced height (adjust as needed)

##############
#K3
##############
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggnewscale)
library(MetBrewer)

hiroshige_palette <- met.brewer("Hiroshige", 8)
k_colors <- hiroshige_palette[c(1, 4, 6)]

# Load your data
test <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/fastSTRUCTURE/canna_k3.csv")

# Create a new column for the dominant K category
test <- test %>%
  mutate(dominant_K = case_when(
    K1 >= K2 & K1 >= K3 ~ "K1",
    K2 >= K1 & K2 >= K3 ~ "K2",
    K3 >= K1 & K3 >= K2 ~ "K3"
  ))

# Sort by the dominant K category and then by the actual values in descending order
test <- test %>%
  arrange(dominant_K, desc(K1), desc(K2), desc(K3))

# Convert to long format after sorting
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3"), names_to = "category", values_to = "value")

# Reorder the 'Name' factor based on the sorted data
long_data$Name <- factor(long_data$Name, levels = unique(test$Name))

type_colors <- c(
  "Aina_cluster_1"      = "darkred",
  "Aina_cluster_2"      = "darkblue",
  "Aina_cluster_3"      = "darkgreen",
  "Aina_cluster_4"      = "purple4",
  "Aina_cluster_5"      = "darkorange4"
)

# Hiroshige palette for K categories
k_colors <- hiroshige_palette[c(1, 4, 6)]

# Plot with the colored bar for Type underneath
final_plot <- ggplot() +
  # Plot the K1 and K2 bars with their own fill colors
  geom_bar(data = long_data, aes(x = Name, y = value, fill = category), stat = "identity") +
  scale_fill_manual(name = "K Categories", values = k_colors) +  # Colors for K1 and K2
  
  # Add the colored bars underneath for Type using a separate geom_tile
  new_scale_fill() +  # Necessary to separate the two scales
  geom_tile(data = test, aes(x = Name, y = -0.05, fill = Type), height = 0.1) +
  scale_fill_manual(name = "Type", values = type_colors) +  # Colors for Type
  
  labs(title = "", x = "", y = "Percent Identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 1),  # Make the text smaller
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 6),       # Smaller legend text
    legend.title = element_text(size = 8),      # Smaller legend title
    legend.key.size = unit(0.4, 'cm'),          # Smaller legend keys
    legend.spacing.y = unit(0.2, 'cm')          # Reduce vertical spacing between legend items
  )

print(final_plot)

# Save the plot as a PDF with reduced height
ggsave(filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/fastSTRUCTURE/fastSTRUCTURE_plot_k3.pdf", 
       plot = final_plot, 
       width = 14,    # Width of the plot (adjust as needed)
       height = 3)   # Reduced height (adjust as needed)

##############
#K4
##############
#load fastSTRUCTURE data

hiroshige_palette <- met.brewer("Hiroshige", 8)
#selected_colors <- hiroshige_palette[c(1, 4,6, 8)]
k_colors <- hiroshige_palette[c(1,4,6,8)]  # Colors for K

test <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/fastSTRUCTURE/canna_k4.csv")

# Create a new column for the dominant K category
test <- test %>%
  mutate(dominant_K = case_when(
    K1 >= K2 & K1 >= K3 & K1 >= K4 ~ "K1",
    K2 >= K1 & K2 >= K3 & K2 >= K4 ~ "K2",
    K3 >= K1 & K3 >= K2 & K3 >= K4 ~ "K3",
    K4 >= K1 & K4 >= K2 & K4 >= K3 ~ "K4"
  ))

# Sort by the dominant K category and then by the actual values in descending order
test <- test %>%
  arrange(dominant_K, desc(K1), desc(K2), desc(K3), desc(K4))

# Convert to long format after sorting
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3", "K4"), names_to = "category", values_to = "value")

# Reorder the 'Name' factor based on the sorted data
long_data$Name <- factor(long_data$Name, levels = unique(test$Name))

# Assign specific colors to each Type and K category
type_colors <- c(
  "Aina_cluster_1"      = "darkred",
  "Aina_cluster_2"      = "darkblue",
  "Aina_cluster_3"      = "darkgreen",
  "Aina_cluster_4"      = "purple4",
  "Aina_cluster_5"      = "darkorange4"
)
# Hiroshige palette for K categories (adjust the palette if needed)
k_colors <- hiroshige_palette[c(1, 4, 6, 8)]  # Choose appropriate colors for K1 to K4

# Plot with the colored bar for Type underneath
final_plot <- ggplot() +
  # Plot the K1 and K2 bars with their own fill colors
  geom_bar(data = long_data, aes(x = Name, y = value, fill = category), stat = "identity") +
  scale_fill_manual(name = "K Categories", values = k_colors) +  # Colors for K1 and K2
  
  # Add the colored bars underneath for Type using a separate geom_tile
  new_scale_fill() +  # Necessary to separate the two scales
  geom_tile(data = test, aes(x = Name, y = -0.05, fill = Type), height = 0.1) +
  scale_fill_manual(name = "Type", values = type_colors) +  # Colors for Type
  
  labs(title = "", x = "", y = "Percent Identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 1),  # Make the text smaller
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 6),       # Smaller legend text
    legend.title = element_text(size = 8),      # Smaller legend title
    legend.key.size = unit(0.4, 'cm'),          # Smaller legend keys
    legend.spacing.y = unit(0.2, 'cm')          # Reduce vertical spacing between legend items
  )

print(final_plot)

# Save the plot as a PDF with reduced height
ggsave(filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/fastSTRUCTURE/fastSTRUCTURE_plot_k4.pdf", 
       plot = final_plot, 
       width = 14,    # Width of the plot (adjust as needed)
       height = 3)   # Reduced height (adjust as needed)

##############
#K5
##############
hiroshige_palette <- met.brewer("Hiroshige", 8)

test <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/fastSTRUCTURE/canna_k5.csv")

# Create a new column for the dominant K category
test <- test %>%
  mutate(dominant_K = case_when(
    K1 >= K2 & K1 >= K3 & K1 >= K4 & K1 >= K5 ~ "K1",
    K2 >= K1 & K2 >= K3 & K2 >= K4 & K2 >= K5 ~ "K2",
    K3 >= K1 & K3 >= K2 & K3 >= K4 & K3 >= K5 ~ "K3",
    K4 >= K1 & K4 >= K2 & K4 >= K3 & K4 >= K5 ~ "K4",
    K5 >= K1 & K5 >= K2 & K5 >= K3 & K5 >= K4 ~ "K5"
  ))

# Sort by the dominant K category and then by the actual values in descending order
test <- test %>%
  arrange(dominant_K, desc(K1), desc(K2), desc(K3), desc(K4), desc(K5))

# Convert to long format after sorting
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3", "K4", "K5"), names_to = "category", values_to = "value")

# Reorder the 'Name' factor based on the sorted data
long_data$Name <- factor(long_data$Name, levels = unique(test$Name))

# Assign specific colors to each Type and K category
type_colors <- c(
  "Aina_cluster_1"      = "darkred",
  "Aina_cluster_2"      = "darkblue",
  "Aina_cluster_3"      = "darkgreen",
  "Aina_cluster_4"      = "purple4",
  "Aina_cluster_5"      = "darkorange4"
)

# Hiroshige palette for K categories (adjust the palette if needed)
k_colors <- hiroshige_palette[c(1, 3, 5, 6, 8)]  # Choose appropriate colors for K1 to K5

# Plot with the colored bar for Type underneath
final_plot <- ggplot() +
  # Plot the K1 and K2 bars with their own fill colors
  geom_bar(data = long_data, aes(x = Name, y = value, fill = category), stat = "identity") +
  scale_fill_manual(name = "K Categories", values = k_colors) +  # Colors for K1 and K2
  
  # Add the colored bars underneath for Type using a separate geom_tile
  new_scale_fill() +  # Necessary to separate the two scales
  geom_tile(data = test, aes(x = Name, y = -0.05, fill = Type), height = 0.1) +
  scale_fill_manual(name = "Type", values = type_colors) +  # Colors for Type
  
  labs(title = "", x = "", y = "Percent Identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 1),  # Make the text smaller
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 6),       # Smaller legend text
    legend.title = element_text(size = 8),      # Smaller legend title
    legend.key.size = unit(0.4, 'cm'),          # Smaller legend keys
    legend.spacing.y = unit(0.2, 'cm')          # Reduce vertical spacing between legend items
  )

print(final_plot)

# Save the plot as a PDF with reduced height
ggsave(filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/fastSTRUCTURE/fastSTRUCTURE_plot_k5.pdf", 
       plot = final_plot, 
       width = 14,    # Width of the plot (adjust as needed)
       height = 3)   # Reduced height (adjust as needed)

##############
#K6
##############
#load fastSTRUCTURE data

# Load your data
test <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/fastSTRUCTURE/canna_k6.csv")

# Create a new column for the dominant K category
test <- test %>%
  mutate(dominant_K = case_when(
    K1 >= K2 & K1 >= K3 & K1 >= K4 & K1 >= K5 & K1 >= K6 ~ "K1",
    K2 >= K1 & K2 >= K3 & K2 >= K4 & K2 >= K5 & K2 >= K6 ~ "K2",
    K3 >= K1 & K3 >= K2 & K3 >= K4 & K3 >= K5 & K3 >= K6 ~ "K3",
    K4 >= K1 & K4 >= K2 & K4 >= K3 & K4 >= K5 & K4 >= K6 ~ "K4",
    K5 >= K1 & K5 >= K2 & K5 >= K3 & K5 >= K4 & K5 >= K6 ~ "K5",
    K6 >= K1 & K6 >= K2 & K6 >= K3 & K6 >= K4 & K6 >= K5 ~ "K6"
  ))

# Sort by the dominant K category and then by the actual values in descending order
test <- test %>%
  arrange(dominant_K, desc(K1), desc(K2), desc(K3), desc(K4), desc(K5), desc(K6))

# Convert to long format after sorting
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3", "K4", "K5", "K6"), names_to = "category", values_to = "value")

# Reorder the 'Name' factor based on the sorted data
long_data$Name <- factor(long_data$Name, levels = unique(test$Name))

# Assign specific colors to each Type and K category
type_colors <- c(
  "Aina_cluster_1"      = "darkred",
  "Aina_cluster_2"      = "darkblue",
  "Aina_cluster_3"      = "darkgreen",
  "Aina_cluster_4"      = "purple4",
  "Aina_cluster_5"      = "darkorange4"
)

# Hiroshige palette for K categories (adjust the palette if needed)
k_colors <- hiroshige_palette[c(1, 3, 5, 6, 7, 8)]  # Choose appropriate colors for K1 to K6

# Plot with the colored bar for Type underneath
final_plot <- ggplot() +
  # Plot the K1 and K2 bars with their own fill colors
  geom_bar(data = long_data, aes(x = Name, y = value, fill = category), stat = "identity") +
  scale_fill_manual(name = "K Categories", values = k_colors) +  # Colors for K1 and K2
  
  # Add the colored bars underneath for Type using a separate geom_tile
  new_scale_fill() +  # Necessary to separate the two scales
  geom_tile(data = test, aes(x = Name, y = -0.05, fill = Type), height = 0.1) +
  scale_fill_manual(name = "Type", values = type_colors) +  # Colors for Type
  
  labs(title = "", x = "", y = "Percent Identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 1),  # Make the text smaller
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 6),       # Smaller legend text
    legend.title = element_text(size = 8),      # Smaller legend title
    legend.key.size = unit(0.4, 'cm'),          # Smaller legend keys
    legend.spacing.y = unit(0.2, 'cm')          # Reduce vertical spacing between legend items
  )

print(final_plot)

# Save the plot as a PDF with reduced height
ggsave(filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/fastSTRUCTURE/fastSTRUCTURE_plot_k6.pdf", 
       plot = final_plot, 
       width = 14,    # Width of the plot (adjust as needed)
       height = 3)   # Reduced height (adjust as needed)


########################################################################################################################################################################
# Figure S9
########################################################################################################################################################################
# PCA by training and test
###########################
library(ggplot2)
library(readr)

# 1) Read the file
df <- read_csv(
  "/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/new_Aina_n760_with_core_flag.csv",
  col_types = cols(
    EV1          = col_double(),
    EV2          = col_double(),
    is_core_310  = col_logical(),
    .default     = col_guess()
  )
)

p<- ggplot(df, aes(x = EV1, y = EV2, color = is_core_310)) +
  geom_point(shape = 16, size = 2, alpha = 0.8) +     # shape 16 = solid circle
  scale_color_manual(
    values = c("TRUE"  = "red", 
               "FALSE" = "grey70"),
    labels = c("TRUE"  = "Training (305)",
               "FALSE" = "Test (455)"), 
    name   = ""
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x = "PC1 (1.95%)",
    y = "PC2 (1.38%)"
  ) +
  theme(
    legend.position    = "right",
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_line(color = "grey90"),
    panel.background   = element_rect(fill = "white", colour = NA),
    plot.background    = element_rect(fill = "white", colour = NA),
    panel.border       = element_blank()
  )
p
ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/FigS10.pdf",
  plot = p,
  width = 8,
  height = 6,
  units = "in"
)




########################################################################################################################################################################
# Figure S10
########################################################################################################################################################################
envdat <- read.csv('/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/new_Aina_n760_with_core_flag.csv', head = TRUE)


# Read training dataset from CSV (only core environmental data)
trainingset <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/new_cannabis_core_n310_greedy.csv", header = TRUE, stringsAsFactors = FALSE)

library(dplyr)

# Ensure ID columns have same type
envdat$sample.id <- as.character(envdat$sample.id)
trainingset$Sample_ID <- as.character(trainingset$Sample_ID)

# Partition envdat into core and non-core
env_core <- envdat %>%
  semi_join(trainingset, by = c("sample.id" = "Sample_ID")) %>%
  filter(is_core_310)   # keep only TRUE flag

env_noncore <- envdat %>%
  anti_join(trainingset, by = c("sample.id" = "Sample_ID"))

write.csv(env_core,
          "/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/new_cannabis_core_n310_greedy_data.csv",
          row.names = FALSE)

######################################################
# Prediction Accuracy
######################################################
library(rrBLUP)  # RR-BLUP package for Ridge Regression Best Linear Unbiased Prediction
library(hibayes) # hibayes package for Bayesian models
library(dplyr)

# Load custom functions for k-fold cross-validation
source("/Users/annamccormick/R/Feral_Cannabis_EGS/EGS/xval_kfold_functions.R")

# Read genotype data in rrBLUP format
#data <- read.delim("~/R/cannabis_GEAV/Inputs/Ren_PRJNA734114_rrBLUP_clean.txt", sep = "\t", header = TRUE)
data <- read.delim("/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/Aina_nn760_rrBLUP_format_hmp.txt", sep = "\t", header = TRUE)

#adjust sample names to match metadata layouts
orig_cols    <- colnames(data)[1:4]
sample_cols  <- colnames(data)[-c(1:4)]
fixed_samples <- gsub("\\.", "-", sample_cols)
colnames(data) <- c(orig_cols, fixed_samples)

# Read environmental data from CSV
#envdat <- read.csv('~/R/cannabis_GEAV/Outputs/ren_n44_extracted_climate_data_all.csv', head = TRUE)  # full environmental dataset
envdat <- read.csv('/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/new_Aina_n760_with_core_flag.csv', head = TRUE)


# Read training dataset from CSV (only core environmental data)
#trainingset <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/new_cannabis_core_n310_greedy.csv", header = TRUE, stringsAsFactors = FALSE)
trainingset <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/new_cannabis_core_n310_greedy_data.csv", header = TRUE, stringsAsFactors = FALSE)

###
head(colnames(data), 20)

#first few sample IDs in your envdat
head(envdat$sample.id, 20)

# overlap
common_ids <- intersect(colnames(data), envdat$sample.id)
length(common_ids)
common_ids


# Set the row names for the environmental and training datasets
row.names(envdat) <- envdat$sample.id
row.names(trainingset) <- trainingset$sample.id

# Filter genotype data to include only columns matching Sample_ID in envdat
gd2 <- data[, c("rs.", "allele", "chrom", "pos", colnames(data)[colnames(data) %in% envdat$sample.id])]

# Set the SNP IDs as the row names for the genotype data
row.names(gd2) <- gd2$rs.

# Remove the first four columns (rs., allele, chrom, pos) to create the genotype matrix
g.in <- gd2[, -c(1:4)]  # Isolate genotype data
g.in <- as.matrix(g.in)
g.in <- t(g.in)  # Transpose to align individuals as rows
cat("Dimensions of g.in", dim(g.in), "\n")

########################
#Remove snps with greater then 20% missingness
#g.in <- g.in[, colSums(is.na(g.in)) == 0]  # Remove SNPs with missing data
g.in <- g.in[, colMeans(is.na(g.in)) <= 0.20]
# Check dimensions of the cleaned genotype matrix
cat("Dimensions of g.in after removing SNPs with missing data:", dim(g.in), "\n")


## where to insert imputation with mean appropriatly - matching the lables used this time
# missingness BEFORE imputation

# overall missing‐data rate, before train/pred split #All data: 871622 / 8090880 missing (10.77%)
total_missing <- sum(is.na(g.in))
total_cells   <- prod(dim(g.in))
cat(sprintf(
  "All data: %d / %d missing (%.2f%%)\n",
  total_missing, total_cells, 100 * total_missing / total_cells
))

# 2) Impute remaining missing values with the SNP (column) mean
marker_means <- colMeans(g.in, na.rm = TRUE)
for (j in seq_len(ncol(g.in))) {
  nas <- is.na(g.in[, j])
  if (any(nas)) {
    g.in[nas, j] <- marker_means[j]
  }
}
cat("Total missing after imputation:", sum(is.na(g.in)), "\n")



### Load Environmental Data ###
# Select unique identifier and environmental data columns for training set
#y.trainset <- trainingset[, c(1, 2:20:(ncol(trainingset)))]  # Select unique identifier and environmental data only
# Select column 1 plus columns 2 through 20
y.trainset <- trainingset[, c(1, 2:20)]
# Move sample.id into row names
rownames(y.trainset) <- y.trainset$sample.id


# Select unique identifier and environmental data columns for the full dataset
#y.in <- envdat[, c(1, 14:ncol(envdat))]
y.in <- envdat[, c(2, 13:ncol(envdat))]
# remove last column



# Convert environmental data to matrix format for RR-BLUP
y.trainset.mat <- as.matrix(y.trainset[, -1])  # Exclude the first column (Sample_ID)
y.in.mat <- as.matrix(y.in[, -1])  # Exclude the first column (Sample_ID)
y.in.mat <- y.in.mat[, -ncol(y.in.mat)]
dim(y.in.mat)


# Ensure row alignment between genotype and environmental data
# Ensure the same for the training set

# Check dimensions before cross-validation
cat("Dimensions of g.in:", dim(g.in), "\n")
cat("Dimensions of y.in.mat:", dim(y.in.mat), "\n")
cat("Dimensions of y.trainset.mat:", dim(y.trainset.mat), "\n")

# RR-BLUP ###
# k.xval function from external script, running 10-fold cross-validation with 50 repetitions
xval_k10_rrblup <- k.xval(g.in = g.in, y.in = y.in.mat, y.trainset = y.trainset.mat, k.fold = 10, reps = 50)
saveRDS(xval_k10_rrblup, "/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/xval_rrblup_kfold_10.RData")


########################
# Gaussian Kernel
K <- A.mat(g.in) #Compute relationship matrix using Gaussian Kernel
k_dist <- dist(K) #Calculate distance matrix from relationship matrix

#Run k-fold cross-validation for Gaussian Kernel
xval_k10_GAUSS <- k.xval.GAUSS(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k_dist = k_dist, k.fold = 10, reps = 50)
saveRDS(xval_k10_GAUSS, "xval_GAUSS_kfold_10.RData")

#Exponential Kernel
#Compute relationship matrix using Exponential Kernel
K.Exp=Kernel_computation(X=g.in, name="exponential", degree=NULL, nL=NULL)

#Set row and column names for the Exponential Kernel matrix
row.names(K.Exp) <- rownames(g.in)
colnames(K.Exp) <- rownames(g.in)

#Calculate distance matrix from Exponential Kernel relationship matrix
exp_dist <- dist(K.Exp) # Calculate Relationship Matrix\

#Run k-fold cross-validation for Exponential Kernel
xval_k10_EXP <- k.xval.EXP(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k_dist = exp_dist, k.fold = 10, reps = 50)
saveRDS(xval_k10_EXP, "xval_EXP_kfold_10.RData")

#BayesCPi
#Run k-fold cross-validation for BayesCPi model
xval_k10_BayesCpi <- k.xval.BayesCpi(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k.fold = 10, reps = 50, niter=3000,nburn=1200)
saveRDS(xval_k10_BayesCpi, "xval_BayesCpi_kfold_10.RData")



#####
#Figure S10 - PLOT
#####
library(ggplot2)
setwd("/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/")

# RR-BLUP - K-Fold Cross-Validation
#Load cross-validation results for rrBLUP
rrblup_kfold10 <- readRDS("xval_rrblup_kfold_10.RData")
rrblup_kfold10 <- rrblup_kfold10$xval.result
rrblup_kfold10$r.mean <- as.numeric(rrblup_kfold10$r.mean)

# Gaussian Kernel - K-fold
#Load cross-validation results for Gaussian Kernel
gauss_kfold_10 <- readRDS('xval_GAUSS_kfold_10.RData')
gauss_kfold_10 <- gauss_kfold_10$xval.result
gauss_kfold_10$r.mean <- as.numeric(gauss_kfold_10$r.mean)

# Exponential Kernel - K-fold
#Load cross-validation results for Exponential Kernel
EXP_kfold_10 <- readRDS('xval_EXP_kfold_10.RData')
EXP_kfold_10 <- EXP_kfold_10$xval.result
EXP_kfold_10$r.mean <- as.numeric(EXP_kfold_10$r.mean)

# BayesCpi - K-fold
#Load cross-validation results for BayesCpi model
bayescpi_kfold_10 <- readRDS('xval_BayesCpi_kfold_10.RData')
bayescpi_kfold_10 <- bayescpi_kfold_10$xval.result
bayescpi_kfold_10$r.mean <- as.numeric(bayescpi_kfold_10$r.mean)


# Organize all dataframes for merging
## Rename model names
rrblup_kfold10$model <- "rrBLUP"
gauss_kfold_10$model <- "Gaussian Kernel"
EXP_kfold_10$model <- "Exponential Kernel"
bayescpi_kfold_10$model <- "BayesCpi"

## Input xval type
rrblup_kfold10$xval <- "Ten-Fold"
gauss_kfold_10$xval <- "Ten-Fold"
EXP_kfold_10$xval <- "Ten-Fold"
bayescpi_kfold_10$xval <- "Ten-Fold"


#Combine all model results into a single list
model_list <- list(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10,bayescpi_kfold_10)


#Remove any NA values from the model results
model_list1 <- lapply(model_list, na.omit)

#Combine all models into a single dataframe
all_models <- do.call("rbind", model_list1)

#Convert standard deviation values to numeric
all_models$r.sd <- as.numeric(all_models$r.sd)

#Ensure cross-validation type is a factor with specified levels
all_models$xval <- factor(all_models$xval, levels = c("Ten-Fold"))

#Filter for specific traits of interest
all_bio <- all_models[all_models$trait %in% c('bio_01', 'bio_02', 'bio_03', 'bio_04', 'bio_05','bio_06',
                                              'bio_07','bio_08','bio_09','bio_010','bio_11',
                                              'bio_12','bio_13','bio_14','bio_15','bio_16','bio_17','bio_18','bio_19'),]

#Plot 
all_bio$trait <- factor(all_bio$trait, levels = paste0("bio", 1:19))

p<- ggplot(all_models, aes(y = r.mean, x = model, color = model)) +
  theme_bw() +
  geom_errorbar(aes(x = model, ymin = r.mean-r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 1) + 
  geom_hline(yintercept = 0.85, color="red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1)
p

ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/FigS11.pdf",
  plot = p,
  width = 15,
  height = 10,
  units = "in"

######################################################
#Extract relevant environmental data from WorldClim for latitude and longitudes 
######################################################
library(raster)
worldclim_path <- "~/R/Data/wc2.0_30s_bio/"
bio_files <- paste0(worldclim_path, "bio_", sprintf("%02d", 1:19), ".tif")
climate_layers <- stack(bio_files)

# Load your latitude and longitude dataset
data <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/PCA_filtered.csv")
head(data)

# Prepare the coordinates (assuming Latitude and Longitude columns are present)
coords <- data.frame(lon = data$Longitude, lat = data$Latitude)

# Extract climate data for these coordinates
climate_values <- extract(climate_layers, coords)

# Combine the extracted climate data with your original dataset
result <- cbind(data, climate_values)

# Save the result to a CSV file
write.csv(result, "/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/Aina_n742_extracted_climate_data.csv", row.names = FALSE)

######################################################
#adding a column for the 310 training in the above file for supplemental tables
library(readr)
library(dplyr)

# 1) Read full metadata
meta <- read_csv(
  "/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/Aina_n742_extracted_climate_data.csv",
  col_types = cols(.default = col_character())
)

# 2) Read core‑310 list (VCF IDs)
core <- read_csv(
  "/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/new_cannabis_core_n310_greedy.csv",
  col_types = cols(Sample_ID = col_character())
)

# 3) Flag using sample.id (the VCF names) instead of sample_final
meta_flagged <- meta %>%
  mutate(
    is_core_310 = sample.id %in% core$Sample_ID
  )

# 4) Verify
cat("Flagged core:", sum(meta_flagged$is_core_310, na.rm=TRUE), "of", nrow(meta_flagged), "\n")

# 5) Save out
write_csv(
  meta_flagged,
  "/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/new_Aina_n760_with_core_flag.csv"
)



######################################################
# EGS
######################################################
# Step 5: Genomic Selection with RR-BLUP
######################################################
library(rrBLUP)
library(dplyr)

# Load genotype and environmental data
gd1 <- read.delim("/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/Aina_nn760_rrBLUP_format_hmp.txt", sep = "\t", header = TRUE)
envdat <- read.csv('/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/new_Aina_n760_with_core_flag.csv', head = TRUE)

# Set row names for genotype data (use 'rs.' column which contains SNP IDs)
row.names(gd1) <- gd1$rs.
head(rownames(gd1))
head(gd1$rs.)

# Remove the non-genotype columns ('rs.', 'allele', 'chrom', 'pos') for analysis
gd3 <- gd1[, -c(1:4)]  # SNP columns only
g.in <- as.matrix(gd3)  # Convert to matrix
g.in <- t(g.in)  # Transpose g.in to align with g.train structure (SNPs as columns, individuals as rows)
print(g.in [1:10, 1:10])  # Visualize 

# Set row names for environmental data using 'Sample_ID'
row.names(envdat) <- envdat$sample.id

# Subset environmental data for analysis (start from the 5th column)
y.in.rr <- envdat[, 13:31]  #14-32
minicore_entries <- which(envdat$is_core_310 == TRUE)  # Subset to core lines where 'Core' is TRUE
y.in.rr <- y.in.rr[minicore_entries, ]  
y.in.mat <- as.matrix(y.in.rr)

print(y.in.mat [1:10, 1:10])  # Visualize 
# Now fix the rownames to match the dots in g.in
rownames(y.in.mat) <- gsub("-", "\\.", rownames(y.in.mat))
print(y.in.mat [1:10, 1:10])  # Visualize 

##########################
# Training set
train <- row.names(y.in.mat)  # Names of training lines

# fix labelling issue # replace every "-" with "."  
train_dots <- gsub("-", "\\.", train)
train <- train_dots

# 2. now you can subset
#g.train <- g.in[ train_dots, ] # Subset rows (individuals) directly
g.train <- g.in[ train, ] # Subset rows (individuals) directly

print(g.train[1:10, 1:10])  # Visualize part of g.train for debugging
dim(g.train)


# Prediction set 
#pred    <- setdiff(row.names(g.in), train_dots) 
pred    <- setdiff(row.names(g.in), train)
g.pred  <- g.in[pred, ]   
dim(g.pred)

########################################################################################################################
# g.train is 310 × M (M = number of SNPs before filtering)
miss_rate <- colMeans(is.na(g.train))
summary(miss_rate)
hist(miss_rate, breaks=50,
     main="Missingness per SNP", xlab="Proportion missing")

# Compute per‑SNP missing rates
miss_rate <- colMeans(is.na(g.train))

# Identify SNPs to keep less then ≤ 20% missing. If greater then 20% missingness discarded
keep_snps <- names(miss_rate)[miss_rate <= 0.20]
cat("After missingness filter, SNP count =", length(keep_snps), "\n")

# Subset both training and prediction matrices to those SNPs
g.train <- g.train[, keep_snps]
g.pred  <- g.pred[,  keep_snps]

dim(g.train)
dim(g.pred)



############################################################
# 2) Impute the remaining missing genotypes with the marker mean
############################################################
# missingness BEFORE imputation

# overall missing‐data rate, before train/pred split
total_missing <- sum(is.na(g.in))
total_cells   <- prod(dim(g.in))
cat(sprintf(
  "All data: %d / %d missing (%.2f%%)\n",
  total_missing, total_cells, 100 * total_missing / total_cells
))


# Function to mean‑impute each column
mean_impute <- function(mat) {
  for (j in seq_len(ncol(mat))) {
    na_idx <- is.na(mat[,j])
    if (any(na_idx)) {
      mu <- mean(mat[,j], na.rm = TRUE)
      # if your algorithm wants integer codes (–1,0,1), you could round(mu)
      mat[na_idx, j] <- mu
    }
  }
  mat
}

g.train <- mean_impute(g.train)
g.pred  <- mean_impute(g.pred)

# sanity check
stopifnot(!any(is.na(g.train)))
stopifnot(!any(is.na(g.pred)))
cat("No missing values remain in g.train or g.pred\n")

############################################################
# RR‑BLUP loop
# List of traits to analyze
traits <- colnames(y.in.mat)

# Initialize object for storing results
gebv_df <- data.frame(matrix(nrow = nrow(g.in), ncol = length(traits)))
colnames(gebv_df) <- traits
row.names(gebv_df) <- row.names(g.in)

# RR-BLUP loop for each trait
for (t in 1:length(traits)) {
  trait <- traits[t]
  
  # Set up training set for the trait
  y.train <- as.matrix(y.in.mat[train, trait])  
  
  # Run RR-BLUP model
  solve.out <- mixed.solve(y = y.train, Z = g.train, SE = FALSE, return.Hinv = FALSE)
  u.hat <- solve.out$u
  
  # Calculate GEBVs for both prediction and training sets
  GEBV <- g.pred %*% u.hat
  GEBV_train <- g.train %*% u.hat
  
  # Store results in the combined dataframe
  pred_rows <- match(row.names(g.pred), row.names(g.in))  # Indices for pred rows
  train_rows <- match(train, row.names(g.in))  # Indices for train rows
  
  gebv_df[pred_rows, t] <- GEBV  # Predictions for test lines
  gebv_df[train_rows, t] <- GEBV_train  # Predictions for training lines
}

# Final output: GEBVs for all individuals
write.csv(gebv_df, '/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/new_rrblup_GEAV_values_feral_cannabis.csv')




######################################################
#Join in metadata
library(dplyr)
library(tibble)    # for rownames_to_column()

# Make the rownames into an explicit column…
gebv2 <- gebv_df %>% 
  rownames_to_column("sample.id")

# Read in your metadata, and join
meta <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/new_Aina_n760_with_core_flag.csv",
                 stringsAsFactors = FALSE)

# make sure meta$sample.id really has the hyphens
# head(meta$sample.id)

final <- meta %>% 
  left_join(gebv2, by = "sample.id")

# Inspect and write out
head(final)
write.csv(final,
          "/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/rrblup_with_metadata.csv",
          row.names = FALSE)


######################################################
# PLOT #ALL 19 PLOT
######################################################
library(dplyr)    
library(tidyr) 
library(tidyverse)
library(ggplot2)

# 1) Read your combined data
df <- read_csv(
  "/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/rrblup_with_metadata_copy.csv",
  guess_max = 10000
)

# 1) define your lookup  
bio_names <- c(
  bio_01 = "Annual Mean Temperature (bio 01)",
  bio_02 = "Mean Diurnal Range (bio 02)",
  bio_03 = "Isothermality (bio 03)",
  bio_04 = "Temperature Seasonality (bio 04)",
  bio_05 = "Max Temperature of Warmest Month (bio 05)",
  bio_06 = "Min Temperature of Coldest Month (bio 06)",
  bio_07 = "Temperature Annual Range (bio 07)",
  bio_08 = "Mean Temperature of Wettest Quarter (bio 08)",
  bio_09 = "Mean Temperature of Driest Quarter (bio 09)",
  bio_10 = "Mean Temperature of Warmest Quarter (bio 10)",
  bio_11 = "Mean Temperature of Coldest Quarter (bio 11)",
  bio_12 = "Annual Precipitation (bio 12)",
  bio_13 = "Precipitation of Wettest Month (bio 13)",
  bio_14 = "Precipitation of Driest Month (bio 14)",
  bio_15 = "Precipitation Seasonality (bio 15)",
  bio_16 = "Precipitation of Wettest Quarter (bio 16)",
  bio_17 = "Precipitation of Driest Quarter (bio 17)",
  bio_18 = "Precipitation of Warmest Quarter (bio 18)",
  bio_19 = "Precipitation of Coldest Quarter (bio 19)"
)

# 2) pivot‐long as before, then re-factor trait
long <- df %>%
  dplyr::select(sample.id, HCPC_cluster, dplyr::starts_with("bio_")) %>%
  tidyr::pivot_longer(
    cols      = dplyr::starts_with("bio_"),
    names_to  = "trait",
    values_to = "GEBV"
  ) %>%
  dplyr::mutate(
    HCPC_cluster = factor(HCPC_cluster),
    trait = gsub("\\.y$", "", trait),   # remove the .x
    trait = factor(trait, levels = paste0("bio_", sprintf("%02d", 1:19)))
  )

# 3) full ggplot
ggplot(long, aes(x = HCPC_cluster, y = GEBV, fill = HCPC_cluster)) +
  # zero line
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.5) +
  # boxplots
  geom_boxplot(outlier.size = 1, alpha = 0.7) +
  # facet by human-readable trait names
  facet_wrap(~ trait, scales = "free_y", ncol = 4) +
  # custom colors & legend labels
  scale_fill_manual(
    name   = "HCPC cluster",
    values = c(
      "1" = "#E41A1C",
      "2" = "#377EB8",
      "3" = "#4DAF4A",
      "4" = "#984EA3",
      "5" = "#FF7F00"
    ),
    labels = c(
      "1" = "Cluster 1 (n=20)",
      "2" = "Cluster 2 (n=342)",
      "3" = "Cluster 3 (n=326)",
      "4" = "Cluster 4 (n=440)",
      "5" = "Cluster 5 (n=14)"
    )
  ) +
  labs(
    x = "HCPC Cluster",
    y = "Predicted GEBV"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position  = "right",
    strip.background = element_rect(fill = "grey90", colour = NA),
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
# 1) pivot-long and clean trait names
long <- df %>%
  select(sample.id, HCPC_cluster, starts_with("bio_")) %>%
  pivot_longer(
    cols      = starts_with("bio_"),
    names_to  = "trait",
    values_to = "GEBV"
  ) %>%
  mutate(
    HCPC_cluster = factor(HCPC_cluster),
    trait = gsub("\\.(x|y)$", "", trait),                      # remove .x or .y
    trait = factor(trait, levels = paste0("bio_", sprintf("%02d", 1:19)))
  )

# 2) plot with a labeller that maps trait -> bio_names
p<- ggplot(long, aes(x = HCPC_cluster, y = GEBV, fill = HCPC_cluster)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.5) +
  geom_boxplot(outlier.size = 1, alpha = 0.7) +
  facet_wrap(
    ~ trait, scales = "free_y", ncol = 4,
    labeller = labeller(trait = as_labeller(bio_names, default = identity))
  ) +
  scale_fill_manual(
    name   = "HCPC cluster",
    values = c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00"),
    labels = c("1"="Cluster 1 (n=20)","2"="Cluster 2 (n=342)","3"="Cluster 3 (n=326)",
               "4"="Cluster 4 (n=440)","5"="Cluster 5 (n=14)")
  ) +
  labs(x = "HCPC Cluster", y = "Predicted GEAV") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        strip.background = element_rect(fill = "grey90", colour = NA),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/FigS12.pdf",
  plot = p,
  width = 15,
  height = 10,
  units = "in"
)





########################################################################################################################################################################
# Figure S11
########################################################################################################################################################################

#  Boxplots of bio variables by Aina cluster
######################################################

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

#-----------------------------
# 1) Read data
#-----------------------------
df <- read_csv(
  "/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/rrblup_with_metadata_copy.csv",
  guess_max = 10000
)

#-----------------------------
# 2) Trait name lookup
#-----------------------------
bio_names <- c(
  bio_01 = "Annual Mean Temperature (bio 01)",
  bio_02 = "Mean Diurnal Range (bio 02)",
  bio_03 = "Isothermality (bio 03)",
  bio_04 = "Temperature Seasonality (bio 04)",
  bio_05 = "Max Temperature of Warmest Month (bio 05)",
  bio_06 = "Min Temperature of Coldest Month (bio 06)",
  bio_07 = "Temperature Annual Range (bio 07)",
  bio_08 = "Mean Temperature of Wettest Quarter (bio 08)",
  bio_09 = "Mean Temperature of Driest Quarter (bio 09)",
  bio_10 = "Mean Temperature of Warmest Quarter (bio 10)",
  bio_11 = "Mean Temperature of Coldest Quarter (bio 11)",
  bio_12 = "Annual Precipitation (bio 12)",
  bio_13 = "Precipitation of Wettest Month (bio 13)",
  bio_14 = "Precipitation of Driest Month (bio 14)",
  bio_15 = "Precipitation Seasonality (bio 15)",
  bio_16 = "Precipitation of Wettest Quarter (bio 16)",
  bio_17 = "Precipitation of Driest Quarter (bio 17)",
  bio_18 = "Precipitation of Warmest Quarter (bio 18)",
  bio_19 = "Precipitation of Coldest Quarter (bio 19)"
)

#-----------------------------
# 3) Aina cluster colours + legend labels
#-----------------------------
aina_cols <- c(
  "1" = "darkred",
  "2" = "darkblue",
  "3" = "darkgreen",
  "4" = "purple4",
  "5" = "darkorange4"
)

# If you want n's in the legend, compute them from THIS dataset
counts <- df %>%
  filter(!is.na(HCPC_cluster)) %>%
  count(HCPC_cluster) %>%
  mutate(
    HCPC_cluster = as.character(HCPC_cluster),
    legend_lab = paste0("Aina Cluster ", HCPC_cluster, " (n=", n, ")")
  )

aina_labels <- setNames(counts$legend_lab, counts$HCPC_cluster)

#-----------------------------
# 4) Long format for plotting
#-----------------------------
long <- df %>%
  select(sample.id, HCPC_cluster, starts_with("bio_")) %>%
  pivot_longer(
    cols = starts_with("bio_"),
    names_to = "trait",
    values_to = "GEBV"
  ) %>%
  mutate(
    HCPC_cluster = factor(as.character(HCPC_cluster), levels = names(aina_cols)),
    trait = gsub("\\.(x|y)$", "", trait),
    trait = factor(trait, levels = paste0("bio_", sprintf("%02d", 1:19)))
  )

#-----------------------------
# 5) Plot
#-----------------------------
p <- ggplot(long, aes(x = HCPC_cluster, y = GEBV, fill = HCPC_cluster)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
  geom_boxplot(outlier.size = 1, alpha = 0.7) +
  facet_wrap(
    ~ trait, scales = "free_y", ncol = 4,
    labeller = labeller(trait = as_labeller(bio_names, default = identity))
  ) +
  scale_fill_manual(
    values = aina_cols,
    labels = aina_labels,
    name   = "Aina cluster",
    drop   = FALSE
  ) +
  labs(
    x = "Aina Cluster",
    y = "Predicted GEBV"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position  = "right",
    strip.background = element_rect(fill = "grey90", colour = NA),
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p)

#-----------------------------
# 6) Save
#-----------------------------
ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/FigS12.pdf",
  plot = p,
  width = 18,
  height = 10,
  units = "in",
  device = "pdf"
)



########################################################################################################################################################################
# Figure S12
########################################################################################################################################################################
#Aina only passes 0.85
library(dplyr)    
library(tidyr) 
library(tidyverse)

# 1) Read your combined data
df <- read_csv(
  "/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/rrblup_with_metadata_copy.csv",
  guess_max = 10000
)


# your human‑readable lookup
bio_names <- c(
  bio_01 = "Annual Mean Temperature (bio 01)",
  bio_02 = "Mean Diurnal Range (bio 02)",
  bio_03 = "Isothermality (bio 03)",
  bio_04 = "Temperature Seasonality (bio 04)",
  bio_05 = "Max Temperature of Warmest Month (bio 05)",
  bio_06 = "Min Temperature of Coldest Month (bio 06)",
  bio_07 = "Temperature Annual Range (bio 07)",
  bio_08 = "Mean Temperature of Wettest Quarter (bio 08)",
  bio_09 = "Mean Temperature of Driest Quarter (bio 09)",
  bio_10 = "Mean Temperature of Warmest Quarter (bio 10)",
  bio_11 = "Mean Temperature of Coldest Quarter (bio 11)",
  bio_12 = "Annual Precipitation (bio 12)",
  bio_13 = "Precipitation of Wettest Month (bio 13)",
  bio_14 = "Precipitation of Driest Month (bio 14)",
  bio_15 = "Precipitation Seasonality (bio 15)",
  bio_16 = "Precipitation of Wettest Quarter (bio 16)",
  bio_17 = "Precipitation of Driest Quarter (bio 17)",
  bio_18 = "Precipitation of Warmest Quarter (bio 18)",
  bio_19 = "Precipitation of Coldest Quarter (bio 19)"
)

# only these codes
keep_codes <- c("bio_01","bio_03","bio_04","bio_05","bio_06","bio_07",
                "bio_08","bio_09","bio_10","bio_11","bio_14",
                "bio_15","bio_17","bio_19")

long_sub <- df %>%
  dplyr::select(sample.id, HCPC_cluster, dplyr::starts_with("bio_")) %>%
  tidyr::pivot_longer(
    cols      = dplyr::starts_with("bio_"),
    names_to  = "trait",
    values_to = "GEBV"
  ) %>%
  # clean trait names and filter
  dplyr::mutate(trait = sub("(\\.?(x|y))$", "", trait)) %>%
  dplyr::filter(trait %in% keep_codes) %>%
  dplyr::mutate(
    HCPC_cluster = factor(HCPC_cluster),
    trait = factor(
      trait,
      levels = keep_codes,
      labels = bio_names[keep_codes]
    )
  )


p<- ggplot(long_sub, aes(x = HCPC_cluster, y = GEBV, fill = HCPC_cluster)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_boxplot(outlier.size = 1, alpha = 0.7) +
  facet_wrap(~ trait, scales = "free_y", ncol = 3) +  # 2 columns to fit 5 panels
  scale_fill_manual(
    name   = "HCPC cluster",
    values = c(
      "1" = "#E41A1C",
      "2" = "#377EB8",
      "3" = "#4DAF4A",
      "4" = "#984EA3",
      "5" = "#FF7F00"
    ),
    labels = c(
      "1" = "Cluster 1 (n=20)",
      "2" = "Cluster 2 (n=342)",
      "3" = "Cluster 3 (n=326)",
      "4" = "Cluster 4 (n=440)",
      "5" = "Cluster 5 (n=14)"
    )
  ) +
  labs(
    x = "HCPC Cluster",
    y = "Predicted GEAV"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position  = "right",
    strip.background = element_rect(fill = "grey90", colour = NA),
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )


p

ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/Fig3.pdf",
  plot = p,
  width = 15,
  height = 10,
  units = "in"
)


#fix colour scheme

#############################
# PLOT (passes 0.85)
# Aina colour scheme + Aina legend labels
#############################

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

#-----------------------------
# 1) Read data
#-----------------------------
df <- read_csv(
  "/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/rrblup_with_metadata_copy.csv",
  guess_max = 10000
)

#-----------------------------
# 2) Human-readable BIO lookup
#-----------------------------
bio_names <- c(
  bio_01 = "Annual Mean Temperature (bio 01)",
  bio_02 = "Mean Diurnal Range (bio 02)",
  bio_03 = "Isothermality (bio 03)",
  bio_04 = "Temperature Seasonality (bio 04)",
  bio_05 = "Max Temperature of Warmest Month (bio 05)",
  bio_06 = "Min Temperature of Coldest Month (bio 06)",
  bio_07 = "Temperature Annual Range (bio 07)",
  bio_08 = "Mean Temperature of Wettest Quarter (bio 08)",
  bio_09 = "Mean Temperature of Driest Quarter (bio 09)",
  bio_10 = "Mean Temperature of Warmest Quarter (bio 10)",
  bio_11 = "Mean Temperature of Coldest Quarter (bio 11)",
  bio_12 = "Annual Precipitation (bio 12)",
  bio_13 = "Precipitation of Wettest Month (bio 13)",
  bio_14 = "Precipitation of Driest Month (bio 14)",
  bio_15 = "Precipitation Seasonality (bio 15)",
  bio_16 = "Precipitation of Wettest Quarter (bio 16)",
  bio_17 = "Precipitation of Driest Quarter (bio 17)",
  bio_18 = "Precipitation of Warmest Quarter (bio 18)",
  bio_19 = "Precipitation of Coldest Quarter (bio 19)"
)

# Only these BIO codes
keep_codes <- c(
  "bio_01","bio_03","bio_04","bio_05","bio_06","bio_07",
  "bio_08","bio_09","bio_10","bio_11","bio_14",
  "bio_15","bio_17","bio_19"
)

#-----------------------------
# 3) Aina colour palette + legend labels
#-----------------------------
aina_cols <- c(
  "1" = "darkred",
  "2" = "darkblue",
  "3" = "darkgreen",
  "4" = "purple4",
  "5" = "darkorange4"
)

# Option A: hard-code legend labels (as you used before)
aina_labels <- c(
  "1" = "Aina Cluster 1 (n=20)",
  "2" = "Aina Cluster 2 (n=342)",
  "3" = "Aina Cluster 3 (n=326)",
  "4" = "Aina Cluster 4 (n=440)",
  "5" = "Aina Cluster 5 (n=14)"
)


long_sub <- df %>%
  select(sample.id, HCPC_cluster, starts_with("bio_")) %>%
  pivot_longer(
    cols      = starts_with("bio_"),
    names_to  = "trait",
    values_to = "GEBV"
  ) %>%
  mutate(
    trait = sub("(\\.?(x|y))$", "", trait)
  ) %>%
  filter(trait %in% keep_codes) %>%
  mutate(
    HCPC_cluster = factor(as.character(HCPC_cluster), levels = names(aina_cols)),
    trait = factor(trait, levels = keep_codes, labels = bio_names[keep_codes])
  )

#-----------------------------
# 5) Plot
#-----------------------------
p <- ggplot(long_sub, aes(x = HCPC_cluster, y = GEBV, fill = HCPC_cluster)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_boxplot(outlier.size = 1, alpha = 0.7) +
  facet_wrap(~ trait, scales = "free_y", ncol = 3) +
  scale_fill_manual(
    values = aina_cols,
    labels = aina_labels,
    name   = "Aina cluster",
    drop   = FALSE
  ) +
  labs(
    x = "Aina Cluster",
    y = "Predicted GEAV"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position  = "right",
    strip.background = element_rect(fill = "grey90", colour = NA),
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p)

#-----------------------------
# 6) Save PDF (edit filename as needed)
#-----------------------------
ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/Figure3.pdf",
  plot = p,
  width = 18,
  height = 10,
  units = "in",
  device = "pdf"
)





########################################################################################################################################################################
# Figure S13 (Supplemental Figure 13)
########################################################################################################################################################################

##########################
# GS for chemistry
##########################

library(readxl)
library(dplyr)
#library(writexl)
library(readr)
# 1. Load data
chem <- read_excel("Aina_chemistry_n526.xlsx")
map  <- read_excel("n760_names.xlsx")

# 2. Clean column names
colnames(chem) <- tolower(gsub(" ", "_", colnames(chem)))
colnames(map)  <- tolower(gsub(" ", "_", colnames(map)))

# 3. Match and filter
matched_ids <- intersect(chem$sample.id, map$sample_final)

# 4. Filter both datasets to matched samples
chem_matched <- chem %>% filter(sample.id %in% matched_ids)
map_matched  <- map  %>% filter(sample_final %in% matched_ids)

# 5. Merge all info by sample ID (join by sample_final ↔ sample.id)
merged_df <- chem_matched %>%
  left_join(map_matched, by = c("sample.id" = "sample_final"))

write_csv(merged_df, "matched_chemistry_with_names.csv")
#n276 samples for training


####################################################################
#1A FROM VCF TO GENOTYPE MATRIX - rows=individuals  columns =SNPs, in 0/1/2 format
####################################################################
#install.packages("vcfR")
#install.packages("adegenet")
library(vcfR)
library(adegenet)
#vcf <- read.vcfR("/Users/annamccormick/R/Feral_Cannabis_EGS/GS/Chemistry/SNPs.NC_only_noBlanks_filtered_GQ20_DP5.recode.vcf.gz")
vcf <- read.vcfR("/Users/annamccormick/R/Feral_Cannabis_EGS/Chemistry_GS_new/Cannabis_sativa_Aina_n760_filtered.vcf.gz")
genotype_matrix <- extract.gt(vcf, element = "GT")

dim(genotype_matrix)  # Should return (SNPs × Samples)
write.csv(genotype_matrix, "/Users/annamccormick/R/Feral_Cannabis_EGS/Chemistry_GS_new/Chem_genotype_matrix_raw.csv", row.names = TRUE)

genotype_matrix_raw <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/Chemistry_GS_new/Chem_genotype_matrix_raw.csv", 
                                row.names = 1, check.names = FALSE)

# Check dimensions
dim(genotype_matrix_raw)  # Should match (SNPs × Samples)

# Preview 
head(genotype_matrix_raw[, 1:5])  #first 5 samples

####################################################################
#1B - Convert genotype matrix to rrBLUP format
####################################################################
# Function to convert genotypes to numeric format for rrBLUP
convert_geno <- function(x) {
  ifelse(x == "0/0", 0, 
         ifelse(x == "0/1" | x == "1/0", 1, 
                ifelse(x == "1/1", 2, NA)))
}

# Apply conversion to the entire matrix
geno_numeric <- apply(genotype_matrix_raw, 2, convert_geno)

# Convert to dataframe
geno_df <- as.data.frame(geno_numeric)

# Check dimensions after conversion
dim(geno_df)  #(340734 × 423)

# Save the converted genotype matrix
write.csv(geno_df, "/Users/annamccormick/R/Feral_Cannabis_EGS/Chemistry_GS_new/genotype_matrix_rrBLUP_format.csv", row.names = TRUE)

####################################################################
#1C - Transpose genotype matrix for downstream use
####################################################################
library(rrBLUP)

# genotype matrix in rrBLUP format
geno_df <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/Chemistry_GS_new/genotype_matrix_rrBLUP_format.csv", 
                    row.names = 1, check.names = FALSE)

# Check dimensions to confirm SNPs x Samples
dim(geno_df)

# Preview first few values
head(geno_df[, 1:5])

geno_df <- t(geno_df)  # Transpose so genotypes (samples) become rows
geno_df <- as.data.frame(geno_df)  # Convert back to a dataframe
head(geno_df[, 1:5])
cat("Corrected genotype matrix dimensions:", dim(geno_df), "\n")


write.csv(geno_df, "/Users/annamccormick/R/Feral_Cannabis_EGS/Chemistry_GS_new/genotype_matrix_transposed.csv", row.names = TRUE)


##################
# 1D - GENOTYPE IMPORT & CLEANING/IMPUTATION
##################

# Import transposed genotype matrix
geno_df <- read.csv(
  "/Users/annamccormick/R/Feral_Cannabis_EGS/Chemistry_GS_new/genotype_matrix_transposed.csv", 
  row.names = 1, check.names = FALSE
)


# Define imputation function (replace missing values with the most common allele)
impute_geno <- function(x) {
  x[is.na(x)] <- as.numeric(names(which.max(table(x, useNA = "no"))))  # Replace NA with mode
  return(x)
}

# Apply the imputation function to each SNP column
cat("Imputing missing values in genotype data...\n")
geno_df_imputed <- apply(geno_df, 2, impute_geno)

# Convert back to dataframe
geno_df_imputed <- as.data.frame(geno_df_imputed)

# Confirm missing values after imputation
num_missing_after <- sum(is.na(geno_df_imputed))
cat("Missing values after imputation:", num_missing_after, "\n") # Should be 0

write.csv(
  geno_df_imputed, 
  "/Users/annamccormick/R/Feral_Cannabis_EGS/Chemistry_GS_new/genotype_matrix_imputed.csv", 
  row.names = TRUE
)

############################################################################################################
# Step 2
############################################################################################################
library(rrBLUP)
library(dplyr)

##################
# 1. GENOTYPE IMPORT (IMPUTED)
##################
geno_df_imputed <- read.csv(
  "/Users/annamccormick/R/Feral_Cannabis_EGS/Chemistry_GS_new/genotype_matrix_imputed.csv",
  row.names = 1, check.names = FALSE
)

##################
# 2. PHENOTYPE IMPORT
##################
# PHENOTYPE IMPORT
pheno <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/Chemistry_GS_new/matched_chemistry_with_names.csv")

# Set proper sample IDs as rownames
rownames(pheno) <- pheno$sample.id.vcf2

##################
# 3. ALIGN SAMPLES AND FILTER TRAITS
##################

# Step 1: Keep only overlapping samples
shared_samples <- intersect(rownames(pheno), rownames(geno_df_imputed))

head(rownames(pheno))
head(rownames(geno_df_imputed))

geno_train <- geno_df_imputed[shared_samples, ]
pheno_train <- pheno[shared_samples, ]

# Step 2: Keep only numeric trait columns (assuming they're the first 9)
pheno_numeric <- pheno_train[, 1:9]

# Step 3: Count non-missing values per trait
trait_valid_counts <- sapply(pheno_numeric, function(x) sum(!is.na(x)))
print(trait_valid_counts)

# Step 4: Keep only traits with at least 50 non-missing entries (adjust if needed)
traits_to_keep <- names(trait_valid_counts[trait_valid_counts >= 50])
traits_to_keep <- setdiff(traits_to_keep, "sample.id")
cat("Traits retained for rrBLUP:\n")
print(traits_to_keep)

# Step 5: Subset and drop samples with missing values in those traits
pheno_train_clean <- pheno_numeric[, traits_to_keep]
pheno_train_clean <- na.omit(pheno_train_clean)

# Step 6: Filter genotype matrix to match phenotype samples
geno_train_clean <- geno_train[rownames(pheno_train_clean), ]

# Step 7: Convert genotype to numeric matrix
geno_train_clean <- as.matrix(geno_train_clean)
mode(geno_train_clean) <- "numeric"

##################
# 4. DEFINE TEST SET (OPTIONAL)
##################
test_genotypes <- setdiff(rownames(geno_df_imputed), rownames(geno_train_clean))
geno_test <- geno_df_imputed[test_genotypes, ]

# Convert to numeric
geno_test <- as.matrix(geno_test)
mode(geno_test) <- "numeric"

# Check dimensions
cat("Training set dimensions (n samples x m SNPs):", dim(geno_train_clean), "\n")
cat("Phenotype training dimensions:", dim(pheno_train_clean), "\n")
cat("Test set dimensions:", dim(geno_test), "\n")

############################################################################################################
library(rrBLUP)

# Function to train rrBLUP and predict GEBVs for multiple traits
run_rrBLUP_all_traits <- function(geno_train, geno_test, pheno_train, output_prefix) {
  
  traits <- colnames(pheno_train)  # Get trait names
  GEBVs_train_list <- list()
  GEBVs_test_list <- list()
  
  for (trait in traits) {
    cat("Training rrBLUP model for", trait, "...\n")
    
    y <- pheno_train[[trait]]  # Extract trait values
    y[is.na(y)] <- mean(y, na.rm = TRUE)  # Impute missing values (just in case)
    
    # Skip if variance is zero (to prevent errors)
    trait_var <- var(y, na.rm = TRUE)
    if (is.na(trait_var) || trait_var == 0) {
      cat("Skipping", trait, "due to NA or zero variance.\n")
      next
    }
    
    # Train rrBLUP model
    rrblup_model <- mixed.solve(y = y, Z = geno_train)
    
    # Predict GEBVs for training and test sets
    GEBVs_train <- geno_train %*% rrblup_model$u
    GEBVs_test <- geno_test %*% rrblup_model$u
    
    # Convert to data frames and store
    GEBVs_train_list[[trait]] <- as.data.frame(GEBVs_train)
    GEBVs_test_list[[trait]] <- as.data.frame(GEBVs_test)
    
    # Rename columns
    colnames(GEBVs_train_list[[trait]]) <- paste0("GEBV_", trait, "_Train")
    colnames(GEBVs_test_list[[trait]]) <- paste0("GEBV_", trait, "_Test")
  }
  
  # Combine GEBVs for all traits
  GEBVs_train_final <- do.call(cbind, GEBVs_train_list)
  GEBVs_test_final <- do.call(cbind, GEBVs_test_list)
  
  # Add row names
  rownames(GEBVs_train_final) <- rownames(geno_train)
  rownames(GEBVs_test_final) <- rownames(geno_test)
  
  # Save results
  write.csv(GEBVs_train_final, paste0(output_prefix, "_Train.csv"), row.names = TRUE)
  write.csv(GEBVs_test_final, paste0(output_prefix, "_Test.csv"), row.names = TRUE)
  
  cat("GEBV predictions saved for all traits.\n")
}

run_rrBLUP_all_traits(
  geno_train = geno_train_clean,
  geno_test  = geno_test,
  pheno_train = pheno_train_clean,
  output_prefix = "/Users/annamccormick/R/Feral_Cannabis_EGS/Chemistry_GS_new/GEBVs_cannabinoids"
)


########
# Load both CSVs
train_df <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/Chemistry_GS_new/GEBVs_cannabinoids_Train.csv", row.names = 1)
test_df  <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/Chemistry_GS_new/GEBVs_cannabinoids_Test.csv",  row.names = 1)

# Standardize column names by removing "_Train" and "_Test"
colnames(train_df) <- gsub("_Train$", "", colnames(train_df))
colnames(test_df)  <- gsub("_Test$",  "", colnames(test_df))

combined_df <- rbind(train_df, test_df)

# Save
write.csv(combined_df, "/Users/annamccormick/R/Feral_Cannabis_EGS/Chemistry_GS_new/GEBVs_cannabinoids_ALL.csv", row.names = TRUE)


########################################################################################################################################
#adding in metadata to plot
########################################################################################################################################
#FIGURE S13B - PLOT
########################################################################################################################################
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# Load your merged GEBV + metadata file
df <- read_csv("/Users/annamccormick/R/Feral_Cannabis_EGS/Chemistry_GS_new/GEBVs_cannabinoids_ALL.csv")

df_long <- df %>%
  pivot_longer(
    cols = starts_with("GEBV_"),  # or explicitly list trait columns
    names_to = "Trait",
    values_to = "GEBV"
  )

ggplot(df_long, aes(x = Population, y = GEBV, fill = Population)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  facet_wrap(~ Trait, scales = "free_y") +
  theme_minimal() +
  labs(title = "",
       y = "GEBV Value",
       x = "Population") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


#plotting subset 
unique(df_long$Trait)

traits_keep <- c(
  "GEBV_total_cannabinoids",
  "GEBV_cbd",
  "GEBV_d9thc",
  "GEBV_cbg"
)

df_long_filt <- df_long %>%
  filter(Trait %in% traits_keep)

df_long_filt <- df_long_filt %>%
  mutate(
    Trait = recode(
      Trait,
      "GEBV_total_cannabinoids" = "Total Cannabinoids",
      "GEBV_cbd"                = "CBD",
      "GEBV_d9thc"              = "Δ9-THC",
      "GEBV_cbg"                = "CBG"
    )
  )

df_long_filt$Trait <- factor(
  df_long_filt$Trait,
  levels = c(
    "Total Cannabinoids",
    "CBG",
    "CBD",
    "Δ9-THC"
  )
)

aina_cols <- c(
  "Cluster_1" = "darkred",
  "Cluster_2" = "darkblue",
  "Cluster_3" = "darkgreen",
  "Cluster_4" = "purple4",
  "Cluster_5" = "darkorange4"
)

p <- ggplot(df_long_filt, aes(x = Population, y = GEBV, fill = Population)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  facet_wrap(~ Trait, scales = "free_y", ncol = 2) +
  scale_fill_manual(
    values = aina_cols,
    drop   = FALSE
  ) +
  theme_minimal() +
  labs(
    title = "",
    y = "GEBV Value",
    x = "Population"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

p

ggsave(filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/gebv_cannabinoids.pdf", 
       plot =p, 
       width = 10,    # Width of the plot (adjust as needed)
       height = 8)   # Reduced height (adjust as needed)

#########################################################################################################
#Step 3: Model Prediction Accuracies via rrBLUP
#########################################################################################################

library(rrBLUP)
library(dplyr)
library(hibayes)
library(ggplot2)

####################
# STEP ONE: Load data
####################

# Load Genotype Data (in rrBLUP format)
geno_data <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/Chemistry_GS_new/genotype_matrix_imputed.csv", row.names = 1, check.names = FALSE)

# Load Phenotypic Data
pheno_data <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/Chemistry_GS_new/matched_chemistry_with_names.csv", row.names = 1)
# Set proper sample IDs as rownames
rownames(pheno_data) <- pheno_data$sample.id.vcf2

# Only keep the traits you want (excluding "sample.id")
pheno_data <- pheno_data[, c("total_cannabinoids", "cbdv", "thcv", "cbd", "cbc", "d8thc", "d9thc", "cbg")]

# Drop lines with all NA
pheno_data <- pheno_data[rowSums(!is.na(pheno_data)) > 0, ]

# Keep only lines present in both genotype and phenotype
common_ids <- intersect(rownames(pheno_data), rownames(geno_data))
geno_data <- geno_data[common_ids, ]
pheno_data <- pheno_data[common_ids, ]

# Convert genotype to numeric matrix
geno_matrix <- as.matrix(geno_data)
mode(geno_matrix) <- "numeric"

####################
# STEP TWO: Prep for CV
####################

g.in <- geno_matrix
y.in <- pheno_data
y.trainset <- y.in  # all rows are used, NA values will be handled by k.xval()

####################
# STEP THREE: Run k-fold CV
####################
source("/Users/annamccormick/R/Feral_Cannabis_EGS/EGS/xval_kfold_functions.R")

# Run k-fold cross-validation for rrBLUP
xval_k10_rrblup <- k.xval(
  g.in = g.in,
  y.in = y.in,
  y.trainset = y.trainset,
  k.fold = 10,
  reps = 50
)

# Save result if needed
saveRDS(xval_k10_rrblup, "/Users/annamccormick/R/Feral_Cannabis_EGS/Chemistry_GS_new/xval_rrblup_kfold_10.RDS")

########################################################################################################################################
#FIGURE S13A - PLOT of PAs
########################################################################################################################################
library(ggplot2)
setwd("/Users/annamccormick/R/Feral_Cannabis_EGS/Chemistry_GS_new/")

# RR-BLUP - K-Fold Cross-Validation
#Load cross-validation results for rrBLUP
rrblup_kfold10 <- readRDS("xval_rrblup_kfold_10.RDS")
rrblup_kfold10 <- rrblup_kfold10$xval.result
rrblup_kfold10$r.mean <- as.numeric(rrblup_kfold10$r.mean)


# Organize all dataframes for merging
## Rename model names
rrblup_kfold10$model <- "rrBLUP"

## Input xval type
rrblup_kfold10$xval <- "Ten-Fold"

#Combine all model results into a single list
model_list <- list(rrblup_kfold10)
#model_list <- list(rrblup_kfold10,gauss_kfold_10,EXP_kfold_10)

#Remove any NA values from the model results
model_list1 <- lapply(model_list, na.omit)

#Combine all models into a single dataframe
all_models <- do.call("rbind", model_list1)

#Convert standard deviation values to numeric
all_models$r.sd <- as.numeric(all_models$r.sd)

#Ensure cross-validation type is a factor with specified levels
all_models$xval <- factor(all_models$xval, levels = c("Ten-Fold"))

#Filter for specific traits of interest

#Plot 
library(ggplot2)

# Make trait a factor to control facet order (optional: order by r.mean or custom)
all_models$trait <- factor(all_models$trait,
                           levels = c("total_cannabinoids", "cbd", "cbg", "cbc", 
                                      "cbdv", "thcv", "d9thc", "d8thc"))

# Plot prediction accuracies
ggplot(all_models, aes(y = r.mean, x = trait, color = model)) +
  theme_bw() +
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd), 
                width = 0.3, position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
  labs(
    y = "Prediction Accuracy (r)",
    x = "Cannabinoid Trait",
    title = ""
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  ylim(0, 1.0)  # <-- Updated to cap y-axis at 1.0


traits_keep <- c(
  "total_cannabinoids",
  "cbg",
  "cbd",
  "d9thc"
)

all_models_filt <- all_models %>%
  dplyr::filter(trait %in% traits_keep)

trait_labels <- c(
  total_cannabinoids = "Total Cannabinoids",
  cbg                = "CBG",
  cbd                = "CBD",
  d9thc              = "Δ9-THC"
)

all_models_filt$trait <- factor(
  all_models_filt$trait,
  levels = c(
    "total_cannabinoids",
    "cbg",
    "cbd",
    "d9thc"
  )
)


p<- ggplot(all_models_filt, aes(y = r.mean, x = trait, color = model)) +
  theme_bw() +
  geom_errorbar(
    aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd),
    width = 0.3,
    position = position_dodge(width = 0.5)
  ) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
  labs(
    y = "Prediction Accuracy (r)",
    x = "Cannabinoid Trait",
    title = ""
  ) +
  scale_x_discrete(labels = trait_labels) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 12),
    axis.title  = element_text(size = 14),
    legend.position = "none"
  ) +
  ylim(0, 1.0)

ggsave(filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/PA_gebv_cannabinoids.pdf", 
       plot =p, 
       width = 8,    # Width of the plot (adjust as needed)
       height = 6)   # Reduced height (adjust as needed)





########################################################################################################################################################################
# Figure S14
########################################################################################################################################################################

############################################################################################################
# HIGHLIGHT: specific SNPs and shade the chr7 CBDAS-cluster region (29–61 Mb)
############################################################################################################
library(readr)
library(dplyr)
library(ggplot2)
library(tibble)

# Load marker effects
eff <- read_csv("/Users/annamccormick/R/Feral_Cannabis_EGS/Chemistry_GS_new/rrBLUP_marker_effects_CBD_TOP1pct.csv")

# Parse chr/pos from SNP name like "NC_044371.1_1234567"
eff2 <- eff %>%
  mutate(
    chr = sub("_[0-9]+$", "", SNP),
    pos = as.numeric(sub("^.*_", "", SNP)),
    pos_Mb = pos / 1e6,
    sign = ifelse(effect >= 0, "positive", "negative")
  ) %>%
  filter(!is.na(pos))

# Map NC accessions -> chromosome number (edit if needed)
chr_map <- tribble(
  ~chr,          ~chr_num,
  "NC_044371.1",  1,
  "NC_044375.1",  2,
  "NC_044372.1",  3,
  "NC_044373.1",  4,
  "NC_044374.1",  5,
  "NC_044377.1",  6,
  "NC_044378.1",  7,
  "NC_044379.1",  8,
  "NC_044376.1",  9,
  "NC_044370.1", 10
)

# Join map + make ordered facet titles: "1 - NC_..."
eff3 <- eff2 %>%
  left_join(chr_map, by = "chr") %>%
  filter(!is.na(chr_num)) %>%
  mutate(
    chr_title = factor(
      paste0(chr_num, " - ", chr),
      levels = paste0(1:10, " - ", chr_map$chr),
      ordered = TRUE
    )
  )

# Plot: |effect| across position, faceted in chromosome-number order with NC accession in title
ggplot(eff3, aes(x = pos_Mb, y = abs_effect)) +
  geom_point(aes(shape = sign), alpha = 0.8, size = 1.4) +
  facet_wrap(~ chr_title, scales = "free_x", ncol = 2) +
  labs(x = "Position (Mb)", y = "|CBD marker effect|", shape = "Effect sign") +
  theme_classic(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))

ggplot(eff3, aes(x = pos_Mb, y = abs_effect)) +
  geom_point(aes(color = sign), alpha = 0.8, size = 1.6) +
  scale_color_manual(
    values = c(
      "negative" = "red3",
      "positive" = "darkgreen"
    )
  ) +
  facet_wrap(~ chr_title, scales = "free_x", ncol = 2) +
  labs(
    x = "Position (Mb)",
    y = "|CBD marker effect|",
    color = "Effect sign"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )


ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/Chemistry_GS_new/CBD_marker_effects_by_chromosome.pdf",
  device = "pdf",
  width = 15,
  height = 10,
  units = "in"
)


#counts
eff_all <- read_csv(
  "/Users/annamccormick/R/Feral_Cannabis_EGS/Chemistry_GS_new/rrBLUP_marker_effects_CBD_ALL.csv"
)

thr <- quantile(eff_all$abs_effect, 0.99, na.rm = TRUE)

eff_top1 <- eff_all %>% filter(abs_effect >= thr)

n_top1 <- nrow(eff_top1)
n_total <- nrow(eff_all)

cat("Top 1% SNPs:", n_top1, "out of", n_total, "\n")


eff_top1 %>%
  mutate(sign = ifelse(effect > 0, "positive", "negative")) %>%
  count(sign)

eff_top1 %>%
  mutate(
    chr = sub("_[0-9]+$", "", SNP),
    sign = ifelse(effect > 0, "positive", "negative")
  ) %>%
  count(chr, sign) %>%
  arrange(chr)
########################
# highlighting the SNPs in the CBDAS window
########################
library(readr)
library(dplyr)
library(ggplot2)
library(tibble)

# Load marker effects
eff <- read_csv("/Users/annamccormick/R/Feral_Cannabis_EGS/Chemistry_GS_new/rrBLUP_marker_effects_CBD_TOP1pct.csv")

# Parse chr/pos from SNP name like "NC_044371.1_1234567"
eff2 <- eff %>%
  mutate(
    chr = sub("_[0-9]+$", "", SNP),
    pos = as.numeric(sub("^.*_", "", SNP)),
    pos_Mb = pos / 1e6,
    sign = ifelse(effect >= 0, "positive", "negative")
  ) %>%
  filter(!is.na(pos))

# Map NC accessions -> chromosome number
chr_map <- tribble(
  ~chr,          ~chr_num,
  "NC_044371.1",  1,
  "NC_044375.1",  2,
  "NC_044372.1",  3,
  "NC_044373.1",  4,
  "NC_044374.1",  5,
  "NC_044377.1",  6,
  "NC_044378.1",  7,
  "NC_044379.1",  8,
  "NC_044376.1",  9,
  "NC_044370.1", 10
)

eff3 <- eff2 %>%
  left_join(chr_map, by = "chr") %>%
  filter(!is.na(chr_num)) %>%
  mutate(
    chr_title = factor(
      paste0(chr_num, " - ", chr),
      levels = paste0(1:10, " - ", chr_map$chr),
      ordered = TRUE
    )
  )


# SNPs to highlight
highlight_snps <- c(
  "NC_044378.1_29008289",
  "NC_044378.1_30301091",
  "NC_044378.1_34851124",
  "NC_044378.1_34851275",
  "NC_044378.1_34992941",
  "NC_044378.1_39816534",
  "NC_044378.1_44182846",
  "NC_044378.1_49206235",
  "NC_044378.1_55588169"
)

eff3 <- eff3 %>%
  mutate(is_highlight = ifelse(SNP %in% highlight_snps, "Highlighted SNPs", "Other SNPs"))

# Data frame for shading region on chromosome 7 facet only
shade_df <- tibble(
  chr_title = factor("7 - NC_044378.1", levels = levels(eff3$chr_title), ordered = TRUE),
  xmin = 29,
  xmax = 61,
  ymin = -Inf,
  ymax = Inf
)

# Plot
p <- ggplot(eff3, aes(x = pos_Mb, y = abs_effect)) +
  # Shade CBDAS cluster region only on chr7 facet
  geom_rect(
    data = shade_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    alpha = 0.12,
    fill = "gold"
  ) +
  # Base points colored by sign
  geom_point(aes(color = sign), alpha = 0.75, size = 1.6) +
  scale_color_manual(values = c("negative" = "red3", "positive" = "darkgreen")) +
  # Overlay highlighted SNPs in a distinct color (and slightly larger)
  geom_point(
    data = eff3 %>% filter(SNP %in% highlight_snps),
    aes(x = pos_Mb, y = abs_effect),
    inherit.aes = FALSE,
    color = "dodgerblue3",
    size = 2.6,
    alpha = 0.95
  ) +
  facet_wrap(~ chr_title, scales = "free_x", ncol = 2) +
  labs(
    x = "Position (Mb)",
    y = "|CBD marker effect|",
    color = "Effect sign"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

print(p)

ggsave(
  plot = p,
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/CBD_marker_effects_by_chromosome_highlight_CBDAScluster.pdf",
  device = "pdf",
  width = 10,
  height = 8,
  units = "in"
)

#table S11
########################################## take 2
library(readr)
library(dplyr)
library(stringr)
library(rtracklayer)
library(GenomicRanges)
library(tibble)

library(readr)
library(dplyr)
library(stringr)

annot <- read_csv("/Users/annamccormick/R/Feral_Cannabis_EGS/Chemistry_GS_new/rrBLUP_marker_effects_CBD_TOP1pct.csv") %>%
  mutate(
    chr = sub("_[0-9]+$", "", SNP),
    pos = as.integer(sub("^.*_", "", SNP))
  )

snps_gr <- GRanges(
  seqnames = annot$chr,
  ranges   = IRanges(start = annot$pos, end = annot$pos),
  SNP      = annot$SNP
)

# 2) Import GFF
gff <- import("~/R/Feral_Cannabis_EGS/GEA/SNP_checks/genomic_main_chromosomes.gff")

# 3) Keep transcript-like features that actually carry 'product'
tx <- gff[gff$type %in% c("mRNA", "lnc_RNA", "transcript")]

# Helper to decode NCBI-escaped strings like %2C and %25
decode_ncbi <- function(x) {
  x %>%
    str_replace_all("%2C", ",") %>%
    str_replace_all("%3B", ";") %>%
    str_replace_all("%25", "%")
}

tx_df <- as.data.frame(mcols(tx)) %>%
  transmute(
    tx_id   = ID,
    gene_id = Parent,  # typically gene-LOC...
    gene    = if ("gene" %in% names(.)) gene else NA_character_,   # LOC...
    Name    = if ("Name" %in% names(.)) Name else NA_character_,   # XM_/XR_
    product = if ("product" %in% names(.)) decode_ncbi(product) else NA_character_,
    Dbxref  = if ("Dbxref" %in% names(.)) Dbxref else NA_character_,
    transcript_id = if ("transcript_id" %in% names(.)) transcript_id else NA_character_
  )

# 4) Overlap SNPs with transcripts; if none, use nearest transcript
hits <- findOverlaps(snps_gr, tx, ignore.strand = TRUE)

overlap_tx <- tibble(
  SNP = mcols(snps_gr)$SNP[queryHits(hits)],
  tx_id = mcols(tx)$ID[subjectHits(hits)]
) %>%
  left_join(tx_df, by = "tx_id") %>%
  mutate(
    hit_type = "overlaps_transcript",
    distance_bp = 0L
  )

# nearest transcript for every SNP
nearest_tx_idx <- nearest(snps_gr, tx, ignore.strand = TRUE)
nearest_tx_dist <- distanceToNearest(snps_gr, tx, ignore.strand = TRUE)@elementMetadata$distance

nearest_tx <- tibble(
  SNP = mcols(snps_gr)$SNP,
  tx_id = mcols(tx)$ID[nearest_tx_idx],
  distance_bp = nearest_tx_dist
) %>%
  left_join(tx_df, by = "tx_id") %>%
  mutate(hit_type = "nearest_transcript")

# 5) Prefer overlap annotations when present
tx_annot <- nearest_tx %>%
  left_join(overlap_tx %>% select(SNP, tx_id, gene_id, gene, Name, product, Dbxref, transcript_id, hit_type, distance_bp),
            by = "SNP",
            suffix = c("_nearest","")) %>%
  mutate(
    tx_id = coalesce(tx_id, tx_id_nearest),
    gene_id = coalesce(gene_id, gene_id_nearest),
    gene = coalesce(gene, gene_nearest),
    Name = coalesce(Name, Name_nearest),
    product = coalesce(product, product_nearest),
    Dbxref = coalesce(Dbxref, Dbxref_nearest),
    transcript_id = coalesce(transcript_id, transcript_id_nearest),
    hit_type = coalesce(hit_type, hit_type_nearest),
    distance_bp = ifelse(hit_type == "overlaps_transcript", 0L, distance_bp_nearest)
  ) %>%
  select(SNP, hit_type, distance_bp, gene_id, gene, Name, transcript_id, product, Dbxref)

# 6) Merge back onto your SNP effect table (keeps abs_effect etc.)
annot2 <- annot %>%
  left_join(tx_annot, by = "SNP")

write_csv(
  annot2,
  "/Users/annamccormick/R/Feral_Cannabis_EGS/Chemistry_GS_new/CBD_top1pct_SNP_annotations_with_transcript_products.csv"
)




########################################################################################################################################################################
# Figure S15
########################################################################################################################################################################

##############
#  PLOT
##############
library(readr)
library(dplyr)
library(purrr)
library(ggplot2)
library(grid)   # for unit()

# 1) Load & stack all WZA outputs
wza_dir <- "/Users/annamccormick/R/Feral_Cannabis_EGS/WZA_all/WZA_out"  # adjust if needed
files   <- list.files(wza_dir, pattern = "_wza\\.csv$", full.names = TRUE)

df <- map_dfr(files, ~{
  tr <- sub("_wza\\.csv$", "", basename(.x))     # e.g., BIO01
  read_csv(.x, show_col_types = FALSE) |> mutate(trait = tr)
})

# 2) Use Z_pVal only → make -log10
stopifnot("Z_pVal" %in% names(df))
df <- df |>
  mutate(
    Z_pVal = pmin(pmax(as.numeric(Z_pVal), .Machine$double.xmin), 1),
    ml10p  = -log10(Z_pVal),
    chr    = as.character(chr),
    pos    = as.numeric(pos)
  ) |>
  filter(is.finite(ml10p))

# 3) Order traits BIO01..BIO19 and chromosomes 1..10 (natural order)
trait_levels <- sprintf("BIO%02d", 1:19)
df$trait <- factor(df$trait,
                   levels = intersect(trait_levels, unique(df$trait)),
                   ordered = TRUE)

chr_levels <- df |>
  distinct(chr) |>
  mutate(n = suppressWarnings(as.numeric(gsub("[^0-9]", "", chr)))) |>
  arrange(n, chr) |>
  pull(chr)
df$chr <- factor(df$chr, levels = chr_levels, ordered = TRUE)


p<- ggplot(df, aes(pos/1e6, ml10p)) +
  geom_point(size=0.35, alpha=0.75) +
  geom_hline(yintercept=-log10(0.01), linetype="dashed", color="red", linewidth=0.4) +
  facet_grid(trait ~ chr, scales = "free_y") +   # <-- free y per panel
  labs(x="Position (Mbp)", y=expression(-log[10](Z[pVal]))) +
  theme_classic(base_size = 9)

p
ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/FigS7.pdf",
  plot = p,
  width = 15,
  height = 10,
  units = "in"
)

############################################
#counting number of snps past threshold

snps_bio1_11 <- df %>%
  filter(
    trait %in% sprintf("BIO%02d", 1:11),
    ml10p >= thr
  ) %>%
  distinct(trait, chr, pos) %>%   # ensures unique SNPs
  nrow()

snps_bio1_11


snps_bio1_11_unique_loci <- df %>%
  filter(trait %in% sprintf("BIO%02d", 1:11), ml10p >= thr) %>%
  distinct(chr, pos) %>%
  nrow()

snps_bio1_11_unique_loci



snps_bio12_19 <- df %>%
  filter(
    trait %in% sprintf("BIO%02d", 12:19),
    ml10p >= thr
  ) %>%
  distinct(trait, chr, pos) %>%
  nrow()

snps_bio12_19


snps_bio12_19_unique_loci <- df %>%
  filter(trait %in% sprintf("BIO%02d", 12:19),
         ml10p >= thr) %>%
  distinct(chr, pos) %>%
  nrow()

snps_bio12_19_unique_loci


##########################################################################################
#red dots for bios appearing 3 times or more and above line for TEMPERATURE
##########################################################################################
thr <- -log10(0.01)  # p ≤ 0.01

# Keep only temperature traits
df_temp <- df %>%
  filter(trait %in% sprintf("BIO%02d", 1:11))

# Count in how many temp traits each locus passes threshold
temp_replication <- df_temp %>%
  filter(ml10p >= thr) %>%
  distinct(chr, pos, trait) %>%
  group_by(chr, pos) %>%
  summarise(
    n_temp_traits = n(),
    traits = paste(trait, collapse = ", "),
    .groups = "drop"
  )

# Loci that replicate in ≥3 temperature variables
rep3 <- temp_replication %>%
  filter(n_temp_traits >= 3) %>%
  select(chr, pos, n_temp_traits)

df_temp2 <- df_temp %>%
  left_join(rep3, by = c("chr", "pos")) %>%
  mutate(
    is_rep3 = !is.na(n_temp_traits),
    highlight = is_rep3 & (ml10p >= thr)
  )

p_temp_rep3 <- ggplot(df_temp2, aes(pos/1e6, ml10p)) +
  geom_point(size = 0.35, alpha = 0.6) +
  geom_hline(yintercept = thr, linetype = "dashed", color = "red", linewidth = 0.4) +
  geom_point(
    data = df_temp2 %>% filter(highlight),
    size = 0.7, alpha = 0.9, color = "red"
  ) +
  facet_grid(trait ~ chr, scales = "free_y") +
  labs(
    x = "Position (Mbp)",
    y = expression(-log[10](Z[pVal])),
    title = "BIO01–BIO11: loci replicated across ≥3 BIO variables (highlighted when significant)"
  ) +
  theme_classic(base_size = 9)

p_temp_rep3

ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/BIO01_BIO11_WZA_manhattan_red.pdf",
  plot     = p_temp_rep3,
  width    = 10,
  height   = 8,
  units    = "in"
)


########################################################################################################################################################################
# Figure S16
########################################################################################################################################################################

########################################################################
# Precipitation WZA: BIO12–BIO19
# Highlight loci that are significant (ml10p >= thr) AND replicated across >= 3 precip traits
##############

library(readr)
library(dplyr)
library(purrr)
library(ggplot2)

# ---- settings ----
wza_dir <- "/Users/annamccormick/R/Feral_Cannabis_EGS/WZA_all/WZA_out"
files   <- list.files(wza_dir, pattern = "_wza\\.csv$", full.names = TRUE)

thr <- -log10(0.01)  # WZA significance threshold: p <= 0.01  (ml10p >= 2)

# ---- 1) load & stack all traits ----
df <- map_dfr(files, ~{
  tr <- sub("_wza\\.csv$", "", basename(.x))  # e.g., "BIO01"
  read_csv(.x, show_col_types = FALSE) |> mutate(trait = tr)
})

stopifnot("Z_pVal" %in% names(df))

df <- df |>
  mutate(
    Z_pVal = pmin(pmax(as.numeric(Z_pVal), .Machine$double.xmin), 1),
    ml10p  = -log10(Z_pVal),
    chr    = as.character(chr),
    pos    = as.numeric(pos)
  ) |>
  filter(is.finite(ml10p), !is.na(chr), !is.na(pos))

# ---- 2) order traits/chromosomes (optional but nice) ----
trait_levels <- sprintf("BIO%02d", 1:19)
df$trait <- factor(df$trait,
                   levels = intersect(trait_levels, unique(df$trait)),
                   ordered = TRUE)

chr_levels <- df |>
  distinct(chr) |>
  mutate(n = suppressWarnings(as.numeric(gsub("[^0-9]", "", chr)))) |>
  arrange(n, chr) |>
  pull(chr)
df$chr <- factor(df$chr, levels = chr_levels, ordered = TRUE)

# ---- 3) subset precipitation traits ----
df_precip <- df %>%
  filter(trait %in% sprintf("BIO%02d", 12:19))

# ---- 4) replication count: in how many precip traits is each locus significant? ----
precip_replication <- df_precip %>%
  filter(ml10p >= thr) %>%
  distinct(chr, pos, trait) %>%       # one hit per trait
  group_by(chr, pos) %>%
  summarise(
    n_precip_traits = n(),
    traits = paste(trait, collapse = ", "),
    .groups = "drop"
  )

# loci replicated across >= 3 precipitation variables
rep3_precip <- precip_replication %>%
  filter(n_precip_traits >= 3) %>%
  select(chr, pos, n_precip_traits)

# ---- 5) flag the points to highlight (ONLY if above threshold) ----
df_precip2 <- df_precip %>%
  left_join(rep3_precip, by = c("chr", "pos")) %>%
  mutate(
    is_rep3 = !is.na(n_precip_traits),
    highlight = is_rep3 & (ml10p >= thr)   # ensures nothing below red line is highlighted
  )

# ---- 6) plot: base points black, replicated significant loci blue ----
p_precip_rep3 <- ggplot(df_precip2, aes(pos/1e6, ml10p)) +
  geom_point(size = 0.35, alpha = 0.60) +
  geom_hline(yintercept = thr, linetype = "dashed", color = "red", linewidth = 0.4) +
  geom_point(
    data = df_precip2 %>% filter(highlight),
    size = 0.75, alpha = 0.95, color = "blue"
  ) +
  facet_grid(trait ~ chr, scales = "free_y") +
  labs(
    x = "Position (Mbp)",
    y = expression(-log[10](Z[pVal])),
    title = "BIO12–BIO19: precipitation-associated loci replicated across ≥3 BIO variables (blue = significant & replicated)"
  ) +
  theme_classic(base_size = 9)

p_precip_rep3


ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/BIO12_BIO19_WZA_manhattan_rep3_blue.pdf",
  plot     = p_precip_rep3,
  width    = 10,
  height   = 8,
  units    = "in"
)


library(readr)

write_csv(
  precip_replication %>% filter(n_precip_traits >= 3),
  "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/BIO12_BIO19_rep3_precipitation_GEA_SNPs.csv"
)



########################################################################################################################################################################
# Figure S17
########################################################################################################################################################################

#################
# Figure S17A
#################

# choose the BIO5 column 
library(dplyr)
library(ggplot2)
library(forcats)

df <- readr::read_csv(
  "/Users/annamccormick/R/Feral_Cannabis_EGS/Aina_new_EGS/rrblup_with_metadata_copy.csv",
  guess_max = 10000
)

# pick the BIO5 column 
bio5_col <- intersect(c("bio_05.x","bio05.x","bio_5.x","bio_05","bio_05.y"), names(df))[1]
stopifnot(length(bio5_col) == 1)

plotdat <- df %>%
  dplyr::select(sample.id, pop_base, dplyr::all_of(bio5_col)) %>%
  dplyr::rename(BIO5 = !!rlang::sym(bio5_col)) %>%
  dplyr::filter(!is.na(pop_base), !is.na(BIO5)) %>%
  dplyr::mutate(pop_base = fct_reorder(pop_base, BIO5, .fun = median, na.rm = TRUE))



p<- ggplot(plotdat, aes(x = pop_base, y = BIO5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red",
             linewidth = 0.5, inherit.aes = FALSE) +
  geom_boxplot(outlier.size = 0.7, width = 0.6, fill = "grey85", color = "black") +
  labs(x = "Population", y = "BIO5 (GEBV)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())

p


ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/FigS15.pdf",
  plot     = p,
  width    = 16,
  height   = 4,
  units    = "in"
)

#install.packages("ggtext")
library(ggtext)

highlight_pops <- c("NE-22-ID-06", "KS-22-PB-04", "NY-23-LG-01")

# base plot without x tick labels
p <- ggplot(plotdat, aes(x = pop_base, y = BIO5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red",
             linewidth = 0.5, inherit.aes = FALSE) +
  geom_boxplot(outlier.size = 0.7, width = 0.6, fill = "grey85", color = "black") +
  labs(x = "Population", y = "BIO5 (GEBV)") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# create custom tick-label data
tick_df <- data.frame(pop_base = levels(plotdat$pop_base)) %>%
  mutate(lbl = ifelse(pop_base %in% highlight_pops,
                      paste0("<span style='color:#d62728; font-weight:700;'>", pop_base, "</span>"),
                      pop_base))
p
# add custom labels (placed slightly below plot area)
p2<- p +
  geom_richtext(
    data = tick_df,
    aes(x = pop_base, y = -Inf, label = lbl),
    inherit.aes = FALSE,
    vjust = -0.5,
    angle = 45,
    size = 3,
    fill = NA,
    label.color = NA
  ) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5.5, 5.5, 30, 5.5))  # add bottom margin for labels

p2
##### highlight specific samples

highlight_pops <- c("NE-22-ID-06", "KS-22-PB-04", "NY-23-LG-01")

# build a colour vector in the *final plotted order*
x_levels <- levels(plotdat$pop_base)
x_cols <- ifelse(x_levels %in% highlight_pops, "#d62728", "black")

p_simple <- ggplot(plotdat, aes(x = pop_base, y = BIO5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
  geom_boxplot(outlier.size = 0.7, width = 0.6, fill = "grey85", color = "black") +
  labs(x = "Population", y = "BIO5 (GEBV)") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = x_cols),
    panel.grid.minor = element_blank()
  )

p_simple
ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/FigS15.pdf",
  plot     = p_simple,
  width    = 16,
  height   = 4,
  units    = "in"
)

###################
# compute the same stat used to order the plot (median BIO5)
pop_order <- plotdat %>%
  group_by(pop_base) %>%
  summarise(
    n = n(),
    median_BIO5 = median(BIO5, na.rm = TRUE),
    mean_BIO5   = mean(BIO5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(median_BIO5)

# bottom/top 10 exactly as in the faceted boxplot order
bottom10_from_plot <- head(pop_order$pop_base, 10)
top10_from_plot    <- tail(pop_order$pop_base, 10)

cat("Bottom 10 (by median, matches plot order):\n",
    paste(bottom10_from_plot, collapse = ", "), "\n\n")
cat("Top 10 (by median, matches plot order):\n",
    paste(top10_from_plot, collapse = ", "), "\n")



####
library(dplyr)

# make sure BIO5 column is there
stopifnot(all(c("sample.id","pop_base","BIO5") %in% names(plotdat)))

library(dplyr)

bottom10_samples <- plotdat %>%
  filter(!is.na(BIO5)) %>%
  dplyr::slice_min(BIO5, n = 10, with_ties = FALSE) %>%
  mutate(rank = row_number())

top10_samples <- plotdat %>%
  filter(!is.na(BIO5)) %>%
  dplyr::slice_max(BIO5, n = 10, with_ties = FALSE) %>%
  mutate(rank = row_number())

cat("Bottom 10 samples (BIO5 GEBV):\n")
print(bottom10_samples %>% dplyr::select(rank, sample.id, pop_base, BIO5))

cat("\nTop 10 samples (BIO5 GEBV):\n")
print(top10_samples %>% dplyr::select(rank, sample.id, pop_base, BIO5))

print(as.data.frame(bottom10_samples)[, c("rank","sample.id","pop_base","BIO5")], row.names = FALSE)
print(as.data.frame(top10_samples   )[, c("rank","sample.id","pop_base","BIO5")], row.names = FALSE)



#################
# Figure S17B
#################
#Field Results
library(readr)
library(dplyr)
library(ggplot2)
library(forcats)
library(gridExtra)

# 1. Load data
df <- read_csv("/Users/annamccormick/R/Feral_Cannabis_EGS/data/inputs/Feral_Hemp_Evaluation_24_nomissing_height.csv")

# 2. Clean + summarize
df_clean <- df %>%
  mutate(
    `Female Plant_Height (cm)` = as.numeric(`Female Plant_Height (cm)`),
    `Female Plt Stem_Diameter (mm)` = as.numeric(`Female Plt Stem_Diameter (mm)`)
  ) %>%
  filter(!is.na(Accession))

df_summary <- df_clean %>%
  group_by(Accession) %>%
  summarise(
    mean_height = mean(`Female Plant_Height (cm)`, na.rm = TRUE),
    se_height   = sd(`Female Plant_Height (cm)`, na.rm = TRUE) /
      sqrt(sum(!is.na(`Female Plant_Height (cm)`))),
    mean_diam   = mean(`Female Plt Stem_Diameter (mm)`, na.rm = TRUE),
    se_diam     = sd(`Female Plt Stem_Diameter (mm)`, na.rm = TRUE) /
      sqrt(sum(!is.na(`Female Plt Stem_Diameter (mm)`)))
  ) %>%
  ungroup()

# 3. Highlight specific accessions
highlight_height <- c("KS-22-PB-04", "NE-22-ID-06", "NY-23-LG-01")
highlight_diam   <- c("KS-22-PB-05")

df_summary <- df_summary %>%
  mutate(
    highlight_height = ifelse(Accession %in% highlight_height, "highlight", "normal"),
    highlight_diam   = ifelse(Accession %in% highlight_diam, "highlight", "normal")
  )

# 4. Plot height (with red bars + red text for highlighted)
p_height <- df_summary %>%
  mutate(Accession = fct_reorder(Accession, mean_height)) %>%
  ggplot(aes(x = Accession, y = mean_height, fill = highlight_height)) +
  geom_col() +
  scale_fill_manual(values = c("highlight" = "red", "normal" = "#6C91BF")) +
  geom_errorbar(aes(ymin = mean_height - se_height,
                    ymax = mean_height + se_height),
                width = 0.3, color = "black") +
  labs(
    title = "",
    x = "",
    y = "Height (cm)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1,
                               colour = ifelse(
                                 levels(fct_reorder(df_summary$Accession, df_summary$mean_height))
                                 %in% highlight_height, "red", "black"
                               ))
  )

p_height

ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/Fig4.pdf",
  plot = p_height,
  width = 16,
  height = 4,
  units = "in"
)



##################
#Figure S17C
##################
library(tidyverse)

# 1) Read data
df <- read.csv("/Users/annamccormick/R/Feral_Cannabis_EGS/chamber_expt/EGS Growth Chamber.csv", stringsAsFactors = FALSE)
colnames(df)
library(tidyverse)

# Long format for heights only
long <- df %>%
  mutate(across(c(X0DPH, X7DPH, X14DPH), ~na_if(as.character(.), "."))) %>% 
  mutate(across(c(X0DPH, X7DPH, X14DPH), as.numeric)) %>%
  pivot_longer(
    cols = c(X0DPH, X7DPH, X14DPH),
    names_to = "Timepoint",
    values_to = "Height"
  ) %>%
  mutate(
    Timepoint = recode(Timepoint,
                       "X0DPH"  = "0 Days",
                       "X7DPH"  = "7 Days",
                       "X14DPH" = "14 Days"),
    Timepoint = factor(Timepoint, levels = c("0 Days","7 Days","14 Days")),
    Temp = factor(Temp)
  )

# Summarise mean + SE
sumdat <- long %>%
  group_by(Temp, Genotype, Timepoint) %>%
  summarise(
    n = sum(!is.na(Height)),
    mean = mean(Height, na.rm = TRUE),
    sd = sd(Height, na.rm = TRUE),
    se = sd / sqrt(n),
    .groups = "drop"
  )

# Plot
p <- ggplot(sumdat, aes(Timepoint, mean, color = Genotype, group = Genotype)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.12, linewidth = 0.6) +
  facet_wrap(~Temp, ncol = 1) +
  labs(
    title = "",
    x = "Timepoint",
    y = "Mean Height (cm)",
    color = "Genotype"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "left"
  )

p

sort(unique(sumdat$Genotype))

geno_cols <- c(
  "G33199 #3"       = "#FF2DA1",  # bright pink
  "G33199 #6"       = "#C40076",  # darker magenta
  
  "KS-22-PB-02-02"  = "#00A9FF",  # light blue/cyan
  "KS-22-PB-02-05"  = "#1F78B4",  # medium blue
  "KS-22-PB-05-03"  = "#005BBB",  # strong blue
  "KS-22-PB-05-07"  = "#08306B",  # navy
  
  "WI-22-AB-07-06"  = "#8A2BE2"   # purple-ish (if yours is more magenta, swap)
)


p <- p +
  scale_color_manual(values = geno_cols, breaks = names(geno_cols))
p

############ FINAL PLOT
library(dplyr)
library(ggplot2)

# 1) mapping: relabel + placement
key <- tibble::tribble(
  ~Genotype,        ~Genotype_label,   ~Placement,
  "KS-22-PB-02-02", "KS-22-PB-02-02",  "top10",
  "KS-22-PB-02-05", "KS-22-PB-02-05",  "top10",
  "KS-22-PB-05-03", "KS-22-PB-05-03",  "top10",
  "KS-22-PB-05-07", "KS-22-PB-05-07",  "top10",
  "WI-22-AB-07-06", "WI-22-AB-07-06",  "bottom10",
  "G33199 #3",      "NY-23-LG-01-03",  "bottom10",
  "G33199 #6",      "NY-23-LG-01-06",  "bottom10"
)

sumdat2 <- sumdat %>%
  left_join(key, by = "Genotype") %>%
  mutate(
    Genotype_label = ifelse(is.na(Genotype_label), Genotype, Genotype_label),
    Placement = factor(Placement, levels = c("top10", "bottom10"))
  ) %>%
  # optional: keep only those 7 genotypes
  filter(!is.na(Placement))

# 2) plot: colour by placement, shape by genotype (label)
p2 <- ggplot(sumdat2, aes(Timepoint, mean,
                          group = Genotype_label,
                          colour = Placement,
                          shape = Genotype_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.12, linewidth = 0.6) +
  facet_wrap(~Temp, ncol = 1) +
  scale_colour_manual(values = c(top10 = "red", bottom10 = "blue")) +
  labs(x = "Timepoint", y = "Mean Height (cm)",
       colour = "Placement", shape = "Genotype") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "left"
  )

p2
ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/FigS17C.pdf",
  plot     = p2,
  width    = 6,
  height   = 8,
  units    = "in"
)



########################################################################################################################################################################
# Figure S18
########################################################################################################################################################################
##############################
# US states Köppen–Geiger change through time
# Consecutive periods + SSP5-8.5
##############################
library(terra)
#packageVersion("terra")
#citation("terra")
library(sf)
library(tigris)
library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(scales)


kg_dir  <- "/Users/annamccormick/R/Feral_Cannabis_EGS/koppen_geiger_tif"
kg_file <- "koppen_geiger_0p1.tif"   # must exist in ALL periods used

# load US states 
options(tigris_use_cache = TRUE)
options(tigris_class = "sf")

states_sf <- tigris::states(cb = TRUE, year = 2022) |>
  filter(!STUSPS %in% c("PR", "VI", "GU", "MP", "AS"))

states_v <- terra::vect(states_sf)

# define time slices 
time_slices <- tibble(
  period   = c("1901_1930","1931_1960","1961_1990","1991_2020","2041_2070","2071_2099"),
  scenario = c(NA, NA, NA, NA, "ssp585", "ssp585"),
  label    = c("1901–1930","1931–1960","1961–1990","1991–2020","2041–2070 (ssp585)","2071–2099 (ssp585)"),
  midyear  = c(1915, 1945, 1975, 2005, 2055, 2085)
)

# build intervals explicitly
n <- nrow(time_slices)

intervals <- tibble(
  from_period    = time_slices$period[1:(n-1)],
  from_scenario  = time_slices$scenario[1:(n-1)],
  to_period      = time_slices$period[2:n],
  to_scenario    = time_slices$scenario[2:n],
  from_label     = time_slices$label[1:(n-1)],
  to_label       = time_slices$label[2:n],
  interval_label = paste0(time_slices$label[1:(n-1)], " → ", time_slices$label[2:n]),
  midyear        = (time_slices$midyear[1:(n-1)] + time_slices$midyear[2:n]) / 2
)

stopifnot(!any(is.na(intervals$from_period)))
stopifnot(!any(is.na(intervals$to_period)))

# raster load
load_kg <- function(period, scenario) {
  f <- if (is.na(scenario)) {
    file.path(kg_dir, period, kg_file)
  } else {
    file.path(kg_dir, period, scenario, kg_file)
  }
  if (!file.exists(f)) stop("Missing raster:\n", f)
  rast(f)
}

# % changed per state between two rasters 
pct_changed_by_state <- function(r1, r2, states_v) {
  
  # ensure same grid
  if (!compareGeom(r1, r2, stopOnError = FALSE)) {
    r2 <- resample(r2, r1, method = "near")
  }
  
  map_dfr(seq_len(nrow(states_v)), function(i) {
    
    st <- states_v[i]
    
    a <- mask(crop(r1, st), st)
    b <- mask(crop(r2, st), st)
    
    va <- values(a, mat = FALSE)
    vb <- values(b, mat = FALSE)
    
    keep <- !is.na(va) & !is.na(vb)
    va <- va[keep]
    vb <- vb[keep]
    
    tibble(
      state = st$NAME,
      pct_changed = ifelse(length(va) == 0, NA_real_, sum(va != vb) / length(va))
    )
  })
}

# loop across intervals 
state_change_ts <- map_dfr(seq_len(nrow(intervals)), function(i) {
  
  row <- intervals[i, ]
  
  r1 <- load_kg(row$from_period, row$from_scenario)
  r2 <- load_kg(row$to_period,   row$to_scenario)
  
  pct_changed_by_state(r1, r2, states_v) |>
    mutate(
      interval = row$interval_label,
      midyear  = row$midyear
    )
})

state_change_ts <- state_change_ts %>%
  mutate(interval = factor(interval, levels = unique(intervals$interval_label)))

# facet plot of % change 
state_change_ts2 <- state_change_ts %>%
  mutate(
    x_period = factor(
      sub("^.*→\\s*", "", interval),
      levels = c("1931–1960","1961–1990","1991–2020","2041–2070 (ssp585)","2071–2099 (ssp585)")
    ),
    x_period_folder = factor(
      c("1931_1960","1961_1990","1991_2020","2041_2070","2071_2099")[as.integer(x_period)],
      levels = c("1931_1960","1961_1990","1991_2020","2041_2070","2071_2099")
    )
  )

p_facet_lines_folderx <- ggplot(
  state_change_ts2,
  aes(x = x_period_folder, y = pct_changed, group = state)
) +
  geom_line(linewidth = 0.6, color = "black") +
  geom_point(size = 2) +
  facet_wrap(~state, ncol = 6) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Time window (folder name; value is change vs previous window)",
    y = "Percent of state area changing Köppen–Geiger class",
    title = ""
  ) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 90, hjust = 1),
    strip.text   = element_text(size = 8)
  )


# colour scheme
kg_colors <- c(
  "1"  = "#0000FF","2"  = "#0078FF","3"  = "#46A6FF",
  "4"  = "#FF0000","5"  = "#FF9696","6"  = "#F5A300","7"  = "#FFDC64",
  "8"  = "#FFFF00","9"  = "#C8FF00","10" = "#96FF00","11" = "#00FF00",
  "12" = "#00C800","13" = "#009600","14" = "#96FF96","15" = "#64C864",
  "16" = "#329632","17" = "#00FFFF","18" = "#37C8FF","19" = "#007DFF",
  "20" = "#0050C8","21" = "#B4E6FF","22" = "#64B4FF","23" = "#3296FF",
  "24" = "#0064C8","25" = "#C8C8FF","26" = "#9696FF","27" = "#6464FF",
  "28" = "#3232C8","29" = "#BEBEBE","30" = "#808080"
)

climate_labels <- c(
  `1` = "Af - Tropical, rainforest",
  `2` = "Am - Tropical, monsoon",
  `3` = "Aw - Tropical, Savanna",
  `4` = "BWh - Arid, desert, hot",
  `5` = "BWk - Arid, desert, cold",
  `6` = "BSh - Arid, steppe, hot",
  `7` = "BSk - Arid, steppe, cold",
  `8` = "Csa - Temperate, dry summer, hot summer",
  `9` = "Csb - Temperate, dry summer, warm summer",
  `10` = "Csc - Temperate, dry summer, cold summer",
  `11` = "Cwa - Temperate, dry winter, hot summer",
  `12` = "Cwb - Temperate, dry winter, warm summer",
  `13` = "Cwc - Temperate, dry winter, cold summer",
  `14` = "Cfa - Temperate, no dry season, hot summer",
  `15` = "Cfb - Temperate, no dry season, warm summer",
  `16` = "Cfc - Temperate, no dry season, cold summer",
  `17` = "Dsa - Cold, dry summer, hot summer",
  `18` = "Dsb - Cold, dry summer, warm summer",
  `19` = "Dsc - Cold, dry summer, cold summer",
  `20` = "Dsd - Cold, dry summer, very cold winter",
  `21` = "Dwa - Cold, dry winter, hot summer",
  `22` = "Dwb - Cold, dry winter, warm summer",
  `23` = "Dwc - Cold, dry winter, cold summer",
  `24` = "Dwd - Cold, dry winter, very cold winter",
  `25` = "Dfa - Cold, no dry season, hot summer",
  `26` = "Dfb - Cold, no dry season, warm summer",
  `27` = "Dfc - Cold, no dry season, cold summer",
  `28` = "Dfd - Cold, no dry season, very cold winter",
  `29` = "ET - Polar, tundra",
  `30` = "EF - Polar, frost"
)

# labels WITH numbers for legend, in the same order as colors
climate_labels_num <- setNames(
  paste0(names(climate_labels), " \u2014 ", unname(climate_labels)),
  names(climate_labels)
)

# KG composition by state function 
kg_composition_by_state <- function(r, states_v, state_name_col = "NAME") {
  
  purrr::map_dfr(seq_len(nrow(states_v)), function(i) {
    
    st <- states_v[i]
    st_name <- as.character(st[[state_name_col]][1])
    
    x <- terra::mask(terra::crop(r, st), st)
    v <- terra::values(x, mat = FALSE)
    v <- v[!is.na(v)]
    
    if (length(v) == 0) {
      return(tibble(
        state = st_name,
        KG_class = NA_character_,
        n_pix = 0L,
        total_pix = 0L,
        prop = NA_real_
      ))
    }
    
    tab <- as.data.frame(table(as.character(v)), stringsAsFactors = FALSE)
    colnames(tab) <- c("KG_class", "n_pix")
    tab$n_pix <- as.integer(tab$n_pix)
    
    tab %>%
      mutate(
        state = st_name,
        total_pix = sum(n_pix),
        prop = n_pix / total_pix
      ) %>%
      select(state, KG_class, n_pix, total_pix, prop)
  })
}

# compute composition across all time windows 
kg_state_time <- map_dfr(seq_len(nrow(time_slices)), function(i) {
  
  row <- time_slices[i, ]
  r <- load_kg(row$period, row$scenario)
  
  kg_composition_by_state(r, states_v) %>%
    mutate(
      period = row$period,
      label  = row$label
    )
})

# enforce factor levels for plotting order + legend stability
kg_state_time <- kg_state_time %>%
  mutate(
    label    = factor(label, levels = time_slices$label),
    KG_class = factor(KG_class, levels = names(kg_colors))
  )

# plot
p_kg_comp <- ggplot(
  kg_state_time,
  aes(x = label, y = prop, fill = KG_class)
) +
  geom_col(width = 0.9) +
  facet_wrap(~state, ncol = 6) +
  scale_fill_manual(
    values = kg_colors,
    limits = names(kg_colors),          # <- prevents blank legend swatches
    breaks = names(kg_colors),
    labels = climate_labels_num,
    drop   = FALSE
  ) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    x = "Time window",
    y = "Percent of state area",
    fill = "Köppen–Geiger climate class",
    title = ""
  ) +
  theme_bw() +
  theme(
    axis.text.x      = element_text(angle = 90, hjust = 1),
    strip.text       = element_text(size = 8),
    legend.position  = "right"
  )
p_kg_comp
# save composition plot
ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/state_koppenA.pdf",
  plot     = p_kg_comp,
  width    = 14,   # a bit wider helps the long legend labels
  height   = 10,
  units    = "in"
)





# colours for labels
kg_colors <- c(
  "1"  = "#0000FF","2"  = "#0078FF","3"  = "#46A6FF",
  "4"  = "#FF0000","5"  = "#FF9696","6"  = "#F5A300","7"  = "#FFDC64",
  "8"  = "#FFFF00","9"  = "#C8FF00","10" = "#96FF00","11" = "#00FF00",
  "12" = "#00C800","13" = "#009600","14" = "#96FF96","15" = "#64C864",
  "16" = "#329632","17" = "#00FFFF","18" = "#37C8FF","19" = "#007DFF",
  "20" = "#0050C8","21" = "#B4E6FF","22" = "#64B4FF","23" = "#3296FF",
  "24" = "#0064C8","25" = "#C8C8FF","26" = "#9696FF","27" = "#6464FF",
  "28" = "#3232C8","29" = "#BEBEBE","30" = "#808080"
)

climate_labels <- c(
  `1` = "Af - Tropical, rainforest",
  `2` = "Am - Tropical, monsoon",
  `3` = "Aw - Tropical, Savanna",
  `4` = "BWh - Arid, desert, hot",
  `5` = "BWk - Arid, desert, cold",
  `6` = "BSh - Arid, steppe, hot",
  `7` = "BSk - Arid, steppe, cold",
  `8` = "Csa - Temperate, dry summer, hot summer",
  `9` = "Csb - Temperate, dry summer, warm summer",
  `10` = "Csc - Temperate, dry summer, cold summer",
  `11` = "Cwa - Temperate, dry winter, hot summer",
  `12` = "Cwb - Temperate, dry winter, warm summer",
  `13` = "Cwc - Temperate, dry winter, cold summer",
  `14` = "Cfa - Temperate, no dry season, hot summer",
  `15` = "Cfb - Temperate, no dry season, warm summer",
  `16` = "Cfc - Temperate, no dry season, cold summer",
  `17` = "Dsa - Cold, dry summer, hot summer",
  `18` = "Dsb - Cold, dry summer, warm summer",
  `19` = "Dsc - Cold, dry summer, cold summer",
  `20` = "Dsd - Cold, dry summer, very cold winter",
  `21` = "Dwa - Cold, dry winter, hot summer",
  `22` = "Dwb - Cold, dry winter, warm summer",
  `23` = "Dwc - Cold, dry winter, cold summer",
  `24` = "Dwd - Cold, dry winter, very cold winter",
  `25` = "Dfa - Cold, no dry season, hot summer",
  `26` = "Dfb - Cold, no dry season, warm summer",
  `27` = "Dfc - Cold, no dry season, cold summer",
  `28` = "Dfd - Cold, no dry season, very cold winter",
  `29` = "ET - Polar, tundra",
  `30` = "EF - Polar, frost"
)

all_classes <- as.character(1:30)

# rebuild labels to include the class number
climate_labels_num <- setNames(
  paste0(all_classes, " — ", unname(climate_labels[all_classes])),
  all_classes
)

# enforce factor levels
kg_state_time <- kg_state_time %>%
  mutate(
    label    = factor(label, levels = time_slices$label),
    KG_class = factor(trimws(KG_class), levels = all_classes)
  )

# add dummy rows for ALL classes so legend keys never go blank 
dummy_legend <- tibble(
  state    = kg_state_time$state[1],          # any real state (must exist for faceting)
  label    = kg_state_time$label[1],          # any real time label
  KG_class = factor(all_classes, levels = all_classes),
  prop     = 0
)

#plot
p_kg_comp <- ggplot(
  kg_state_time,
  aes(x = label, y = prop, fill = KG_class)
) +
  geom_col(width = 0.9) +
  
  # this layer has zero-height bars but forces ggplot to train legend fills for ALL classes
  geom_col(
    data = dummy_legend,
    inherit.aes = FALSE,
    aes(x = label, y = prop, fill = KG_class),
    width = 0.9
  ) +
  
  facet_wrap(~state, ncol = 6) +
  scale_fill_manual(
    values = kg_colors[all_classes],
    breaks = all_classes,
    labels = climate_labels_num,
    drop   = FALSE,
    na.translate = FALSE
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Time window",
    y = "Percent of state area",
    fill = "Köppen–Geiger climate class",
    title = ""
  ) +
  theme_bw() +
  theme(
    axis.text.x      = element_text(angle = 90, hjust = 1),
    strip.text       = element_text(size = 8),
    legend.position  = "right"
  )

p_kg_comp

# save 
ggsave(
  filename = "/Users/annamccormick/R/Feral_Cannabis_EGS/RESULTS PPT/Figure_pannels_pdf/state_koppenA.pdf",
  plot     = p_kg_comp,
  width    = 14,   # a bit wider helps the long legend labels
  height   = 10,
  units    = "in"
)

########################################################################################################################################################################
########################################################################## FIN #########################################################################################
########################################################################################################################################################################


