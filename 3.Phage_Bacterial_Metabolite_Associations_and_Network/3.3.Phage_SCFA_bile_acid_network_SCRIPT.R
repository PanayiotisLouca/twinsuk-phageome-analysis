## Author: 
#  Panayiotis Louca 

## Clear environment 
rm(list = ls()) 

## Set seed 
set.seed(1234)

## load up packages: 

### core 
library(tidyverse)

library(Hmisc)
library(igraph)
library(corrplot)
library(ggraph)
library(RColorBrewer)

# -------------------------------------------------------------------------- # 

# ************************* # 
#   IMPORT & PREP DATA   ---- 
# ************************* # 

df <- read_rds(file.path("/Users/panayiotislouca/Documents/KCL_files/Phages_metagenome_characterisation/Data/Revised_phages_microbes_metabs_DATASET.rds")) %>%
  as.data.frame()

# import phage - SCFA/BA results 
phage_metab_res <- read.csv(file.path("/Users/panayiotislouca/Documents/KCL_files/Phages_metagenome_characterisation/Analysis/Phages_SCFA_and_Bile_acids/Combined_phage_bile_acid_and_phage_SCFA_results.csv"))

phage_metab_res_sig <- phage_metab_res %>%
  filter(fdr_combined < 0.01)

# phage_assoc metab_names 
phage_assoc_metab_names = phage_metab_res_sig %>%
  pull(depend_var) %>%
  unique()
  
metab_assoc_phage_names = phage_metab_res_sig %>%
  pull(independ_var) %>%
  unique()

# ---------------------------------------------------------------------------- #  

df_plot <- df %>%
  select(all_of(phage_assoc_metab_names),
         all_of(metab_assoc_phage_names))

# -------------------------------------------------------------------------- #  

# Correlation Network Visualization 
phage_data <- df[ ,metab_assoc_phage_names]
metab_data <- df[ ,c(phage_assoc_metab_names)]

# -------------------------------------------------------------------------- #  

# Compute Spearman correlations between phages and bacteria 
cor_matrix <- cor(phage_data, metab_data, method = "spearman", use = "pairwise.complete.obs")

# Convert to long format and filter for strong correlations
cor_df <- as.data.frame(as.table(cor_matrix)) %>%
  rename(Phage = Var1, Metabolite = Var2, Correlation = Freq) %>%
  filter(abs(Correlation) > 0.2) %>% # reduce edges 
  mutate(EdgeColour = ifelse(Correlation > 0, "Positive", "Negative"))

#  Compute p-values 
cor_df <- cor_df %>%
  mutate(p_value = pmap_dbl(list(Phage, Metabolite), function(phg, metab) {
    x <- phage_data[[phg]]
    y <- metab_data[[metab]]
    suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE)$p.value)
  }))


# -------------------------------------------------------------------------- #  

# Get unique combinations 
unique_pairs <- unique(phage_metab_res[c("depend_var", "Metabolite_label")])
metab_mapping <- setNames(unique_pairs$depend_var, unique_pairs$Metabolite_label)

# Reverse the mapping so names match cor_df$Metabolite 
metab_lookup <- setNames(names(metab_mapping), metab_mapping)

# Create a data frame for the mapping 
metab_lookup_df <- tibble(
  Metabolite = names(metab_lookup),
  Metabolite_label = unname(metab_lookup)
)

# Apply mapping 
cor_df <- cor_df %>%
  left_join(metab_lookup_df, by = "Metabolite") %>%
  mutate(Metabolite = Metabolite_label) %>%
  select(-Metabolite_label)


# -------------------------------------------------------------------------- #  

# Create network and set attributes 
net <- graph_from_data_frame(cor_df[, c("Phage", "Metabolite", "Correlation")], directed = FALSE)

# Set edge weights and colors 
E(net)$weight <- abs(cor_df$Correlation)
E(net)$color <- cor_df$EdgeColour

# Set node attributes 
V(net)$type <- ifelse(grepl("phage", V(net)$name), "Phage", "Metabolite")

# Labels 
degree_vals <- degree(net)

V(net)$label <- gsub("_", " ", V(net)$name) 

# -------------------------------------------------------------------------- #  

n_connections <- 100  # Threshold for labeling nodes 
Metabolite_nodes <- V(net)$name[V(net)$type == "Metabolite"]

# -------------------------------------------------------------------------- #  

# set seed 
set.seed(12345)

# save plot 
png(file.path("/Users/panayiotislouca/Documents/KCL_files/Phages_metagenome_characterisation/Analysis/Phages_SCFA_and_Bile_acids/Phages_SCFA_bile_acid_network/Phage_SCFA_bile_acid_network_PLOT.png"), 
    width = 16, height = 14, units = 'in', res = 300)

# Plot network 
ggraph(net, layout = "fr", niter = 5000) +
  geom_edge_link0(aes(edge_alpha = weight, edge_colour = color, edge_width = weight), 
                  show.legend = FALSE) +
  scale_edge_alpha_continuous(range = c(0.2, 1)) +
  scale_edge_width_continuous(range = c(0.1, 1)) +
  geom_node_point(aes(color = type, size = degree_vals, 
                      shape = ifelse(degree_vals > quantile(degree_vals, 0.95), "triangle", "circle")), 
                  alpha = 0.9) +
  geom_node_text(aes(label = 
                       ifelse(degree_vals >= n_connections, label, "")),
  size = 4,
repel = TRUE,
force = 8,
box.padding = 1.5,
point.padding = 0,
max.overlaps = 4500,
                  color = "black", fontface = "bold", segment.color = "grey80") +
  scale_edge_colour_manual(values = c("Positive" = "#E41A1C", "Negative" = "#377EB8"),  
                           name = "Correlation") +
  scale_color_manual(values = c("Phage" = "#6A3D9A", "Metabolite" = "#33A02C"),
                     name = "Node Type") +
  scale_shape_manual(values = c("circle" = 16, "triangle" = 17), guide = "none") +
  scale_size_area(max_size = 6) +
  guides(
    size = "none",
    color = guide_legend(
      title = "Node Type",
      override.aes = list(size = 6, shape = 16),
      keywidth = 1.5,
      keyheight = 1.2
    )
  ) +
  theme_void() +
  theme(
    # Panel and plot appearance 
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    
    # Legend styling 
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    panel.grid.major = element_line(color = "grey90", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 0.5)        ) +
  labs(title = "Phage-SCFA/Bile Acid Co-abundance Network",
       subtitle = "Rho > 0.2; Red = Positive, Blue = Negative" 
       )


dev.off()

############################################################################################### 
##################################### END ##################################################### 
############################################################################################### 