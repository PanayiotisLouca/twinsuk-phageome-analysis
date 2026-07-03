## Author: 
#  Panayiotis Louca 

## Clear environment 
rm(list = ls()) 

## Set seed 
set.seed(1234)

## set working directory 
setwd("/Users/panayiotislouca/Documents")

## load up packages: 

### core 
library(tidyverse)
library(skimr)

library(Hmisc)
library(igraph)
library(corrplot)
library(ggraph)
library(RColorBrewer)

# -------------------------------------------------------------------------- # 

# ************************* # 
#   IMPORT & PREP DATA   ---- 
# ************************* # 

##   Dataset ---- 
df <- read_rds(file.path("/Users/panayiotislouca/Documents/KCL_files/Phages_metagenome_characterisation/Data/Revised_phages_microbes_metabs_DATASET.rds")) %>% 
  as.data.frame() %>%
  tibble::remove_rownames() 

# vectorise columns
phage_columns <- grep("^phage_[0-9]+", names(df), value = TRUE)
taxa_columns = grep("species", names(df), value = TRUE)

# -------------------------------------------------------------------------- # 

# CLR transform phages and taxa 

# function to CLR normalise 
CLRnorm <- function(features) {
  if (!is.data.frame(features) && !is.matrix(features)) {
    stop("Input must be a data frame or matrix.")
  }
  features_norm <- as.matrix(features)
  features_norm[features_norm == 0] <- 1e-6
  features_CLR <- chemometrics::clr(features_norm)
  as.data.frame(features_CLR, colnames = colnames(features))
}

# Apply CLR normalisation 
df[, phage_columns] = CLRnorm(df[, phage_columns])
df[, taxa_columns] = CLRnorm(df[, taxa_columns])

# -------------------------------------------------------------------------- #  

phage_data <- df[ ,phage_columns]
taxa_data <- df[ ,taxa_columns]

# modify names with annotated names 
path <- file.path("/Users/panayiotislouca/Documents/KCL_files/Data/TwinsUK/Microbiome/Metagenome/GenomeScan/GenomeScan_taxa_annotation_PL.csv")
gm_annot <- read.csv(path) %>% as.data.frame()

tmp <- setNames(gm_annot$clade_name_clean, gm_annot$annotation)

names(taxa_data) <- tmp[names(taxa_data)] 

# -------------------------------------------------------------------------- #  

# Compute Spearman correlations between phages and bacteria 
cor_matrix <- cor(phage_data, taxa_data, method = "spearman")

# Convert to long format and filter for strong correlations
cor_df <- as.data.frame(as.table(cor_matrix)) %>%
  rename(Phage = Var1, Bacteria = Var2, Correlation = Freq) %>%
  filter(Correlation < -0.3) %>% # to reduce edges 
  mutate(EdgeColour = ifelse(Correlation > 0, "Positive", "Negative"))

#  Compute p-values 
cor_df <- cor_df %>%
  mutate(p_value = pmap_dbl(list(Phage, Bacteria), function(phg, bac) {
    x <- phage_data[[phg]]
    y <- taxa_data[[bac]]
    suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE)$p.value)
  }))

# -------------------------------------------------------------------------- #  

# Create network and set attributes 
net <- graph_from_data_frame(cor_df[, c("Phage", "Bacteria", "Correlation")], directed = FALSE)

# Set edge weights and colours 
E(net)$weight <- abs(cor_df$Correlation)
E(net)$color <- cor_df$EdgeColour

# Set node attributes 
V(net)$type <- ifelse(grepl("phage", V(net)$name), "Phage", "Bacteria")

# Labels 
degree_vals <- degree(net)
V(net)$label <- gsub("_", " ", V(net)$name) # Store all names 

n_connections <- 100  # Threshold for labeling nodes 

bacteria_nodes <- V(net)$name[V(net)$type == "Bacteria"]

# -------------------------------------------------------------------------- # 

# save plot 
set.seed(123456)

png(file.path("/Users/panayiotislouca/Documents/KCL_files/Phages_metagenome_characterisation/Analysis/Phage_species_networks/Phage_species_ggraph_network_NEGATIVE_PLOT.png"),
    width = 16, height = 11, units = 'in', res = 300)

# Plot network 
ggraph(net, layout = "fr", niter = 5000) +
  geom_edge_link0(aes(edge_alpha = weight, edge_colour = color, edge_width = weight), 
                  show.legend = FALSE) +
  scale_edge_alpha_continuous(range = c(0.2, 1)) +
  scale_edge_width_continuous(range = c(0.1, 1)) +
  geom_node_point(aes(color = type, size = degree_vals, 
                      shape = ifelse(degree_vals > quantile(degree_vals, 0.95), "triangle", "circle")), 
                  alpha = 0.9) +
   geom_node_text(aes(label = ifelse(degree_vals >= n_connections, label, "")),
                  size = 5,
                  repel = TRUE,
                  force = 2,
                  box.padding = 0.5,
                  point.padding = 0,
                  max.overlaps = 4500,
                  color = "black", fontface = "bold", segment.color = "grey80") +
  scale_edge_colour_manual(values = c("Positive" = "#E41A1C", "Negative" = "#377EB8"),  
                           name = "Correlation") +
  scale_color_manual(values = c("Phage" = "#6A3D9A", "Bacteria" = "#33A02C"),
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
  labs(title = "Phage-Bacteria Correlation Network",
       subtitle = "Rho < -0.3; Red = Positive, Blue = Negative" ,
        caption = paste0("Nodes labeled for  bacteria with at least 200 phage correlations.")
       )


dev.off()

############################################################################################### 
##################################### END ##################################################### 
############################################################################################### 