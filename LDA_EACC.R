###############################################################
#   Linear Discriminant Analysis (LDA) for Species            #
#   Delimitation with Morphometric Characters in Macromycetes #               #
###############################################################
# Author: Ezequiel Alberto Cruz Campuzano                     #
# Citation: Cruz-Campuzano et al., 2025: _DOI_                #
###############################################################

## NOTE: This script performs LDA using standardized (Z-score) data.
## It produces a biplot of species separation and variable contribution,
## and a barplot of variable importance across LD axes.

# ===============================================================
# STEP 1: Install and Load Required Packages
# ===============================================================

required_packages <- c("ggplot2", "MASS", "dplyr", "tidyr", "gridExtra", "showtext", "sysfonts", "grid")
invisible(lapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}))

font_add_google("Roboto", "roboto")
showtext_auto()

# ===============================================================
# STEP 2: Load and Transform Data (Z-scores)
# ===============================================================

rawdata <- read.csv("Micromorphology_Measurements_Example_NonParametric.csv")

# Convert numeric columns to Z-scores, excluding 'Species'
data_z <- rawdata %>%
  mutate(across(where(is.numeric), ~ (.-mean(.)) / sd(.)))

write.csv(data_z, "Micromorphology_Measurements_Example_Z-Scores.csv", row.names = FALSE)

# ===============================================================
# STEP 3: Perform LDA and Extract Results
# ===============================================================
data <- read.csv("Micromorphology_Measurements_Example_Z-Scores.csv")
data$Species <- as.factor(data$Species)

lda_model <- lda(Species ~ ., data = data)
lda_values <- predict(lda_model)

# Variance explained by each axis
lda_var_explained <- lda_model$svd^2 / sum(lda_model$svd^2) * 100

# Projected coordinates
lda_df <- data.frame(Species = data$Species, LD1 = lda_values$x[,1], LD2 = lda_values$x[,2])

# Contribution of each variable
lda_contrib <- lda_model$scaling^2
lda_contrib <- apply(lda_contrib, 2, function(x) x / sum(x) * 100)
lda_contrib_df <- data.frame(
  Variable = rownames(lda_contrib),
  LD1_Contribution = round(lda_contrib[,1], 2),
  LD2_Contribution = round(lda_contrib[,2], 2)
)

write.csv(lda_contrib_df, "VariableContribution_LDA_Example.csv", row.names = FALSE)

# ===============================================================
# STEP 3B: Calculate and Export Mahalanobis Distances
# ===============================================================

# I. Extract projections of the data in LDA space
lda_proj <- lda_values$x

# II. Calculate centroids per species in LDA space
species_levels <- levels(data$Species)
centroids <- aggregate(lda_proj, by = list(Species = data$Species), FUN = mean)
rownames(centroids) <- centroids$Species
centroids <- centroids[, -1]  # Remove the 'Species' column
centroids <- as.matrix(centroids)  # Ensure it's numeric matrix

# III. Covariance matrix of LDA projections
cov_matrix <- cov(lda_proj)

# IV. Mahalanobis distance of each sample to its own group centroid
mahal_dist <- numeric(nrow(data))
for (i in 1:nrow(data)) {
  sp <- as.character(data$Species[i])
  centroid_vec <- centroids[sp, , drop = FALSE]  # drop = FALSE keeps the result as a matrix
  mahal_dist[i] <- mahalanobis(lda_proj[i, ], center = centroid_vec, cov = cov_matrix)
}

# V. Save individual distances
mahal_df <- data.frame(Sample_ID = seq_len(nrow(data)), Species = data$Species, Mahalanobis_Distance = mahal_dist)
write.csv(mahal_df, "Mahalanobis_Distances_LDA.csv", row.names = FALSE)

# VI. Mahalanobis distances between group centroids
mahal_centroid_matrix <- matrix(NA, nrow = length(species_levels), ncol = length(species_levels),
                                dimnames = list(species_levels, species_levels))

for (i in 1:length(species_levels)) {
  for (j in 1:length(species_levels)) {
    mahal_centroid_matrix[i, j] <- mahalanobis(centroids[i, ], center = centroids[j, ], cov = cov_matrix)
  }
}

write.csv(mahal_centroid_matrix, "Mahalanobis_CentroidDistances.csv", row.names = TRUE)

# ===============================================================
# STEP 3C: Optional Filtering of Top Contributing Variables
# ===============================================================

## If desired, you can filter to show only the variables with the highest contribution
## to LD1 and/or LD2. This is useful when a few variables dominate the plot,
## making others visually negligible.

## Example: Filter variables in the top 25% of contribution to LD1 or LD2
top_vars <- lda_contrib_df %>%
  filter(LD1_Contribution > quantile(LD1_Contribution, 0.75) | 
           LD2_Contribution > quantile(LD2_Contribution, 0.75))

## You may replace 'lda_contrib_df' with 'top_vars' in the plot section to visualize only those arrows:
##   - Replace geom_segment(data = ...) and geom_text(data = ...) accordingly.

# ===============================================================
# STEP 4: Create LDA Biplot with Arrows
# ===============================================================

scale_factor <- 0.065  # Scaling factor for arrow length in the biplot.
                       # Increase this value if arrows are too short, decrease if they are too long.
                       # To use raw (unscaled) contributions, set to 1 or remove the scaling from the plot code.
species_colors <- c("Species A" = "#9933FF", "Species B" = "#66FFFF", "Species C" = "#FF0000")

lda_plot <- ggplot(lda_df, aes(x = LD1, y = LD2, color = Species)) +
  geom_point(size = 3, alpha = 0.6) +
  stat_ellipse(aes(fill = Species), alpha = 0.2, geom = "polygon", color = NA) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  labs(
    x = paste("LD1 (", round(lda_var_explained[1]), "%)", sep = ""),
    y = paste("LD2 (", round(lda_var_explained[2]), "%)", sep = "")
  ) +
  geom_segment(data = lda_contrib_df,
               aes(x = 0, y = 0,
                   xend = LD1_Contribution * scale_factor,
                   yend = LD2_Contribution * scale_factor),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", size = 0.8) +
  geom_text(data = lda_contrib_df,
            aes(x = LD1_Contribution * scale_factor,
                y = LD2_Contribution * scale_factor,
                label = Variable),
            color = "black", size = 3, hjust = 0.5, vjust = 1.5) +
  theme_minimal(base_family = "Roboto") +
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "grey50"),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linetype = "dashed"),
    legend.position = c(0.925, 0.915),
    legend.title = element_blank()
  )

# ===============================================================
# STEP 5: Create Contribution Barplot
# ===============================================================

lda_contrib_long <- lda_contrib_df %>%
  pivot_longer(cols = starts_with("LD"), names_to = "Axis", values_to = "Contribution")

bar_plot <- ggplot(lda_contrib_long, aes(x = reorder(Variable, -Contribution), y = Contribution, fill = Axis)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("LD1_Contribution" = "lightblue", "LD2_Contribution" = "lightcoral")) +
  labs(x = "Variable", y = "Contribution (%)", fill = "Axis") +
  theme_minimal(base_family = "Roboto") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "grey50"),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linetype = "dashed"),
    legend.position = c(0.90, 0.915),
    legend.title = element_blank()
  )

# ===============================================================
# STEP 6: Combine Both Plots and Label Panels
# ===============================================================

label_a <- textGrob("A", x = unit(0.02, "npc"), y = unit(0.5, "npc"),
                    just = c("left", "top"), gp = gpar(fontsize = 16, fontface = "bold"))
label_b <- textGrob("B", x = unit(0.02, "npc"), y = unit(0.5, "npc"),
                    just = c("left", "top"), gp = gpar(fontsize = 16, fontface = "bold"))

lda_plot_labeled <- arrangeGrob(lda_plot, top = label_a)
bar_plot_labeled <- arrangeGrob(bar_plot, top = label_b)

grid.arrange(lda_plot_labeled, bar_plot_labeled, ncol = 2)

# NOTE:
# After generating the final figure, you can customize sizes using vector graphic editing software
##################################################################################################