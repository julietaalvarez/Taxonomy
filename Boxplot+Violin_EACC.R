##################################################################
#    Micromorphological Character Comparison in Macromycetes     #
#    Statistical & Visualization Workflow in R. Supports         #
#    parametric (ANOVA/Tukey) and non-parametric (Kruskal/Dunn)  #
##################################################################
# Author: Ezequiel Alberto Cruz Campuzano                     #
# Citation: Cruz-Campuzano et al., 2025: _DOI_                #
###############################################################

## NOTE: This script is configured to run Non-Parametric analyses by default.
## If your data meet parametric assumptions, you can switch to the Parametric pipeline 
## (ANOVA + Tukey) by following the comments provided throughout the script (see STEP 4 and STEP 5).

# ===============================================================
# STEP 1: Install and Load Required Packages
# ===============================================================

required_packages <- c("ggplot2", "dplyr", "rstatix", "showtext", "sysfonts", "cowplot", "car")
invisible(lapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}))

font_add_google("Roboto", "roboto")
showtext_auto()

# ===============================================================
# STEP 2: Load Dataset
# ===============================================================

data <- read.csv("Micromorphology_Measurements_Example_NonParametric.csv")

# ===============================================================
# STEP 3: Normality and Variance Diagnostics
# ===============================================================

characters_to_check <- c("Spore.length", "Spore.width", "Spore.Q", "Basidia.length", "Basidia.width")

diagnostic_summary <- data.frame(
  Variable = character(), Normality_All_Groups = character(),
  Homogeneity = character(), Suggested_Analysis = character(), stringsAsFactors = FALSE
)

for (var in characters_to_check) {
  normality_results <- data %>%
    group_by(Species) %>%
    summarise(p_value = shapiro.test(.data[[var]])$p.value)
  
  normal_pass <- all(normality_results$p_value > 0.05)
  
  levene_test <- car::leveneTest(reformulate("Species", var), data = data)
  homogeneity_pass <- levene_test$`Pr(>F)`[1] > 0.05
  
  suggested <- if (normal_pass & homogeneity_pass) "Parametric (ANOVA)" else "Non-parametric (Kruskal-Wallis)"
  
  diagnostic_summary <- rbind(diagnostic_summary, data.frame(
    Variable = var,
    Normality_All_Groups = ifelse(normal_pass, "PASS", "FAIL"),
    Homogeneity = ifelse(homogeneity_pass, "PASS", "FAIL"),
    Suggested_Analysis = suggested
  ))
}

cat("\nNormality and Variance Diagnostics:\n")
print(diagnostic_summary)

# ===============================================================
# STEP 4: Define Analysis and Plot Function (Flexible)
# ===============================================================

# If diagnostic tests suggest Normality and Homogeneity of data, determine 'parametric = TRUE'
# If the opposite, determine 'parametric = False'.
# If you wish to add more species, adjust the 'Comparisons <- list' command, and add
# as many as needed.
create_violin_plot <- function(data, response, letter, parametric = FALSE) {
  
  species_levels <- levels(factor(data$Species))
  group_positions <- function(g1, g2) {
    return((match(g1, species_levels) + match(g2, species_levels)) / 2)
  }

# If you wish to add more species, adjust the 'Comparisons <- list' 
# command, and add as many as needed.    
  comparisons <- list(
    c("Species A", "Species B"),
    c("Species A", "Species C"),
    c("Species B", "Species C")
  )
  
  y_min_raw <- min(data[[response]], na.rm = TRUE)
  y_max_raw <- max(data[[response]], na.rm = TRUE)
  y_range <- y_max_raw - y_min_raw
  y_min <- y_min_raw - 0.2 * y_range
  y_max <- y_max_raw + 0.1 * y_range
  
  if (parametric) {
    anova_res <- aov(reformulate("Species", response), data = data)
    tukey_res <- TukeyHSD(anova_res)
    tukey_df <- as.data.frame(tukey_res[[1]])
    tukey_df$Comparison <- rownames(tukey_df)
    tukey_df <- tukey_df %>%
      mutate(
        group1_raw = sub("^(.*)-.*$", "\\1", Comparison),
        group2_raw = sub("^.*-(.*)$", "\\1", Comparison),
        group1 = pmin(group1_raw, group2_raw),
        group2 = pmax(group1_raw, group2_raw)
      )
    
    # Establecer las comparaciones manuales (también con orden estandarizado)
    comparisons_df <- data.frame(
      group1 = sapply(comparisons, function(x) min(x)),
      group2 = sapply(comparisons, function(x) max(x))
    )
    
    dunn_labels <- comparisons_df %>%
      left_join(tukey_df, by = c("group1", "group2")) %>%
      mutate(
        x = mapply(group_positions, group1, group2),
        y = c(y_max * 0.98, y_max * 1.009, y_max * 1.035), # Heights ordered/Adjust here to modify position of Tukey labels
        label = paste0("P Tukey = ", signif(`p adj`, digits = 2)) 
      )
    
    kruskal_or_anova_p <- summary(anova_res)[[1]]$`Pr(>F)`[1]

  } else {
    kruskal_res <- kruskal.test(reformulate("Species", response), data = data)
    dunn_res <- data %>% dunn_test(reformulate("Species", response), p.adjust.method = "holm")
    
    dunn_labels <- data.frame(
      group1 = sapply(comparisons, `[[`, 1),
      group2 = sapply(comparisons, `[[`, 2)
    ) %>%
      left_join(dunn_res, by = c("group1", "group2")) %>%
      mutate(
        x = mapply(group_positions, group1, group2),
        y = c(y_max * 0.95, y_max * 1.0, y_max * 1.05), # Heights ordered/Adjust here to modify position of Dunn labels
        label = paste0("P Holm = ", signif(p.adj, digits = 2))
      )
    kruskal_or_anova_p <- kruskal_res$p.value
  }
  
  p <- ggplot(data, aes(x = Species, y = .data[[response]], fill = Species)) +
    geom_violin(trim = FALSE, alpha = 0.15) +
    geom_boxplot(width = 0.2, alpha = 0.3, color = "black", outlier.shape = NA) +
    stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "red", color = "black") +
    stat_summary(fun = mean, geom = "label", aes(label = sprintf("x = %.2f", ..y..)),
                 vjust = -0.6, hjust = -0.4, fill = "white", color = "black", family = "Roboto",
                 label.padding = unit(0.1, "lines"), size = 2.5) +
    geom_text(data = dunn_labels, aes(x = x, y = y, label = label), inherit.aes = FALSE,
              vjust = -0.5, size = 2.5, family = "Roboto") +
    geom_segment(data = dunn_labels, aes(x = match(group1, species_levels), xend = match(group2, species_levels),
                                         y = y * 0.995, yend = y * 0.995), inherit.aes = FALSE,
                 color = "black", linetype = "dashed") +
    annotate("text", x = 0.5, y = y_min + 0.05 * y_range, 
             label = paste0(ifelse(parametric, "ANOVA P-value = ", "Kruskal-Wallis P-value = "),
                            format(kruskal_or_anova_p, scientific = TRUE, digits = 3)),
             hjust = 0, size = 3, family = "Roboto") +
    annotate("text", x = -Inf, y = Inf, label = letter, hjust = -0.7, vjust = 1.5, size = 4) +
    labs(y = paste(gsub("\\.", " ", response), "(μm)"), x = NULL) +
    scale_fill_manual(values = c(
      "Species A" = "cyan", "Species B" = "purple", "Species C" = "salmon"
    )) +
# If you wish to adjust the Y axis manually to reach a certain number, you can edit
# 'limits = c(y_min, y_max)' for any specific numbers need as follows: 'limits = c(0, 20)'
# If you wish your Y axis to be divided into a certain number of breaks, 
# modify this code with the following: 
## scale_y_continuous(labels = scales::number_format(accuracy = 0.1), ##
## limits = c(y_min, y_max), breaks = pretty(c(y_min, y_max), n = 6)) ## <-- n = 6 sets the number of desired breaks
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1),
                       limits = c(y_min, y_max*1.055)) +  ### <- Modify your maximum Y axis value to meet requirements of your data
    theme_minimal(base_family = "Roboto") +
    theme(
      axis.title.y = element_text(face = "bold", size = 10),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(linetype = "dashed", color = "gray"),
      axis.line.x.bottom = element_line(color = "gray40"),
      axis.line.y.left = element_line(color = "gray40"),
      legend.position = "none",
      panel.background = element_blank(),
      plot.background = element_blank()
    )
  
  return(p)
}

# ===============================================================
# STEP 5: Generate Plots
# ===============================================================
# If diagnostic tests suggest Normality and Homogeneity of data, determine 'parametric = TRUE'
# if the opposite, determine 'parametric = False'

char_list <- list(
  list(var = "Spore.length", label = "A"),
  list(var = "Spore.width", label = "B"),
  list(var = "Spore.Q", label = "C"),
  list(var = "Basidia.length", label = "D"),
  list(var = "Basidia.width", label = "E")
)

plot_list <- lapply(char_list, function(x) create_violin_plot(data, x$var, x$label, parametric = FALSE))

# ===============================================================
# STEP 6: Combine All Plots
# ===============================================================

combined <- cowplot::plot_grid(plotlist = plot_list, ncol = 2, align = 'hv', rel_heights = rep(1, length(plot_list)))

# ===============================================================
# STEP 7: Export Final Figure
# ===============================================================

# Show on screen
print(combined)

# ===============================================================
# STEP 8: Export Statistical Results
# ===============================================================

# Save diagnostic summary table (normality and variance check)
write.csv(diagnostic_summary, file = "Example_Diagnostic_Summary.csv", row.names = FALSE)

# Define type of analysis (FALSE = Kruskal-Wallis, TRUE = ANOVA)
parametric <- FALSE

# Containers for test results
global_tests <- data.frame(
  Variable = character(),
  Test = character(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

comparison_results <- list()

# Loop through each character
for (x in char_list) {
  var <- x$var
  
  if (parametric) {
    aov_res <- aov(reformulate("Species", var), data = data)
    p_val <- summary(aov_res)[[1]][["Pr(>F)"]][1]
    
    global_tests <- rbind(global_tests, data.frame(
      Variable = var,
      Test = "ANOVA",
      P_value = p_val
    ))
    
    tukey_res <- TukeyHSD(aov_res)
    tukey_df <- as.data.frame(tukey_res[[1]])
    tukey_df$Comparison <- rownames(tukey_df)
    tukey_df$Variable <- var
    comparison_results[[var]] <- tukey_df
    
  } else {
    kruskal_res <- kruskal.test(reformulate("Species", var), data = data)
    p_val <- kruskal_res$p.value
    
    global_tests <- rbind(global_tests, data.frame(
      Variable = var,
      Test = "Kruskal-Wallis",
      P_value = p_val
    ))
    
    dunn_res <- data %>%
      dunn_test(reformulate("Species", var), p.adjust.method = "holm")
    dunn_res$Variable <- var
    comparison_results[[var]] <- dunn_res
  }
}

# Save global test results (if Parametric = ANOVA results / if Non-Parametric = Kruskal-Wallis results)
write.csv(global_tests, file = "Example_Global_Tests_NonParametric(KW).csv", row.names = FALSE)

# Save pairwise comparisons
all_comparisons <- do.call(rbind, comparison_results)
write.csv(all_comparisons, file = "Example_Multiple_Comparisons_NonParametric(Dunn).csv", row.names = FALSE)

  # NOTE:
  # After generating the final figure, you can optionally add a customized unified legend
  # or any additional labels manually using vector graphic editing software if desired.
  #######################################################################################