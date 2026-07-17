############################################################
# Dudleya leaf morphology analysis — no DULA
# Revised 15 July 2026 for the Madroño submission analysis
#
# Plant-level means, plots, Kruskal-Wallis tests, and
# Dunn post hoc comparisons with Holm adjustment.
#
# DULA is retained in the original input file but excluded
# explicitly below. Boulder Ridge plants coded DUCYxDUSE in
# the input file are relabeled BOU for neutral presentation.
############################################################

# ---- 1. Packages ----

packages <- c("tidyverse", "readr", "FSA", "multcompView", "patchwork")

installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) install.packages(p)
}

library(tidyverse)
library(readr)
library(FSA)
library(multcompView)
library(patchwork)

# ---- 2. Read data ----

df <- read_csv(file.choose(), na = c("", "NA", "N/A"))

species_order <- c("DUSE", "DUCY", "DUAB", "BOU")

# Colors align as closely as possible with Figs. 4 and 5:
# focal DUSE = black, DUCY = orange, DUAB = green.
# BOU is gray to avoid implying a taxonomic or hybrid assignment.
species_colors <- c(
  "DUSE" = "#000000",
  "DUCY" = "#F28E2B",
  "DUAB" = "#009E73",
  "BOU"  = "#8C8C8C"
)

species_shapes <- c(
  "DUSE" = 16,
  "DUCY" = 15,
  "DUAB" = 4,
  "BOU"  = 17
)

# ---- 3. Clean species column ----

df <- df %>%
  mutate(
    Species_original = str_trim(as.character(`Actual Species (JBW)`)),
    Species_original = na_if(Species_original, ""),
    Species = recode(Species_original, "DUCYxDUSE" = "BOU")
  ) %>%
  filter(!is.na(Species)) %>%
  filter(Species_original != "DULA") %>%
  filter(Species %in% species_order) %>%
  mutate(Species = factor(Species, levels = species_order)) %>%
  droplevels()

print(table(df$Species, useNA = "ifany"))

# Reproducibility checks: analysis must contain only the four intended groups.
stopifnot(!any(df$Species_original == "DULA", na.rm = TRUE))
stopifnot(setequal(as.character(unique(df$Species)), species_order))

# ---- 4. Measurement columns ----

measure_cols <- c(
  "Length (mm)",
  "Width at Midpoint (mm)",
  "Width at Widest (mm)",
  "Width at Base (mm)"
)

# ---- 5. Plant-level means ----

plant_means <- df %>%
  group_by(`Sample #`, Species) %>%
  summarize(
    across(all_of(measure_cols), ~ mean(.x, na.rm = TRUE)),
    n_replicates = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Ratio_Length_Midpoint = `Length (mm)` / `Width at Midpoint (mm)`,
    Ratio_Length_Widest   = `Length (mm)` / `Width at Widest (mm)`,
    Ratio_Length_Base     = `Length (mm)` / `Width at Base (mm)`
  ) %>%
  droplevels()

write_csv(plant_means, "Dudleya_plant_level_means_with_ratios.csv")

# ---- 6. Long format ----

traits_direct <- measure_cols

traits_ratios <- c(
  "Ratio_Length_Midpoint",
  "Ratio_Length_Widest",
  "Ratio_Length_Base"
)

all_traits <- c(traits_direct, traits_ratios)

plant_long <- plant_means %>%
  pivot_longer(
    cols = all_of(all_traits),
    names_to = "Trait",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value)) %>%
  droplevels()

# ---- 7. Labels ----

trait_labels <- c(
  "Length (mm)" = "Length (mm)",
  "Width at Midpoint (mm)" = "Width at Midpoint (mm)",
  "Width at Widest (mm)" = "Width at Widest (mm)",
  "Width at Base (mm)" = "Width at Base (mm)",
  "Ratio_Length_Midpoint" = "Ratio Length / Width at Midpoint",
  "Ratio_Length_Widest" = "Ratio Length / Width at Widest",
  "Ratio_Length_Base" = "Ratio Length / Width at Base"
)

# ---- 8. Kruskal-Wallis + Dunn letters ----

get_dunn_letters <- function(data, trait_name) {
  
  sub <- data %>%
    filter(Trait == trait_name) %>%
    filter(!is.na(Value)) %>%
    filter(!is.na(Species)) %>%
    droplevels()
  
  kw <- kruskal.test(Value ~ Species, data = sub)
  
  if (kw$p.value >= 0.05) {
    letters_df <- sub %>%
      distinct(Species) %>%
      mutate(letter = "a")
    
    return(list(
      kruskal = kw,
      dunn = NULL,
      letters = letters_df
    ))
  }
  
  dunn <- FSA::dunnTest(
    Value ~ Species,
    data = sub,
    method = "holm"
  )
  
  dunn_res <- dunn$res
  
  pvals <- dunn_res$P.adj
  names(pvals) <- gsub(" ", "", dunn_res$Comparison)
  
  letters <- multcompLetters(pvals, threshold = 0.05)$Letters
  
  letters_df <- tibble(
    Species = names(letters),
    letter = letters
  ) %>%
    mutate(Species = factor(Species, levels = species_order)) %>%
    filter(!is.na(Species))
  
  return(list(
    kruskal = kw,
    dunn = dunn_res,
    letters = letters_df
  ))
}

# ---- 9. Run and save stats ----

stats_summary <- lapply(all_traits, function(tr) {
  out <- get_dunn_letters(plant_long, tr)
  
  tibble(
    Trait = tr,
    Kruskal_statistic = as.numeric(out$kruskal$statistic),
    Kruskal_p = out$kruskal$p.value
  )
}) %>%
  bind_rows()

write_csv(stats_summary, "Dudleya_Kruskal_summary_all_traits.csv")
print(stats_summary)

for (tr in all_traits) {
  out <- get_dunn_letters(plant_long, tr)
  
  if (!is.null(out$dunn)) {
    write_csv(
      out$dunn,
      paste0("Dunn_posthoc_", make.names(tr), ".csv")
    )
  }
}

# ---- 10. Plot function ----

plot_trait <- function(trait_name, plot_type = "boxplot") {
  
  sub <- plant_long %>%
    filter(Trait == trait_name) %>%
    filter(!is.na(Value)) %>%
    filter(!is.na(Species)) %>%
    droplevels()
  
  stat_out <- get_dunn_letters(plant_long, trait_name)
  
  # Place each significance letter the same vertical distance above the
  # highest observed value for that group. This keeps letters clear of
  # outliers and jittered points while preserving group-specific placement.
  value_range <- diff(range(sub$Value, na.rm = TRUE))
  if (!is.finite(value_range) || value_range == 0) value_range <- max(sub$Value, na.rm = TRUE) * 0.10
  letter_offset <- 0.07 * value_range

  letters_df <- stat_out$letters %>%
    left_join(
      sub %>%
        group_by(Species) %>%
        summarize(
          group_max = max(Value, na.rm = TRUE),
          .groups = "drop"
        ),
      by = "Species"
    ) %>%
    mutate(y_pos = group_max + letter_offset)

  plot_upper <- max(letters_df$y_pos, na.rm = TRUE) + 0.06 * value_range
  
  base <- ggplot(
    sub,
    aes(
      x = Species,
      y = Value,
      fill = Species
    )
  )
  
  if (plot_type == "boxplot") {
    p <- base +
      geom_boxplot(
        width = 0.65,
        outlier.shape = 21,
        alpha = 0.75,
        color = "black"
      )
  }
  
  if (plot_type == "violin") {
    p <- base +
      geom_violin(
        trim = FALSE,
        alpha = 0.55,
        color = "black"
      ) +
      geom_boxplot(
        width = 0.14,
        outlier.shape = NA,
        alpha = 0.85,
        color = "black"
      )
  }
  
  if (plot_type == "violin_points") {
    p <- base +
      geom_violin(
        trim = FALSE,
        alpha = 0.45,
        color = "black"
      ) +
      geom_jitter(
        width = 0.08,
        size = 1.8,
        alpha = 0.55,
        show.legend = FALSE
      ) +
      geom_boxplot(
        width = 0.12,
        outlier.shape = NA,
        alpha = 0.70,
        color = "black"
      )
  }
  
  p <- p +
    geom_text(
      data = letters_df,
      aes(x = Species, y = y_pos, label = letter),
      inherit.aes = FALSE,
      size = 5
    ) +
    scale_x_discrete(limits = species_order, drop = TRUE) +
    scale_y_continuous(
      limits = c(0, plot_upper),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_fill_manual(values = species_colors, drop = TRUE) +
    labs(
      title = NULL,
      x = NULL,
      y = trait_labels[[trait_name]]
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      aspect.ratio = 0.75
    )
  
  return(p)
}

# ---- 11. Generate all boxplots and violin plots ----

for (tr in all_traits) {
  
  ggsave(
    filename = paste0("Boxplot_", make.names(tr), ".png"),
    plot = plot_trait(tr, "boxplot"),
    width = 7,
    height = 5,
    dpi = 300
  )
  
  ggsave(
    filename = paste0("Violin_", make.names(tr), ".png"),
    plot = plot_trait(tr, "violin"),
    width = 7,
    height = 5,
    dpi = 300
  )
  
  ggsave(
    filename = paste0("Violin_points_", make.names(tr), ".png"),
    plot = plot_trait(tr, "violin_points"),
    width = 7,
    height = 5,
    dpi = 300
  )
}

# ---- 12. Scatterplots ----

width_traits <- c(
  "Width at Midpoint (mm)",
  "Width at Widest (mm)",
  "Width at Base (mm)"
)

for (w in width_traits) {
  
  p <- ggplot(
    plant_means,
    aes(
      x = `Length (mm)`,
      y = .data[[w]],
      color = Species,
      shape = Species
    )
  ) +
    geom_point(size = 3, alpha = 0.85) +
    scale_color_manual(values = species_colors, breaks = species_order, drop = TRUE) +
    scale_shape_manual(values = species_shapes, breaks = species_order, drop = TRUE) +
    scale_x_continuous(
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_y_continuous(
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      title = paste("Length vs", w),
      x = "Length (mm)",
      y = w,
      color = "Species"
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      aspect.ratio = 0.75
    )
  
  ggsave(
    filename = paste0("Scatter_Length_vs_", make.names(w), ".png"),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}

# ---- 13. Four-panel figure ----

# Panel order follows the manuscript figure deck:
# A = leaf length, B = maximum leaf width,
# C = length / maximum width, D = length vs. maximum width.
p_length <- plot_trait("Length (mm)", "boxplot") +
  labs(y = "Length (mm)")

p_width <- plot_trait("Width at Widest (mm)", "boxplot") +
  labs(y = "Maximum Width (mm)")

p_ratio <- plot_trait("Ratio_Length_Widest", "boxplot") +
  labs(y = "Length / Maximum Width")

p_scatter <- ggplot(
  plant_means,
  aes(
    x = `Length (mm)`,
    y = `Width at Widest (mm)`,
    color = Species,
    shape = Species
  )
) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_manual(
    values = species_colors,
    breaks = species_order,
    drop = TRUE
  ) +
  scale_shape_manual(
    values = species_shapes,
    breaks = species_order,
    drop = TRUE
  ) +
  scale_x_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Length vs. Maximum Width",
    x = "Length (mm)",
    y = "Maximum Width (mm)",
    color = NULL,
    shape = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    aspect.ratio = 0.75
  )

four_panel <- (p_length | p_width) / (p_ratio | p_scatter) +
  plot_layout(widths = c(1, 1), heights = c(1, 1)) +
  plot_annotation(tag_levels = "A")

ggsave(
  filename = "Figure6_Dudleya_leaf_morphology_no_DULA.png",
  plot = four_panel,
  width = 11,
  height = 8.5,
  dpi = 600,
  bg = "white"
)

# ---- 14. QC checks ----

duse_check <- plant_means %>%
  filter(Species == "DUSE") %>%
  mutate(Ratio_Length_Widest = `Length (mm)` / `Width at Widest (mm)`) %>%
  arrange(Ratio_Length_Widest) %>%
  select(
    `Sample #`,
    Species,
    `Length (mm)`,
    `Width at Widest (mm)`,
    Ratio_Length_Widest
  )

write_csv(duse_check, "DUSE_broad_leaf_QC_check.csv")

print("Lowest DUSE length/widest-width ratios:")
print(head(duse_check, 10))

ratio_summary <- plant_means %>%
  group_by(Species) %>%
  summarize(
    n = n(),
    mean_ratio = mean(Ratio_Length_Widest, na.rm = TRUE),
    sd_ratio = sd(Ratio_Length_Widest, na.rm = TRUE),
    se_ratio = sd_ratio / sqrt(n),
    .groups = "drop"
  )

write_csv(ratio_summary, "Species_summary_Ratio_Length_Widest.csv")
print(ratio_summary)

############################################################
# End
############################################################