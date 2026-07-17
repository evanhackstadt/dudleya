#### This is to create Mantel tests using Dudleya genetic & geographic distance matrices vibe coded with GPT

install.packages("vegan")
library(vegan)

###update workign directory and input filenames as necessary

setwd("~/Desktop/Research/Dudleya/Madrono Paper.260327/IBD/Mantel Tests and IBD Graphs with Stu's curves in R.260716/within QUIN")

gen <- read.csv("QUIN_genetic_matrix_aligned.csv", row.names = 1, check.names = FALSE)
geo_cm <- read.csv("QUIN_geographic_matrix_aligned.csv", row.names = 1, check.names = FALSE)

# Confirm alignment
stopifnot(identical(rownames(gen), rownames(geo_cm)))
stopifnot(identical(colnames(gen), colnames(geo_cm)))

# Remove extra members of field-detected clone groups
remove_ids <- c("QUIN_LP_43_", "QUIN_LP_23_", "QUIN_LP_28_", "QUIN_LP_55_")

gen_dc <- gen[!rownames(gen) %in% remove_ids, !colnames(gen) %in% remove_ids]
geo_dc_cm <- geo_cm[!rownames(geo_cm) %in% remove_ids, !colnames(geo_cm) %in% remove_ids]

# Convert cm to m
geo_dc_m <- geo_dc_cm / 100

# Mantel tests
set.seed(1)

mantel_linear <- mantel(
  as.dist(gen_dc),
  as.dist(geo_dc_m),
  method = "pearson",
  permutations = 9999
)

mantel_log <- mantel(
  as.dist(gen_dc),
  as.dist(log10(geo_dc_m + 1)),
  method = "pearson",
  permutations = 9999
)

mantel_linear
mantel_log

# Export decloned matrices
write.csv(gen_dc, "QUIN_genetic_matrix_decloned.csv")
write.csv(geo_dc_m, "QUIN_geographic_matrix_decloned_m.csv")
gen <- read.csv("QUIN_genetic_matrix_aligned.csv", row.names = 1, check.names = FALSE)
geo_cm <- read.csv("QUIN_geographic_matrix_aligned.csv", row.names = 1, check.names = FALSE)





###### Plot convention used below
# Dashed black line = ordinary linear regression on geographic distance.
# Solid black curve = regression on log10(geographic distance + 1),
# transformed back for display against the original geographic-distance axis.

######Making the X-Y IBD Plot for Within Pops (Quicksilver North)

# ---- Make XY points and plots from decloned matrices ----

library(ggplot2)

# Function to convert two matched distance matrices into long-format XY data
make_xy <- function(gen_mat, geo_mat) {
  stopifnot(identical(rownames(gen_mat), rownames(geo_mat)))
  stopifnot(identical(colnames(gen_mat), colnames(geo_mat)))
  
  pair_index <- which(upper.tri(gen_mat), arr.ind = TRUE)
  
  data.frame(
    sample1 = rownames(gen_mat)[pair_index[, 1]],
    sample2 = colnames(gen_mat)[pair_index[, 2]],
    genetic_distance = gen_mat[pair_index],
    geographic_distance_m = geo_mat[pair_index]
  )
}

xy_dc <- make_xy(gen_dc, geo_dc_m)

xy_dc$log10_geographic_distance_m_plus1 <- log10(xy_dc$geographic_distance_m + 1)

# Save XY points
write.csv(xy_dc, "QUIN_decloned_mantel_xy_points.csv", row.names = FALSE)


# ---- Linear geographic distance plot ----

p_linear <- ggplot(
  xy_dc,
  aes(x = geographic_distance_m, y = genetic_distance)
) +
#  geom_point(alpha = 0.65, size = 2.5) + #this creates gray filled points that turn black with overlap
  geom_point(
    shape = 1,
    color = "black",
    size = 2.6,
    stroke = 0.6
  ) +
  
  
  #  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) + #this makes a blue solid line, instead use below for black dashed line
  # Dashed line: ordinary linear fit on untransformed geographic distance
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    color = "black",
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  # Solid curve: linear model fitted to log10(distance + 1), plotted on the original distance axis
  geom_smooth(
    method = "lm",
    formula = y ~ log10(x + 1),
    se = FALSE,
    color = "black",
    linewidth = 1.2,
    linetype = "solid"
  ) +
#    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.02)) #this makes the y-axis go to zero (too much white space)

  
  theme_classic(base_size = 18) +
  labs(
    x = "Geographic distance (m)",
    y = "Genetic distance",
    title = "Within-population IBD: Quicksilver North, decloned"
  )

p_linear

ggsave(
  "QUIN_decloned_linear_mantel_scatterplot.png",
  p_linear,
  width = 8,
  height = 6,
  dpi = 300
)


# ---- Log10 geographic distance plot ---- UNUSED, but would need updating (hollow black circles, black hashed line)

p_log <- ggplot(
  xy_dc,
  aes(x = log10_geographic_distance_m_plus1, y = genetic_distance)
) +
  geom_point(alpha = 0.65, size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
  theme_classic(base_size = 18) +
  labs(
    x = "log10 Geographic distance (m + 1)",
    y = "Genetic distance",
    title = "Within-population IBD: Quicksilver North, decloned"
  )

p_log

ggsave(
  "QUIN_decloned_log_mantel_scatterplot.png",
  p_log,
  width = 8,
  height = 6,
  dpi = 300
)







###### Among-population Mantel tests and IBD graph

library(vegan)
library(ggplot2)

# ---- Read files ----

gen_multi <- read.csv(
  "multi-pop-DUSE-only.genetic distance matrix.csv",
  row.names = 1,
  check.names = FALSE
)

geo_multi_m <- read.csv(
  "multi-pop.geographic distance matrix.DUSE only.csv",
  row.names = 1,
  check.names = FALSE
)


# ---- Collapse sample-level genetic matrix to population-level matrix ----

get_pop <- function(x) {
  pop <- sub("_LP.*", "", x)
  pop[grepl("^QUIN", x)] <- "QUIN"
  pop
}

gen_multi_mat <- as.matrix(gen_multi)

row_pops <- get_pop(rownames(gen_multi_mat))
col_pops <- get_pop(colnames(gen_multi_mat))

pop_order <- rownames(geo_multi_m)

gen_pop <- matrix(
  NA,
  nrow = length(pop_order),
  ncol = length(pop_order),
  dimnames = list(pop_order, pop_order)
)

for (i in pop_order) {
  for (j in pop_order) {
    if (i == j) {
      gen_pop[i, j] <- 0
    } else {
      vals <- gen_multi_mat[row_pops == i, col_pops == j]
      gen_pop[i, j] <- mean(vals, na.rm = TRUE)
    }
  }
}

gen_pop <- as.data.frame(gen_pop)

# Check for problems
sum(is.na(gen_pop))
gen_pop["QUIN", ]
gen_pop[, "QUIN"]


# ---- Geographic matrix: convert meters to kilometers ----

geo_multi_km <- geo_multi_m / 1000

# Confirm alignment
stopifnot(identical(rownames(gen_pop), rownames(geo_multi_km)))
stopifnot(identical(colnames(gen_pop), colnames(geo_multi_km)))


# ---- Mantel tests ----

set.seed(1)

mantel_among_linear <- mantel(
  as.dist(gen_pop),
  as.dist(geo_multi_km),
  method = "pearson",
  permutations = 9999
)

mantel_among_log <- mantel(
  as.dist(gen_pop),
  as.dist(log10(geo_multi_km + 1)),
  method = "pearson",
  permutations = 9999
)

mantel_among_linear
mantel_among_log


# ---- Save final population-level matrices ----

write.csv(
  gen_pop,
  "DUSE_among_population_genetic_matrix_collapsed.csv"
)

write.csv(
  geo_multi_km,
  "DUSE_among_population_geographic_matrix_km.csv"
)


# ---- Make XY table ----

make_xy <- function(gen_mat, geo_mat) {
  gen_mat <- as.matrix(gen_mat)
  geo_mat <- as.matrix(geo_mat)
  
  stopifnot(identical(rownames(gen_mat), rownames(geo_mat)))
  stopifnot(identical(colnames(gen_mat), colnames(geo_mat)))
  
  pair_index <- which(upper.tri(gen_mat), arr.ind = TRUE)
  
  data.frame(
    pop1 = rownames(gen_mat)[pair_index[, 1]],
    pop2 = colnames(gen_mat)[pair_index[, 2]],
    genetic_distance = gen_mat[pair_index],
    geographic_distance_km = geo_mat[pair_index]
  )
}

xy_among <- make_xy(gen_pop, geo_multi_km)

xy_among$log10_geographic_distance_km_plus1 <-
  log10(xy_among$geographic_distance_km + 1)

# Keep only finite rows
xy_among_clean <- xy_among[
  is.finite(xy_among$genetic_distance) &
    is.finite(xy_among$geographic_distance_km) &
    is.finite(xy_among$log10_geographic_distance_km_plus1),
]

nrow(xy_among_clean)
# should be 78

write.csv(
  xy_among_clean,
  "DUSE_among_population_mantel_xy_points.csv",
  row.names = FALSE
)


# ---- Ordinary lm comparison: linear geographic distance vs log-transformed ----
# Note: this is exploratory only because pairwise distances are non-independent.

lm_among_linear <- lm(
  genetic_distance ~ geographic_distance_km,
  data = xy_among_clean
)

lm_among_log <- lm(
  genetic_distance ~ log10_geographic_distance_km_plus1,
  data = xy_among_clean
)

summary(lm_among_linear)
summary(lm_among_log)

AIC(lm_among_linear, lm_among_log)

# Optional: compare R-squared values directly
summary(lm_among_linear)$r.squared
summary(lm_among_log)$r.squared


# ---- Among-population linear IBD plot ----

p_among_linear <- ggplot(
  xy_among_clean,
  aes(x = geographic_distance_km, y = genetic_distance)
) +
  geom_point(
    shape = 1,
    color = "black",
    size = 2.6,
    stroke = 0.6
  ) +
  # Dashed line: ordinary linear fit on untransformed geographic distance
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    color = "black",
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  # Solid curve: linear model fitted to log10(distance + 1), plotted on the original distance axis
  geom_smooth(
    method = "lm",
    formula = y ~ log10(x + 1),
    se = FALSE,
    color = "black",
    linewidth = 1.2,
    linetype = "solid"
  ) +
  theme_classic(base_size = 18) +
  labs(
    x = "Geographic distance (km)",
    y = "Genetic distance",
    title = "Among-population IBD"
  )

p_among_linear

ggsave(
  "DUSE_among_population_linear_mantel_scatterplot.png",
  p_among_linear,
  width = 8,
  height = 6,
  dpi = 300
)



###### Exploratory east/west structure in among-population IBD

library(ggplot2)

# ---- Define western and eastern populations ----

western_pops <- c("STI", "RSV", "QUI1", "QUIN")

all_pops <- sort(unique(c(xy_among_clean$pop1, xy_among_clean$pop2)))
eastern_pops <- setdiff(all_pops, western_pops)

eastern_pops
western_pops


# ---- Classify each pairwise comparison ----

xy_among_clean$region1 <- ifelse(xy_among_clean$pop1 %in% western_pops, "West", "East")
xy_among_clean$region2 <- ifelse(xy_among_clean$pop2 %in% western_pops, "West", "East")

xy_among_clean$comparison_type <- ifelse(
  xy_among_clean$region1 == "East" & xy_among_clean$region2 == "East",
  "East-East",
  ifelse(
    xy_among_clean$region1 == "West" & xy_among_clean$region2 == "West",
    "West-West",
    "East-West"
  )
)

xy_among_clean$comparison_type <- factor(
  xy_among_clean$comparison_type,
  levels = c("East-East", "West-West", "East-West")
)

table(xy_among_clean$comparison_type)


# ---- Color-coded IBD graph ----

p_among_region <- ggplot(
  xy_among_clean,
  aes(x = geographic_distance_km, y = genetic_distance)
) +
  geom_point(
    aes(color = comparison_type),
    shape = 16,
    size = 3,
    alpha = 0.85
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "black",
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  theme_classic(base_size = 18) +
  labs(
    x = "Geographic distance (km)",
    y = "Genetic distance",
    color = "Comparison",
    title = "Among-population IBD by geographic group"
  )

p_among_region

ggsave(
  "DUSE_among_population_IBD_by_east_west_group.png",
  p_among_region,
  width = 8,
  height = 6,
  dpi = 300
)


# ---- Exploratory lm: does comparison type matter after distance? ----
# Pairwise comparisons are non-independent, so treat this as exploratory.

lm_distance_only <- lm(
  genetic_distance ~ geographic_distance_km,
  data = xy_among_clean
)

lm_distance_group <- lm(
  genetic_distance ~ geographic_distance_km + comparison_type,
  data = xy_among_clean
)

lm_distance_group_interaction <- lm(
  genetic_distance ~ geographic_distance_km * comparison_type,
  data = xy_among_clean
)

summary(lm_distance_only)
summary(lm_distance_group)
summary(lm_distance_group_interaction)

anova(lm_distance_only, lm_distance_group)
anova(lm_distance_group, lm_distance_group_interaction)

AIC(lm_distance_only, lm_distance_group, lm_distance_group_interaction)


# ---- Examine residuals from distance-only model ----
# Positive residual = more genetically differentiated than expected from distance alone.

xy_among_clean$resid_distance_only <- residuals(lm_distance_only)

aggregate(
  resid_distance_only ~ comparison_type,
  data = xy_among_clean,
  FUN = mean
)

aggregate(
  resid_distance_only ~ comparison_type,
  data = xy_among_clean,
  FUN = median
)

boxplot(
  resid_distance_only ~ comparison_type,
  data = xy_among_clean,
  ylab = "Residual genetic distance after accounting for geography",
  xlab = "Comparison type",
  main = "Residual genetic distance by comparison type"
)

###### Population-label permutation test
# Tests whether East-West comparisons have higher residual genetic distance
# than expected by random assignment of 4 populations to "West".

set.seed(1)

obs_stat <- mean(
  xy_among_clean$resid_distance_only[
    xy_among_clean$comparison_type == "East-West"
  ]
)

nperm <- 9999

perm_stats <- numeric(nperm)

for (k in 1:nperm) {
  
  perm_west <- sample(all_pops, length(western_pops), replace = FALSE)
  
  r1 <- ifelse(xy_among_clean$pop1 %in% perm_west, "West", "East")
  r2 <- ifelse(xy_among_clean$pop2 %in% perm_west, "West", "East")
  
  perm_type <- ifelse(
    r1 == "West" & r2 == "West",
    "West-West",
    ifelse(r1 == "East" & r2 == "East", "East-East", "East-West")
  )
  
  perm_stats[k] <- mean(
    xy_among_clean$resid_distance_only[perm_type == "East-West"]
  )
}

# One-sided p-value:
# Are observed East-West comparisons more genetically differentiated
# than expected after accounting for distance?

p_perm_one_sided <- (sum(perm_stats >= obs_stat) + 1) / (nperm + 1)

obs_stat
p_perm_one_sided

hist(
  perm_stats,
  breaks = 40,
  main = "Permutation test: East-West residual genetic distance",
  xlab = "Mean residual genetic distance for randomized East-West comparisons"
)

abline(v = obs_stat, col = "red", lwd = 3)


#to make labels appear on the IBD colored graph east-west
xy_among_clean$pair <-
  paste(xy_among_clean$pop1,
        xy_among_clean$pop2,
        sep = "-")

install.packages("ggrepel")   # first time only
library(ggrepel)

p_among_labels <- ggplot(
  xy_among_clean,
  aes(x = geographic_distance_km,
      y = genetic_distance)
) +
  
  geom_point(
    aes(color = comparison_type),
    size = 3
  ) +
  
  geom_text_repel(
    aes(label = pair,
        color = comparison_type),
    size = 2.8,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.15,
    segment.size = 0.2,
    show.legend = FALSE
  ) +
  
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "black",
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  
  theme_classic(base_size = 18)

p_among_labels