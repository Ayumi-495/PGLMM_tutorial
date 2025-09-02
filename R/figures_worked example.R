pacman::p_load("coda", "tidyverse", "here", 
               "MCMCglmm", "brms", "MASS", "patchwork", 
               "phytools", "patchwork", "bayesplot", "tidybayes")

select <- dplyr::select
##### ordinal data ####
mcmcglmm_mo3 <- readRDS(here("Rdata_tutorial", "ordinal", "m6.rds"))
brm_mo3 <- readRDS(here("Rdata_tutorial", "ordinal", "brm_m3.rds"))

summary(mcmcglmm_mo3)
summary(brm_mo3)

# mcmcglmm
## fixed effects
fixed_effects_samples <- as.data.frame(mcmcglmm_mo3$Sol) %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")

mcmcglmm_p1 <- ggplot(fixed_effects_samples, aes(y = parameter, x = value)) +
    stat_halfeye(
      normalize = "xy", 
      point_interval = "mean_qi", 
      fill = "lightcyan3", 
      color = "lightcyan4"
               ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
  scale_x_continuous(breaks = seq(-3.0, 3.0, 1), limits = c(-3.0, 4.0)) + 
  labs(
    title = "Posterior distributions of fixed effects - MCMCglmm",
    y = "Fixed effects"
    ) +
  theme_classic()
# 
# mcmcglmm_p1_1 <- ggplot(fixed_effects_samples, aes(y = parameter, x = value)) +
#   stat_halfeye(
#     normalize = "xy", 
#     point_interval = "median_qi", 
#     fill = "lightcyan3", 
#     color = "lightcyan4"
#   ) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
#   scale_x_continuous(breaks = seq(-5.0, 5.0, 2), limits = c(-5.0, 7.0)) + 
#   labs(
#     title = "Posterior distributions of fixed effects - MCMCglmm",
#     y = "Fixed effects"
#   ) +
#   theme_classic()

# samples <- as.matrix(mcmcglmm_mo3$Sol) 
# mcmc_areas_ridges(samples, prob = 0.95, prob_outer = 1) +
#   ggtitle("Posterior Distributions of Fixed Effects")

## random effects
random_samples <- as.data.frame(mcmcglmm_mo3$VCV[, "Phylo"]) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")

mcmcglmm_p2 <- ggplot(random_samples, aes(y = parameter, x = value)) +
    stat_halfeye(
      normalize = "xy",
      point_interval = "mean_qi", 
      fill = "olivedrab3", 
      color = "olivedrab4"
               ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
  scale_x_continuous(breaks = seq(0, 5.0, 2), limits = c(0, 7.0)) + 
  labs(
    title = "Posterior distributions of random effects- MCMCglmm",
    y = "Random effects"
  ) +
  theme_classic() 


# samples <- as.matrix(mcmcglmm_mo3$VCV) 
# mcmc_areas_ridges(samples, prob = 0.95, prob_outer = 1) +
#   ggtitle("Posterior Distributions of random effects")

## cutpoints
colnames(mcmcglmm_mo3$CP)
cutpoint_samples <- as.data.frame(mcmcglmm_mo3$CP) %>%
  pivot_longer(cols = everything(), names_to = "cutpoint", values_to = "value")

mcmcglmm_p3 <- ggplot(cutpoint_samples, aes(y = cutpoint, x = value)) +
    stat_halfeye(
      normalize = "xy", 
      point_interval = "mean_qi", 
      fill = "pink2", 
      color = "pink4"
               ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
  scale_x_continuous(breaks = seq(-3.5, 3.5, 1), limits = c(-3.5, 5.0)) + 
  labs(
    title = "Posterior distributions of cutpoints - MCMCglmm",
     y = "Cutpoints"
  ) +
  theme_classic()

# mcmcglmm_plot <- mcmcglmm_p1 + mcmcglmm_p2 + mcmcglmm_p3 + plot_layout(ncol = 1)
# mcmcglmm_plot

# brms
get_variables(brm_mo3)

## fixed effects
fixed_effects_samples_brms <- brm_mo3 %>%
  spread_draws(b_logMass, b_Habitat.Densityopen, b_Habitat.DensitysemiMopen)
fixed_effects_samples_brms <- fixed_effects_samples_brms %>%
  pivot_longer(cols = starts_with("b_"), 
               names_to = ".variable", 
               values_to = ".value")

brms_p1 <- ggplot(fixed_effects_samples_brms, aes(x = .value, y = .variable)) +
   stat_halfeye(
     normalize = "xy", 
     point_interval = "mean_qi", 
     fill = "lightcyan3", 
     color = "lightcyan4"
   ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
  scale_x_continuous(breaks = seq(-3.0, 3.0, 1), limits = c(-3.0, 4.0)) + 
  labs(title = "Posterior distributions of fixed effects - brms",
       y = "Fixed effects"
       ) +
  theme_classic()
head(fixed_effects_samples_brms)

## random effects
random_effects_samples_brms <- brm_mo3 %>%
  spread_draws(sd_Phylo__Intercept)
random_effects_samples_brms <- random_effects_samples_brms %>%
  pivot_longer(cols = starts_with("sd_"), 
               names_to = ".variable", 
               values_to = ".value")　%>% 
  mutate(.value = .value^2)

head(random_effects_samples_brms)

brms_p2 <- ggplot(random_effects_samples_brms, aes(x = .value, y = "sd(Intercept)")) +
   stat_halfeye(
     normalize = "xy",
     point_interval = "mean_qi", 
     fill = "olivedrab3", 
     color = "olivedrab4"
   ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
  scale_x_continuous(breaks = seq(0, 5.0, 2), limits = c(0, 7.0)) + 
  labs(
    title = "Posterior distributions of random effects - brms",
    y = "Random effects"
  ) +
  theme_classic() 


## cutpoint
cutpoint_samples_brms <- brm_mo3 %>%
  spread_draws(b_Intercept[cutpoint]) 

brms_p3 <- ggplot(cutpoint_samples_brms, aes(x = b_Intercept, y = factor(cutpoint)))+  
  stat_halfeye(
    normalize = "xy", 
    point_interval = "mean_qi", 
    fill = "pink2", 
    color = "pink4"
   ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
  scale_x_continuous(breaks = seq(-5.0, 5.0, 2), limits = c(-5.0, 7.0)) + 
  labs(
    title = "Posterior distributions of cutpoints - brms",
    y = "Cutpoints"
  ) +
  theme_classic()

# brms_plot <- brms_p1 + brms_p2 + brms_p3 + plot_layout(ncol = 1)
# brms_plot

# compare
p_ordinal <- mcmcglmm_p1 + mcmcglmm_p2 + mcmcglmm_p3 + brms_p1 + brms_p2 + brms_p3 


##### nominal data ####
mcmcglmm_mn3 <- readRDS(here("Rdata_tutorial", "multinomial", "mcmcglmm_m3_PL_prior3_1.rds"))
brm_mn3 <- readRDS(here("R", "worked_ex_multinominal_model", "AM", "brms_m_AM_piror8.rds"))
c2 <- (16 * sqrt(3) / (15 * pi))^2
c2c <- 1 + c2
summary(mcmcglmm_mn3$Sol/sqrt(c2c))
summary(mcmcglmm_mn3$VCV/c2c)
summary(brm_mn3)


# MCMCglmm
fixed_effects_samples <- as.data.frame(mcmcglmm_mn3$Sol/sqrt(c2c)) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")

mcmcglmm_p1 <- ggplot(fixed_effects_samples, aes(y = parameter, x = value)) +
  stat_halfeye(
    normalize = "xy", 
    point_interval = "mean_qi", 
    fill = "lightcyan3", 
    color = "lightcyan4"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
  scale_x_continuous(breaks = seq(-10.0, 10.0, 2), limits = c(-10.0, 10.0)) + 
  labs(
    title = "Posterior distributions of fixed effects - MCMCglmm",
    y = "Fixed effects"
  ) +
  theme_classic()
mcmcglmm_p1
# samples <- as.matrix(mcmcglmm_mo3$Sol) 
# mcmc_areas_ridges(samples, prob = 0.95, prob_outer = 1) +
#   ggtitle("Posterior Distributions of Fixed Effects")

## random effects
random_samples <- as.data.frame(mcmcglmm_mn3$VCV/c2c) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")

random_samples_filtered <- random_samples %>%
  filter(parameter %in% c("traitPrimary.Lifestyle.Insessorial:traitPrimary.Lifestyle.Insessorial.Phylo", 
                          "traitPrimary.Lifestyle.Terrestrial:traitPrimary.Lifestyle.Terrestrial.Phylo"))

mcmcglmm_p2 <- ggplot(random_samples_filtered, aes(y = parameter, x = value)) +
  stat_halfeye(
    normalize = "xy",
    point_interval = "mean_qi", 
    fill = "olivedrab3", 
    color = "olivedrab4"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
  scale_x_continuous(breaks = seq(0.0, 20.0, 2), limits = c(0.0, 20.0)) + 
  labs(
    title = "Posterior distributions of selected random effects - MCMCglmm",
    y = "Random effects"
  ) +
  theme_classic() 

mcmcglmm_p2
# samples <- as.matrix(mcmcglmm_mo3$VCV) 
# mcmc_areas_ridges(samples, prob = 0.95, prob_outer = 1) +
#   ggtitle("Posterior Distributions of random effects")

## correlation
cor_df <- as.data.frame(as.matrix(cor))
random_samples <- cor_df %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")

mcmcglmm_p3 <- ggplot(random_samples, aes(y = parameter, x = value)) +
  stat_halfeye(
    normalize = "xy",
    point_interval = "mean_qi", 
    fill = "#FF6347", 
    color = "#8B3626"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
  scale_x_continuous(breaks = seq(-1.0, 1.0, 0.2), limits = c(-1.0, 1.0)) + 
  labs(
    title = "Posterior distributions of selected random effects - MCMCglmm",
    y = "Correlation"
  ) +
  theme_classic() 

mcmcglmm_p3

##brms

get_variables(brm_mn3)

## fixed effects
fixed_effects_samples_brms <- brm_mn3 %>%
  spread_draws(b_muInsessorial_Intercept, 
               b_muTerrestrial_Intercept, 
               b_muInsessorial_log_Tail_Length_centered, 
               b_muInsessorial_IsOmnivore, 
               b_muTerrestrial_log_Tail_Length_centered, 
               b_muTerrestrial_IsOmnivore)

fixed_effects_samples_brms <- fixed_effects_samples_brms %>%
  pivot_longer(cols = starts_with("b_"), 
               names_to = ".variable", 
               values_to = ".value")

brms_p1 <- ggplot(fixed_effects_samples_brms, aes(x = .value, y = .variable)) +
  stat_halfeye(
    normalize = "xy", 
    point_interval = "mean_qi", 
    fill = "lightcyan3", 
    color = "lightcyan4"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
  scale_x_continuous(breaks = seq(-10.0, 10.0, 2), limits = c(-10.0, 10.0)) + 
  labs(title = "Posterior distributions of fixed effects - brms",
       y = "Fixed effects"
  ) +
  theme_classic()
brms_p1
head(fixed_effects_samples_brms)

## random effects
random_effects_samples_brms <- brm_mn3 %>%
  spread_draws(sd_Phylo__muInsessorial_Intercept, 
               sd_Phylo__muTerrestrial_Intercept
  )

random_effects_samples_brms_long <- random_effects_samples_brms %>%
  pivot_longer(cols = contains("phylo"), 
               names_to = ".variable", 
               values_to = ".value") %>% 
  mutate(.value = .value^2)

head(random_effects_samples_brms_long)

brms_p2 <- ggplot(random_effects_samples_brms_long, aes(x = .value, y = .variable)) +
  stat_halfeye(
    normalize = "xy",
    point_interval = "mean_qi", 
    fill = "olivedrab3", 
    color = "olivedrab4"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
  scale_x_continuous(breaks = seq(0.0, 80.0, 20), limits = c(0.0, 80.0)) + 
  labs(
    title = "Posterior distributions of random effects - brms",
    y = "Random effects"
  ) +
  theme_classic()

brms_p2

## correlation

correlation_effects_samples_brms <- brm_mn3 %>%
  spread_draws(cor_Phylo__muInsessorial_Intercept__muTerrestrial_Intercept)

correlation_effects_samples_brms_long <- correlation_effects_samples_brms %>%
  pivot_longer(cols = contains("phylo"), 
               names_to = ".variable", 
               values_to = ".value")

head(correlation_effects_samples_brms_long)

brms_p3 <- ggplot(correlation_effects_samples_brms_long, aes(x = .value, y = .variable)) +
  stat_halfeye(
    normalize = "xy",
    fill = "#FF6347", 
    color = "#8B3626"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
  scale_x_continuous(breaks = seq(-1.0, 1.0, 0.2), limits = c(-1.0, 1.0)) + 
  labs(
    title = "Posterior distributions of phylogenetic correlation - brms",
    y = "Correlation"
  ) +
  theme_classic()

brms_p3
# compare
p_nominal <- mcmcglmm_p1 + mcmcglmm_p2 + mcmcglmm_p3 + brms_p1 + brms_p2 + brms_p3 


#### additional - figure 5 ####
## ordered ----

dat <- read.csv(here("data", "bird body mass", "accipitridae_sampled.csv"))
dat <- dat %>% 
  mutate(across(c(Trophic.Level, Trophic.Niche, Primary.Lifestyle, Migration, Habitat, Species.Status), as.factor))

# rename: numbers to descriptive name
dat$Migration_ordered <- factor(dat$Migration,  levels = c("sedentary", "partially_migratory", "migratory"), ordered = TRUE)
dat <- dat %>%
  mutate(Habitat.Density = case_when(
    Habitat.Density == 1 ~ "dense",
    Habitat.Density == 2 ~ "semi-open",
    Habitat.Density == 3 ~ "open",
    TRUE ~ as.character(Habitat.Density)
  ))

table(dat$Habitat.Density)
dat$logMass <- log(dat$Mass)

brm_mo3 <- readRDS(here("Rdata_tutorial", "ordinal", "brm_m3.rds"))

# --- STEP 0: Load your original data frame 'dat' ---
# Ensure 'dat' is loaded and has 136 rows.
brm_mo3 <- readRDS(here("Rdata_tutorial", "ordinal", "brm_m3.rds"))


# --- STEP 1: Get raw posterior samples of the latent scale (linear predictor eta) ---
# posterior_linpred with summary=FALSE (default) gives the raw samples.
# The output is a matrix: MCMC samples (rows) x observations (columns)
latent_eta_raw <- posterior_linpred(brm_mo3)

# --- STEP 2: Calculate summary statistics (mean, lower, upper CIs) from raw samples ---
# Use apply to calculate statistics for each column (observation)
latent_eta_summary_calculated <- as.data.frame(t(apply(latent_eta_raw, 2, function(x) {
  c(Estimate = mean(x),
    Q2.5 = quantile(x, 0.025),
    Q97.5 = quantile(x, 0.975))
})))

# Rename the columns to match what you expect for combining
colnames(latent_eta_summary_calculated) <- c("Estimate", "Q2.5", "Q97.5")

dat_with_latent <- bind_cols(dat,
                             latent_scale_mean = latent_eta_summary_calculated$Estimate,
                             latent_scale_lower = latent_eta_summary_calculated$Q2.5,
                             latent_scale_upper = latent_eta_summary_calculated$Q97.5)

message("\nHead of dat_with_latent (after binding):")
print(head(dat_with_latent))

# --- STEP 4: Get thresholds and fixed effect coefficients ---
threshold0 <- summary(brm_mo3)$fixed["Intercept[1]", "Estimate"]
threshold1 <- summary(brm_mo3)$fixed["Intercept[2]", "Estimate"]
log_mass_coef <- summary(brm_mo3)$fixed["logMass", "Estimate"]

p1 <- ggplot(dat_with_latent, aes(x = logMass, y = latent_scale_mean)) +
  # Plot data points. Color them by the actual response category (Migration_ordered).
  # Setting alpha (transparency) helps visualize overlapping points.
  geom_point(aes(color = Migration_ordered), alpha = 0.6) +
  
  # Manually set colors corresponding to the Migration_ordered categories.
  # Adjust colors to match those in Figure 3's A, B, C.
  # Example: "Sedentary" = "darkgreen", "Partially migratory" = "darkblue", "Migratory" = "darkred"
  scale_color_manual(values = c("sedentary" = "#A2CD5A",
                                "partially_migratory" = "#00688B",
                                "migratory" = "#CD2626")) +
  
  # Add the regression line for the fixed effect.
  # geom_smooth(method = "lm") fits a linear model to the data and displays its prediction line.
  # se = TRUE also displays the standard error (confidence interval).
  geom_smooth(method = "lm", se = FALSE, color = "#8B8989", linetype = "solid", size = 1) +
  
  # Optionally, you can also display the credible interval of the latent scale variable as a ribbon.
  # geom_ribbon(aes(ymin = latent_scale_lower, ymax = latent_scale_upper), alpha = 0.2, fill = "lightblue") +
  
  # Add horizontal lines for the thresholds.
  # threshold0 and threshold1 are values obtained from the brm_mo3 summary.
  geom_hline(yintercept = threshold0, linetype = "dashed", color = "#CDC9C9", size = 0.8) +
  geom_hline(yintercept = threshold1, linetype = "dashed", color = "#CDC9C9", size = 0.8) +
  
  # Add labels for the thresholds.
  # Adjust x and y coordinates to fit within the plot's range.
  annotate("text", x = max(dat_with_latent$logMass), y = threshold0 + 0.1,
           label = paste0("Threshold0 (c0): ", round(threshold0, 2)),
           hjust = 1, vjust = 0, color = "#8B8989", size = 3) +
  annotate("text", x = max(dat_with_latent$logMass), y = threshold1 + 0.1,
           label = paste0("Threshold1 (c1): ", round(threshold1, 2)),
           hjust = 1, vjust = 0, color = "#8B8989", size = 3) +
  
  # Add labels indicating the range of response categories (refer to Figure 3a).
  # Adjust these Y-coordinates based on the min/max of latent_scale_mean and the thresholds.
  # X-coordinates should be placed towards the left side of the graph (e.g., min(logMass) + adjustment).
  # Labels should correspond to the categories A, B, C in the figure.
  annotate("text", x = min(dat_with_latent$logMass) + (max(dat_with_latent$logMass) - min(dat_with_latent$logMass)) * 0.05,
           y = (threshold0 + min(dat_with_latent$latent_scale_mean)) / 2, # Center of the region below threshold0
           label = "Sedentary", color = "#A2CD5A", size = 5, fontface = "bold", hjust = 0) +
  annotate("text", x = min(dat_with_latent$logMass) + (max(dat_with_latent$logMass) - min(dat_with_latent$logMass)) * 0.05,
           y = (threshold1 + threshold0) / 2, # Center of the region between threshold0 and threshold1
           label = "Partially\nmigratory", color = "#00688B", size = 5, fontface = "bold", hjust = 0) +
  annotate("text", x = min(dat_with_latent$logMass) + (max(dat_with_latent$logMass) - min(dat_with_latent$logMass)) * 0.05,
           y = (max(dat_with_latent$latent_scale_mean) + threshold1) / 2, # Center of the region above threshold1
           label = "Migratory", color = "#CD2626", size = 5, fontface = "bold", hjust = 0) +
  
  # Set axis labels and title.
  labs(
    title = "Latent variable vs. Log(body mass)",
    x = "Log(body mass)",
    y = "Latent scale variable"
  ) +
  # Set a theme (theme_minimal provides a clean look).
  theme_classic() +
  # Set legend position.
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"), # Center and bold title
        axis.title = element_text(face = "bold"), # Bold axis labels
        legend.title = element_blank()) # Hide legend title (category names serve as title)

# To save the plot
# ggsave("logMass_latent_plot.png", width = 8, height = 6, dpi = 300)


# --- STEP 5: Plotting the latent variable by Habitat Density ---
# Check the categories of Habitat.Density and set their order if necessary.
# Check the actual levels of dat_with_latent$Habitat.Density and set them in the desired order.
# For example, if you want the order 'Dense', 'Semi-open', 'Open':
# dat_with_latent$Habitat.Density <- factor(dat_with_latent$Habitat.Density,
#                                            levels = c("Dense", "Semi-open", "Open"))

p2 <- ggplot(dat_with_latent, aes(x = Habitat.Density, y = latent_scale_mean)) +
  # Display the data distribution for each category using violin plots.
  # Use fill to color by category and alpha for transparency.
  geom_violin(aes(fill = Habitat.Density), alpha = 0.7, trim = FALSE) +
  
  # Overlay individual data points, colored by the actual response category (Migration_ordered).
  # position_jitter slightly scatters the points to reduce overlap.
  geom_point(position = position_jitter(width = 0.1), alpha = 0.4, aes(color = Migration_ordered)) +
  
  # Manually set colors for Migration_ordered categories (same as above).
  scale_color_manual(values = c("sedentary" = "#A2CD5A",
                                "partially_migratory" = "#00688B",
                                "migratory" = "#CD2626")) +
  
  # Adjust fill colors for the violin plots (optional)
  scale_fill_manual(values = c("dense" = "#EEE9E9", "semi-open" = "#EEE9E9", "open" = "#EEE9E9")) +
  
  # Add horizontal lines for the thresholds.
  geom_hline(yintercept = threshold0, linetype = "dashed", color = "#CDC9C9", size = 0.8) +
  geom_hline(yintercept = threshold1, linetype = "dashed", color = "#CDC9C9", size = 0.8) +
  
  # Add labels for the thresholds.
  # Adjust x-coordinate based on the number of categories. Here, it's at the position of the rightmost category.
  annotate("text", x = max(as.numeric(factor(dat_with_latent$Habitat.Density))), y = threshold0 + 0.1,
           label = paste0("Threshold0 (c0): ", round(threshold0, 2)),
           hjust = 0.3, vjust = 0, color = "#8B8989", size = 3) +
  annotate("text", x = max(as.numeric(factor(dat_with_latent$Habitat.Density))), y = threshold1 + 0.1,
           label = paste0("Threshold1 (c1): ", round(threshold1, 2)),
           hjust = 0.3, vjust = 0, color = "#8B8989", size = 3) +
  
  # Add labels indicating the range of response categories (refer to Figure 3b).
  # X-coordinates should be placed towards the left side of the plot space.
  annotate("text", x = 0.7, # Near the left edge of the plot (position 0.5 for categories)
           y = (threshold0 + min(dat_with_latent$latent_scale_mean)) / 2,
           label = "Sedentary", color = "#A2CD5A", size = 4, fontface = "bold") +
  annotate("text", x = 0.7,
           y = (threshold1 + threshold0) / 2,
           label = "Partially\nmigratory", color = "#00688B", size = 4, fontface = "bold") +
  annotate("text", x = 0.7,
           y = (max(dat_with_latent$latent_scale_mean) + threshold1) / 2,
           label = "Migratory", color = "#CD2626", size = 4, fontface = "bold") +
  
  # Set axis labels and title.
  labs(
    title = "Latent variable by Habitat Density",
    x = "Habitat Density",
    y = "Latent scale variable"
  ) +
  theme_classic() +
  # For this plot, the legend is typically not shown as fill indicates categories.
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"))

# To save the plot
# ggsave("habitat_density_latent_plot.png", width = 7, height = 6, dpi = 300)

p1/p2


# nominal ----
# ----------------------------------------------------------------------------
# 3. Effect of IsOmnivore (Insessorial vs Generalist)
# ----------------------------------------------------------------------------
# brm_mn3 <- readRDS(here("R", "worked_ex_multinominal_model", "AM", "brms_m_AM_piror8.rds"))
# fe <- brms::fixef(brm_mn3)[, "Estimate"]
# 
# muInsessorial_Intercept                 <- unname(fe["muInsessorial_Intercept"])
# muInsessorial_log_Tail_Length_centered  <- unname(fe["muInsessorial_log_Tail_Length_centered"])
# muInsessorial_IsOmnivore_effect         <- unname(fe["muInsessorial_IsOmnivore"])
# 
# muTerrestrial_Intercept                 <- unname(fe["muTerrestrial_Intercept"])
# muTerrestrial_log_Tail_Length_centered  <- unname(fe["muTerrestrial_log_Tail_Length_centered"])
# muTerrestrial_IsOmnivore_effect         <- unname(fe["muTerrestrial_IsOmnivore"])
# 
# plot_data_omn_ins <- dat %>%
#   filter(Primary.Lifestyle %in% c("Insessorial", "Generalist"))
# 
# plot_data_omn_ins_predicted <- plot_data_omn_ins %>%
#   mutate(
#     predicted_logit_odds = case_when(
#       Primary.Lifestyle == "Insessorial" & IsOmnivore_factor == "No" ~ muInsessorial_Intercept + muInsessorial_log_Tail_Length_centered * 0,
#       Primary.Lifestyle == "Insessorial" & IsOmnivore_factor == "Yes" ~ muInsessorial_Intercept + muInsessorial_log_Tail_Length_centered * 0 + muInsessorial_IsOmnivore_effect,
#       Primary.Lifestyle == "Generalist" ~ 0,
#       TRUE ~ NA_real_
#     )
#   )
# 
# effect_omn_ins <- data.frame(
#   IsOmnivore_factor = factor(c("No", "Yes"), levels = c("No", "Yes")),
#   logit_odds = c(muInsessorial_Intercept, muInsessorial_Intercept + muInsessorial_IsOmnivore_effect)
# )
# 
# p_omni_ins <- ggplot(plot_data_omn_ins_predicted, aes(x = IsOmnivore_factor, y = predicted_logit_odds)) +
#   geom_violin(aes(fill = Primary.Lifestyle), alpha = 0.3, trim = FALSE) +
#   geom_jitter(
#   aes(color = Primary.Lifestyle),
#   position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9, jitter.height = 0.5),
#   alpha = 0.6,
#   size = 1.5
# ) +
#   # geom_point(data = effect_omn_ins, aes(x = IsOmnivore_factor, y = logit_odds),
#   #            size = 5, color = "black", shape = 18) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
#   scale_color_manual(values = c("Generalist" = "gray50", "Insessorial" = "#d95f02")) +
#   scale_fill_manual(values = c("Generalist" = "gray50", "Insessorial" = "#d95f02")) +
#   labs(
#     # title = "Insessorial vs. Generalist",
#     subtitle = "Effect of omnivory",
#     x = "Is omnivore?",
#     y = "Logit(P(Insessorial) / P(Generalist))"
#   ) +
#   theme_classic() +
#   theme(legend.position = "none") +
#   ylim(-3, 3)
# 
# # ----------------------------------------------------------------------------
# # 4. Effect of IsOmnivore (Terrestrial vs Generalist)
# # ----------------------------------------------------------------------------
# plot_data_omn_ter <- dat %>%
#   filter(Primary.Lifestyle %in% c("Terrestrial", "Generalist"))
# 
# plot_data_omn_ter_predicted <- plot_data_omn_ter %>%
#   mutate(
#     predicted_logit_odds = case_when(
#       Primary.Lifestyle == "Terrestrial" & IsOmnivore_factor == "No" ~ muTerrestrial_Intercept + muTerrestrial_log_Tail_Length_centered * 0,
#       Primary.Lifestyle == "Terrestrial" & IsOmnivore_factor == "Yes" ~ muTerrestrial_Intercept + muTerrestrial_log_Tail_Length_centered * 0 + muTerrestrial_IsOmnivore_effect,
#       Primary.Lifestyle == "Generalist" ~ 0,
#       TRUE ~ NA_real_
#     )
#   )
# 
# effect_omn_ter <- data.frame(
#   IsOmnivore_factor = factor(c("No", "Yes"), levels = c("No", "Yes")),
#   logit_odds = c(muTerrestrial_Intercept, muTerrestrial_Intercept + muTerrestrial_IsOmnivore_effect)
# )
# 
# p_omni_ter <- ggplot(plot_data_omn_ter_predicted, aes(x = IsOmnivore_factor, y = predicted_logit_odds)) +
#   geom_violin(aes(fill = Primary.Lifestyle), alpha = 0.3, trim = FALSE) +
#   geom_jitter(
#   aes(color = Primary.Lifestyle),
#   position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9, jitter.height = 0.5),
#   alpha = 0.6,
#   size = 1.5
# ) +
#   # geom_point(data = effect_omn_ter, aes(x = IsOmnivore_factor, y = logit_odds),
#   #            size = 5, color = "black", shape = 18) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
#   scale_color_manual(values = c("Generalist" = "gray50", "Terrestrial" = "#1b9e77")) +
#   scale_fill_manual(values = c("Generalist" = "gray50", "Terrestrial" = "#1b9e77")) +
#   labs(
#     # title = "Terrestrial vs. Generalist",
#     # subtitle = "Effect of omnivory",
#     x = "Is Omnivore?",
#     y = "Logit(P(Terrestrial) / P(Generalist))"
#   ) +
#   theme_classic() +
#   theme(legend.position = "none") +
#   ylim(-3, 3)
# 
# # Combine all plots
# (p_tail_ins | p_tail_ter) / (p_omni_ins | p_omni_ter)


# nominal model with brms ----

brm_mn3 <- readRDS(here("R", "worked_ex_multinominal_model", "AM", "brms_m_AM_piror8.rds"))
dat <- read.csv(here("data", "potential", "avonet", "turdidae.csv"))

# Centring for continuous variables
dat <- dat %>%
  mutate(
    log_Mass_centered = scale(log(Mass), center = TRUE, scale = FALSE),
    log_Tail_Length_centered = scale(log(Tail.Length), center = TRUE, scale = FALSE)
  )

dat <- dat %>% 
  mutate(across(c(Trophic.Level, Trophic.Niche, Primary.Lifestyle, Migration, Habitat, Species.Status), as.factor))

dat$IsOmnivore <- ifelse(dat$Trophic.Level == "Omnivore", 1, 0)

# Convert IsOmnivore into a factor with labels
dat <- dat %>%
  mutate(
    IsOmnivore_factor = factor(IsOmnivore, levels = c(0, 1), labels = c("No", "Yes"))
  )


## Visualising brms prediction with random effect and marginalising Omnivore  ----

# Sequence of tail lengths to predict over
tail_seq <- seq(min(dat$log_Tail_Length_centered),
                max(dat$log_Tail_Length_centered),
                length.out = 100)
phylo_levels <- unique(dat$Phylo)

newdata0 <- dat %>%
  mutate(IsOmnivore = 0) %>%
  dplyr::select(log_Tail_Length_centered, IsOmnivore, Phylo)

newdata1 <- dat %>%
  mutate(IsOmnivore = 1) %>%
  dplyr::select(log_Tail_Length_centered, IsOmnivore, Phylo)

# Get expected predictions (posterior_epred includes random effects)
epred0 <- posterior_epred(brm_mn3, newdata = newdata0, re_formula = NULL)
epred1 <- posterior_epred(brm_mn3, newdata = newdata1, re_formula = NULL)

# Check category order
cat_levels <- levels(dat$Primary.Lifestyle)  # should be "Generalist", "Insessorial", "Terrestrial"
print(cat_levels)

# Calculate weights for marginalisation
w0 <- mean(dat$IsOmnivore == 0)
w1 <- mean(dat$IsOmnivore == 1)

# Compute marginalised class probabilities
pG_marg <- w0 * epred0[, , "Generalist"]   + w1 * epred1[, , "Generalist"]
pI_marg <- w0 * epred0[, , "Insessorial"]  + w1 * epred1[, , "Insessorial"]
pT_marg <- w0 * epred0[, , "Terrestrial"]  + w1 * epred1[, , "Terrestrial"]

# logit
logitI <- log(pI_marg / pG_marg)
logitT <- log(pT_marg / pG_marg)

logitI_mean <- apply(logitI, 2, mean)
logitT_mean <- apply(logitT, 2, mean)

# Summary (mean and 95% CI)
plot_points <- dat %>%
  mutate(
    logitI = logitI_mean,
    logitT = logitT_mean
  )

# Insessorial vs Generalist
p_ins <- ggplot(plot_points %>% filter(Primary.Lifestyle %in% c("Insessorial", "Generalist")),
                aes(x = log_Tail_Length_centered, y = logitI, color = Primary.Lifestyle)) +
  geom_jitter(width = 0.02, height = 0.05, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "#ec6a07") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Insessorial" = "#d95f02", "Generalist" = "gray50")) +
  labs(title = "Insessorial vs. Generalist (marginalised Omnivore)",
       x = "Tail length (centered-log)",
       y = "logit(P(Insessorial) / P(Generalist))") +
       ylim(-3, 3) +
       xlim(-0.5, 0.5) +
       theme(legend.position = "bottom") +
  theme_classic()

# Terrestrial vs Generalist
p_ter <- ggplot(plot_points %>% filter(Primary.Lifestyle %in% c("Terrestrial", "Generalist")),
                aes(x = log_Tail_Length_centered, y = logitT, color = Primary.Lifestyle)) +
  geom_jitter(width = 0.02, height = 0.05, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "#13600f") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Terrestrial" = "#1b9e77", "Generalist" = "gray50")) +
  labs(title = "Terrestrial vs. Generalist (marginalised Omnivore)",
       x = "Tail length (centered-log)",
       y = "logit(P(Terrestrial) / P(Generalist))") +
       ylim(-3, 3) +
       xlim(-0.5, 0.5) +
       theme(legend.position = "bottom") +
  theme_classic()

# Display plots side-by-side
p_ins + p_ter +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
# =============================================
# Visualising the effect of Omnivore (IsOmnivore)
# using brms model with random effect (Phylo)
# =============================================

# Fix Tail length at the centred mean (i.e., 0)
tail_fixed <- 0

# Create new data for Omnivore = 0 and Omnivore = 1
newdata0 <- data.frame(
  log_Tail_Length_centered = tail_fixed,
  IsOmnivore = 0,
  Phylo = dat$Phylo
)

newdata1 <- data.frame(
  log_Tail_Length_centered = tail_fixed,
  IsOmnivore = 1,
  Phylo = dat$Phylo
)
posterior_predict(brm_mn3, newdata = newdata0, re_formula = NULL)
# Generate expected predictions including random effects
epred0 <- posterior_epred(brm_mn3, newdata = newdata0, re_formula = NULL)
epred1 <- posterior_epred(brm_mn3, newdata = newdata1, re_formula = NULL)

# Calculate logit values
logitI_omn0 <- log(epred0[, , "Insessorial"] / epred0[, , "Generalist"])
logitI_omn1 <- log(epred1[, , "Insessorial"] / epred1[, , "Generalist"])

logitT_omn0 <- log(epred0[, , "Terrestrial"] / epred0[, , "Generalist"])
logitT_omn1 <- log(epred1[, , "Terrestrial"] / epred1[, , "Generalist"])

# Compute mean logits per observation
plot_marg <- data.frame(
  Phylo = dat$Phylo,
  logitI_omn0 = apply(logitI_omn0, 2, mean),
  logitI_omn1 = apply(logitI_omn1, 2, mean),
  logitT_omn0 = apply(logitT_omn0, 2, mean),
  logitT_omn1 = apply(logitT_omn1, 2, mean)
) %>%
  pivot_longer(cols = -Phylo, names_to = "Condition", values_to = "logit") %>%
  mutate(
    Lifestyle = ifelse(str_detect(Condition, "^logitI"), "Insessorial", "Terrestrial"),
    IsOmnivore_factor = ifelse(str_detect(Condition, "omn1"), "Yes", "No")
  )

# Plot: Insessorial vs Generalist
plot_points_ins <- plot_points %>%
  filter(Primary.Lifestyle %in% c("Insessorial", "Generalist"))

plot_points_ter <- plot_points %>%
  filter(Primary.Lifestyle %in% c("Terrestrial", "Generalist"))
p_ins_om <- ggplot(plot_points_ins, aes(x = IsOmnivore_factor, y = logitI)) +
  geom_violin(aes(fill = IsOmnivore_factor), trim = FALSE, alpha = 0.4) +
　geom_jitter(aes(color = Primary.Lifestyle),
  alpha = 0.6, size = 2,
  position = position_jitter(width = 0.1, seed = 42)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("No" = "#d8d7fb", "Yes" = "#f8b6b4")) +
  scale_color_manual(values = c(
    "Insessorial" = "#d95f02",
    "Generalist" = "#2c2e2b"
  )) +
  labs(title = "Omnivore (mean tail length)",
       x = "Is Omnivore?", y = "logit(P(Insessorial) / P(Generalist)") +
  ylim(-4, 4) +
  theme(legend.position = "bottom") +
  theme_classic()

# Plot: Terrestrial vs Generalist
p_ter_om <- ggplot(plot_points_ter, aes(x = IsOmnivore_factor, y = logitT)) +
  geom_violin(aes(fill = IsOmnivore_factor), trim = FALSE, alpha = 0.4) +
  geom_jitter(aes(color = Primary.Lifestyle),
  alpha = 0.6, size = 2,
  position = position_jitter(width = 0.1, seed = 42)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("No" = "#d8d7fb", "Yes" = "#f8b6b4")) +
  scale_color_manual(values = c(
    "Terrestrial" = "#1b9e77",
    "Generalist" = "#2c2e2b"
  )) +
  labs(title = "Omnivore (mean tail length)",
       x = "Is Omnivore?", y = "logit(P(Terrestrial) / P(Generalist)") +
 ylim(-4, 4) +
 theme(legend.position = "bottom") +
theme_classic()

# Combine plots side by side
p_ins_om + p_ter_om +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
