#-------------------------------------------------------------------------------
# PACKAGES
#-------------------------------------------------------------------------------

library(tidyverse)
library(patchwork)
library(ggpattern)
library(ggtext)
library(rstan)
library(tidybayes)
library(shinystan)

#-------------------------------------------------------------------------------
# DATA READING
#-------------------------------------------------------------------------------

data <- read.csv("https://raw.githubusercontent.com/MFajgenblat/hgp-chlorpyrifos/main/HGP_Chlorpyrifos_Data.csv") %>%
  drop_na()
# The drop_na function removes 8 data points
# These correspond to one day where no zooplankton biomass was available for two mesocosms

#-------------------------------------------------------------------------------
# METADATA FOR PLOTTING PURPOSES
#-------------------------------------------------------------------------------

species_metadata <- data.frame(speciesid = 1:6, 
                               speciesname = factor(c("Daphnia magna", "Daphnia galeata", "Daphnia pulicaria", "Scapholeberis mucronata", "Chlorophyll-a", "Phycocyanin"),
                                                    levels = c("Daphnia magna", "Daphnia pulicaria", "Daphnia galeata", "Scapholeberis mucronata", "Chlorophyll-a", "Phycocyanin")),
                               speciesname_br = factor(c("Daphnia\nmagna", "Daphnia\ngaleata", "Daphnia\npulicaria", "Scapholeberis\nmucronata", "Chlorophyll-a", "Phycocyanin"),
                                                       levels = c("Daphnia\nmagna", "Daphnia\npulicaria", "Daphnia\ngaleata", "Scapholeberis\nmucronata", "Chlorophyll-a", "Phycocyanin")),
                               speciesformatted = factor(c("<i>Daphnia magna</i>", "<i>Daphnia galeata</i>", "<i>Daphnia pulicaria</i>", "<i>Scapholeberis mucronata</i>", "Chlorophyll-a", "Phycocyanin"),
                                                         levels = c("<i>Daphnia magna</i>", "<i>Daphnia pulicaria</i>", "<i>Daphnia galeata</i>", "<i>Scapholeberis mucronata</i>", "Chlorophyll-a", "Phycocyanin")))

#-------------------------------------------------------------------------------
# EXPLORATORY DATA VISUALISATION
#-------------------------------------------------------------------------------

data %>%
  mutate(Mix = paste0("Mix ",as.numeric(substr(Replicate, 6, 6))),
         Treatment = c("Control", "Pesticide")[as.numeric(factor(substr(Replicate, 1, 3), levels = c("CTL", "CPF")))],
         Species = factor(Species, c("Daphnia magna", "Daphnia galeata", "Daphnia pulicaria", "Scapholeberis mucronata", "Chlorophyll-a", "Phycocyanin"))) %>%
  left_join(species_metadata, by = c("Species" = "speciesname")) %>%
  ggplot() +
  geom_line(aes(x = Day, y = Biomass, color = Treatment, group = Replicate), alpha = 0.75, linewidth = 0.15) +
  geom_vline(data = data.frame(x = c(10, 24, 38), type = "dashed"), aes(linetype = type, xintercept = x), color = "grey50", size = 0.3) +
  scale_color_manual("Treatment", values = c("#0090C8", "#CE178C"), guide = guide_legend(override.aes = list(alpha = 1, size = 1))) +
  scale_fill_manual("Treatment", values = c("#0090C8", "#CE178C"), guide = "none") +
  scale_alpha_manual("", values = 0.1, labels = "> 95% probability\nof a treatment\neffect") +
  facet_wrap(~ speciesformatted, scales = "free", ncol = 2) +
  scale_x_continuous("Time (days)", breaks = c(0, seq(10, 51, by = 7)), labels = c("<span style = 'font-size:6pt'>Inoculation<br>moment</span>", seq(10, 51, by = 7)), expand = c(0,0)) +
  scale_y_log10("Estimated biomass (µg/L)", labels = scales::label_number()) +
  scale_linetype_manual("", breaks = "dashed", values = "dashed", labels = "Pesticide pulse") +
  coord_cartesian(xlim = c(0,51)) +
  theme(panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "grey95", size = 0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        strip.text.x = element_markdown(size = 9, face = "bold"),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black", size = 0.3),
        axis.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_markdown(size = 9),
        axis.text.y = element_text(size = 9),
        axis.ticks = element_line(size = 0.3),
        legend.position = "right",
        legend.key = element_blank(),
        legend.key.width = unit(5, "mm"),
        legend.key.height = unit(5, "mm"),
        legend.spacing.y = unit(0.5, "mm"),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 9))
ggsave("HGP_Chlorpyrifos_Figure_0.pdf", width = 17.7, height = 10, units = "cm")

#-------------------------------------------------------------------------------
# FITTING THE HIERARCHICAL GAUSSIAN PROCESS MODEL
#-------------------------------------------------------------------------------

options(mc.cores = parallel::detectCores())
fit_HGP <- stan("HGP_Chlorpyrifos_Model.stan",
                data = list(N = nrow(data),
                            N_observed = sum(!is.na(data$Biomass)),
                            N_missing = sum(is.na(data$Biomass)),
                            N_times = diff(range(data$Day)) + 1,
                            N_times_observed = length(unique(data$Day)),
                            N_times_common = 10,
                            N_mix = length(unique(substr(data$Replicate, 6, 6))),
                            N_rep = length(unique(data$Replicate)),
                            N_treat = length(unique(substr(data$Replicate, 1, 3))),
                            N_species = length(unique(data$Species)),
                            y_obs = log(1+data$Biomass[!is.na(data$Biomass)]),
                            observed = which(!is.na(data$Biomass)),
                            missing = which(is.na(data$Biomass)),
                            species = as.numeric(factor(data$Species, levels = c("Daphnia magna", "Daphnia galeata", "Daphnia pulicaria", "Scapholeberis mucronata", "Chlorophyll-a", "Phycocyanin"))),
                            treat = as.numeric(factor(substr(data$Replicate, 1, 3), levels = c("CTL", "CPF"))),
                            mix = as.numeric(substr(data$Replicate, 6, 6)),
                            rep = as.numeric(factor(data$Replicate)),
                            time = data$Day,
                            time_observed = as.numeric(factor(data$Day)),
                            times_scaled = seq(0, 1, length.out = diff(range(data$Day)) + 1),
                            times_observed = sort(unique(data$Day))),
                iter = 5000, chains = 4)
saveRDS(fit_HGP, "HGP_Chlorpyrifos_Fit.rds")

# Or alternatively, retrieve an existing fit_HGP object
#fit_HGP <- readRDS("HGP_Chlorpyrifos_Fit.rds")

#-------------------------------------------------------------------------------
# CHECK CONVERGENCE AND MODEL DIAGNOSTICS
#-------------------------------------------------------------------------------

launch_shinystan(fit_HGP)

#-------------------------------------------------------------------------------
# PLOTTING THE HIERARCHICAL GAUSSIAN PROCESS FITS
#-------------------------------------------------------------------------------

fit_HGP %>%
  spread_draws(f_main[speciesid,Treatment,x]) %>%
  filter(!(speciesid > 4 & x < 11)) %>%
  left_join(species_metadata) %>%
  mutate(f_main = exp(f_main)-1,
         Treatment = factor(case_when(Treatment == 1 ~ "Control",
                                      Treatment == 2 ~ "Pesticide"),
                            levels = c("Control", "Pesticide")),
         x = x - 1) -> plotdata

importance <- plotdata %>%
  pivot_wider(names_from = Treatment, values_from = f_main) %>%
  mutate(difference = Control - Pesticide) %>%
  group_by(speciesformatted, x) %>%
  summarise(importance = mean(difference > 0)) %>%
  filter(importance > 0.95 | importance < 0.05)

ggplot(plotdata) +
  geom_rect(data = importance, aes(xmin = x - 0.5, xmax = x + 0.5, ymin = -Inf, ymax = +Inf, alpha = "0.1"), fill = "black") +
  stat_lineribbon(aes(x = x, y = f_main, color = Treatment, fill = Treatment), .width = c(0.5, 0.8, 0.95, 0.99), alpha = 1/5, size = 0.5) +
  #geom_line(data = data, aes(x = Week, y = density, color = Treatment, group = Replicate), alpha = 0.6, size = 0.12) +
  geom_vline(data = data.frame(x = c(10, 24, 38), type = "dashed"), aes(linetype = type, xintercept = x), color = "grey50", size = 0.3) +
  scale_color_manual("Treatment", values = c("#0090C8", "#CE178C"), guide = guide_legend(override.aes = list(alpha = 1, size = 1, fill = NA))) +
  scale_fill_manual("Treatment", values = c("#0090C8", "#CE178C"), guide = "none") +
  scale_alpha_manual("", values = 0.15, labels = "> 95% probability\nof a treatment\neffect") +
  facet_wrap(~ speciesformatted, scales = "free", ncol = 2) +
  scale_x_continuous("Time (days)", breaks = c(0, seq(10, 51, by = 7)), labels = c("<span style = 'font-size:6pt'>Inoculation<br>moment</span>", seq(10, 51, by = 7)), expand = c(0,0)) +
  scale_y_continuous("Estimated biomass (µg/L)") +
  scale_linetype_manual("", breaks = "dashed", values = "dashed", labels = "Pesticide pulse") +
  coord_cartesian(xlim = c(0,60)) +
  theme(panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "grey95", size = 0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        strip.text.x = element_markdown(size = 9, face = "bold"),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black", size = 0.3),
        axis.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_markdown(size = 9),
        axis.text.y = element_text(size = 9),
        axis.ticks = element_line(size = 0.3),
        legend.position = "right",
        legend.key = element_blank(),
        legend.key.width = unit(5, "mm"),
        legend.key.height = unit(5, "mm"),
        legend.spacing.y = unit(0.5, "mm"),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 9))
ggsave("HGP_Chlorpyrifos_Figure_2.pdf", width = 17.7, height = 10, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# VARIATION PARTITIONING
#-------------------------------------------------------------------------------

# Bar plot to include in Figure 2
fit_HGP %>%
  spread_draws(alpha_main[speciesid,Treatment], alpha_mix[speciesid,Treatment], alpha_rep[speciesid,Treatment]) %>%
  pivot_longer(c(alpha_main, alpha_mix, alpha_rep)) %>%
  left_join(species_metadata) %>%
  mutate(Treatment = factor(case_when(Treatment == 1 ~ "Control",
                                      Treatment == 2 ~ "Pesticide"),
                            levels = c("Control", "Pesticide")),
         name = factor(case_when(name == "alpha_main" ~ "Main",
                                 name == "alpha_mix" ~ "Mix",
                                 name == "alpha_rep" ~ "Replicate"),
                       levels = rev(c("Main", "Mix", "Replicate")))) %>%
  group_by(speciesformatted, Treatment, name) %>%
  summarize(value = mean(value)) %>%
  ggplot() +
  geom_bar_pattern(aes(x = Treatment, y = value, fill = Treatment, pattern_angle = name, pattern_fill = Treatment), pattern_key_scale_factor = 0.25, pattern_spacing = 0.075, pattern_density = 0.5, pattern_color = NA, color = "white", width = 0.8, stat = "identity", position = "fill") +
  scale_fill_manual("Treatment", values = c("#0090C8", "#CE178C"), guide = "none") +
  scale_pattern_fill_manual("Treatment", values = c("#004f6e", "#700449"), guide = "none") +
  scale_pattern_fill2_manual("Treatment", values = c("#0090C8", "#CE178C"), guide = "none") +
  scale_pattern_angle_manual("name", values = c(90, 45, 0)) +
  facet_wrap(~ speciesformatted, ncol = 2) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.text.x = element_markdown(size = 9, face = "bold"),
        strip.background = element_blank())
ggsave("HGP_Chlorpyrifos_Figure_2_VP.pdf", width = 6, height = 10, units = "cm", dpi = 600)

# Dedicated figure
fit_HGP %>%
  spread_draws(alpha_main[speciesid,Treatment], alpha_mix[speciesid,Treatment], alpha_rep[speciesid,Treatment]) %>%
  pivot_longer(c(alpha_main, alpha_mix, alpha_rep)) %>%
  left_join(species_metadata) %>%
  mutate(Treatment = factor(case_when(Treatment == 1 ~ "Control",
                                      Treatment == 2 ~ "Pesticide"),
                            levels = c("Control", "Pesticide")),
         name = factor(case_when(name == "alpha_main" ~ "Main",
                                 name == "alpha_mix" ~ "Mix",
                                 name == "alpha_rep" ~ "Replicate"),
                       levels = c("Main", "Mix", "Replicate"))) %>%
  ggplot() +
  stat_interval(aes(x = name, y = value, color = Treatment), position = position_dodge(), .width = c(0.5, 0.8, 0.95, 0.99), alpha = 1/2, size = 7) +
  scale_color_manual("Treatment", values = c("#0090C8", "#CE178C"), guide = guide_legend(override.aes = list(alpha = 0.9, size = 7))) +
  scale_x_discrete("Hierarchical level") +
  scale_y_continuous("GP marginal standard deviation", expand = c(0,0), limits = c(0, 5.2)) +
  facet_wrap(~ speciesformatted, scales = "free_x", ncol = 2) +
  theme(panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "grey95", size = 0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        strip.text.x = element_markdown(size = 9, face = "bold"),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black", size = 0.3),
        axis.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_markdown(size = 9),
        axis.text.y = element_text(size = 9),
        axis.ticks = element_line(size = 0.3),
        legend.position = "right",
        legend.key = element_blank(),
        legend.key.width = unit(4, "mm"),
        legend.key.height = unit(4, "mm"),
        legend.spacing.y = unit(0.5, "mm"),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 9))
ggsave("HGP_Chlorpyrifos_Figure_S2.pdf", width = 17.7, height = 12, units = "cm", dpi = 600)

# Calculating the credible intervals for the overall variances of the three levels (averaged over species) 
fit_HGP %>%
  spread_draws(alpha_main[speciesid,Treatment], alpha_mix[speciesid,Treatment], alpha_rep[speciesid,Treatment]) %>%
  ungroup() %>%
  mutate(main = (alpha_main)/(alpha_main+alpha_mix+alpha_rep),
         mix = (alpha_mix)/(alpha_main+alpha_mix+alpha_rep),
         rep = (alpha_rep)/(alpha_main+alpha_mix+alpha_rep)) -> VP
mean(VP$main)
quantile(VP$main, c(0.025, 0.975))
mean(VP$mix)
quantile(VP$mix, c(0.025, 0.975))
mean(VP$rep)
quantile(VP$rep, c(0.025, 0.975))

#-------------------------------------------------------------------------------
# PLOTTING THE HIERARCHICAL GAUSSIAN PROCESS FITS - ALL SPECIES TOGETHER
#-------------------------------------------------------------------------------

A <- data %>%
  mutate(Mix = paste0("Mix ",as.numeric(substr(Replicate, 6, 6))),
         Treatment = c("Control", "Pesticide")[as.numeric(factor(substr(Replicate, 1, 3), levels = c("CTL", "CPF")))],
         Species = factor(Species, c("Daphnia magna", "Daphnia pulicaria", "Daphnia galeata", "Scapholeberis mucronata", "Chlorophyll-a", "Phycocyanin"))) %>%
  filter(!(Species %in% c("Chlorophyll-a", "Phycocyanin"))) %>%
  group_by(Day, Species, Treatment) %>%
  summarise(Median = median(Biomass),
            Lbound = quantile(Biomass, 0.25),
            Ubound = quantile(Biomass, 0.75)) %>%
  left_join(species_metadata, by = c("Species" = "speciesname")) %>%
  ggplot() +
  geom_ribbon(aes(x = Day, ymin = Lbound, ymax = Ubound, fill = speciesformatted), alpha = 0.25) +
  geom_line(aes(x = Day, y = Median, color = speciesformatted), alpha = 0.75, linewidth = 0.15) +
  geom_vline(data = data.frame(x = c(10, 24, 38), type = "dashed"), aes(linetype = type, xintercept = x), color = "grey50", size = 0.3) +
  scale_color_manual("Species", values = c("#EB4235", "#FCBC05", "#32A953", "#3E85F4"), limits = levels(species_metadata$speciesformatted)[1:4], guide = guide_legend(override.aes = list(alpha = 1, linewidth = 1, fill = NA))) +
  scale_fill_manual("Species", values = c("#EB4235", "#FCBC05", "#32A953", "#3E85F4"), limits = levels(species_metadata$speciesformatted)[1:4], guide = "none") +
  facet_wrap(~ Treatment, ncol = 2) +
  scale_x_continuous("Time (days)", breaks = c(0, seq(10, 51, by = 7)), labels = c("<span style = 'font-size:6pt'>Inoculation<br>moment</span>", seq(10, 51, by = 7)), expand = c(0,0)) +
  scale_y_log10("Estimated biomass (µg/L)", labels = scales::label_number()) +
  scale_linetype_manual("", breaks = "dashed", values = "dashed", labels = "Pesticide pulse") +
  coord_cartesian(xlim = c(0,51)) +
  theme(panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "grey95", size = 0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        strip.text.x = element_markdown(size = 9, face = "bold"),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black", size = 0.3),
        axis.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_markdown(size = 9),
        axis.text.y = element_text(size = 9),
        axis.ticks = element_line(size = 0.3),
        legend.position = "right",
        legend.key = element_blank(),
        legend.key.width = unit(5, "mm"),
        legend.key.height = unit(5, "mm"),
        legend.spacing.y = unit(0.5, "mm"),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_markdown(size = 9))

fit_HGP %>%
  spread_draws(f_main[speciesid,Treatment,x]) %>%
  filter(!(speciesid > 4 & x < 11)) %>%
  left_join(species_metadata) %>%
  mutate(f_main = f_main,
         Treatment = factor(case_when(Treatment == 1 ~ "Control",
                                      Treatment == 2 ~ "Pesticide"),
                            levels = c("Control", "Pesticide")),
         x = x - 1,
         speciesformatted = droplevels(speciesformatted)) -> plotdata
B <- plotdata %>%
  filter(speciesid < 5) %>%
  ggplot() +
  stat_lineribbon(aes(x = x, y = f_main, color = speciesformatted, fill = speciesformatted), .width = c(0.5, 0.8, 0.95, 0.99), alpha = 1/5, size = 0.1) +
  #geom_line(data = data, aes(x = Week, y = density, color = Treatment, group = Replicate), alpha = 0.6, size = 0.12) +
  geom_vline(data = data.frame(x = c(10, 24, 38), type = "dashed"), aes(linetype = type, xintercept = x), color = "grey50", size = 0.3) +
  scale_color_manual("Species", values = c("#EB4235", "#FCBC05", "#32A953", "#3E85F4"), limits = levels(species_metadata$speciesformatted)[1:4], guide = guide_legend(override.aes = list(alpha = 1, linewidth = 1, fill = NA))) +
  scale_fill_manual("Species", values = c("#EB4235", "#FCBC05", "#32A953", "#3E85F4"), limits = levels(species_metadata$speciesformatted)[1:4], guide = "none") +
  scale_x_continuous("Time (days)", breaks = c(0, seq(10, 51, by = 7)), labels = c("<span style = 'font-size:6pt'>Inoculation<br>moment</span>", seq(10, 51, by = 7)), expand = c(0,0)) +
  scale_y_continuous("Estimated biomass (µg/L)", breaks = log(1+round(exp(0:8) - 1)), labels = round(exp(0:8) - 1), expand = c(0,0)) +
  scale_linetype_manual("", breaks = "dashed", values = "dashed", labels = "Pesticide pulse") +
  coord_cartesian(ylim = c(0, NA)) +
  facet_wrap(~ Treatment, ncol = 2) +
  theme(panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "grey95", size = 0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        strip.text.x = element_markdown(size = 9, face = "bold"),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black", size = 0.3),
        axis.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_markdown(size = 9),
        axis.text.y = element_text(size = 9),
        axis.ticks = element_line(size = 0.3),
        legend.position = "right",
        legend.key = element_blank(),
        legend.key.width = unit(5, "mm"),
        legend.key.height = unit(5, "mm"),
        legend.spacing.y = unit(0.5, "mm"),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_markdown(size = 9))

A / B + plot_layout(guides = 'collect') + plot_annotation(tag_levels = c("a"), tag_prefix = '(', tag_suffix = ')') & 
  theme(plot.tag = element_text(size = 9, face = "bold", hjust = 0, vjust = 0))
ggsave("HGP_Chlorpyrifos_Figure_speciestogether.pdf", width = 17.7, height = 10, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# DERIVED QUANTITIES
#-------------------------------------------------------------------------------

peakdata <- data.frame(speciesid = NA, treatmentid = NA, .draw = NA, maxbiomass = NA, peak = NA, AUC = NA)[-1,]
for (i in 1:6) {
  for (j in 1:2) {
    fit_HGP %>%
      spread_draws(f_main[speciesid,Treatment,x]) %>%
      filter(speciesid == i, Treatment == j) %>%
      select(.draw, x, f_main) %>%
      pivot_wider(id_cols = .draw, names_from = x, values_from = f_main) %>%
      select(-.draw) -> temp
    peak <- sapply(1:nrow(temp), function(i) seq(0, 51, by = 0.001)[which.max(spline(1:52 - 1, as.numeric((temp[i,])), xout = seq(0, 51, by = 0.001))$y)])
    peakdata <- rbind(peakdata,
                      data.frame(speciesid = i,
                                 treatmentid = j,
                                 .draw = 1:nrow(temp),
                                 maxbiomass = sapply(1:nrow(temp), function(i) max(spline(1:52 - 1, exp(as.numeric((temp[i,])))-1, xout = seq(0, 51, by = 0.001))$y)),
                                 peak = sapply(1:nrow(temp), function(i) seq(0, 51, by = 0.001)[which.max(spline(1:52 - 1, as.numeric((temp[i,])), xout = seq(0, 51, by = 0.001))$y)]),
                                 AUC = rowSums(temp)))
  }
}

peakdata %>%
  pivot_longer(c(maxbiomass, peak, AUC)) %>%
  pivot_wider(id_cols = c(speciesid, .draw, name), names_from = treatmentid, values_from = value) %>%
  mutate(difference = `1` - `2`) %>%
  group_by(speciesid, name) %>%
  summarise(importance = mean(difference > 0),
            posteriormean = mean(difference),
            lowerCrI = quantile(difference, 0.025),
            upperCrI = quantile(difference, 0.975)) %>%
  left_join(species_metadata) %>%
  arrange(name)

# Day of maximum biomass subplot
peakday <- peakdata %>%
  left_join(species_metadata) %>%
  group_by(speciesid, treatmentid) %>%
  filter(peak > quantile(peak, 0.001) & peak < quantile(peak, 0.999),
         speciesid < 5) %>%
  ggplot() +
  stat_slab(aes(x = peak, y = 0, fill = factor(treatmentid)), alpha = 0.5, adjust = 2, normalize = "panels") +
  geom_hline(aes(yintercept = 0)) +
  stat_pointinterval(aes(x = peak, y = -0.15, color = factor(treatmentid)), point_alpha = 1, .width = 0.95, interval_alpha = 1, point_size = 0.75, interval_size = 0.5, position = position_dodge(width = 0.15)) +
  scale_color_manual("Treatment", labels = c("Control", "Pesticide"), values = c("#0090C8", "#CE178C")) +
  scale_fill_manual("Treatment", labels = c("Control", "Pesticide"), values = c("#0090C8", "#CE178C")) +
  scale_x_continuous("Day of maximum biomass") +
  scale_y_continuous("Posterior density", breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  facet_wrap(~ speciesformatted, scales = "free_x", ncol = 4) +
  coord_cartesian(ylim = c(-0.25, 1)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(color = "grey93", size = 0.3),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        strip.text.x = element_markdown(size = 9, face = "bold"),
        strip.background = element_blank(),
        axis.title.x = element_text(size = 9, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_markdown(size = 9),
        axis.text.y = element_text(size = 9),
        axis.ticks = element_line(size = 0.3),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.key.width = unit(5, "mm"),
        legend.key.height = unit(3, "mm"),
        legend.spacing.y = unit(0.5, "mm"),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 9))

# Maximum biomass subplot
maxbiomass <- peakdata %>%
  left_join(species_metadata) %>%
  group_by(speciesid, treatmentid) %>%
  filter(maxbiomass > quantile(maxbiomass, 0.005) & maxbiomass < quantile(maxbiomass, 0.995),
         speciesid < 5) %>%
  ggplot() +
  stat_slab(aes(x = maxbiomass, y = 0, fill = factor(treatmentid)), alpha = 0.5, adjust = 2, normalize = "panels") +
  geom_hline(aes(yintercept = 0)) +
  stat_pointinterval(aes(x = maxbiomass, y = -0.15, color = factor(treatmentid)), point_alpha = 1, .width = 0.95, interval_alpha = 1, point_size = 0.75, interval_size = 0.5, position = position_dodge(width = 0.15)) +
  scale_color_manual("Treatment", labels = c("Control", "Pesticide"), values = c("#0090C8", "#CE178C")) +
  scale_fill_manual("Treatment", labels = c("Control", "Pesticide"), values = c("#0090C8", "#CE178C")) +
  scale_x_continuous("Maximum biomass (µg/L)") +
  scale_y_continuous("Posterior density", breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  facet_wrap(~ speciesformatted, scales = "free_x", ncol = 4) +
  coord_cartesian(ylim = c(-0.25, 1)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(color = "grey93", size = 0.3),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        axis.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_markdown(size = 9),
        axis.text.y = element_text(size = 9),
        axis.ticks = element_line(size = 0.3),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.key.width = unit(5, "mm"),
        legend.key.height = unit(3, "mm"),
        legend.spacing.y = unit(0.5, "mm"),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 9))

# Cumulated biomass subplot
AUC <- peakdata %>%
  left_join(species_metadata) %>%
  group_by(speciesid, treatmentid) %>%
  filter(AUC > quantile(AUC, 0.005) & AUC < quantile(AUC, 0.995),
         speciesid < 5) %>%
  ggplot() +
  stat_slab(aes(x = AUC, y = 0, fill = factor(treatmentid)), alpha = 0.5, adjust = 2, normalize = "panels") +
  geom_hline(aes(yintercept = 0)) +
  stat_pointinterval(aes(x = AUC, y = -0.15, color = factor(treatmentid)), point_alpha = 1, .width = 0.95, interval_alpha = 1, point_size = 0.75, interval_size = 0.5, position = position_dodge(width = 0.15)) +
  scale_color_manual("Treatment", labels = c("Control", "Pesticide"), values = c("#0090C8", "#CE178C")) +
  scale_fill_manual("Treatment", labels = c("Control", "Pesticide"), values = c("#0090C8", "#CE178C")) +
  scale_x_continuous("Cumulated biomass (µg·day/L)") +
  scale_y_continuous("Posterior density", breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  facet_wrap(~ speciesformatted, scales = "free_x", ncol = 4) +
  coord_cartesian(ylim = c(-0.25, 1)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(color = "grey93", size = 0.3),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        axis.title.x = element_text(size = 9, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_markdown(size = 9),
        axis.text.y = element_text(size = 9),
        axis.ticks = element_line(size = 0.3),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.key.width = unit(5, "mm"),
        legend.key.height = unit(3, "mm"),
        legend.spacing.y = unit(0.5, "mm"),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 9))

# Combining the three subplots
peakday / maxbiomass / AUC + plot_layout(guides = 'collect') & theme(legend.position = "bottom")
ggsave("HGP_Chlorpyrifos_Figure_3.pdf", width = 17.7, height = 12, units = "cm")

#-------------------------------------------------------------------------------
# DAPHNIID PHASE PORTRAITS
#-------------------------------------------------------------------------------

fit_HGP %>%
  spread_draws(f_main[speciesid,treat,time]) %>%
  mutate(biomass = exp(f_main)-1,
         time = time - 1) %>%
  select(treat, speciesid, time, biomass, .draw) -> basedata

times <- seq(0, 51, by = 0.2)
dynamics <- data.frame(treat = NA, time = NA, .draw = NA, from = NA, N = NA, to = NA, r = NA)[-1,]
for (treat in 1:2) {
  for (from in 1:4) {
    for (to in 1:4) {
      for (.draw in seq(1, 10000, by = 20)) {
        biomass_from <- basedata$biomass[which(basedata$treat == treat & basedata$speciesid == from & basedata$.draw == .draw)]
        biomass_to <- basedata$biomass[which(basedata$treat == treat & basedata$speciesid == to & basedata$.draw == .draw)]
        splinetemp <- spline(0:51, biomass_to, xout = sort(c(times, times + 0.001)))$y
        dynamics <- rbind(dynamics,
                          cbind(treat = treat,
                                time = times,
                                .draw = .draw,
                                from = from,
                                N = spline(0:51, biomass_from, xout = times)$y,
                                to = to,
                                r = 1000*(splinetemp[(1:length(times))*2] - splinetemp[(1:length(times))*2-1])))
        print(paste0("Treatment: ", treat, "; from species: ", from, "; to species: ", to, "; draw: ", .draw, "."))
      }
    }
  }
}

dynamics %>%
  filter(from <= 4, to <= 4) %>%
  left_join(species_metadata, by = c("from" = "speciesid")) %>%
  left_join(species_metadata, by = c("to" = "speciesid")) %>%
  mutate(from = speciesname_br.x,
         to = speciesname_br.y,
         Treatment = factor(case_when(treat == 1 ~ "Control",
                                      treat == 2 ~ "Pesticide"),
                            levels = c("Control", "Pesticide"))) -> dynamicsplotdata

dynamicsplotdata %>%
  group_by(treat, Treatment, from, to, time) %>%
  summarise(r = mean(r),
            N = mean(N)) %>%
  ggplot() +
  geom_path(data = dynamicsplotdata, aes(x = r, y = N, group = paste(treat, .draw), color = Treatment), alpha = 0.05, size = 0.2) +
  geom_path(aes(x = r, y = N, color = Treatment), alpha = 1) +
  scale_color_manual("Treatment", values = c("#0090C8", "#CE178C"),
                     guide = guide_legend(override.aes = list(alpha = 1, size = 1))) +
  scale_x_continuous("Growth rate (µg/(L·day))") +
  scale_y_continuous("Biomass (µg/L)") +
  facet_grid(from ~ to, scales = "free") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid = element_line(color = "grey93"),
        strip.text.x = element_text(size = 9, face = "bold.italic"),
        strip.text.y = element_text(size = 9, face = "bold.italic"),
        strip.background = element_blank(),
        axis.title = element_text(size = 9, face = "bold"),
        axis.text = element_text(size = 9),
        axis.ticks = element_line(size = 0.3),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.key.width = unit(5, "mm"),
        legend.key.height = unit(5, "mm"),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_markdown(size = 9))
ggsave("HGP_Chlorpyrifos_Figure_4.png", width = 17.7, height = 16, units = "cm", dpi = 600)
ggsave("HGP_Chlorpyrifos_Figure_4.pdf", width = 17.7, height = 16, units = "cm")
