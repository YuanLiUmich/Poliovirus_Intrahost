
### Project: Poliovirus sequencing from Matlab
### Purpose: Longitudinal analysis of iSNVs
### Working directory: Poliovirus_Intrahost

# ================== Import packages, read in data =================

library(tidyverse)
library(wesanderson)
library(ggbeeswarm)
library(reshape2)
library(lubridate)
library(cowplot)
library(ggridges)
library(patchwork)

metadata <- read.csv("data/metadata/specimens_Sabin2Positive_CTbelow37_processedFinal_deID_runInfo_covInfo_week_final.csv")
variants <- read.csv("data/processed/all.concat.varQuality.cleaned.variants.csv")

palette <- wes_palette("Darjeeling1")
darjeeling2_palette <- wes_palette("Darjeeling2")

# ======================= Support Functions ======================

# Get "zeros" in a new dataframe
GetCrossSection <- function(var_df, position, mut, var_base, meta)
{
  
  var_df_position <- filter(var_df, pos == position)
  ids_position <- unique(var_df_position$anonSpecId)
  var_df_position <- filter(var_df_position, var == var_base)
  new_df <- data.frame(anonSpecId = ids_position, collectionWeekRelative = NA, mutation = mut)
  new_df <- mutate(new_df, collectionWeekRelative = meta$collectionWeekRelative[match(anonSpecId, meta$anonSpecId)])
  new_df <- mutate(new_df, frequency = var_df_position$freq.var[match(anonSpecId, var_df_position$anonSpecId)])
  new_df <- mutate(new_df, sid = variants$sid[match(anonSpecId, variants$anonSpecId)])
  new_df$frequency[is.na(new_df$frequency)] <- 0
  
  return(new_df)
}

# Cross-sectional plot
PlotCrossSection <- function(data, data_box, mut, colors)
{
  
  plot <- ggplot(filter(data, mutation == mut), aes(x = as.factor(collectionWeekRelative), y = frequency, fill = mutation)) + 
    geom_boxplot(data = filter(data_box, mutation == mut), aes(as.factor(collectionWeekRelative)), notch = FALSE, outlier.shape = NA, alpha = 0.9) + 
    theme_bw() + 
    geom_jitter(width = 0.2) +
    xlab("Week Post Vaccination") + 
    ylab("Frequency") + 
    scale_fill_manual(name = "Mutation", values = colors) +
    theme(legend.position = "none") +
    ylim(c(0, 1))
  
  return(plot)
}

# Longitudinal plot
PlotLongitudinal <- function(data, mut, color)
{
  plot <- ggplot(filter(data, mutation == mut), aes(x = collectionWeekRelative, y = frequency, color = as.factor(sid))) + 
    theme_bw() + 
    geom_point(color = "black", size = 2) + 
    geom_line(aes(color = as.factor(sid)), alpha = 0.7) +
    xlab("Week Post Vaccination") + 
    ylab("Frequency") + 
    theme(legend.position = "none") +
    scale_color_manual(values = c(rep(color, 100))) 
  
  return(plot)
}


# ==================== Get data for gatekeeper mutations in mOPV2 recipients ====================

variants_mOPV2_recipients <- filter(variants, mOPV2 == 1 & transmissionAcquired == 0 & collectionWeekRelative > 0 & S2copyNumberPerGram > 900000)

complete_data_481 <- GetCrossSection(variants_mOPV2_recipients, 481, "A481G", "G", metadata)
complete_data_2909 <- GetCrossSection(variants_mOPV2_recipients, 2909, "U2909C", "C", metadata)
complete_data_398 <- GetCrossSection(variants_mOPV2_recipients, 398, "U398C", "C", metadata)

complete_data_bigthree <- rbind(complete_data_481, complete_data_2909, complete_data_398)
#write.csv(complete_data_bigthree, "data/processed/gatekeepers_cross_section.csv")

complete_data_VP1_143 <- read.csv("data/processed/VP1_143_RevertantCodonFrequencies.csv", stringsAsFactors = FALSE) %>% select(-X)
complete_data_VP1_143$frequency[complete_data_VP1_143$frequency < 0.05] <- 0 # Apply the 5% cutoff here to match other data

complete_data_all <- rbind(complete_data_481, complete_data_VP1_143, complete_data_398)


# ================================== Cross-sectional and longitudinal plots ==========================

# Cross-Sectional
complete_data_all %>%
  group_by(mutation, collectionWeekRelative) %>%
  filter(n() > 4) -> complete_data_all_box

cross.sectional.481.plot <- PlotCrossSection(complete_data_all, complete_data_all_box, "A481G", c(palette[4]))
cross.sectional.143.plot <- PlotCrossSection(complete_data_all, complete_data_all_box, "VP1-I143X", c(darjeeling2_palette[2]))
cross.sectional.398.plot <- PlotCrossSection(complete_data_all, complete_data_all_box, "U398C", c(palette[5]))

cross.section.plot <- plot_grid(cross.sectional.481.plot, cross.sectional.143.plot, cross.sectional.398.plot, ncol = 3) # 10 by 3.5

# Longitudinal
complete_data_all %>%
  group_by(mutation, sid) %>%
  filter(n() > 1) %>%
  ungroup() -> complete_data_all_longitudinal

long.481.plot <- PlotLongitudinal(complete_data_all_longitudinal, "A481G", palette[4])
long.143.plot <- PlotLongitudinal(complete_data_all_longitudinal, "VP1-I143X", darjeeling2_palette[2])
long.398.plot <- PlotLongitudinal(complete_data_all_longitudinal, "U398C", palette[5])

long.plot <- plot_grid(long.481.plot, long.143.plot, long.398.plot, ncol = 3) # 10 by 3.5

# ============================== Difference in fixation rates by trial arm? ===============================

unique(variants_mOPV2_recipients$arm)

variants_mOPV2_recipients_IPV <- filter(variants_mOPV2_recipients, arm %in% c("IPVx1", "IPVx2"))
variants_mOPV2_recipients_OPV <- filter(variants_mOPV2_recipients, arm %in% c("tOPV"))

# 481
var_IPV_481 <- GetCrossSection(variants_mOPV2_recipients_IPV, 481, "A481G", "G", metadata)
var_IPV_481 %>%
  group_by(mutation, collectionWeekRelative) %>%
  filter(n() > 4) -> var_IPV_481_box
var_OPV_481 <- GetCrossSection(variants_mOPV2_recipients_OPV, 481, "A481G", "G", metadata)
var_OPV_481 %>%
  group_by(mutation, collectionWeekRelative) %>%
  filter(n() > 4) -> var_OPV_481_box
IPV.481.plot <- PlotCrossSection(var_IPV_481, var_IPV_481_box, "A481G", c(palette[4])) + ggtitle("IPVx1 and IPVx2: A481G")
OPV.481.plot <- PlotCrossSection(var_OPV_481, var_OPV_481_box, "A481G", c(palette[4])) + ggtitle("tOPV: A481G")

# 2909
var_IPV_2909 <- GetCrossSection(variants_mOPV2_recipients_IPV, 2909, "U2909C", "C", metadata)
var_IPV_2909 %>%
  group_by(mutation, collectionWeekRelative) %>%
  filter(n() > 4) -> var_IPV_2909_box
var_OPV_2909 <- GetCrossSection(variants_mOPV2_recipients_OPV, 2909, "U2909C", "C", metadata)
var_OPV_2909 %>%
  group_by(mutation, collectionWeekRelative) %>%
  filter(n() > 4) -> var_OPV_2909_box
IPV.2909.plot <- PlotCrossSection(var_IPV_2909, var_IPV_2909_box, "U2909C", c(darjeeling2_palette[2])) + ggtitle("IPVx1 and IPVx2: U2909C")
OPV.2909.plot <- PlotCrossSection(var_OPV_2909, var_OPV_2909_box, "U2909C", c(darjeeling2_palette[2])) + ggtitle("tOPV: U2909C")

# 398
var_IPV_398 <- GetCrossSection(variants_mOPV2_recipients_IPV, 398, "U398C", "C", metadata)
var_IPV_398 %>%
  group_by(mutation, collectionWeekRelative) %>%
  filter(n() > 4) -> var_IPV_398_box
var_OPV_398 <- GetCrossSection(variants_mOPV2_recipients_OPV, 398, "U398C", "C", metadata)
var_OPV_398 %>%
  group_by(mutation, collectionWeekRelative) %>%
  filter(n() > 4) -> var_OPV_398_box
IPV.398.plot <- PlotCrossSection(var_IPV_398, var_IPV_398_box, "U398C", c(palette[5])) + ggtitle("IPVx1 and IPVx2: U398C")
OPV.398.plot <- PlotCrossSection(var_OPV_398, var_OPV_398_box, "U398C", c(palette[5])) + ggtitle("tOPV: U398C")

compare.arms.plot <- (IPV.481.plot | IPV.2909.plot | IPV.398.plot) / (OPV.481.plot | OPV.2909.plot | OPV.398.plot) # 10 by 7

### Test whether there's a difference at week 1
var_IPV_481_week1 <- filter(var_IPV_481, collectionWeekRelative == 1)
var_OPV_481_week1 <- filter(var_OPV_481, collectionWeekRelative == 1)
wilcox.test(var_IPV_481_week1$frequency, var_OPV_481_week1$frequency) # p = 0.66

var_IPV_2909_week1 <- filter(var_IPV_2909, collectionWeekRelative == 1)
var_OPV_2909_week1 <- filter(var_OPV_2909, collectionWeekRelative == 1)
wilcox.test(var_IPV_2909_week1$frequency, var_OPV_2909_week1$frequency)

var_IPV_398_week1 <- filter(var_IPV_398, collectionWeekRelative == 1)
var_OPV_398_week1 <- filter(var_OPV_398, collectionWeekRelative == 1)
wilcox.test(var_IPV_398_week1$frequency, var_OPV_398_week1$frequency)

# ================================== Change in frequency per week ==========================

complete_data_all_longitudinal %>%
  group_by(sid, mutation) %>%
  mutate(delta = ((frequency - lag(frequency)) / (collectionWeekRelative - lag(collectionWeekRelative)))) -> complete_data_all_longitudinal_delta
complete_data_all_longitudinal_delta <- filter(complete_data_all_longitudinal_delta, !is.na(delta) & !(frequency == 1 & delta == 0))

delta.plot <- ggplot(complete_data_all_longitudinal_delta, aes(x = as.factor(mutation), y = delta)) +
  geom_jitter(width = 0.2, aes(color = mutation)) +
  geom_boxplot(notch = FALSE, outlier.shape = NA, alpha = 0) +
  theme_bw() +
  scale_color_manual(name = "Mutation", values = c(palette[4], darjeeling2_palette[2], palette[5])) +
  theme(legend.position = "none") +
  xlab("Mutation") +
  ylab("Change in Frequency Per Week")

aov_delta <- aov(delta ~ mutation, data = complete_data_all_longitudinal_delta)
summary(aov_delta)
TukeyHSD(aov_delta)

# ================================== Probability of transmission from, given the bottleneck estimate ==========================

bottleneck_info <- read.csv("data/processed/bottleneck_estimates.csv", stringsAsFactors = FALSE)
bottleneck_PA <- filter(bottleneck_info, model == "PA")
Pt_PA <- function(x, l, max_Nb)
{
  s <- 0
  for(i in 1:max_Nb)
  {
    c <- ((1-x)^i) *  (l^i)/((exp(l)-1)*factorial(i))
    s <- s + c
  }
  return(1 - s)
}

complete_data_all <- mutate(complete_data_all, prob_transmission = Pt_PA(frequency, bottleneck_PA$lambda, 100),
                                 prob_transmission_lower = Pt_PA(frequency, bottleneck_PA$lower, 100),
                                 prob_transmission_upper = Pt_PA(frequency, bottleneck_PA$upper, 100))
complete_data_all <- filter(complete_data_all, collectionWeekRelative < 7)

complete_data_bigthree <- mutate(complete_data_bigthree, prob_transmission = Pt_PA(frequency, bottleneck_PA$lambda, 100),
                            prob_transmission_lower = Pt_PA(frequency, bottleneck_PA$lower, 100),
                            prob_transmission_upper = Pt_PA(frequency, bottleneck_PA$upper, 100))
complete_data_bigthree <- filter(complete_data_bigthree, collectionWeekRelative < 7)

complete_data_bigthree %>%
  group_by(mutation, collectionWeekRelative) %>%
  mutate(median_PT = median(prob_transmission), 
         lower_PT = median(prob_transmission_lower), 
         upper_PT = median(prob_transmission_upper)) -> transmission_prob_data

transmission.model.plot <- ggplot() +
  geom_line(data = transmission_prob_data, aes(x = collectionWeekRelative, y = median_PT, color = mutation)) +
  scale_color_manual(name = "Mutation", values = c(palette[4], darjeeling2_palette[2], palette[5])) +
  scale_fill_manual(name = "Mutation", values = c(palette[4], darjeeling2_palette[2], palette[5])) +
  geom_ribbon(data = transmission_prob_data, aes(x = collectionWeekRelative, ymin = lower_PT, ymax = upper_PT, fill = mutation), alpha = 0.5) +
  xlab("Week Post Vaccination") +
  ylab("Probability of Transmission") +
  theme_bw() +
  theme(legend.position = "right", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(seq(1, 6, by = 1)))

# =========================== Compare to transmission incidence by week ============================

transmission_counts <- read.csv("data/processed/householdcontact_transmission_incidence.csv", stringsAsFactors = FALSE)

transmission.count.plot <- ggplot(transmission_counts, aes(x = first, y = count)) + 
  geom_line() +
  geom_point() +
  xlab("Week Post Vaccination") +
  ylab("Household Transmission Events") +
  scale_x_continuous(limits = c(1, 6), breaks = c(seq(1, 6, 1))) +
  scale_y_continuous(limits = c(0, 15), breaks = c(seq(0, 15, 1)))

# =========================== Compare to what actually transmits ======================

transmission_data <- read.csv("data/processed/bigthree_fractionTransmitted.csv", stringsAsFactors = FALSE) 
#transmission_data$mutation[transmission_data$mutation == "U2909C"] <- "VP1-I143X" # would need to get VP1-143 haplotype info for transmission samples

fraction.transmitted.plot <- ggplot() +
  geom_line(data = transmission_data, aes(x = collectionWeekRelative, y = fractionTransmitted, color = mutation)) +
  geom_point(data = transmission_data, aes(x = collectionWeekRelative, y = fractionTransmitted)) +
  scale_color_manual(name = "Mutation", values = c(palette[4], darjeeling2_palette[2], palette[5])) +
  xlab("Week Post Vaccination") +
  ylab("Fraction of Samples with Mutation Present") +
  theme_bw() +
  theme(legend.position = "right", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(seq(1, 6, by = 1)))

model.plot <- plot_grid(transmission.model.plot, fraction.transmitted.plot, ncol = 3) # save as 12 by 4.5
