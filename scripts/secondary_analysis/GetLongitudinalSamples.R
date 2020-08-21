
### Project: Poliovirus Intra-host
### Purpose: Get number of samples we have on each individual, per time point, etc. to set up longitudinal analyses.

# ================== Import packages, read in data =================

library(tidyverse)
library(wesanderson)
library(gplots)
library(reshape2)
library(lubridate)

metadata <- read.csv("data/metadata/specimens_Sabin2Positive_CTbelow37_processedFinal_deID_runInfo_covInfo_week_final.csv.csv")

# ============================ Samples per individual ======================

metadata_any <- filter(metadata, any_data == TRUE)
as.data.frame(table(metadata_any$sid)) %>% rename(sid = Var1, numSamples = Freq) -> anySab2DataPerIndividual

metadata_partial <- filter(metadata, partialGenome == TRUE)
as.data.frame(table(metadata_partial$sid)) %>% rename(sid = Var1, numSamples = Freq) -> partialGenomesPerIndividual

metadata_variant <- filter(metadata, Sab2_variant_qualified == TRUE)
as.data.frame(table(metadata_variant$sid)) %>% rename(sid = Var1, numSamples = Freq) -> variantQualityGenomesPerIndividual

# ================================ Samples by time ========================

metadata_campaign <- filter(metadata, studyPeriod == "mOPV2_campaign", any_data == TRUE)

metadata_campaign %>% 
  select(sid, anonSpecId, collectionWeekRelative) %>%
  spread(collectionWeekRelative, anonSpecId) -> metadata_campaign_longitudinal_weeks

metadata_campaign %>% 
  select(sid, anonSpecId, SPECDT) %>%
  spread(SPECDT, anonSpecId) -> metadata_campaign_longitudinal_days

# Can we graph this?
metadata_graph <- filter(metadata, collectionDateRelativeExact >= 0 & any_data == TRUE)
data.frame(table(metadata_graph$sid)) %>% mutate(sid = as.integer(as.character(Var1)), num = Freq) %>% select(-Var1, -Freq) -> num_samples_per_sid
data.frame(sid = unique(metadata_graph$sid)) -> new_id_df
new_id_df$sid_plot = 1:nrow(new_id_df)
metadata_graph <- left_join(metadata_graph, new_id_df, by = "sid")
metadata_graph <- left_join(metadata_graph, num_samples_per_sid, by = "sid")
metadata_graph <- metadata_graph[order(metadata_graph$num, decreasing = TRUE),]
metadata_graph$sort_order <- 1:nrow(metadata_graph)

metadata_graph %>% group_by(sid_plot) %>% mutate(week_min = min(collectionWeekRelative), week_max = max(collectionWeekRelative)) %>% ungroup() -> metadata_graph

ggplot(metadata_graph, aes(x = week_min, xend = week_max, y = sid_plot, yend = sid_plot, group = as.factor(num))) + 
  geom_segment(color = "black", size = 0.1) +
  ylab("Individual") +
  xlab("Weeks Post Vaccination") +
  geom_point(aes(y = sid_plot, x = collectionWeekRelative, color = as.factor(num)), size = 1.2) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) + 
  theme_classic() 

# Samples by titer; plot successes with different shapes/saturation etc.
metadata_graph <- mutate(metadata_graph, titer_group = as.factor(floor(log10(S2_copiesPerUl_inTNA))))
metadata_graph <- mutate(metadata_graph, titer_group = ifelse(S2_copiesPerUl_inTNA < 900, "< 9E2", titer_group))
metadata_graph <- mutate(metadata_graph, titer_group = ifelse(S2_copiesPerUl_inTNA > 900 & S2_copiesPerUl_inTNA < 45000, "9E2 - 4.5E4", titer_group))
metadata_graph <- mutate(metadata_graph, titer_group = ifelse(S2_copiesPerUl_inTNA > 45000 & S2_copiesPerUl_inTNA < 450000, "4.5E4 - 4.5E5", titer_group))
metadata_graph <- mutate(metadata_graph, titer_group = ifelse(S2_copiesPerUl_inTNA > 450000, "> 4.5E5", titer_group))
metadata_graph$titer_group <- factor(metadata_graph$titer_group, levels = c("> 4.5E5", "4.5E4 - 4.5E5", "9E2 - 4.5E4", "< 9E2"))

ggplot(metadata_graph, aes(x = week_min, xend = week_max, y = sid_plot, yend = sid_plot, group = titer_group)) + 
  geom_segment(color = "black", size = 0.1) +
  ylab("Individual") +
  xlab("Weeks Post Vaccination") +
  geom_point(aes(y = sid_plot, x = collectionWeekRelative, color = titer_group), size = 1, alpha = 0.2) +
  geom_point(data = filter(metadata_graph, any_data == TRUE), aes(y = sid_plot, x = collectionWeekRelative, color = titer_group), size = 1, alpha = 1) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) + 
  theme_classic()

# ================================ Which individuals were positive for OPV2 before the vaccination campaign? ===========================

metadata_full <- read.csv("data/metadata/specimens_Sabin2Positive_CTbelow37_all_deID.csv")
metadata_full <- filter(metadata_full, sid %in% metadata$sid & studyPeriod == "mOPV2_campaign")
metadata_full %>% mutate(collectionDateFloor = floor_date(mdy(SPECDT), "week")) %>%
  mutate(collectionDateRelative = collectionDateFloor - ymd("2016-1-24")) %>%
  mutate(collectionWeekRelative = round(collectionDateRelative / 7, 0)) -> metadata_full
metadata_full <- filter(metadata_full, collectionWeekRelative == 0 & S2 == 1)
sid_positiveAtWeekZero <- metadata_full$sid

metadata_any_positiveAtWeekZero <- filter(metadata_any, sid %in% sid_positiveAtWeekZero) # We should follow up with these samples to see if they have mixed infections of Sabin2.

# =================================== Sampling graph per individual, color by coverage type ======================

metadata_graph <- mutate(metadata_graph, cov_category = "NA")
metadata_graph <- mutate(metadata_graph, cov_category = ifelse(any_data == TRUE, "Partial", cov_category))
metadata_graph <- mutate(metadata_graph, cov_category = ifelse(Sab2_consensus_qualified == TRUE, "Consensus", cov_category))
metadata_graph <- mutate(metadata_graph, cov_category = ifelse(Sab2_variant_qualified == TRUE, "Variant-quality", cov_category))

ggplot(metadata_graph, aes(x = week_min, xend = week_max, y = sid_plot, yend = sid_plot, group = cov_category)) + 
  geom_segment(color = "black", size = 0.1) +
  ylab("Individual") +
  xlab("Weeks Post Vaccination") +
  geom_point(aes(y = sid_plot, x = collectionWeekRelative, color = cov_category), size = 1.2) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) + 
  theme_classic() +
  scale_color_manual(name = "", values = c(palette[2], palette[5], palette[1]))
