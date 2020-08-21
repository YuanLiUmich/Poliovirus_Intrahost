

### Project: Poliovirus sequencing from Matlab
### Purpose: Analyze codon frequencies.
### Working directory: Poliovirus_Intrahost

# ================================= Read in data, import packages ============================

library(tidyverse)
library(wesanderson)
library(ggridges)
library(LaCroixColoR)

palette <- wes_palette("Darjeeling1")
darjeeling2_palette <- wes_palette("Darjeeling2")

codon_2908 <- read.csv("data/processed/CodonFrequencies/codon_2908_frequencies.csv", stringsAsFactors = FALSE)
codon_2782 <- read.csv("data/processed/CodonFrequencies/codon_2782_frequencies.csv", stringsAsFactors = FALSE)
codon_2608 <- read.csv("data/processed/CodonFrequencies/codon_2608_frequencies.csv", stringsAsFactors = FALSE)

metadata <- read.csv("data/metadata/specimens_Sabin2Positive_CTbelow37_processedFinal_deID_runInfo_covInfo_week_final.csv")

# =================================== Look at the data =====================================

codon_2908 <- filter(codon_2908, Depth > 150 & ID %in% metadata$sequencingID_1)
codon_2908 <- mutate(codon_2908, sequencingID_1 = as.integer(ID))
codon_2908 <- left_join(codon_2908, metadata, by = "sequencingID_1")

codon_2908_mOPV2 <- filter(codon_2908, mOPV2 == 1 & transmissionAcquired == 0 & collectionWeekRelative > 0)

aa2909.cross.sectional.plot <- ggplot(codon_2908_mOPV2, aes(x = as.factor(collectionWeekRelative), y = Frequency, fill = Residue)) + 
  geom_boxplot(aes(as.factor(collectionWeekRelative)), notch = FALSE, outlier.shape = NA) + 
  theme_bw() + 
  geom_jitter(width = 0.2) +
  xlab("Week Post mOPV2 Administration") + 
  ylab("Frequency") + 
  facet_wrap(~Residue)

codon_2908_mOPV2_major <- filter(codon_2908_mOPV2, Frequency > 0.05)
codon_2908_mOPV2_major$Residue <- factor(codon_2908_mOPV2_major$Residue, levels = c("I", "T", "V", "N", "F", "S", "H", "Q", "W", "Y"))
aminoacid.choice.plot <- ggplot(codon_2908_mOPV2_major, aes(x = as.factor(Residue), y = ..count.., fill = as.factor(Residue))) + 
  geom_bar() +
  xlab("Residue") +
  ylab("Number of Samples Above 5% Frequency") +
  theme_bw() +
  scale_fill_manual(name = "", values = rev(c(palette, darjeeling2_palette[2], "black"))) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # 5 by 4
aminoacid.choice.plot

# Example plot
codon_2908_mOPV2 %>%
  filter(Frequency > 0.03) %>%
  group_by(sid) %>%
  filter(length(unique(collectionWeekRelative)) > 1 & length(unique(Residue)) > 2) -> codon_2908_mOPV2_barplot
aa2909.bar.plot <- ggplot(codon_2908_mOPV2_barplot, aes(x = as.factor(collectionWeekRelative), y = Frequency, fill = Residue)) + 
  geom_bar(stat = "identity", position = position_fill()) + 
  theme_bw() + 
  xlab("Week Post mOPV2 Administration") + 
  ylab("Frequency") + 
  facet_wrap(~sid)

# now need to sum frequencies of revertant group?
codon_2908_mOPV2_group <- mutate(codon_2908_mOPV2, ResidueGroup = NA)
codon_2908_mOPV2_group <- mutate(codon_2908_mOPV2_group, ResidueGroup = ifelse(Residue != "I", "Revertant", ResidueGroup))
codon_2908_mOPV2_group <- mutate(codon_2908_mOPV2_group, ResidueGroup = ifelse(Residue == "I", "Ancestral", ResidueGroup))

codon_2908_mOPV2_group %>%
  group_by(ID, ResidueGroup) %>%
  mutate(FrequencySum = sum(Frequency)) %>%
  ungroup() %>%
  group_by(sid) %>%
  filter(length(unique(collectionWeekRelative)) > 0) -> codon_2908_mOPV2_group_plot

ggplot(codon_2908_mOPV2_group_plot, aes(x = collectionWeekRelative, y = FrequencySum)) +
  geom_point() + 
  geom_line(aes(color = as.factor(sid))) +
  facet_wrap(~ResidueGroup) +
  theme(legend.position = "none")

# Add zeros
IDs_2909 <- unique(codon_2908_mOPV2_group_plot$anonSpecId)
codon_2908_mOPV2_group_plot_I143X <- filter(codon_2908_mOPV2_group_plot, ResidueGroup == "Revertant")
complete_data_2909 <- data.frame(anonSpecId = IDs_2909, collectionWeekRelative = NA, mutation = "VP1-I143X")
complete_data_2909 <- mutate(complete_data_2909, collectionWeekRelative = metadata$collectionWeekRelative[match(anonSpecId, metadata$anonSpecId)])
complete_data_2909 <- mutate(complete_data_2909, frequency = codon_2908_mOPV2_group_plot_I143X$FrequencySum[match(anonSpecId, codon_2908_mOPV2_group_plot_I143X$anonSpecId)])
complete_data_2909 <- mutate(complete_data_2909, sid = codon_2908_mOPV2_group_plot$sid[match(anonSpecId, codon_2908_mOPV2_group_plot$anonSpecId)])
complete_data_2909$frequency[is.na(complete_data_2909$frequency)] <- 0

write.csv(complete_data_2909, "data/processed/VP1_143_RevertantCodonFrequencies.csv")

VP1.143.cross.sectional.plot <- ggplot(complete_data_2909, aes(x = collectionWeekRelative, y = frequency, fill = mutation)) + 
  #geom_boxplot(aes(as.factor(collectionWeekRelative)), notch = FALSE, outlier.shape = NA) + 
  theme_bw() + 
  geom_jitter(width = 0.2) +
  xlab("Week Post mOPV2 Administration") + 
  ylab("Probability of Transmission") + 
  geom_smooth(span = 1, method = "loess", formula = y ~ x)
  #facet_wrap(~mutation) +
  #scale_fill_manual(name = "Mutation", values = c(palette[3], palette[2], palette[5]))

# ======================= Need to add in zeros for the residue groups. That way the lines connect better. ===================
# Can we shade the area under the curve? Muller-like plot?
# https://stackoverflow.com/questions/50065471/invert-colored-area-geom-area
# maybe here, take out I143, and fill in any spot where I143X is effectively 0. Then fill in above/below the I143X line.

codon_2908_mOPV2_test <- mutate(codon_2908_mOPV2, ResidueGroup = NA)
codon_2908_mOPV2_test <- mutate(codon_2908_mOPV2_test, ResidueGroup = ifelse(Residue == "I", "I143", "I143X"))

codon_2908_mOPV2_test %>%
  group_by(ID, ResidueGroup) %>%
  mutate(FrequencySum = sum(Frequency)) %>%
  ungroup() %>%
  group_by(sid) %>%
  filter(length(unique(collectionWeekRelative)) > 0) -> codon_2908_mOPV2_test_plot

# Eh just to try it out
codon_2908_mOPV2_ridge <- filter(codon_2908_mOPV2_test_plot, ResidueGroup == "I143X") 
ggplot(data = codon_2908_mOPV2_ridge, aes(x = FrequencySum, y = as.factor(collectionWeekRelative))) +
  geom_density_ridges(jittered_points = TRUE) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges()

ggplot(codon_2908_mOPV2_test_plot, aes(x = collectionWeekRelative, y = FrequencySum)) +
  geom_point(aes(color = as.factor(ResidueGroup))) + 
  geom_line(aes(color = as.factor(ResidueGroup))) +
  facet_wrap(~sid) +
  theme(legend.position = "right") +
  #geom_area(aes(fill = as.factor(ResidueGroup))) +
  ylim(c(0, 1)) +
  theme_classic() +
  xlab("Week Post Vaccination") +
  ylab("Frequency") +
  scale_fill_manual(name = "Mutation", values = c(palette[2], palette[3], palette[1])) +
  scale_color_manual(name = "Mutation", values = c(palette[2], palette[3], palette[1]))
  
# Add in zeros and plot
# This is ugly but it works
codon_2908_mOPV2_test_plot %>%
  filter(ResidueGroup == "I143X") %>%
  select(sid, anonSpecId, ResidueGroup, FrequencySum, collectionWeekRelative) %>%
  distinct() -> codon_2908_mOPV2_test_df

new_data_df = data.frame()
for(sid_index in unique(codon_2908_mOPV2_test_df$sid))
{
  sub_df <- filter(codon_2908_mOPV2_test_df, sid == sid_index)
  
  # what weeks in the primary data do we have data on? If they are not here, make a row where FreqSum is 0.
  weeks <- unique(filter(codon_2908_mOPV2_test_plot, sid == sid_index)$collectionWeekRelative)
  weeks_absent <- weeks[!weeks %in% unique(sub_df$collectionWeekRelative)]
  if(length(weeks_absent) == 0)
  {
    next
  }

  
  weeks_absent_df <- data.frame()
  for(week in weeks_absent)
  {
    base_row <- as.data.frame(sub_df[1,])
    base_row <- mutate(base_row, collectionWeekRelative = week)
    base_row <- mutate(base_row, FrequencySum = 0)
    #print(base_row)
    weeks_absent_df <- rbind(weeks_absent_df, base_row)
  }
  
  #print(sid_index)
  #print(sub_df)
  #print(weeks_absent_df)
  #all <- rbind(sub_df, weeks_absent_df)
  new_data_df <- rbind(new_data_df, as.data.frame(sub_df), as.data.frame(weeks_absent_df))
}


muller.plot <- ggplot(new_data_df, aes(x = collectionWeekRelative, y = FrequencySum)) +
  #geom_point(aes(color = as.factor(ResidueGroup))) + 
  geom_line(aes(color = as.factor(ResidueGroup))) +
  facet_wrap(~sid) +
  theme(legend.position = "right") +
  #geom_area(aes(fill = as.factor(ResidueGroup))) +
  ylim(c(0, 1)) +
  theme_bw() +
  xlab("Week Post Vaccination") +
  ylab("Frequency") +
  scale_fill_manual(name = "Mutation", values = c(palette[2], palette[1], palette[1])) +
  #scale_color_manual(name = "Mutation", values = c(palette[2], palette[3], palette[1])) +
  theme(legend.position = "none") +
  geom_ribbon(aes(x = collectionWeekRelative, ymin = 0, ymax = FrequencySum, fill = palette[1])) +
  geom_ribbon(aes(x = collectionWeekRelative, ymin = FrequencySum, ymax = 1, fill = palette[2]))

# =================================== Codon 2782 =====================================

codon_2782 <- filter(codon_2782, Depth > 150 & ID %in% metadata$sequencingID_1)
codon_2782 <- mutate(codon_2782, sequencingID_1 = as.integer(ID))
codon_2782 <- left_join(codon_2782, metadata, by = "sequencingID_1")

codon_2782_mOPV2 <- filter(codon_2782, mOPV2 == 1 & transmissionAcquired == 0 & collectionWeekRelative > 0)

aa2782.cross.sectional.plot <- ggplot(codon_2782_mOPV2, aes(x = as.factor(collectionWeekRelative), y = Frequency, fill = Residue)) + 
  geom_boxplot(aes(as.factor(collectionWeekRelative)), notch = FALSE, outlier.shape = NA) + 
  theme_bw() + 
  geom_jitter(width = 0.2) +
  xlab("Week Post mOPV2 Administration") + 
  ylab("Frequency") + 
  facet_wrap(~Residue)

codon_2782_mOPV2 %>%
  filter(Frequency > 0.03) %>%
  group_by(sid) %>%
  filter(length(unique(collectionWeekRelative)) > 1) -> codon_2782_mOPV2_barplot
aa2782.bar.plot <- ggplot(codon_2782_mOPV2_barplot, aes(x = as.factor(collectionWeekRelative), y = Frequency, fill = Residue)) + 
  geom_bar(stat = "identity", position = position_fill()) + 
  theme_bw() + 
  xlab("Week Post mOPV2 Administration") + 
  ylab("Frequency") + 
  facet_wrap(~sid)

# =================================== Codon 2608 =====================================

codon_2608 <- filter(codon_2608, Depth > 150 & ID %in% metadata$sequencingID_1)
codon_2608 <- mutate(codon_2608, sequencingID_1 = as.integer(ID))
codon_2608 <- left_join(codon_2608, metadata, by = "sequencingID_1")

codon_2608_mOPV2 <- filter(codon_2608, mOPV2 == 1 & transmissionAcquired == 0 & collectionWeekRelative > 0)

aa2608.cross.sectional.plot <- ggplot(codon_2608_mOPV2, aes(x = as.factor(collectionWeekRelative), y = Frequency, fill = Residue)) + 
  geom_boxplot(aes(as.factor(collectionWeekRelative)), notch = FALSE, outlier.shape = NA) + 
  theme_bw() + 
  geom_jitter(width = 0.2) +
  xlab("Week Post mOPV2 Administration") + 
  ylab("Frequency") + 
  facet_wrap(~Residue)

codon_2608_mOPV2 %>%
  filter(Frequency > 0.03) %>%
  group_by(sid) %>%
  filter(length(unique(collectionWeekRelative)) > 1) -> codon_2608_mOPV2_barplot
aa2608.bar.plot <- ggplot(codon_2608_mOPV2_barplot, aes(x = as.factor(collectionWeekRelative), y = Frequency, fill = Residue)) + 
  geom_bar(stat = "identity", position = position_fill()) + 
  theme_bw() + 
  xlab("Week Post mOPV2 Administration") + 
  ylab("Frequency") + 
  facet_wrap(~sid)

