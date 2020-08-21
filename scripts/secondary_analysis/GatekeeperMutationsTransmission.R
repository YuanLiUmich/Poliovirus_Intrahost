
### Project: Poliovirus sequencing from Matlab
### Purpose: Are the gatekeeper mutations found in transmission recipients?
### Working directory: Poliovirus_Intrahost

# ================== Import packages, read in data =================

library(tidyverse)
library(wesanderson)

palette <- wes_palette("Darjeeling1")
darjeeling2_palette <- wes_palette("Darjeeling2")
aquatic_palette <- wes_palette("Zissou1")

metadata <- read.csv("data/metadata/specimens_Sabin2Positive_CTbelow37_processedFinal_deID_runInfo_covInfo_week_final.csv")
variants <- read.csv("data/processed/all.concat.varQuality.cleaned.variants.csv")

give.n <- function(x)
{
  return(c(y = 1.05, label = length(x))) # experiment with the multiplier to find the perfect position
}

# ======================= Function to get "zeros" in a new dataframe ======================

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

# ==================== Examine gatekeeper mutations ====================

variants %>% 
  filter(mOPV2 == 0 & transmissionAcquired == 1 & collectionWeekRelative > 0) %>% 
  mutate(coverage = cov.tst.fw + cov.tst.bw) %>% 
  filter(coverage > 150) %>%
  filter(!sid %in% c(430, 211, 507)) %>% # our sample is not the first PCR positive one for these individuals
  group_by(sid) %>%
  filter(collectionWeekRelative == min(collectionWeekRelative)) -> variants_transmission

variants_transmission_481 <- filter(variants_transmission, pos == 481)
variants_transmission_2909 <- filter(variants_transmission, pos == 2909)
variants_transmission_398 <- filter(variants_transmission, pos == 398)

length(unique(variants_transmission$anonSpecId)) # 23 transmission-acquired individuals with some data, and we have the first RT-qPCR positive sample.
length(unique(variants_transmission$sid)) # 23, since we got the first one if they had more

length(unique(variants_transmission_481$anonSpecId))
length(unique(variants_transmission_481$sid))
length(unique(variants_transmission_2909$anonSpecId)) 
length(unique(variants_transmission_2909$sid))
length(unique(variants_transmission_398$anonSpecId))
length(unique(variants_transmission_398$sid))

  
transmission_481 <- GetCrossSection(variants_transmission, 481, "A481G", "G", metadata)
transmission_2909 <- GetCrossSection(variants_transmission, 2909, "U2909C", "C", metadata)
transmission_398 <- GetCrossSection(variants_transmission, 398, "U398C", "C", metadata)

transmission_data_bigthree <- rbind(transmission_481, transmission_2909, transmission_398)

bigthree.transmitted.plot <- ggplot(transmission_data_bigthree, aes(x = as.factor(collectionWeekRelative), y = frequency, fill = mutation)) + 
  geom_boxplot(aes(as.factor(collectionWeekRelative)), notch = FALSE, outlier.shape = NA) + 
  theme_bw() + 
  geom_jitter(width = 0.2) +
  xlab("Week Post mOPV2 Administration") + 
  ylab("Frequency") + 
  facet_wrap(~mutation) +
  scale_fill_manual(name = "Mutation", values = c(palette[3], palette[2], palette[5])) + 
  stat_summary(fun.data = give.n, geom = "text", fun.y = median)

transmission_data_bigthree %>%
  group_by(mutation, collectionWeekRelative) %>%
  summarize(fractionTransmitted = sum(frequency > 0) / sum(frequency >= 0), count = n()) -> transmission_data

write.csv(transmission_data, "data/processed/bigthree_fractionTransmitted.csv")

