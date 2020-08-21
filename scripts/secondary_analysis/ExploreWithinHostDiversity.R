

### Project: Poliovirus sequencing from Matlab
### Purpose: Explore within-host diversity by iSNVs
### Working directory: Poliovirus_Intrahost

# ================== Import packages, read in data =================

library(tidyverse)
library(wesanderson)
library(ggsci)
library(LaCroixColoR)
library(patchwork)

palette <- wes_palette("Darjeeling1")
darjeeling2_palette <- wes_palette("Darjeeling2")
aquatic_palette <- wes_palette("Zissou1")
palette_viridis <- c("#FDE725FF", "#73D055FF", "#29AF7FFF", "#238A8DFF", "#39568CFF", "black")
palette_lime <- lacroix_palette("Lime", type = "discrete")

metadata <- read.csv("data/metadata/specimens_Sabin2Positive_CTbelow37_processedFinal_deID_runInfo_covInfo_week_final.csv")
variants <- read.csv("data/processed/all.concat.varQuality.cleaned.variants.csv")

plot.median <- function(x) 
{
  m <- median(x)
  c(y = m, ymin = m, ymax = m)
}

give.n <- function(x)
{
  return(c(y = 1.05, label = length(x))) # experiment with the multiplier to find the perfect position
}

# ================= Get minor SNV per sample ====================

# Only from mOPV2 recipients from campaign period
metadata_quality <- filter(metadata, collectionWeekRelative > 0 & Sab2_variant_qualified == TRUE & mOPV2 == 1 & transmissionAcquired == 0 & S2copyNumberPerGram > 900000)
variants_quality <- filter(variants, anonSpecId %in% metadata_quality$anonSpecId)
minor <- filter(variants_quality, freq.var < 0.5)

minor %>% 
  group_by(anonSpecId) %>% 
  dplyr::summarize(iSNV = length(unique(mutation)), nonsyn_iSNV = length(which(class_factor == "Nonsynonymous")), syn_iSNV = length(which(class_factor == "Synonymous"))) -> minor_counts
metadata_quality_isnv <- merge(metadata_quality, minor_counts, by = "anonSpecId", all.x = TRUE)
metadata_quality_isnv$iSNV[is.na(metadata_quality_isnv$iSNV)] <- 0 # NA means there were no variants; change to 0
metadata_quality_isnv$nonsyn_iSNV[is.na(metadata_quality_isnv$nonsyn_iSNV)] <- 0
metadata_quality_isnv$syn_iSNV[is.na(metadata_quality_isnv$syn_iSNV)] <- 0

# ============================= Raw NS/S ratio over time ==================

minor %>%
  group_by(collectionWeekRelative) %>%
  summarize(ratio = length(which(class_factor == "Nonsynonymous")) / length(which(class_factor == "Synonymous"))) -> ratios

# ========================== Plot different coverage groups =========================

metadata_campaign <- filter(metadata, collectionWeekRelative > 0 & collectionWeekRelative < 11 & xor(mOPV2 == 1, transmissionAcquired == 1))
metadata_campaign <- mutate(metadata_campaign, cov_category = ifelse(any_data == FALSE, "No data", NA))
metadata_campaign <- mutate(metadata_campaign, cov_category = ifelse(Sab2_variant_qualified == TRUE, "Variant-quality", cov_category))
metadata_campaign <- mutate(metadata_campaign, cov_category = ifelse(Sab2_consensus_qualified == TRUE & !Sab2_variant_qualified == TRUE, "Consensus", cov_category))
metadata_campaign <- mutate(metadata_campaign, cov_category = ifelse(any_data == TRUE & !Sab2_consensus_qualified == TRUE & !Sab2_variant_qualified == TRUE, "Partial", cov_category))
metadata_campaign$cov_category <- factor(metadata_campaign$cov_category, levels = c("Variant-quality", "Consensus", "Partial", "No data"))
titer_palette <- c(palette[1], darjeeling2_palette[2], palette[5], "grey")
titer.by.time <- ggplot(metadata_campaign, aes(x = collectionWeekRelative, y = S2copyNumberPerGram, group = cov_category)) +
  geom_jitter(aes(color = cov_category), width = 0.25) +
  ylab("OPV2 Copies Per Gram of Stool") +
  xlab("Week post vaccination") +
  scale_y_log10() +
  theme_bw() + 
  geom_hline(yintercept = 900000, linetype = "dashed", color = "black", size = 0.5) +
  geom_hline(yintercept = 45000000, linetype = "dashed", color = "black", size = 0.5) +
  scale_color_manual(name = "" , values = titer_palette) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) +
  theme(panel.grid.minor = element_blank()) # save as 9 by 6

cov.category.pie <- ggplot(metadata_campaign, aes(x = factor(1), fill = cov_category)) +
  geom_bar(width = 1, position = position_fill()) +
  coord_polar("y", start = 0) +
  scale_fill_manual(name = "" , values = titer_palette) +
  facet_wrap(~collectionWeekRelative) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid  = element_blank()) # save as 6 by 6

# How many samples per individual in each category?
metadata_any <- filter(metadata_campaign, any_data == TRUE)
as.data.frame(table(metadata_any$sid)) %>% rename(sid = Var1, numSamples = Freq) %>% mutate(category = "Partial") -> anySab2DataPerIndividual
metadata_consensus <- filter(metadata_campaign, Sab2_consensus_qualified == TRUE)
as.data.frame(table(metadata_consensus$sid)) %>% rename(sid = Var1, numSamples = Freq) %>% mutate(category = "Consensus") -> consensusGenomesPerIndividual
metadata_variant <- filter(metadata_campaign, cov_category == "Variant-quality")
as.data.frame(table(metadata_variant$sid)) %>% rename(sid = Var1, numSamples = Freq) %>% mutate(category = "Variant-quality") -> variantQualityGenomesPerIndividual
samples_per_indiv <- rbind(anySab2DataPerIndividual, consensusGenomesPerIndividual, variantQualityGenomesPerIndividual)
samples_per_indiv$facet = factor(samples_per_indiv$category, levels = c("Partial", "Consensus", "Variant-quality"))
cov.category.hist <- ggplot(samples_per_indiv, aes(numSamples, fill = facet)) + 
  geom_histogram(binwidth = 1, position = "identity") + 
  theme_classic() + 
  scale_fill_manual(name = "", values = c(palette[5], darjeeling2_palette[2], palette[1])) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8)) +
  xlab("Number of Samples") + 
  ylab("Number of Individuals") # save as 4 by 5

# ===================== iSNV and viral load ========================

quantile(metadata_quality_isnv$iSNV)

snv.by.copynum.reps <- ggplot(metadata_quality_isnv, aes(x = S2copyNumberPerGram, y = iSNV, fill = as.factor(reps))) + 
  geom_point(shape = 21, size = 3) + 
  ylab("iSNV Per Sample") + 
  theme_bw() +
  xlab("Log (base 10) of genomes per gram of stool") +
  scale_x_log10() +
  scale_fill_manual(name = "Replicates" , values = c(darjeeling2_palette[2], palette[2])) # save as 5 by 4

snv.by.copynum.week <- ggplot(metadata_quality_isnv, aes(x = S2copyNumberPerGram, y = iSNV, fill = as.factor(collectionWeekRelative))) + 
  geom_point(shape = 21, size = 3) + 
  ylab("iSNV Per Sample") + 
  theme_bw() +
  xlab("Log (base 10) of genomes per gram of stool") +
  scale_x_log10() +
  scale_fill_manual(name = "Week", values = rev(c(palette, darjeeling2_palette[2])))

snv.by.copynum.lm <- ggplot(metadata_quality_isnv, aes(x = log(S2copyNumberPerGram, 10), y = iSNV)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  ylab("iSNV Per Sample") + 
  xlab("Log (base 10) of genomes per gram of stool") +
  theme_bw()

snv.by.week.plot <- ggplot(metadata_quality_isnv, aes(x = collectionWeekRelative, y = iSNV)) + 
  geom_jitter(width = 0.2) +
  ylab("iSNV Per Sample") + 
  xlab("Week Post Vaccination") +
  theme_bw() +
  geom_smooth(method = "lm") +
  scale_x_continuous(breaks = c(seq(1, 6)))
  
snv_by_gc <- lm(data = metadata_quality_isnv, formula = iSNV ~ log(S2copyNumberPerGram, 10))
summary(snv_by_gc)
AIC(snv_by_gc)

snv_by_week <- lm(data = metadata_quality_isnv, formula = iSNV ~ collectionWeekRelative)
summary(snv_by_week)
AIC(snv_by_week)

snv_lm <- lm(data = metadata_quality_isnv, formula = iSNV ~ log(S2copyNumberPerGram, 10) + collectionWeekRelative)
summary(snv_lm)
AIC(snv_lm)

anova(snv_lm, snv_by_week) # week explains data well; gc does not add to the model.


# ===================== iSNV plots ========================

minor_plot <- filter(minor, !is.na(class_factor))
freq.by.pos.palette <- c(palette[5], palette[3], darjeeling2_palette[2])
freq.by.pos <- ggplot(minor_plot, aes(x = pos, y = freq.var, fill = class_factor)) +
  geom_point(size = 2, shape = 21) +
  theme_bw() + 
  xlab("Genome Position") + 
  ylab("Frequency") +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(seq(0, 7000, by = 1000))) +
  scale_fill_manual(name = "" , values = freq.by.pos.palette) # 10 by 4

variants_quality_compare_NS_S <- filter(minor, class_factor %in% c("Nonsynonymous", "Synonymous"))
freq_histogram <- ggplot(variants_quality_compare_NS_S, aes(x = freq.var, fill = class_factor)) + 
  geom_histogram(color = "white", binwidth = 0.05, position = "dodge", boundary = 0.02) +
  xlab("Frequency") + 
  ylab("Number of iSNV") +
  scale_fill_manual(name = "" , values = c(palette[3], darjeeling2_palette[2])) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # 5 by 4

variants_quality_compare_NS_S$gene <- factor(variants_quality_compare_NS_S$gene, levels = c("VP4", "VP2", "VP3", "VP1", "2A", "2B", "2C", "3A", "3B", "3C", "3D"))

mutation_type_by_gene <- ggplot(variants_quality_compare_NS_S, aes(x = gene, fill = class_factor)) + 
  geom_bar(position = "dodge") +
  xlab("Gene in coding region") + 
  ylab("Number of iSNV") +
  scale_fill_manual(name = "" , values = c(palette[3], darjeeling2_palette[2])) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # 5 by 4

freq_histogram / mutation_type_by_gene

fig2.plot <- freq.by.pos / (freq_histogram | mutation_type_by_gene | snv.by.copynum.week)


# ========================= Add information on antigenic sites / receptor binding sites ==================

antigenic_sites <- read.csv("data/reference/OPV2_AntigenicSites.csv") # amino acid sites. From Shaw, J. et al., JVI 2018

# general order is VP4, VP2, VP3, VP1. VP1-143 is polyprotien 721.
# VP4 is 69 amino acids
# VP2 is 271 amino acids
# VP3 is 238 amino acids
# VP1 is 301 amino acids
GetAbsolutePosition <- function(gene, position)
{
  if(gene == "VP2")
  {
    x = position + 69
  } 
  else if (gene == "VP3")
  {
    x = position + 69 + 271
  } 
  else if (gene == "VP1")
  {
    x = position + 69 + 271 + 238
  }
}

GetAbsolutePosition.v <- Vectorize(GetAbsolutePosition)
antigenic_sites_absolute <- mutate(antigenic_sites, polyprotein_position = GetAbsolutePosition.v(gene_antigenic, aa_position))
#write.csv(antigenic_sites_absolute, "data/reference/OPV2_AntigenicSites_Absolute.csv")

antigenic_sites_absolute_hypervariable <- filter(antigenic_sites_absolute, name == "hypervariable")
antigenic_sites_absolute_receptor <- filter(antigenic_sites_absolute, name == "receptorBinding")
antigenic_sites_absolute_antigenic <- filter(antigenic_sites_absolute, name %in% c("NAg1", "NAg2", "NAg3a", "NAg3b"))

# get amino acid sites as numbers
minor_withantigenic <- filter(minor, !is.na(class_factor))
minor_withantigenic <- mutate(minor_withantigenic, pos_string = as.character(AA_pos))
minor_withantigenic <- mutate(minor_withantigenic, pos_string = as.character(AA_pos))
minor_withantigenic$pos_string <- gsub("\\[", "", minor_withantigenic$pos_string)
minor_withantigenic$pos_string <- gsub("\\]", "", minor_withantigenic$pos_string)
minor_withantigenic %>% mutate(polyprotein_position = as.numeric(pos_string)) %>% select(-pos_string) -> minor_withantigenic

# attach antigenic site information
minor_withantigenic <- mutate(minor_withantigenic, in_hypervariable = ifelse(polyprotein_position %in% antigenic_sites_absolute_hypervariable$polyprotein_position, "In HVR", "Not in HVR."))
minor_withantigenic <- mutate(minor_withantigenic, in_receptorBinding = ifelse(polyprotein_position %in% antigenic_sites_absolute_receptor$polyprotein_position, "In RBS", "Not in RBS"))
minor_withantigenic <- mutate(minor_withantigenic, in_antigenic_site = ifelse(polyprotein_position %in% antigenic_sites_absolute_antigenic$polyprotein_position, "Antigenic", "Non-Antigenic"))
minor_withantigenic <- mutate(minor_withantigenic, in_special_site = ifelse(polyprotein_position %in% antigenic_sites_absolute$polyprotein_position, "Special", "Non-Special"))

minor_withantigenic_capsid <- filter(minor_withantigenic, gene %in% c("VP1", "VP2", "VP3", "VP4"))
freq.by.pos.antigenic <- ggplot(minor_withantigenic_capsid, aes(x = pos, y = freq.var, fill = class_factor, alpha = as.factor(in_special_site))) +
  geom_point(size = 2, shape = 21) +
  theme_bw() + 
  xlab("Genome Position (Capsid)") + 
  ylab("Frequency") +
  theme(panel.grid.minor = element_blank()) +
  #annotate("rect", xmin = 2470, xmax = 2577, ymin = 0, ymax = Inf, alpha = 0.5, fill = "grey") + # hypervariable area
  scale_fill_manual(name = "" , values = palette[c(4, 5)]) +
  scale_alpha_discrete(range = c(0.2, 1))

minor_withantigenic_capsid_no721 <- filter(minor_withantigenic_capsid, polyprotein_position != 721)
mutation_type_by_antigenic <- ggplot(minor_withantigenic_capsid_no721, aes(x = as.factor(in_antigenic_site), fill = class_factor)) + 
  geom_bar(position = "dodge") +
  xlab("Capsid Minor iSNV (excluding VP1-143)") + 
  ylab("Number of iSNV") +
  scale_fill_manual(name = "" , values = c(palette[3], darjeeling2_palette[2])) +
  theme(legend.position = c(0.5, 0.5)) + 
  theme_classic()

