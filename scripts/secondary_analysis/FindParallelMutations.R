
### Project: Poliovirus sequencing from Matlab
### Purpose: Find mutations shared across multiple people
### Working directory: Poliovirus_Intrahost

# ================== Import packages, read in data =================

library(tidyverse)
library(wesanderson)
library(viridis)

palette <- wes_palette("Darjeeling1")
darjeeling2_palette <- wes_palette("Darjeeling2")
aquatic_palette <- wes_palette("Zissou1")

metadata <- read.csv("data/metadata/specimens_Sabin2Positive_CTbelow37_processedFinal_deID_runInfo_covInfo_week_final.csv")
variants <- read.csv("data/processed/all.concat.varQuality.cleaned.variants.csv")
#variants <- read.csv("data/processed/all.concat.varQuality.cleaned.variants.ctl.csv") # to see coding changes relative to OPV2 reference

# =================================== Find mutations ==========================

variants_mOPV2_var <- filter(variants, ref != var & in_primer_site == FALSE & mOPV2 == 1 & transmissionAcquired == 0 & collectionWeekRelative > 0 & Sab2_variant_qualified == TRUE)

# Attach antigenic site information
antigenic_sites_absolute <- read.csv("data/reference/OPV2_AntigenicSites_Absolute.csv")
antigenic_sites_absolute_antigenic <- filter(antigenic_sites_absolute, name %in% c("NAg1", "NAg2", "NAg3a", "NAg3b"))
variants_mOPV2_var <- filter(variants_mOPV2_var, !is.na(class_factor))
variants_mOPV2_var <- mutate(variants_mOPV2_var, pos_string = as.character(AA_pos))
variants_mOPV2_var <- mutate(variants_mOPV2_var, pos_string = as.character(AA_pos))
variants_mOPV2_var$pos_string <- gsub("\\[", "", variants_mOPV2_var$pos_string)
variants_mOPV2_var$pos_string <- gsub("\\]", "", variants_mOPV2_var$pos_string)
variants_mOPV2_var %>% mutate(polyprotein_position = as.numeric(pos_string)) %>% select(-pos_string) -> variants_mOPV2_var
variants_mOPV2_var <- mutate(variants_mOPV2_var, in_antigenic_site = ifelse(polyprotein_position %in% antigenic_sites_absolute_antigenic$polyprotein_position, "Antigenic", "Non-Antigenic"))

# Find shared mutations
variants_mOPV2_var %>% group_by(mutation) %>% dplyr::summarize(counts = length(unique(sid))) -> var_count
var_count_multiple <- filter(var_count, counts > 1)

# Assign group
var_count_multiple <- mutate(var_count_multiple, status = "P2/P3")
var_count_multiple <- mutate(var_count_multiple, status = ifelse(mutation %in% filter(variants_mOPV2_var, gene %in% c("VP1", "VP2", "VP3", "VP4"))$mutation, "Capsid", status))
var_count_multiple <- mutate(var_count_multiple, status = ifelse(mutation %in% filter(variants_mOPV2_var, in_antigenic_site == "Antigenic")$mutation, "Capsid - Antigenic", status)) # get whether it is noncoding, antigenic, or gatekeeper
var_count_multiple <- mutate(var_count_multiple, status = ifelse(mutation %in% filter(variants_mOPV2_var, class_factor == "Noncoding")$mutation, "Noncoding", status))
var_count_multiple <- mutate(var_count_multiple, status = ifelse(mutation %in% c("Sabin2ref_A481G", "Sabin2ref_T398C", "Sabin2ref_T2909C", "Sabin2ref_T2909A", "Sabin2ref_A2908G", "Sabin2ref_A2908T"), "Gatekeeper", status))

# ================================== Plot ======================================

var_count_five <- filter(var_count_multiple, counts >= 5)

parallel_palette <- c(darjeeling2_palette[2], palette[4], palette[1], palette[5], "black")

multiple.hist <- ggplot(var_count_five, aes(x = reorder(mutation, -counts), y = counts)) + 
  geom_bar(stat = "identity", aes(fill = status)) + 
  #scale_x_discrete(breaks = NULL) + 
  theme_classic() + 
  xlab("Mutation") + 
  ylab("Number of Individuals") + 
  scale_fill_manual(name = "Mutation Type", values = parallel_palette) + 
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) + 
  ggtitle("Mutations Found in Five or More Individuals")


var_count_multiple <- mutate(var_count_multiple, group = as.character(counts))
var_count_multiple <- mutate(var_count_multiple, group = ifelse(counts > 4, ">4", group))
var_count_multiple$group <- factor(var_count_multiple$group, levels = c("2", "3", "4", ">4"))

shared.mutations.bar <- ggplot(var_count_multiple, aes(group, y = ..count.., fill = status)) + 
  geom_bar() +
  scale_fill_manual(name = "Mutation Type", values = parallel_palette) + 
  ylab("Number of Mutations") + 
  xlab("Number of Individuals") + 
  theme_bw() # 5 by 4

write.csv(var_count_multiple, "data/processed/var_count_multiple.csv")

# =================================== For the significant parallel mutations, which arms are the individuals in? ============================

variants_parallel <- filter(variants_mOPV2_var, mutation %in% var_count_five$mutation)

variants_parallel_OPV <- filter(variants_parallel, arm == "tOPV")
length(unique(variants_parallel_OPV$sid)) # 11
variants_parallel_IPV <- filter(variants_parallel, arm != "tOPV")
length(unique(variants_parallel_IPV$sid)) # 67

getNumIndividuals <- function(df, mut, arms)
{
  df_arm <- filter(df, arm %in% arms & mutation == mut)
  numIndivs <- length(unique(df_arm$sid))
  return(numIndivs)
}

arm_data <- data.frame(mutation = var_count_five$mutation, count_IPV = rep(NA, nrow(var_count_five)), count_OPV = rep(NA, nrow(var_count_five)))
for(mut in arm_data$mutation)
{
  count_IPV = getNumIndividuals(variants_parallel, mut, c("IPVx1", "IPVx2"))
  #print(count_IPV)
  count_OPV = getNumIndividuals(variants_parallel, mut, c("tOPV"))
  #print(count_OPV)
  arm_data$count_IPV[match(mut, arm_data$mutation)] <- count_IPV
  arm_data$count_OPV[match(mut, arm_data$mutation)] <- count_OPV
}

arm_data <- left_join(arm_data, select(var_count_multiple, mutation, status), by = "mutation")
arm_data <- mutate(arm_data, counts = count_IPV + count_OPV) %>% filter(mutation != "Sabin2ref_T7401G")

write.csv(arm_data, "data/processed/var_count_multiple_arms.csv") # No major difference. Can do chi-squared for some.


# Do the synonymous mutations disrupt CpG or UpA?

variants_mOPV2 <- filter(variants, in_primer_site == FALSE & mOPV2 == 1 & transmissionAcquired == 0 & collectionWeekRelative > 0 & Sab2_variant_qualified == TRUE)
samples <- filter(variants_mOPV2_var, mutation == "Sabin2ref_T6693A")$anonSpecId
variants_mOPV2_mut <- filter(variants_mOPV2, anonSpecId %in% samples & pos %in% c(6692, 6693, 6694))

# G6084U. No.
# U882C. No.
# U4374C. Yes. Disrupts a UpA.
# C2580U. Yes. Disrupts a CpG.
# U1641C. Yes. Disrupts a UpA.
# U5811A. Dunno. Disrupts a UpA, but creates another one.
# U6693A. No.
