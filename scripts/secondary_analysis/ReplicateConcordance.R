
### Project: Poliovirus sequencing from Matlab
### Purpose: Explore within-host diversity by iSNVs
### Working directory: Poliovirus_Intrahost

# ================== Import packages, read in data =================

library(tidyverse)
library(wesanderson)

palette <- wes_palette("Darjeeling1")
palette2 <- wes_palette("Darjeeling2")

metadata <- read.csv("data/metadata/specimens_Sabin2Positive_CTbelow37_processedFinal_deID_runInfo_covInfo_week_final.csv")
variants_notcollapsed <- read.csv("data/processed/all.concat.varQuality.variants.csv")

# ============================= Function to find primer binding sites =======================

IsPositionInPrimerSite <- function(pos)
{
  # Primer binding sites
  seg1_fwd_start <- 98
  seg1_fwd_stop  <- 114
  seg1_rev_start <- 2419
  seg1_rev_stop  <- 2436
  seg2_fwd_start <- 1472
  seg2_fwd_stop  <- 1491
  seg2_rev_start <- 4129
  seg2_rev_stop  <- 4151
  seg3_fwd_start <- 3214
  seg3_fwd_stop  <- 3230
  seg3_rev_start <- 5842
  seg3_rev_stop  <- 5861
  seg4_fwd_start <- 4966
  seg4_fwd_stop  <- 4982
  seg4_rev_start <- 7384
  seg4_rev_stop  <- 7400
  
  if((pos >= seg1_fwd_start & pos <= seg1_fwd_stop) | 
     (pos >= seg2_fwd_start & pos <= seg2_fwd_stop) | 
     (pos >= seg3_fwd_start & pos <= seg3_fwd_stop) | 
     (pos >= seg4_fwd_start & pos <= seg4_fwd_stop) |
     (pos >= seg1_rev_start & pos <= seg1_rev_stop) |
     (pos >= seg2_rev_start & pos <= seg2_rev_stop) |
     (pos >= seg3_rev_start & pos <= seg3_rev_stop) |
     (pos >= seg4_rev_start & pos <= seg4_rev_stop))
  {
    x = TRUE
  }
  else
  {
    x = FALSE
  }
}

IsPositionInPrimerSite.v <- Vectorize(IsPositionInPrimerSite)

# ======================== Do analysis ===================

metadata_quality_dups <- filter(metadata, Sab2_variant_qualified == TRUE & reps == 2)
variants_notcollapsed_quality <- filter(variants_notcollapsed, anonSpecId %in% metadata_quality_dups$anonSpecId)
variants_notcollapsed_minor <- filter(variants_notcollapsed_quality, freq.var > 0.05 & freq.var < 0.5)
isnv_not_collapsed_rep1 <- filter(variants_notcollapsed_minor, Id == sequencingID_1)
isnv_not_collapsed_rep2 <- filter(variants_notcollapsed_minor, Id == sequencingID_2)
isnv_not_collapsed_rep1 <- mutate(isnv_not_collapsed_rep1, mut_id = paste0(mutation, "_", anonSpecId))
isnv_not_collapsed_rep2 <- mutate(isnv_not_collapsed_rep2, mut_id = paste0(mutation, "_", anonSpecId))
replicates_merged <- merge(isnv_not_collapsed_rep1, isnv_not_collapsed_rep2, by = "mut_id")

replicate_concordance <- ggplot(data = replicates_merged, aes(x = freq.var.x, y = freq.var.y)) +
  geom_point(aes(fill = as.factor(floor(log10(S2copyNumberPerGram.x/1000)))),  shape = 21, size = 2.5) + 
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  scale_y_continuous(limits = c(0, 0.5)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  xlab("Frequency in replicate 1") + 
  ylab("Frequency in replicate 2") +
  theme_bw() +
  theme(legend.position = "right") +
  scale_fill_manual(name = "Genomes per microliter" , values = palette[c(5,3)])

replicate_lm <- lm(data = replicates_merged, freq.var.y ~ freq.var.x)
summary(replicate_lm)
cor.test(replicates_merged$freq.var.y, replicates_merged$freq.var.x)

# Look at effect of PBS mutations:
variants_notcollapsed_quality <- mutate(variants_notcollapsed_quality, in_primer_site = IsPositionInPrimerSite.v(pos))
variants_pbs <- filter(variants_notcollapsed_quality, in_primer_site == TRUE)
variants_pbs <- filter(variants_pbs, freq.var > 0.05)
variants_pbs <- monomorphic(variants_pbs, Id, pos)
variants_pbs <- filter(variants_pbs, ref != var) # Everything that is different than the reference. THIS may need to change. Get anything that is different than the PRIMER sequence.
variants_pbs <- mutate(variants_pbs, coverage = cov.tst.fw + cov.tst.bw) %>% filter(coverage > 10)

variants_pbs_seg1 <- filter(variants_pbs, pos < 1000)
variants_pbs_seg3 <- filter(variants_pbs, pos > 1000)

replicates_merged <- mutate(replicates_merged, pbs_affected = "No PBS mutations")
replicates_merged <- mutate(replicates_merged, pbs_affected = ifelse(anonSpecId.x %in% variants_pbs_seg3$anonSpecId & pos.x %in% seq(3230, 5842, by = 1), "PBS mutation found", pbs_affected))
replicates_merged <- mutate(replicates_merged, pbs_affected = ifelse(anonSpecId.x %in% variants_pbs_seg1$anonSpecId & pos.x %in% seq(114, 2419, by = 1), "PBS mutation found", pbs_affected))

replicate_concordance_pbs <- ggplot(data = replicates_merged, aes(x = freq.var.x, y = freq.var.y)) +
  geom_point(aes(fill = as.factor(pbs_affected)),  shape = 21, size = 2.5) + 
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  scale_y_continuous(limits = c(0, 0.5)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  xlab("Frequency in replicate 1") + 
  ylab("Frequency in replicate 2") +
  theme_bw() +
  theme(legend.position = "right") +
  scale_fill_manual(name = "" , values = c(palette2[2], palette[3])) #+
  #facet_wrap(~pbs_affected) +
  #geom_smooth(method = 'lm', formula = y ~ x) # pdf, 6 by 4

# Adj. R-sqaured for each group
replicates_merged_pbs <- filter(replicates_merged, pbs_affected == "PBS mutation found")
mutations_lm <- lm(data = replicates_merged_pbs, freq.var.y ~ freq.var.x)
summary(mutations_lm)

replicates_merged_nopbs <- filter(replicates_merged, pbs_affected == "No PBS mutations")
no_mutations_lm <- lm(data = replicates_merged_nopbs, freq.var.y ~ freq.var.x)
summary(no_mutations_lm)
