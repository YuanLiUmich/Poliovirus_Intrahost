

### Project: Poliovirus sequencing from Matlab
### Purpose: Find data from transmission pairs and compare the genetic diversity across pairs
### Working directory: Poliovirus_Intrahost

# ================================= Read in data, import packages ============================

library(tidyverse)
library(lubridate)
library(wesanderson)

palette <- wes_palette("Darjeeling1")
darjeeling2_palette <- wes_palette("Darjeeling2")

metadata <- read.csv("data/metadata/specimens_Sabin2Positive_CTbelow37_processedFinal_deID_runInfo_covInfo_week_final.csv", stringsAsFactors = FALSE)
variants <- read.csv("data/processed/all.concat.varQuality.cleaned.variants.csv")

# ================== Helper functions =========================

# Yields mutations that are polymorphic in the donor.
polish_freq <- function(df, relative, min_freq = 0)
{
  rel_quo <- enquo(relative)
  min_freq_quo <- enquo(min_freq)
  m = paste0("Requiring polymophic at least in ", quo_name(rel_quo))
  message(m)
  
  # We apply a freqeuncy cut off. it will be at least 0
  frequency_cut <- rlang::expr(UQ(rel_quo) > UQ(min_freq_quo))
  df <- dplyr::filter(df, UQ(frequency_cut))
  
  # We want polymorphic sites in the given sample. Due to frequency oddities the major allels may not
  # be exactly 1-min_freq i.e. we don't set any frequencies by hand. Minor frequencies are set by deepSNV
  # and this includes an estimated error rate. Major frequencies are set by allele counts.
  
  frequency_cut_minor <- rlang::expr(UQ(rel_quo) < 1)
  
  # How many alleles have frequencies in the choosen sample above the min and below 1 at each loci
  df %>% dplyr::group_by(chr, pos, anonSpecId1, anonSpecId2) %>% dplyr::summarize(alleles = length(which(UQ(frequency_cut) & UQ(frequency_cut_minor)))) -> counts
  
  df <- merge(df, counts)
  return(subset(df, alleles == 2, select = -c(alleles)))
}

equal_compare <- function(position)
{ 
  # Take in all the variants found in a sample pair at a given position and correct for differences in inference. 
  # Meaning the reference base was called here but only in one sample because in the other sample there wasn't a minor allele. 
  # Each position with a minor variant should have at least 2 variants, in both samples after this is done. A minor variant and the reference. 
  # If a position has no variants in either sample then it is excluded from this analysis. In this case only the reference base was present in both samples.  
  # If only major variants are present then the major allele is fixed and there are no other minor variants to infer.
  
  pos_sum_1 <- sum(position$freq1) # What is the sum of all variants at this site in the first sample
  pos_sum_2 <- sum(position$freq2) # What is the sum of all variants at this site in the second sample
  stopifnot(pos_sum_1 + pos_sum_2 > 0) # No variant at this position - something terrible happened somewhere
  if((pos_sum_1 > 0 & min(position$freq1) < 0.95 & pos_sum_1 == min(position$freq1))) warning(paste0("There is only 1 allele in sample", position$anonSpecId, ", its frequency is : ", position$freq1), ", it may be removed.\n")
  if((pos_sum_2 > 0 & min(position$freq2) < 0.95 & pos_sum_2 == min(position$freq2))) warning(paste0("There is only 1 allele in sample", position$anonSpecId, ", its frequency is : ", position$freq2), ", it may be removed.\n") # only a minor variant called
  
  if(pos_sum_1 == 0)
  { # there isn't a variant in the first sample. The reference is fixed there. 
    
    x <- which(position$ref == position$var) # Where are the infered calls
    stopifnot(length(x) < 2) # should be 1 or 0.
    if(length(x) == 1)
    { # The reference has been infered in the second sample but there are no varaints in the first so it was not infered
      position$freq1[x] <- 1 # The reference base is fixed here
    } else if(length(x) == 0)
    { # the reference was not infered in the second sample and no variants here. Add the reference now. There must have been another base fixed in the second sample.
      extra <- dplyr::mutate(position, var = ref, mutation = paste0(chr, "_", ref, pos, var), freq1 = 1, freq2 = 0) # adding a line with the reference base fixed in sample 1 but not sample 2
      position <- rbind(position, extra)
    }
    # Do it again if there are no varaints in the second sample.
  } else if(pos_sum_2 == 0)
  { # there isn't a variant in the second sample. The reference is fixed.
    x <- which(position$ref == position$var)
    stopifnot(length(x) < 2) # should be 1 or 0
    if(length(x) == 1)
    { # The reference has been infered in the first sample but there are no varaints in the second so it was not infered
      position$freq2[x] <- 1
    } else if(length(x) == 0)
    { # the reference was not infered in the first sample and no variants here. Add the reference now
      extra <- dplyr::mutate(position, var = ref, mutation = paste0(chr, "_", ref, pos, var), freq1 = 0, freq2 = 1) # adding a line with the reference base fixed in sample 2 but not sample 1
      position <- rbind(position, extra)
    }
  }
  
  if(nrow(position) > 1)
  { 
    return(position)
  }else
  {
    return(position[F,])
  }
}

get_freqs <- function(pairs, snv)  
{
  # Input: List of SNV calls and a dataframe of one ID pair.
  # Output: Dataframe comparing each SNV's frequency in both samples. Only those where the frequency changes from one sample to another.
  
  snv <- subset(snv, anonSpecId %in% pairs, select = c(anonSpecId, mutation, chr, pos, ref, var, freq.var))
  
  if(nrow(snv) > 0)
  { # There are mutations.
    # We only compare samples from the same season and strain so in each case the reference base is the same
    
    mut_table <- tidyr::spread(snv, anonSpecId, freq.var, fill = 0) # a data frame with mutation down the first row and then frequency in either sample in the next 2.
    
    mut_table$anonSpecId1 <- pairs[1] # add column with first sample ID
    mut_table$anonSpecId2 <- pairs[2] # add column with second sample ID
    names(mut_table)[which(names(mut_table) == as.character(pairs[1]))] <- 'freq1' # rename this column as the frequency in the first sample
    names(mut_table)[which(names(mut_table) == as.character(pairs[2]))] <- 'freq2'
    
    # This function can only be run on samples that qualified for variant identification.
    # If no variants were found in the sample then the SPECID will be missing from mut_table column and so
    # freq1 or freq2 will be missing since nothing was found we set that to 0 here and add the column.
    # equal compare will replace these cases with the reference at 1.
    if(!('freq1' %in% names(mut_table)))
    {
      mut_table$freq1 <- 0
    }
    if(!('freq2' %in% names(mut_table)))
    {
      mut_table$freq2 <- 0
    }
    mut_table <- dplyr::select(mut_table, mutation, chr, pos, ref, var, freq1, freq2, anonSpecId1, anonSpecId2)
    mut_table <- mut_table[order(mut_table$chr, mut_table$pos),]
    
    # Fill in differences based on inferred. In this case we are left with only sites that are polymorphic in one or between the 2 samples
    mut_table %>% dplyr::group_by(chr, pos) %>% dplyr::do(equal_compare(.)) -> all.freq
    
    if(nrow(all.freq) > 0)
    { # Some differences exist
      return(all.freq)
    }
    else
    { # some snv were initially present (before frequency filter) but after taking into account differences in inference the samples are the same.
      x <- tibble(mutation = NA, freq1 = NA, freq2 = NA, anonSpecId1 = NA, anonSpecId2 = NA, chr = NA, ref = NA, "pos" = NA, var = NA)
      return(x[F,])
    }
  }
  else
  { # No variants found in either sample
    x <- tibble(mutation = NA, freq1 = NA, freq2 = NA, anonSpecId1 = NA, anonSpecId2 = NA, chr = NA, ref = NA, "pos" = NA, var = NA)
    return(x[F,])
  }
}

# ====================================== Find household pairs ================================

metadata %>%
  filter(any_data == TRUE & collectionWeekRelative > 0) %>% 
  group_by(houseId) %>%
  summarize(house_pairs = length(unique(sid))) %>%
  filter(house_pairs > 1) -> meta_households

metadata_anydata <- filter(metadata, any_data == TRUE, houseId %in% meta_households$houseId & collectionWeekRelative > 0)
metadata_anydata_donor <- filter(metadata_anydata, transmissionAcquired == 0 & mOPV2 == 1)
metadata_anydata_recipient <- filter(metadata_anydata, transmissionAcquired == 1 & mOPV2 == 0)

meta_households <- mutate(meta_households, donor = metadata_anydata_donor$sid[match(houseId, metadata_anydata_donor$houseId)])
meta_households <- mutate(meta_households, recipient = metadata_anydata_recipient$sid[match(houseId, metadata_anydata_recipient$houseId)])
meta_households <- filter(meta_households, !is.na(recipient))
metadata_inpairs <- filter(metadata_anydata, houseId %in% meta_households$houseId) %>% select(houseId, sid, SPECDT, anonSpecId, any_data, Sab2_variant_qualified) %>% arrange(-houseId)

meta_households <- mutate(meta_households, recipient_sample = metadata_inpairs$anonSpecId[match(recipient, metadata_inpairs$sid)])
meta_households <- mutate(meta_households, recipient_sampledate = metadata_inpairs$SPECDT[match(recipient, metadata_inpairs$sid)])

metadata_inpairs <- mutate(metadata_inpairs, recipient_sampledate = meta_households$recipient_sampledate[match(houseId, meta_households$houseId)])
metadata_inpairs %>%
  filter(nchar(sid) < 4) %>%
  mutate(time_between_pairs = ymd(recipient_sampledate) - ymd(SPECDT)) %>%
  group_by(sid) %>%
  filter(time_between_pairs > -8) %>%
  filter(time_between_pairs == min(time_between_pairs)) -> metadata_inpairs_donorsamples

meta_households <- mutate(meta_households, donor_sample = metadata_inpairs_donorsamples$anonSpecId[match(donor, metadata_inpairs_donorsamples$sid)])
meta_households <- mutate(meta_households, donor_sampledate = metadata_inpairs_donorsamples$SPECDT[match(donor, metadata_inpairs_donorsamples$sid)])
meta_households <- filter(meta_households, !is.na(donor_sample))
meta_households <- mutate(meta_households, time_difference = ymd(recipient_sampledate) - ymd(donor_sampledate)) # six household pairs with sampling times within a week.
meta_households <- mutate(meta_households, donor_sample = as.character(donor_sample), recipient_sample = as.character(recipient_sample))
meta_households <- mutate(meta_households, donor_mOPV2_date = metadata$vaccination_date[match(donor, metadata$sid)])

meta_households <- filter(meta_households, houseId %in% c(115, 171, 702, 927))

write.csv(meta_households, "data/processed/meta_households.csv")

# ================================== Get diversity across pairs =================================

transmission_variants <- plyr::adply(meta_households, 1, function(x) {get_freqs(c(x$donor_sample, x$recipient_sample), variants)})

# Plot the diversity across transmission pairs
transmission_variants <- filter(transmission_variants, pos != 2909)
tv.plot <- ggplot(transmission_variants, aes(x = freq1, y = freq2)) + 
  geom_point() + 
  xlab("Frequency in donor") + 
  ylab("Frequency in recipient") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # PDF 5 by 4

transmission_variants_donorPolymorphic <- polish_freq(transmission_variants, freq1, 0.05)
tv.donor.plot <- ggplot(transmission_variants_donorPolymorphic, aes(x = freq1, y = freq2)) + 
  geom_point(size = 1.5, aes(color = as.factor(houseId))) + 
  xlab("Frequency in donor") + 
  ylab("Frequency in recipient") + 
  theme_classic() +
  facet_wrap(~houseId)

write.csv(transmission_variants_donorPolymorphic, "data/processed/variants_donor_polymorphic.csv")

# Get minor variants from donor 702
transmission_variants_donorPolymorphic_702 <- filter(transmission_variants_donorPolymorphic, houseId == 702)
unique(transmission_variants_donorPolymorphic_702$pos) # sites to investigate
transmission_variants_donorPolymorphic_702_minor <- filter(transmission_variants_donorPolymorphic_702, freq1 < 0.5)
transmission_variants_donorPolymorphic_702_minor <- mutate(transmission_variants_donorPolymorphic_702_minor, sequencingID_1 = metadata$sequencingID_1[match(anonSpecId1, metadata$anonSpecId)])
write.csv(transmission_variants_donorPolymorphic_702_minor, "data/processed/variants_donor_polymorphic_702_minor.csv")
