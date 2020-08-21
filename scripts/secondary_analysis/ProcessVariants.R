

### Project: Poliovirus sequencing from Matlab
### Purpose: Processing for iSNVs
### Working directory: Poliovirus_Intrahost

# ================================= Read in data, import packages ============================

library(tidyverse)
library(lubridate)

metadata <- read.csv("data/metadata/specimens_Sabin2Positive_CTbelow37_processedFinal_deID_runInfo_covInfo_week.csv")

# Read in OPV2 variants
variants_all_concat <- read.csv("data/processed/all.concat.variants.csv")
variants_all_concat_Sab2 <- filter(variants_all_concat, chr == "Sabin2ref")

# ==================================== Some metadata processing ============================

# Link house ID and vaccination date into metadata.
read.csv("data/metadata/Matlab_HouseID_VaxDate.csv", header = TRUE, stringsAsFactors = FALSE, na.strings = c("","NA")) %>% 
  select(-X, -X.1, -X.2, -X.3, -X.4, -X.5, -X.6, -X.7) -> house_and_vaxdate
house_and_vaxdate_00 <- filter(house_and_vaxdate, !is.na(date_CHALL_DAY00)) %>% mutate(vaccination_date = mdy(date_CHALL_DAY00), estimated_vaxdate = FALSE)
house_and_vaxdate_07 <- filter(house_and_vaxdate, !is.na(date_CHALL_DAY07) & is.na(date_CHALL_DAY00)) %>% mutate(vaccination_date = mdy(date_CHALL_DAY07) - 7, estimated_vaxdate = TRUE)
house_and_vaxdate_14 <- filter(house_and_vaxdate, !is.na(date_CHALL_DAY14) & is.na(date_CHALL_DAY00) & is.na(date_CHALL_DAY07)) %>% mutate(vaccination_date = mdy(date_CHALL_DAY14) - 14, estimated_vaxdate = TRUE)
house_and_vaxdate_21 <- filter(house_and_vaxdate, !is.na(date_CHALL_DAY21) & is.na(date_CHALL_DAY00) & is.na(date_CHALL_DAY07) & is.na(date_CHALL_DAY14)) %>% mutate(vaccination_date = mdy(date_CHALL_DAY21) - 21, estimated_vaxdate = TRUE)
rbind(house_and_vaxdate_00, house_and_vaxdate_07, house_and_vaxdate_14, house_and_vaxdate_21) %>% 
  select(-date_CHALL_DAY00, -date_CHALL_DAY07, -date_CHALL_DAY14, -date_CHALL_DAY21, -date_CHALL_DAY28) -> house_and_vaxdate_final
metadata <- left_join(metadata, house_and_vaxdate_final, by = "sid")
metadata <- mutate(metadata, days_since_vaccination = ymd(SPECDT) - vaccination_date)
write.csv(metadata, "data/metadata/specimens_Sabin2Positive_CTbelow37_processedFinal_deID_runInfo_covInfo_week_final.csv")


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

# =================================== Filter by coverage: variant-level bins ===================================

# Get all the positions that have variant-quality coverage.

variant_quality_regions <- read.csv("data/processed/coverage_variant_quality_regions.csv")
ids <- list(unique(variants_all_concat_Sab2$Id))
quality_positions <- list()
for(i in 1:length(ids[[1]]))
{
  id <- ids[[1]][i]
  group <- filter(variant_quality_regions, Id == id)
  if(nrow(group) == 0)
  {
    next
  }
  positions <- c()
  for(r in 1:nrow(group))
  {
    row <- group[r,]
    new_positions <- seq(row$start, row$end, 1)
    positions <- c(positions, new_positions)
  }
  quality_positions[[id]] <- positions
}

# Now, use those variant-quality positions to filter the variant dataframe. This loop takes a while.
variants_all_concat_Sab2_forVar <- mutate(variants_all_concat_Sab2, Id = as.character(Id))
variants_all_concat_Sab2_forVar <- filter(variants_all_concat_Sab2_forVar, Id %in% ids[[1]])
variants_all_concat_Sab2_forVar_variantQuality <- data.frame()
for(i in 1:length(ids[[1]]))
{
  id <- ids[[1]][i]
  print(paste0("Working on ID: ", id, ". ", as.character(i), " out of ", as.character(length(ids[[1]]))))
  group <- filter(variants_all_concat_Sab2_forVar, Id == id)
  good_positions <- quality_positions[[id]]
  group_quality <- filter(group, pos %in% good_positions)
  variants_all_concat_Sab2_forVar_variantQuality <- rbind(variants_all_concat_Sab2_forVar_variantQuality, group_quality)
}

length(unique(variants_all_concat_Sab2_forVar_variantQuality$Id))

# =========================== Join variant data and metadata ===============================

seq_locations <- read_csv("data/metadata/SeqID_to_location.csv")
seq_locations %>% select(SampleNumber, SampleName) %>% rename(Id = SampleNumber, location = SampleName) -> seq_locations
metadata %>% select(anonSpecId, location1, location2) %>% gather(v, value, location1:location2) %>% filter(!is.na(value)) %>% select(-v) %>% rename(location = value) -> meta_samples
merge(meta_samples, seq_locations, by = "location") %>% select(-location) -> seq_info
variants_withID <- merge(variants_all_concat_Sab2_forVar_variantQuality, seq_info, by = "Id") # variant dataframe from above
variants_with_meta <- merge(variants_withID, metadata, by = "anonSpecId")

# =========================== Add in gene info ================================

# coordinates from data/reference/AY184220.1_genomeFeatures.csv

GetGeneInfo <- function(pos)
{
  if(pos >= 1 & pos <= 747){
    x="5UTR"} else if (pos >= 748 & pos <= 954){
      x="VP4"} else if (pos >= 955 & pos <= 1767){
        x="VP2"} else if (pos >= 1768 & pos <= 2481){
          x="VP3"} else if (pos >= 2482 & pos <= 3384){
            x="VP1"} else if (pos >= 3385 & pos <= 3831){
              x="2A"} else if (pos >= 3832 & pos <= 4122){
                x="2B"} else if (pos >= 4123 & pos <= 5109){
                  x="2C"} else if (pos >= 5110 & pos <= 5370){
                    x="3A"} else if (pos >= 5371 & pos <= 5436){
                      x="3B"} else if (pos >= 5437 & pos <= 5985){
                        x="3C"} else if (pos >= 5986 & pos <= 7368){
                          x="3D"} else if (pos >= 7369 & pos <= 7440){
                            x="3UTR"}
}

GetGeneInfo.v <- Vectorize(GetGeneInfo)

mutate(variants_with_meta, gene = GetGeneInfo.v(pos)) -> variants_with_meta_genes

write.csv(variants_with_meta_genes, "data/processed/all.concat.varQuality.variants.csv")


# =========================== Function for getting iSNV found in both sequencing replicates ===============================

sift_duplicates <- function(df)
{
  if(nrow(df) > 2) stop("Too many mutations here!")
  
  df <- dplyr::mutate(df, coverage = cov.tst.bw + cov.tst.fw)
  higher_qual <- subset(df, coverage == max(df$coverage))
  if(nrow(higher_qual) > 1)
  { 
    higher_qual <- higher_qual[1,]
  }
  return(higher_qual)
}

collapse_localcov <- function(df)
{
  stopifnot(length(unique(df$anonSpecId)) == 1)
  replicates <- unique(df$reps)
  
  if(replicates == 2)
  {
    df %>% dplyr::group_by(mutation) %>% dplyr::summarize(found = length(mutation)) -> count_mutations
    stopifnot(max(count_mutations$found) < 3)
    count_mutations %>% dplyr::filter(found == 2) -> good_mut
    df %>% dplyr::filter(mutation %in% good_mut$mutation) -> good_var
    good_var %>% dplyr::group_by(mutation) %>% dplyr::do(sift_duplicates(.)) -> dups_good
    return(dplyr::ungroup(dups_good))
    
  } else if(replicates == 1)
  {
    return(df)
  }
}

# =========================== Filtering: helper functions ===============================

# After filtering false positives, there could be a monomorphic site that has frequency <1. This sets those sites to 1.
monomorphic <- function(df,...)
{
  group_var <- rlang::quos(...)
  df %>% dplyr::group_by(!!!group_var) %>% dplyr::summarize(alleles = length(unique(var))) -> df_alleles
  df <- dplyr::left_join(df, df_alleles)
  stopifnot(any(!is.na(df$alleles)))
  df$freq.var[df$alleles == 1] <- 1
  df %>% dplyr::select(-alleles) -> df
  return(df)
}

# ========================= Processing ==========================

#variants_with_meta_genes <- read.csv("data/processed/all.concat.varQuality.variants.csv")

# Current variant filtering:
variants_with_meta_genes %>% group_by(anonSpecId) %>% do(collapse_localcov(.)) -> quality_var
quality_var$class_factor = NA
quality_var$class_factor[grep("Noncoding", quality_var$Class)] <- "Noncoding"
quality_var$class_factor[grep('Syn', quality_var$Class)] <- "Synonymous"
quality_var$class_factor[grep('Nonsyn', quality_var$Class)] <- "Nonsynonymous"
quality_var$class_factor <- as.factor(quality_var$class_factor)
quality_var <- filter(quality_var, freq.var > 0.05)
quality_var_monomorphic <- monomorphic(quality_var, Id, pos)
quality_var_monomorphic <- mutate(quality_var_monomorphic, in_primer_site = IsPositionInPrimerSite.v(pos))

write.csv(quality_var_monomorphic, "data/processed/all.concat.varQuality.cleaned.variants.csv")

