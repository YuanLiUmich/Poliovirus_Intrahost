

### Project: Polio Mixing Experiment 2.0
### Purpose: Collapse output variant CSV by sequencing duplicate.

# ======================== Import packages, read in data =============================

library(tidyverse)

sampleInfo <- read.csv("./data/reference/MixingStudySamples2.csv")

variants_1000x <- read.csv("./data/raw/shiny2.verysens.onesided.1000.variants.csv")
variants_500x <- read.csv("./data/raw/shiny2.verysens.onesided.500.variants.csv")
variants_200x <- read.csv("./data/raw/shiny2.verysens.onesided.200.variants.csv")


# ======================== Collapse functions ============================

sift_dups <- function(df)
{
  if(nrow(df) > 2) stop("Too many mutations here")
  
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
  stopifnot(length(unique(df$group)) == 1)
  
  df %>% dplyr::group_by(mutation) %>% dplyr::summarize(found = length(mutation)) -> count_mutations
    
  stopifnot(max(count_mutations$found) < 3)
    
  count_mutations %>% dplyr::filter(found == 2) -> good_mut
    
  df %>% dplyr::filter(mutation %in% good_mut$mutation) -> good_var
    
  good_var %>% dplyr::group_by(mutation) %>% dplyr::do(sift_dups(.)) -> dups_good
    
  return(dplyr::ungroup(dups_good))
}

# ========================== Main ==============================

expectedVariants <- read.csv("./data/reference/MixingStudyExpectedVariants.csv")
#expectedVariants <- mutate(expectedVariants, mutation_old = mutation)
#expectedVariants <- mutate(expectedVariants, mutation = paste0("Sabin1_ref_", Sabin_nt, Position, WT_nt))
#write.csv(expectedVariants, "./data/reference/MixingStudyExpectedVariants.csv")
sampleInfo <- mutate(sampleInfo, Id = SampleNumber)

expectedVariantsLevel <- read.csv("./data/reference/MixingStudyExpectedVariants_ByLevel.csv")
#expectedVariantsLevel <- mutate(expectedVariantsLevel, mutation_old = mutation)
#expectedVariantsLevel <- mutate(expectedVariantsLevel, mutation = paste0("Sabin1_ref_", Sabin_nt, Position, WT_nt))
#write.csv(expectedVariantsLevel, "./data/reference/MixingStudyExpectedVariants_ByLevel.csv")

# Annotate variant files for each coverage group

# 1000x
variants_with_meta <- left_join(variants_1000x, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
notcollapsed_1000x <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
collapsed_1000x <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write.csv(notcollapsed_1000x, "./data/processed/shiny2.verysens.onesided.1000.notcollapsed.variants.csv")
write.csv(collapsed_1000x, "./data/processed/shiny2.verysens.onesided.1000.collapsed.variants.csv")

# 500x
variants_with_meta <- left_join(variants_500x, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
notcollapsed_500x <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
collapsed_500x <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write.csv(notcollapsed_500x, "./data/processed/shiny2.verysens.onesided.500.notcollapsed.variants.csv")
write.csv(collapsed_500x, "./data/processed/shiny2.verysens.onesided.500.collapsed.variants.csv")

# 200x
variants_with_meta <- left_join(variants_200x, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
notcollapsed_200x <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
collapsed_200x <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write.csv(notcollapsed_200x, "./data/processed/shiny2.verysens.onesided.200.notcollapsed.variants.csv")
write.csv(collapsed_200x, "./data/processed/shiny2.verysens.onesided.200.collapsed.variants.csv")

### Prepare coverage data ###

sampleInfo <- mutate(sampleInfo, Id = as.character(Id))

# downsampled to coverage level
coverage_1000x <- read.csv("./data/raw/shiny2.verysens.onesided.1000.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))
coverage_500x <- read.csv("./data/raw/shiny2.verysens.onesided.500.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))
coverage_200x <- read.csv("./data/raw/shiny2.verysens.onesided.200.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))

coverage_1000x_wmeta <- left_join(coverage_1000x, sampleInfo, by = "Id")
coverage_500x_wmeta <- left_join(coverage_500x, sampleInfo, by = "Id")
coverage_200x_wmeta <- left_join(coverage_200x, sampleInfo, by = "Id")

write_csv(coverage_1000x_wmeta, "./data/processed/shiny2.verysens.onesided.1000.coverage.csv")
write_csv(coverage_500x_wmeta, "./data/processed/shiny2.verysens.onesided.500.coverage.csv")
write_csv(coverage_200x_wmeta, "./data/processed/shiny2.verysens.onesided.200.coverage.csv")
