
### Project: Poliovirus sequencing from Matlab
### Purpose: Prepare the variant file for processing. Read in the CSV for each sample for which we have OPV2 data, make sure the columns are ordered correctly, and write.
### Working directory: Poliovirus_Intrahost

# ============================= Read in metadata =========================

library(tidyverse)

samples_to_use <- read.csv("data/processed/variant_filenames_with_data.csv", stringsAsFactors = FALSE) # These are the samples that have some level of data

# ============================ Parse variant data files ==================

variants_all <- data.frame()
read.csv(as.character(samples_to_use[1,]$filename_original)) %>% select(-X, -'Unnamed..0.1') -> first_file
variants_all <- rbind(variants_all, first_file)
correctColumns <- colnames(variants_all)

for(index in seq(2, nrow(samples_to_use), by = 1))
{
  message <- paste0(as.character(index), " of ", as.character(nrow(samples_to_use)))
  print(message)
  filename <- as.character(samples_to_use[index,]$filename_original)
  read.csv(filename) %>% select(-X, -'Unnamed..0.1') -> next_file
  if(is.null(next_file))
  {
    next
  }
  if(!all(colnames(next_file) == correctColumns))
  {
    next_file_rearr <- select(next_file, Id, MapQ, Phred, Read_pos, everything())
    variants_all <- rbind(variants_all, next_file_rearr)
  }
  else
  {
    variants_all <- rbind(variants_all, next_file)
  }
}

write.csv(variants_all, "data/processed/all.concat.variants.csv")

# ============================ Same thing, but with coding info relative to OPV2 reference seq, not sample consensus ==================

samples_to_use <- mutate(samples_to_use, filename_control = paste0("data/raw/", platen, "ctl/Final_variants/", seq_id, ".removed.ref.sum.filtered.AA.csv"))

variants_all <- data.frame()
read.csv(as.character(samples_to_use[1,]$filename_control)) %>% select(-X, -'Unnamed..0.1') -> first_file
variants_all <- rbind(variants_all, first_file)
correctColumns <- colnames(variants_all)

for(index in seq(2, nrow(samples_to_use), by = 1))
{
  message <- paste0(as.character(index), " of ", as.character(nrow(samples_to_use)))
  print(message)
  filename <- as.character(samples_to_use[index,]$filename_control)
  read.csv(filename) %>% select(-X, -'Unnamed..0.1') -> next_file
  if(is.null(next_file))
  {
    next
  }
  if(!all(colnames(next_file) == correctColumns))
  {
    next_file_rearr <- select(next_file, Id, MapQ, Phred, Read_pos, everything())
    variants_all <- rbind(variants_all, next_file_rearr)
  }
  else
  {
    variants_all <- rbind(variants_all, next_file)
  }
}

write.csv(variants_all, "data/processed/all.concat.variants.ctl.csv")
  