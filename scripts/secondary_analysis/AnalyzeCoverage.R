
### Project: Poliovirus sequencing from Matlab
### Purpose: Assess coverage on samples; depth and evenness across the genome.
### Working directory: Poliovirus_Intrahost


# ================== Read in data, import packages =========================

library(tidyverse)
library(wesanderson)
library(lubridate)

coverage_files <- c("data/raw/plate1/all.coverage.csv",
                    "data/raw/plate2/all.coverage.csv",
                    "data/raw/plate3/all.coverage.csv",
                    "data/raw/plate4/all.coverage.csv",
                    "data/raw/plate5/all.coverage.csv",
                    "data/raw/plate6/all.coverage.csv",
                    "data/raw/plate7/all.coverage.csv",
                    "data/raw/plate8/all.coverage.csv")
coverage_df_all <- data.frame(chr = NA, chr.pos = NA, coverage = NA, concat.pos = NA, Id = NA, plate = NA)
for(file in coverage_files)
{
  data <- read.csv(file)
  plate_num <- str_split(file, "/")[[1]][3]
  data <- mutate(data, plate = plate_num)
  coverage_df_all <- rbind(coverage_df_all, data)
}

# ==================== Get metadata ==========================

metadata <- read_csv("data/metadata/specimens_Sabin2Positive_CTbelow37_processedFinal_deID.csv")
sequencing_IDs <- read_csv("data/metadata/SeqID_to_location.csv")
sequencing_IDs <- select(sequencing_IDs, -PoolNumber)
sequencing_IDs <- rename(sequencing_IDs, location1 = SampleName)
sequencing_IDs <- rename(sequencing_IDs, SequencingID = SampleNumber)
sequencing_IDs <- mutate(sequencing_IDs, location2 = location1)
sequencing_IDs_1 <- select(sequencing_IDs, -location2)
sequencing_IDs_2 <- select(sequencing_IDs, -location1)

metadata <- left_join(metadata, sequencing_IDs_1, by = "location1", all.x = TRUE)
metadata <- left_join(metadata, sequencing_IDs_2, by = "location2", all.x = TRUE)
metadata <- rename(metadata, sequencingID_1 = SequencingID.x)
metadata <- rename(metadata, sequencingID_2 = SequencingID.y)

metadata$plate1 <- sapply(strsplit(metadata$location1, "_"), `[`, 1)
metadata$plate2 <- sapply(strsplit(metadata$location2, "_"), `[`, 1)
write.csv(metadata, "data/metadata/specimens_Sabin2Positive_CTbelow37_processedFinal_deID_runInfo.csv")

# ==================== Functions for plotting coverage with sliding window ==========================

slide <- function(cov.df, setup.df)
{
  coverage = rep(NA, nrow(setup.df))
  for(i in 1:nrow(setup.df))
  {
    s = setup.df$starts[i]
    e = setup.df$ends[i]
    subset(cov.df, concat.pos >= s & concat.pos < e, select = c(coverage)) -> position
    mean(position$coverage) -> coverage[i]
  }
  out <- data.frame(mean = coverage, concat.pos = setup.df$concat.pos, chr = setup.df$chr)
  out$Id = unique(cov.df$Id)
  return(out)
}

cov_plot <- function(cov.df, title)
{
  cov.df %>% group_by(chr) %>% summarize(first = min(concat.pos), last = max(concat.pos)) %>% plyr::adply(1,function(x) data.frame(starts = seq(x$first,x$last,by=100))) %>% mutate(ends = ifelse(starts + 200 < last, starts + 200, last)) -> setup
  setup %>% select(starts, ends) -> setup_means
  setup$concat.pos <- apply(setup_means, 1, function(x) mean(x))
  plyr::ddply(cov.df, ~Id, slide, setup) -> cov.slid.df
  
  x.labels <- plyr::ddply(cov.slid.df,~chr,plyr::summarize,concat.pos=concat.pos[which(abs(concat.pos-mean(concat.pos))==(min(abs(concat.pos-mean(concat.pos)))))])
  x.labels <- plyr::ddply(x.labels, ~chr, function(x) return(x[1,]))
  
  cov.plot <- ggplot(cov.slid.df, mapping = aes(x = as.factor(concat.pos), y = mean)) + geom_boxplot(fill="white")
  cov.plot <- cov.plot + ggtitle(title) + ylab("Read Depth") + scale_x_discrete(labels = x.labels$chr, breaks = x.labels$concat.pos) + xlab("Genome Position")
  cov.plot <- cov.plot + theme(axis.title.y = element_text(vjust=1.2))
  cov.plot <- cov.plot + theme(legend.position = "none") + theme_classic()
  return(cov.plot)
}

# ======================== Functions for filtering by coverage level ====================

# Pass in the dataframe, the cutoff value, and method ("site" or "bin"). Default windows are 50 bp with no overlap.

filter_by_coverage <- function(df, cutoff, method, step_size = 50, window_size = 50)
{
  stopifnot(cutoff > 0) # Error: cutoff must be greater than zero
  stopifnot(nrow(df) > 0) # Error: the dataframe is empty
  
  if(method == "site")
  {
    num_positions <- length(unique(df$chr.pos))
    df %>% filter(coverage > cutoff) %>% group_by(Id) %>% filter(n() == num_positions) -> quality
    df_return <- filter(df, Id %in% quality$Id)
    
  } else if(method == "bin")
  {
    df %>% 
      group_by(chr) %>%
      summarize(first = min(concat.pos), last = max(concat.pos)) %>% 
      plyr::adply(1,function(x) data.frame(starts = seq(x$first, x$last, by = step_size))) %>% 
      mutate(ends = ifelse(starts + window_size < last, starts + window_size, last)) -> setup
    setup %>% select(starts, ends) %>% mutate(window_length = ends - starts) %>% filter(window_length == window_size) %>% select(-window_length) -> setup_means
    setup %>% mutate(window_length = ends - starts) %>% filter(window_length == window_size) %>% select(-window_length) -> setup
    setup$concat.pos <- apply(setup_means, 1, function(x) mean(x))
    plyr::ddply(df, ~Id, slide, setup) %>% rename(bin = concat.pos) -> coverage_by_bins
    num_bins <- length(unique(coverage_by_bins$bin))
    coverage_by_bins %>% filter(mean > cutoff) %>% group_by(Id) %>% filter(n() == num_bins) -> quality
    df_return <- filter(df, Id %in% quality$Id)
    
  } else 
  {
    print("Error: improper method; choose from site or bin")
    return()
  }
  
  return(df_return)
}


# Get regions of each sample that have variant-quality coverage.
get_variant_qualified_regions <- function(df, cutoff, step_size = 50, window_size = 50)
{
  stopifnot(cutoff > 0) # Error: cutoff must be greater than zero
  stopifnot(nrow(df) > 0) # Error: the dataframe is empty
  
  slide_helper <- function(cov.df, setup.df)
  {
    coverage = rep(NA, nrow(setup.df))
    start = rep(NA, nrow(setup.df))
    end = rep(NA, nrow(setup.df))
    out <- data.frame(concat.pos = setup.df$concat.pos, chr = setup.df$chr, start = start, end =  end)
    for(i in 1:nrow(setup.df))
    {
      s = setup.df$starts[i]
      e = setup.df$ends[i]
      subset(cov.df, concat.pos >= s & concat.pos < e, select = c(coverage)) -> position
      mean(position$coverage) -> coverage[i]
      
      out[i,]$start <- s
      out[i,]$end <- e
    }
    
    out$Id = unique(cov.df$Id)
    out$mean = coverage
    return(out)
  }
  
  df %>% 
    group_by(chr) %>%
    summarize(first = min(concat.pos), last = max(concat.pos)) %>% 
    plyr::adply(1,function(x) data.frame(starts = seq(x$first, x$last, by = step_size))) %>% 
    mutate(ends = ifelse(starts + window_size < last, starts + window_size, last)) -> setup
  setup %>% select(starts, ends) %>% mutate(window_length = ends - starts) %>% filter(window_length == window_size) %>% select(-window_length) -> setup_means
  setup %>% mutate(window_length = ends - starts) %>% filter(window_length == window_size) %>% select(-window_length) -> setup
  setup$concat.pos <- apply(setup_means, 1, function(x) mean(x))
  plyr::ddply(df, ~Id, slide_helper, setup) %>% rename(bin = concat.pos) -> coverage_by_bins
  
  coverage_by_bins %>% filter(mean > cutoff) -> quality

  return(quality)
}


# ====================== Filter by coverage ========================

coverage_df_all_Sab2_NoTrim <- filter(coverage_df_all, chr == "Sabin2ref")
variant_quality_regions <- get_variant_qualified_regions(coverage_df_all_Sab2_NoTrim, 200)
write.csv(variant_quality_regions, file = "data/processed/coverage_variant_quality_regions.csv")

# Filter out positions we didn't expect to amplify
Nextera_end_loss <- 50
forward_primer_position <- 98 + Nextera_end_loss
reverse_primer_position <- 7400 - Nextera_end_loss
coverage_df_all <- filter(coverage_df_all, chr.pos > forward_primer_position & chr.pos < reverse_primer_position)

### Sabin2 Consensus Filters ###

# Filter out any Id that has a coverage less than consensus_cutoff at any relative position on Sabin2. True recombinants will be filtered out here.
consensus_cutoff <- 10

coverage_df_all_Sab2 <- filter(coverage_df_all, chr == "Sabin2ref")

coverage_df_success_Sab2_consensus <- filter_by_coverage(coverage_df_all_Sab2, consensus_cutoff, "site")
length(unique(coverage_df_success_Sab2_consensus$Id))

### Sabin2 Variant Filters ###

variant_cutoff <- 200

coverage_df_success_Sab2_variant <- filter_by_coverage(coverage_df_all_Sab2, variant_cutoff, "site")
length(unique(coverage_df_success_Sab2_variant$Id))

coverage_df_success_Sab2_variant_byBin <- filter_by_coverage(coverage_df_all_Sab2, variant_cutoff, "bin")
length(unique(coverage_df_success_Sab2_variant_byBin$Id))


# ============================== Find partial genomes ========================

# Get those with partial genomes by percent coverage
cutoff <- 30
coverage_df_all_Sab2 <- filter(coverage_df_all, chr == "Sabin2ref")
positions <- length(unique(coverage_df_all_Sab2$chr.pos))
coverage_df_all_Sab2 %>% filter(coverage > cutoff) %>% group_by(Id) %>% summarize(percent_covered = (n()/positions)*100) -> Sab2_percentCovered
Sab2_percentCovered_partial <- filter(Sab2_percentCovered, percent_covered < 100 & percent_covered > 10)
coverage_df_partial <- filter(coverage_df_all, Id %in% Sab2_percentCovered_partial$Id)

# Get those with all but segment 1
coverage_allButSeg1 <- filter(coverage_df_all_Sab2, chr.pos > 1900 & chr.pos < 7300 & !Id %in% coverage_df_success_Sab2_consensus$Id)
coverage_allButSeg1_consensus <- filter_by_coverage(coverage_allButSeg1, 30, "site")
length(unique(coverage_allButSeg1_consensus$Id)) # 76 sequencing samples that have everything but segment 1!

# ============================== Update metadata based on coverage groups ========================

metadata <- mutate(metadata, partialGenome = ifelse(sequencingID_1 %in% coverage_df_partial$Id | sequencingID_2 %in% coverage_df_partial$Id, TRUE, FALSE))
metadata <- mutate(metadata, allButSeg1 = ifelse(sequencingID_1 %in% coverage_allButSeg1_consensus$Id | sequencingID_2 %in% coverage_allButSeg1_consensus$Id, TRUE, FALSE))
metadata <- mutate(metadata, Sab2_consensus_qualified = ifelse(sequencingID_1 %in% coverage_df_success_Sab2_consensus$Id | sequencingID_2 %in% coverage_df_success_Sab2_consensus$Id, TRUE, FALSE))
metadata <- mutate(metadata, is_putative_recombinant = ifelse(sequencingID_1 %in% putative_recom_ids$Id | sequencingID_2 %in% putative_recom_ids$Id, TRUE, FALSE))

# Variant-calling qualified
metadata <- mutate(metadata, Sab2_variant_qualified = NA)
metadata <- mutate(metadata, Sab2_variant_qualified = ifelse(reps == 1 & sequencingID_1 %in% coverage_df_success_Sab2_variant_byBin$Id, TRUE, NA))
metadata <- mutate(metadata, Sab2_variant_qualified = ifelse(reps == 2 & (sequencingID_1 %in% coverage_df_success_Sab2_variant_byBin$Id & sequencingID_2 %in% coverage_df_success_Sab2_variant_byBin$Id), TRUE, Sab2_variant_qualified))
metadata <- mutate(metadata, Sab2_variant_qualified = ifelse(is.na(Sab2_variant_qualified), FALSE, Sab2_variant_qualified))
metadata <- mutate(metadata, any_data = ifelse(partialGenome == TRUE | allButSeg1 == TRUE | Sab2_consensus_qualified == TRUE | Sab2_variant_qualified == TRUE, TRUE, FALSE))
 
metadata %>% mutate(collectionDateFloor = floor_date(ymd(SPECDT), "week")) %>%
  mutate(collectionDateRelative = collectionDateFloor - ymd("2016-1-24")) %>%
  mutate(collectionDateRelativeExact = ymd(SPECDT) - ymd("2016-1-24")) %>%
  mutate(collectionWeekRelative = round(collectionDateRelative / 7, 0)) -> metadata

write.csv(metadata, "data/metadata/specimens_Sabin2Positive_CTbelow37_processedFinal_deID_runInfo_covInfo.csv")

# Link copy number per microliter and collection week into metadata.
metadata <- mutate(metadata, S2_copiesPerUl_inTNA = S2copyNumberPerGram/1000) # extracted from 0.2 g of stool and used 5 uL in RT-qPCR
metadata %>% mutate(collectionDateFloor = floor_date(ymd(SPECDT), "week")) %>%
  mutate(collectionDateRelative = collectionDateFloor - ymd("2016-1-24")) %>%
  mutate(collectionDateRelativeExact = ymd(SPECDT) - ymd("2016-1-24")) %>%
  mutate(collectionWeekRelative = round(collectionDateRelative / 7, 0)) -> metadata
write.csv(metadata, "data/metadata/specimens_Sabin2Positive_CTbelow37_processedFinal_deID_runInfo_covInfo_week.csv")

# ========================== Get all the samples for running the variant pipeline =================

cutoff <- 30
coverage_df_all_Sab2 <- filter(coverage_df_all, chr == "Sabin2ref")
positions <- length(unique(coverage_df_all_Sab2$chr.pos))
coverage_df_all_Sab2 %>% filter(coverage > cutoff) %>% group_by(Id) %>% summarize(percent_covered = (n()/positions)*100) -> Sab2_percentCovered
Sab2_percentCovered_partial <- filter(Sab2_percentCovered, percent_covered > 10)
coverage_df_partial <- filter(coverage_df_all, Id %in% Sab2_percentCovered_partial$Id)
any_data_IDs <- unique(coverage_df_partial$Id)

all_files <- Sys.glob("data/raw/plate*/Final_variants/*.csv")
all_files_df <- data.frame(filename = all_files)
all_files_df %>% 
  mutate(filename_original = filename) %>% 
  separate(filename, c("data", "raw", "platen", "Final_variants", "file"), sep = "/") %>% 
  separate(file, c("seq_id", "removed", "ref", "sum", "filtered", "AA", "csv"), sep = "\\.") %>%
  filter(seq_id %in% any_data_IDs) -> files_use_df

write.csv(files_use_df, "data/processed/variant_filenames_with_data.csv") # This goes into PrepareVariantData.R

# ========================== Get all filenames with consensus genomes for phylogenetic analysis ======================

metadata <- read.csv("data/metadata/specimens_Sabin2Positive_CTbelow37_processedFinal_deID_runInfo_covInfo_week.csv")
metadata %>% 
  filter(Sab2_consensus_qualified == TRUE, collectionWeekRelative > 0) %>%
  group_by(sid) %>%
  filter(collectionWeekRelative == max(collectionWeekRelative)) %>% 
  ungroup() -> metadata_consensus

metadata_treetime <- select(metadata_consensus, anonSpecId, SPECDT) %>% rename(name = anonSpecId, date = SPECDT) %>% mutate(name = as.character(name), date = as.character(date))
write.csv(metadata_treetime, "data/metadata/metadata_treetime_single.csv", row.names = FALSE, quote = FALSE)

# by sequencing ID
consensus_files <- data.frame(Id = unique(coverage_df_success_Sab2_consensus$anon))
consensus_files <- mutate(consensus_files, filename = paste0("data/processed/consensus_filtered/", Id, ".removed.parsed.filtered.fasta"))
consensus_files <- mutate(consensus_files, filename_new = paste0("data/processed/GenerateTree_Campaign/consensus/", Id, ".fasta"))

# by anonSpecId
consensus_files <- data.frame(anonSpecId = unique(metadata_consensus$anonSpecId))
consensus_files <- mutate(consensus_files, filename = paste0("data/processed/consensus_resolved/", anonSpecId, ".removed.parsed.filtered.resolved.fasta"))
consensus_files <- mutate(consensus_files, filename_new = paste0("data/processed/GenerateTree_Campaign/consensus/", anonSpecId, ".fasta"))

for(i in 1:nrow(consensus_files))
{
  row <- consensus_files[i,]
  filename = row$filename
  filename_new <- row$filename_new
  command <- paste0("cp ", filename, " ", filename_new)
  #system(command) # uncomment to run
}

# concatenate into single file: change_fasta_names.py
# then add Sabin 2 reference
# then muscle/RAxML

#################################################################################################################
# =================================== The rest is just for data visualization purposes ===========================
#################################################################################################################

# =================================== Composition of partial genomes ===========================

# For those that are only partial, what is the percent covered by variant-quality bins?

# Get those with partial genomes by percent coverage
cutoff <- 30
coverage_df_all_Sab2 <- filter(coverage_df_all, chr == "Sabin2ref")
positions <- length(unique(coverage_df_all_Sab2$chr.pos))
coverage_df_all_Sab2 %>% filter(coverage > cutoff) %>% group_by(Id) %>% summarize(percent_covered = (n()/positions)*100) -> Sab2_percentCovered
Sab2_percentCovered_partial <- filter(Sab2_percentCovered, percent_covered < 100 & percent_covered > 10)
coverage_df_partial <- filter(coverage_df_all_Sab2, Id %in% Sab2_percentCovered_partial$Id & Id %in% metadata$sequencingID_1)

percent_quality <- function(df, cutoff, method, step_size = 50, window_size = 50)
{
  stopifnot(cutoff > 0) # Error: cutoff must be greater than zero
  stopifnot(nrow(df) > 0) # Error: the dataframe is empty
  
  if(method == "site")
  {
    num_positions <- length(unique(df$chr.pos))
    df %>% filter(coverage > cutoff) %>% group_by(Id) %>% summarize(percent_covered = (n()/num_positions)*100) -> percentCovered
    df_return <- percentCovered
    
  } else if(method == "bin")
  {
    num_positions <- length(unique(df$chr.pos))
    df %>% 
      group_by(chr) %>%
      summarize(first = min(concat.pos), last = max(concat.pos)) %>% 
      plyr::adply(1,function(x) data.frame(starts = seq(x$first, x$last, by = step_size))) %>% 
      mutate(ends = ifelse(starts + window_size < last, starts + window_size, last)) -> setup
    setup %>% select(starts, ends) %>% mutate(window_length = ends - starts) %>% filter(window_length == window_size) %>% select(-window_length) -> setup_means
    setup %>% mutate(window_length = ends - starts) %>% filter(window_length == window_size) %>% select(-window_length) -> setup
    setup$concat.pos <- apply(setup_means, 1, function(x) mean(x))
    plyr::ddply(df, ~Id, slide, setup) %>% rename(bin = concat.pos) -> coverage_by_bins
    num_bins <- length(unique(coverage_by_bins$bin))
    coverage_by_bins %>% filter(mean > cutoff) %>% group_by(Id) %>% summarize(percent_covered = ((n() * window_size)/num_positions)*100) -> percentCovered
    df_return <- percentCovered
    
  } else 
  {
    print("Error: improper method; choose from site or bin")
    return()
  }
  
  return(df_return)
}

partial_percent_quality <- percent_quality(coverage_df_partial, 200, "bin")

palette <- wes_palette("Darjeeling1")
partial.composition.plot <- ggplot(partial_percent_quality, aes(percent_covered)) +
  geom_histogram(binwidth = 2.5, fill = palette[5], color = "white") +
  xlab("Percent of Genome Covered") +
  ylab("Number of Samples") +
  theme_classic() # save as 5 by 4

# ======================================= Make "representative" plot =============================

# Partial: 131746, 131889
# Consensus: 131924
# Variant: 131620

coverage_df_all_Sab2 <- filter(coverage_df_all, chr == "Sabin2ref")
examples <- c("131746", "131889", "131924", "131620")
example_cov <- filter(coverage_df_all_Sab2, Id %in% examples)

base = 15000
interval = 1500
linesize = 0.7

palette1 <- wes_palette("Darjeeling1")
palette2 <- wes_palette("Darjeeling2")
palette_viridis <- c("#FDE725FF", "#73D055FF", "#29AF7FFF", "#238A8DFF", "#39568CFF", "black")

example.plot <- ggplot(example_cov, aes(x = chr.pos, y = coverage, group = Id)) + 
  theme_light() +
  geom_line(aes(color = Id), alpha = 1) + 
  xlab("Genome Position") + 
  ylab("Read Depth") + 
  theme(legend.position = "none") +
  scale_y_log10(limits = c(10, 24000)) +
  scale_color_manual(name = "", values = c(palette1[1], palette1[5], palette1[5], palette2[2])) +
  geom_hline(yintercept = 200, linetype = "dotted", color = "black", size = 0.5) +
  geom_hline(yintercept = 10, linetype = "dotted", color = "black", size = 0.5) +
  geom_segment(size = linesize, color = "black", aes(x = 98, y = base + interval*0, xend = 2436, yend = base + interval*0)) +
  geom_segment(size = linesize, color = "black", aes(x = 1472, y = base + interval*2, xend = 4151, yend = base + interval*2)) +
  geom_segment(size = linesize, color = "black", aes(x = 3214, y = base + interval*4, xend = 5861, yend = base + interval*4)) +
  geom_segment(size = linesize, color = "black", aes(x = 4966, y = base + interval*6, xend = 7400, yend = base + interval*6)) # save as 6 by 4

# ========================== Test different variant cutoffs =================

metadata_final <- read.csv("data/metadata/specimens_Sabin2Positive_CTbelow37_processedFinal_deID_runInfo_covInfo_week_final.csv")
metadata_final <- filter(metadata_final, collectionWeekRelative > 0 & xor(mOPV2 == 1, transmissionAcquired == 1) & S2copyNumberPerGram > 900000)

coverage_df_success_Sab2_variant_200 <- filter_by_coverage(coverage_df_all_Sab2, 200, "bin") # 108 variant-quality samples. 111 including those in "hail mary" group
coverage_df_success_Sab2_variant_500 <- filter_by_coverage(coverage_df_all_Sab2, 500, "bin") # 81
coverage_df_success_Sab2_variant_1000 <- filter_by_coverage(coverage_df_all_Sab2, 1000, "bin") # 48

coverage_df_success_Sab2_variant_test <- coverage_df_success_Sab2_variant_1000
metadata_final <- mutate(metadata_final, Sab2_variant_qualified = NA)
metadata_final <- mutate(metadata_final, Sab2_variant_qualified = ifelse(reps == 1 & sequencingID_1 %in% coverage_df_success_Sab2_variant_test$Id, TRUE, NA))
metadata_final <- mutate(metadata_final, Sab2_variant_qualified = ifelse(reps == 2 & (sequencingID_1 %in% coverage_df_success_Sab2_variant_test$Id & sequencingID_2 %in% coverage_df_success_Sab2_variant_test$Id), TRUE, Sab2_variant_qualified))
metadata_final <- mutate(metadata_final, Sab2_variant_qualified = ifelse(is.na(Sab2_variant_qualified), FALSE, Sab2_variant_qualified))
nrow(filter(metadata_final, Sab2_variant_qualified == TRUE))
