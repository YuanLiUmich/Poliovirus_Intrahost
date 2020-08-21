

### Project: Polio Mixing Experiment 2.0

# ======================== Import packages, read in metadata =============================

library(tidyverse)
source("./supportFunctions.R")

sampleInfo <- read.csv("./data/reference/MixingStudySamples2.csv")
primer_fwd_outer <- 95 # 95 to 115
primer_rev_outer <- 7440 # 7415 to 7440
expectedVariants <- read_csv("./data/reference/MixingStudyExpectedVariants.csv")
expectedTruePos <- nrow(expectedVariants)
possible_vars <- (primer_rev_outer - primer_fwd_outer - 1)*3 # Positions in primer range, times 3.

# ============================= Function to find primer binding sites =======================

IsPositionInPrimerSiteOPV1 <- function(pos)
{
  # Primer binding sites
  seg1_fwd_start <- 98
  seg1_fwd_stop  <- 114
  seg1_rev_start <- 2416
  seg1_rev_stop  <- 2434
  seg2_fwd_start <- 1470
  seg2_fwd_stop  <- 1489
  seg2_rev_start <- 4130
  seg2_rev_stop  <- 4152
  seg3_fwd_start <- 3212
  seg3_fwd_stop  <- 3228
  seg3_rev_start <- 5843
  seg3_rev_stop  <- 5862
  seg4_fwd_start <- 4967
  seg4_fwd_stop  <- 4983
  seg4_rev_start <- 7384
  seg4_rev_stop  <- 7401
  
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

IsPositionInPrimerSiteOPV1.v <- Vectorize(IsPositionInPrimerSiteOPV1)


# ======================= Options for reading in variant data ====================

dups <- "both" # options: both, first, second
subset <- "200" # options: 100, 50, 25, 10
disp <- "onesided" # options: onesided or binomial
inputLevel <- "9_E2" # options: 4_E4, 9_E3, 9_E2, 9_E1
tna <- "yes" # options: yes, no

p.val.cutoff <- 0.1
MapQ.cutoff <- 20
Phred.cutoff <- 35
freq.var.cutoff <- 0
Read_pos_cutoff_1 <- 31
Read_pos_cutoff_2 <- 219

# ============================== Read in variant data ===============================

replicate = ifelse(dups == "first" | dups == "second", "notcollapsed", "collapsed")
file <- paste0("./data/processed/shiny2.verysens.", disp, ".", subset, ".", replicate, ".variants.csv")
variants_raw <- read.csv(file)

if(dups == "first") {
  data <- filter(variants_raw, Rep == 1)
} else if(dups == "second") {
  data <- filter(variants_raw, Rep == 2)
} else {
  data <- variants_raw
}

data <- filter(data, InputLevel == inputLevel, inTNA == tna)
data <- mutate(data, SampleNumber = Id)
data <- mutate(data, exp.freq = PercWT/100)
data <- mutate(data, Id = as.factor(as.character(PercWT)))
data <- filter(data, pos > primer_fwd_outer & pos < primer_rev_outer)
variant_data <- filter(data, p.val <= p.val.cutoff & 
                              MapQ >= MapQ.cutoff & 
                              Phred >= Phred.cutoff & 
                              freq.var >= freq.var.cutoff & 
                              Read_pos <= Read_pos_cutoff_2 & 
                              Read_pos >= Read_pos_cutoff_1)
variant_data <- mutate(variant_data, level = ifelse(PercWT == 1, "1_percent", 
                                                    ifelse(PercWT == 2, "2_percent", 
                                                           ifelse(PercWT == 5, "5_percent", 
                                                                  ifelse(PercWT == 10, "10_percent", 
                                                                         ifelse(PercWT == 100, "100_percent", NA))))))
variant_data <- mutate(variant_data, mutation_level = paste0(mutation, "_", level))

# ===================================== Read in coverage data ================================

file <- paste0("./data/processed/shiny2.verysens.", disp, ".", subset, ".coverage.csv")
coverage_raw <- read.csv(file)

# ============================ Plot coverage ==============================

coverage <- filter(coverage_raw, Id %in% variant_data$SampleNumber) # get the samples currently being analyzed in the variant data

base = 1000
interval = 300
linesize = 1.2
cov.plot <- ggplot(coverage, aes(x = chr.pos, y = coverage, group = Id)) + 
  geom_line(aes(color = Id)) + 
  theme_bw() + 
  xlab("Genome Position") + 
  ylab("Coverage (Read Depth)") + 
  theme(legend.position = "none") +
  xlim(c(95, 7440)) +
  geom_segment(data = coverage, size = linesize, color = "orange1", aes(x = 98, y = base + interval*0, xend = 2436, yend = base + interval*0)) +
  geom_segment(data = coverage, size = linesize, color = "mediumseagreen", aes(x = 1472, y = base + interval*1, xend = 4151, yend = base + interval*1)) +
  geom_segment(data = coverage, size = linesize, color = "indianred4", aes(x = 3214, y = base + interval*2, xend = 5861, yend = base + interval*2)) +
  geom_segment(data = coverage, size = linesize, color = "slateblue3", aes(x = 4966, y = base + interval*3, xend = 7400, yend = base + interval*3)) +
  geom_abline(intercept = 200, slope = 0, linetype = 3, size = 0.5, color = "black") + 
  geom_abline(intercept = 1000, slope = 0, linetype = 3, size = 0.5, color = "black")


# ================================== Make table ===========================

exSeg2 <- "yes"

if(exSeg2 == "yes")
{
  variant_data <- filter(variant_data, pos < 1500 | pos > 4100)
  expectedVariants <- filter(expectedVariants, Position < 1500 | Position > 4100)
  expectedTruePos <- nrow(expectedVariants)
  possible_vars <- (primer_rev_outer - primer_fwd_outer - 1 - (4100-1500))*3 # Positions in primer range, times 3.
} else {
  expectedTruePos <- nrow(expectedVariants)
  possible_vars <- (primer_rev_outer - primer_fwd_outer - 1)*3 # Positions in primer range, times 3.
}

dd = 4
m.roc.table <- miseq.roc.table(variant_data, 1, expectedTruePos, possible_vars, ">")
m.roc.table <- rename(m.roc.table, c("exp.freq"="Frequency","adj.sensitivity"="Sensitivity","TP"="True\nPositives","adj.specificity"="Specificity","FP"="False\nPositives"))
m.roc.table$Frequency <- c("100%", "10%", "5%", "2%", "1%")
m.roc.table$Sensitivity <- round(m.roc.table$Sensitivity, digits = dd)
m.roc.table$Specificity <- round(m.roc.table$Specificity, digits = dd)
m.roc.table

# ============================ Frequency by Position =======================


variant_data_TP <- filter(variant_data, category == TRUE & !is.na(level))

levels <- c("100_percent", "10_percent", "5_percent", "2_percent", "1_percent")
chrs <- data.frame("level" = levels)
chrs <- mutate(chrs, start = 0, stop = 7440)
palette <- wes_palette("Darjeeling1")
variant_data_TP$level <- factor(variant_data_TP$level, levels = rev(c("100_percent","10_percent","5_percent","2_percent","1_percent")))
chrs$level <- factor(chrs$level, levels = levels(variant_data_TP$level))

base = 0.6
interval = 0.02
linesize = 1.2

ggplot(variant_data_TP)  +
  geom_point(aes(x = pos, y = freq.var, color = level)) +
  ylab("Frequency") +
  xlab("Genome Position") +
  theme_minimal() +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  theme(panel.grid.major = element_line(colour = "gray96"), panel.grid.minor = element_line(colour = "white")) + 
  theme(legend.position = "right") +
  ylim(c(0, 0.7)) +
  scale_x_continuous(breaks = c(0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 7500)) +
  scale_color_manual(values = palette, name = "Expected Frequency") +
  geom_abline(intercept = 0.1, slope = 0, linetype = 3, size = 0.5, color = palette[4]) + 
  geom_abline(intercept = 0.05, slope = 0, linetype = 3, size = 0.5, color = palette[3]) + 
  geom_abline(intercept = 0.02, slope = 0, linetype = 3, size = 0.5, color = palette[2]) + 
  geom_abline(intercept = 0.01, slope = 0, linetype = 3, size = 0.5, color = palette[1]) +
  geom_segment(data = data, size = linesize, color = "orange1", aes(x = 98, y = base + interval*0, xend = 2436, yend = base + interval*0)) +
  geom_segment(data = data, size = linesize, color = "mediumseagreen", aes(x = 1472, y = base + interval*1, xend = 4151, yend = base + interval*1)) +
  geom_segment(data = data, size = linesize, color = "indianred4", aes(x = 3214, y = base + interval*2, xend = 5861, yend = base + interval*2)) +
  geom_segment(data = data, size = linesize, color = "slateblue3", aes(x = 4966, y = base + interval*3, xend = 7400, yend = base + interval*3))

# =================================== Plot true positives and false negatives by position =====================

expectedVariants <- read_csv("./data/reference/MixingStudyExpectedVariants_ByLevel.csv")
expectedVariants <- mutate(expectedVariants, mutation_level = paste0(mutation, "_", level))
expectedVariants <- mutate(expectedVariants, found = ifelse(mutation_level %in% variant_data$mutation_level, "True Positive", "False Negative"))

levels <- c("100_percent", "10_percent", "5_percent", "2_percent", "1_percent")
chrs <- data.frame("level" = levels)
chrs <- mutate(chrs, start = 0, stop = 7440)
palette <- wes_palette("Darjeeling1")
expectedVariants$level <- factor(expectedVariants$level, levels = rev(c("100_percent","10_percent","5_percent","2_percent","1_percent")))
chrs$level <- factor(chrs$level, levels = levels(expectedVariants$level))

ggplot(expectedVariants, aes(x = Position, y = level)) +
  geom_point(aes(color = found), size = 5, shape = 108) +
  geom_segment(data = chrs, aes(x = start, y = level, xend = stop, yend = level)) +
  ylab("Expected Frequency Group") +
  xlab("") +
  theme_minimal() +
  scale_color_manual(name = "", values = palette[c(1,2)]) +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  scale_x_continuous(breaks = c()) +
  theme(panel.grid.major = element_line(colour = "gray96"), panel.grid.minor = element_line(colour = "white")) + theme(legend.position = "bottom") +
  scale_x_continuous(breaks = c(0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 7500))
