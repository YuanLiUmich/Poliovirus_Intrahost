
### Project: Poliovirus sequencing from Matlab
### Purpose: Look at simple correlation between number of individuals with a mutation and proportion of cVDPV genomes with that mutation.
### Working directory: Poliovirus_Intrahost

# ================== Import packages, read in data =================

library(tidyverse)
library(wesanderson)

palette <- wes_palette("Darjeeling1")

stern_data <- read.csv("data/reference/VDPV_alignments/LocalGlobalMutationResults.csv")

vdpv.genome.plot.stern <- ggplot(data, aes(x = Count, y = Proportion, color = Type)) +
  geom_point() +
  scale_color_manual(name = "Mutation Type", values = c(palette[5], palette[3], palette[1], palette[2], "black")) + 
  xlab("Number of Individuals with Mutation in Matlab cohort") +
  ylab("Proportion of cVDPV Genomes with Mutation") +
  theme_bw() +
  facet_wrap(~Location)

idm_data <- read.csv("data/reference/VDPV_alignments/LocalGlobalMutationResults_IDMalignments.csv")

vdpv.genome.plot.idm <- ggplot(idm_data, aes(x = Count, y = Proportion, color = Type)) +
  geom_point() +
  scale_color_manual(name = "Mutation Type", values = c(palette[5], palette[3], palette[1], palette[2], "black")) + 
  xlab("Number of Individuals with Mutation in Matlab cohort") +
  ylab("Proportion of cVDPV Genomes with Mutation") +
  theme_bw() +
  facet_wrap(~Location) +
  ylim(c(0, 1))
