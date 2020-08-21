

### Project: Poliovirus sequencing from Matlab
### Purpose: Narrow the poliovirus alignments from IDM by those that are from AFP cases and not from iVDPV, etc.
### Working directory: Poliovirus_Intrahost

# ================== Import packages, read in data =================

library(tidyverse)

categories <- c("cVDPV")

attributes_5NCR <- read.table("data/reference/IDM_Sabin2_Alignments/Sabin2.nt5NCR.alignment.final.attributes.txt", stringsAsFactors = FALSE, sep = '\t', header = TRUE)
attributes_5NCR_AFP <- filter(attributes_5NCR, samplingFrame == "AFP" & classification %in% categories)

attributes_P1 <- read.table("data/reference/IDM_Sabin2_Alignments/Sabin2.ntP1.alignment.final.attributes.txt", stringsAsFactors = FALSE, sep = '\t', header = TRUE)
attributes_P1_AFP <- filter(attributes_P1, samplingFrame == "AFP" & classification %in% categories)

attributes_P2 <- read.table("data/reference/IDM_Sabin2_Alignments/Sabin2.ntP2.alignment.final.attributes.txt", stringsAsFactors = FALSE, sep = '\t', header = TRUE)
attributes_P2_AFP <- filter(attributes_P2, samplingFrame == "AFP" & classification %in% categories)

attributes_P3 <- read.table("data/reference/IDM_Sabin2_Alignments/Sabin2.ntP3.alignment.final.attributes.txt", stringsAsFactors = FALSE, sep = '\t', header = TRUE)
attributes_P3_AFP <- filter(attributes_P3, samplingFrame == "AFP" & classification %in% categories)

accessions <- unique(c(attributes_5NCR_AFP$accession, attributes_P1_AFP$accession, attributes_P2_AFP$accession, attributes_P3_AFP$accession))

write.table(accessions, "data/reference/IDM_Sabin2_Alignments/accessionsToUse.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)
