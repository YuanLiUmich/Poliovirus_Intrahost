

### Project: Poliovirus sequencing from Matlab
### Purpose: Investigate linkage in donor 702
### Working directory: Poliovirus_Intrahost

# ================================= Read in data, import packages ============================

library(tidyverse)
library(cowplot)

data <- read.csv("data/processed/LinkageInfo_Donor702.csv", stringsAsFactors = FALSE)
transmission_variants_donorPolymorphic_702_minor <- read.csv("data/processed/variants_donor_polymorphic_702_minor.csv", stringsAsFactors = FALSE)

data <- mutate(data, fractionLinked_1 = FreqLinked / Freq1, fractionLinked_2 = FreqLinked / Freq2)
data <- mutate(data, fractionLinked_1 = FreqLinked / Freq1, fractionLinked_2 = FreqLinked / Freq2)
data <- mutate(data, freq1 = transmission_variants_donorPolymorphic_702_minor$freq1[match(Position1, transmission_variants_donorPolymorphic_702_minor$pos)],
               freq2 = transmission_variants_donorPolymorphic_702_minor$freq1[match(Position2, transmission_variants_donorPolymorphic_702_minor$pos)])

data_singlecomparison <- filter(data, Position1 != Position2 & Position1 < Position2 & TotalReads > 0)

# Compare frequency measured here to deepSNV to make sure we aren't waaay off
freq1.plot <- ggplot(data_singlecomparison, aes(x = freq1, y = Freq1)) + 
  geom_point(aes(color = cut(TotalReads, c(-Inf, 1000, Inf)))) + 
  geom_smooth(method = "lm", se = TRUE) +
  xlab("Frequency from deepSNV") +
  ylab("Frequency from Haplotype Script") +
  xlim(c(0, 0.25)) +
  ylim(c(0, 0.25)) +
  scale_color_manual(name = "Number of Total Reads", values = c("red", "black")) +
  theme(legend.position = "none")

freq2.plot <- ggplot(data_singlecomparison, aes(x = freq2, y = Freq2)) + # important to note that the measured frequency here is not all the reads, just those that overlap with other site
  geom_point(aes(color = cut(TotalReads, c(-Inf, 1000, Inf)))) + 
  geom_smooth(method = "lm", se = TRUE) +
  xlab("Frequency from deepSNV") +
  ylab("Frequency from Haplotype Script") +
  xlim(c(0, 0.25)) +
  ylim(c(0, 0.25)) +
  scale_color_manual(name = "Number of Total Reads", values = c("red", "black")) +
  theme(legend.position = "none")

plot_grid(freq1.plot, freq2.plot)

# Ok, now let's look at amount of linkage.
data_singlecomparison <- mutate(data_singlecomparison, pair = paste0(as.character(Position1), " and ", as.character(Position2)))

frac.link.1.plot <- ggplot(data_singlecomparison, aes(x = pair, y = fractionLinked_1)) + 
  geom_bar(stat = "identity") +
  ylim(c(0, 1)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle= 90, hjust = 1, vjust = 0.5)) +
  ylab("Fraction of Variant 1 Linked with Variant 2")

frac.link.2.plot <- ggplot(data_singlecomparison, aes(x = pair, y = fractionLinked_2)) + 
  geom_bar(stat = "identity") +
  ylim(c(0, 1)) +
  theme(axis.text.x=element_text(angle= 90, hjust = 1, vjust = 0.5)) +
  ylab("Fraction of Variant 2 Linked with Variant 1")

ggplot(data_singlecomparison, aes(x = freq1, fractionLinked_1)) + geom_point() + xlim(c(0, 1)) + ylim(c(0, 1))

# Let's make a matrix. 20 variants by 20 variants. In the cells, put 1) total reads, and 2) fraction linkage
data_lowerlefttri <- filter(data, Position1 <= Position2)

data_totalreads <- select(data_lowerlefttri, Position1, Position2, TotalReads)
matrix_totalreads <- spread(data_totalreads, Position1, TotalReads)
write.csv(matrix_totalreads, "data/processed/Linkage_TotalReads.csv")

data_fractionlinked <- select(data_lowerlefttri, Position1, Position2, fractionLinked_1)
matrix_fractionlinked <- spread(data_fractionlinked, Position1, fractionLinked_1)
write.csv(matrix_fractionlinked, "data/processed/Linkage_FractionLinked.csv")

