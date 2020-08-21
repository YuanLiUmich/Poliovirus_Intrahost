

### Project: Poliovirus sequencing from Matlab
### Purpose: Get the expected number of shared mutations based on a simple null distribution
### Working directory: Poliovirus_Intrahost

# ======================= Idea =====================

# An extremely basic permutation test: Get the number of individuals to be compared. 
# For each of those individuals, draw a number of variants based on the distribution of within-host variants 
# Assign random locations within the genome. 
# Tally up the number of positions seen in 2, 3, 4, or >4 individuals. 
# Run that 1000x and get the median number of expected shared variants across 2, 3, 4, or >4 individuals.

# ================== Import packages, read in data =================

library(tidyverse)
library(wesanderson)
set.seed(42)

palette <- wes_palette("Darjeeling1")

variants <- read.csv("data/processed/all.concat.varQuality.cleaned.variants.csv")

# ========================== How many mutations to draw per individual? ====================

variants_mOPV2_var <- filter(variants, ref != var & in_primer_site == FALSE & mOPV2 == 1 & transmissionAcquired == 0 & collectionWeekRelative > 0 & Sab2_variant_qualified == TRUE)

variants_mOPV2_var %>%
  group_by(sid) %>%
  summarize(numMutations = length(unique(mutation))) -> mutationsPerIndividual # distribution to draw from

numIndividuals <- length(unique(variants_mOPV2_var$sid)) # number of individuals to simulate

#ggplot(mutationsPerIndividual, aes(numMutations)) + geom_histogram(binwidth = 1) # not a normal distribution
#medianMutationsPerIndividual <- round(mean(mutationsPerIndividual$numMutations), 0)

# =========================== Support functions ===========================

single_permutation <- function(sites)
{
  sampled_sites <- c()
  for(index in 1:numIndividuals)
  {
    mutationsPerIndividual_drawn <- sample(mutationsPerIndividual$numMutations, size = 1)
    sites_one_individual <- sample(sites, size = mutationsPerIndividual_drawn, replace = FALSE)
    sampled_sites <- c(sampled_sites, sites_one_individual)
  }
  
  frequencies <- data.frame(table(sampled_sites)) %>% filter(Freq > 1)
  
  shared_2 <- nrow(filter(frequencies, Freq == 2))
  shared_3 <- nrow(filter(frequencies, Freq == 3))
  shared_4 <- nrow(filter(frequencies, Freq == 4))
  shared_high <- nrow(filter(frequencies, Freq > 4))
  
  mutations_shared <- data.frame(shared2 = shared_2, shared3 = shared_3, shared4 = shared_4, sharedHigh = shared_high)
  
  return(mutations_shared)
}

# ======================= Function to run one test =====================

run_permutation_test <- function(permutations = 1000, fraction_mutable)
{
  num_primer_sites <- (2436-2419) + (1491-1472) + (4151-4129) + (3230-3214) + (5861-5842) + (4982-4966)
  possible_sites <- seq(1, (max(variants$pos) - num_primer_sites - min(variants$pos) + 1) * fraction_mutable, by = 1)
  
  results <- data.frame()
  for(index in 1:permutations)
  {
    #print(paste0("On permutation ", as.character(index)))
    permutation_result <- single_permutation(possible_sites)
    results <- rbind(results, permutation_result)
  }
  
  return(results)
}

# ============================ Plot observed mutations vs. null distribution ==========================

results_1000 <- run_permutation_test(fraction_mutable = 0.6)

results_plot <- results_1000

num_shared_2 <- median(results_plot$shared2)
num_shared_3 <- median(results_plot$shared3)
num_shared_4 <- median(results_plot$shared4)
num_shared_high <- median(results_plot$sharedHigh)

var_count_multiple <- read.csv("data/processed/var_count_multiple.csv", stringsAsFactors = FALSE)
var_count_multiple$group <- factor(var_count_multiple$group, levels = c("2", "3", "4", ">4"))

shared.mutations.bar.withbox <- ggplot() + 
  geom_bar(data = var_count_multiple, aes(group, y = ..count.., fill = status)) +
  scale_fill_manual(name = "Mutation Type", values = c(palette[5], palette[3], palette[1], palette[2], "black")) + 
  ylab("Number of Mutations") + 
  xlab("Number of Individuals") + 
  theme_bw()

# Plot observed and simulated next to each other on the same axis.
var_count_multiple <- mutate(var_count_multiple, type = "Observed")
results_plot <- mutate(results_plot, type = "Simulated")
shared.mutations.null <- ggplot() + 
  geom_bar(data = var_count_multiple, aes(group, y = ..count.., fill = status), alpha =  1) +
  scale_fill_manual(name = "", values = c(palette[5], palette[3], palette[1], palette[2], "black")) + 
  ylab("Number of Mutations") + 
  xlab("Number of Individuals") + 
  theme_bw() +
  geom_boxplot(data = results_plot, aes(x = "2", y = shared2), alpha = 1) +
  geom_boxplot(data = results_plot, aes(x = "3", y = shared3), alpha = 1) +
  geom_boxplot(data = results_plot, aes(x = "4", y = shared4), alpha = 1) +
  geom_boxplot(data = results_plot, aes(x = ">4", y = sharedHigh), alpha = 1) +
  theme(legend.position = "left") +
  facet_wrap(~type)

# Observed numbers of shared mutations
var_count_multiple %>% group_by(group) %>% summarize(num = n()) -> observed_counts
observed_counts_num <- observed_counts$num
num_observed_2 <- observed_counts_num[1]
num_observed_3 <- observed_counts_num[2]
num_observed_4 <- observed_counts_num[3]
num_observed_high <- observed_counts_num[4]

plot2 <- ggplot(results_plot) + 
  geom_histogram(aes(shared2), color = "white", fill = "black") +
  geom_vline(xintercept = num_observed_2, lty = 2, color = "red") +
  ylab("Count") + 
  xlab("Mutations Shared")

plot3 <- ggplot(results_plot) + 
  geom_histogram(aes(shared3), color = "white", fill = "black") +
  geom_vline(xintercept = num_observed_3, lty = 2, color = "red") +
  ylab("Count") + 
  xlab("Mutations Shared")

plot4 <- ggplot(results_plot) + 
  geom_histogram(aes(shared4), color = "white", fill = "black") +
  geom_vline(xintercept = num_observed_4, lty = 2, color = "red") +
  ylab("Count") + 
  xlab("Mutations Shared")

plot4plus <- ggplot(results_plot) + 
  geom_histogram(aes(sharedHigh), color = "white", fill = "black") +
  geom_vline(xintercept = num_observed_high, lty = 2, color = "red") +
  ylab("Count") + 
  xlab("Mutations Shared")

null.dist.plot <- plot_grid(plot2, plot3, plot4, plot4plus, labels = c('A', 'B', 'C', 'D'), label_size = 20)

# Try a violin plot with points for observed data.
shared.mutations.null.box <- ggplot() + 
  ylab("Shared Mutations") + 
  xlab("Number of Individuals") + 
  theme_bw() +
  geom_point(data = results_plot, aes(x = " 2 ", y = num_observed_2), alpha = 1, shape = 19, color = "red", size = 2) +
  geom_point(data = results_plot, aes(x = " 3 ", y = num_observed_3), alpha = 1, shape = 19, color = "red", size = 2) +
  geom_point(data = results_plot, aes(x = " 4 ", y = num_observed_4), alpha = 1, shape = 19, color = "red", size = 2) +
  geom_point(data = results_plot, aes(x = ">4", y = num_observed_high), alpha = 1, shape = 19, color = "red", size = 2) +
  geom_boxplot(data = results_plot, aes(x = " 2 ", y = shared2), alpha = 0.3, outlier.shape = NA) +
  geom_boxplot(data = results_plot, aes(x = " 3 ", y = shared3), alpha = 0.3, outlier.shape = NA) +
  geom_boxplot(data = results_plot, aes(x = " 4 ", y = shared4), alpha = 0.3, outlier.shape = NA) +
  geom_boxplot(data = results_plot, aes(x = ">4", y = sharedHigh), alpha = 0.3, outlier.shape = NA) +
  theme(legend.position = "left")

# ========================= Get p-values at each fraction =====================

# p-value is the fraction of permutations with a number of shared mutations at or above what is observed for that group.

fractions_to_test <- seq(0.1, 1, by = 0.1)
pvalues_all <- data.frame()
for(fraction in fractions_to_test)
{
  print(paste0("On fraction ", as.character(fraction)))
  test_result <- run_permutation_test(fraction_mutable = fraction)
  
  pvalues <- c(nrow(filter(test_result, shared2 >= num_observed_2)) / nrow(test_result), 
               nrow(filter(test_result, shared3 >= num_observed_3)) / nrow(test_result), 
               nrow(filter(test_result, shared4 >= num_observed_4)) / nrow(test_result), 
               nrow(filter(test_result, sharedHigh >= num_observed_high)) / nrow(test_result) )
  groups <- c("2", "3", "4", ">4")
  pvalues_singlefraction <- data.frame(frac = fraction, group = groups, pvalue = pvalues)
  pvalues_all <- rbind(pvalues_all, pvalues_singlefraction)
}

pval_palette <- wes_palette("Zissou1")
pvalues_all$group <- factor(pvalues_all$group, levels = c("2", "3", "4", ">4"))
pval.plot <- ggplot(pvalues_all) +
  geom_point(aes(x = frac, y = pvalue, color = group)) +
  geom_line(aes(x = frac, y = pvalue, color = group)) +
  xlab("Fraction of Sites Available") +
  ylab("p-value") +
  theme_classic() + 
  geom_abline(slope = 0, intercept = 0.05, lty = 2) +
  scale_color_manual(name = "Group", values = c(pval_palette[1], pval_palette[2], pval_palette[3], pval_palette[5], pval_palette[5])) +
  scale_x_continuous(labels = c(seq(0, 1, by = 0.1)), breaks = c(seq(0, 1, by = 0.1))) +
  scale_y_continuous(labels = c(seq(0, 1, by = 0.1)), breaks = c(seq(0, 1, by = 0.1)))

