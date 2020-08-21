
### Project: Poliovirus sequencing from Matlab
### Purpose: Get non-synonymous capsid mutations for structural mapping
### Working directory: Poliovirus_Intrahost

# ================== Import packages, read in data =================

library(tidyverse)

metadata <- read.csv("data/metadata/specimens_Sabin2Positive_CTbelow37_processedFinal_deID_runInfo_covInfo_week_final.csv")
variants <- read.csv("data/processed/all.concat.varQuality.cleaned.variants.ctl.csv")

# =============================== Get variants from the mOPV2 recipients from the campaign period ============================

variants_mOPV2 <- filter(variants, mOPV2 == 1 & transmissionAcquired == 0 & collectionWeekRelative > 0 & Sab2_variant_qualified == TRUE)
variants_mOPV2_capsid <- filter(variants_mOPV2, gene %in% c("VP1", "VP2", "VP3", "VP4"))
variants_mOPV2_capsid_ns <- filter(variants_mOPV2_capsid, class_factor == "Nonsynonymous")
variants_mOPV2_capsid_ns <- mutate(variants_mOPV2_capsid_ns, Var_AA_clean = substr(as.character(Var_AA), 3, nchar(as.character(Var_AA)) - 2))

# Number of individuals with a non-synonymous variant at each residue, regardless of what the residues are

variants_mOPV2_capsid_ns %>%
  group_by(AA_pos) %>%
  summarize(numIndivs = length(unique(sid))) %>%
  ungroup() -> capsid_ns_counts

# Prepare pymol commands

capsid_ns_counts <- mutate(capsid_ns_counts, gene = NA)
capsid_ns_counts <- mutate(capsid_ns_counts, gene = variants_mOPV2$gene[match(AA_pos, variants_mOPV2$AA_pos)])
capsid_ns_counts <- mutate(capsid_ns_counts, AA_pos_clean = as.numeric(substr(as.character(AA_pos), 2, nchar(as.character(AA_pos)) - 1)))
capsid_ns_counts <- mutate(capsid_ns_counts, residue = 0)
capsid_ns_counts <- mutate(capsid_ns_counts, residue = ifelse(gene == "VP4", AA_pos_clean, residue))
capsid_ns_counts <- mutate(capsid_ns_counts, residue = ifelse(gene == "VP2", AA_pos_clean - 69, residue))
capsid_ns_counts <- mutate(capsid_ns_counts, residue = ifelse(gene == "VP3", AA_pos_clean - 340, residue))
capsid_ns_counts <- mutate(capsid_ns_counts, residue = ifelse(gene == "VP1", AA_pos_clean - 578, residue))

capsid_ns_counts <- mutate(capsid_ns_counts, chain = "0")
capsid_ns_counts <- mutate(capsid_ns_counts, chain = ifelse(gene == "VP4", "4", chain))
capsid_ns_counts <- mutate(capsid_ns_counts, chain = ifelse(gene == "VP2", "2", chain))
capsid_ns_counts <- mutate(capsid_ns_counts, chain = ifelse(gene == "VP3", "3", chain))
capsid_ns_counts <- mutate(capsid_ns_counts, chain = ifelse(gene == "VP1", "1", chain))

capsid_ns_counts <- mutate(capsid_ns_counts, residue_str = as.character(residue))

capsid_ns_counts <- mutate(capsid_ns_counts, color = "0")
capsid_ns_counts <- mutate(capsid_ns_counts, color = ifelse(numIndivs == 3, "alv1", color))
capsid_ns_counts <- mutate(capsid_ns_counts, color = ifelse(numIndivs == 4, "alv2", color))
capsid_ns_counts <- mutate(capsid_ns_counts, color = ifelse(numIndivs == 5, "alv3", color))
capsid_ns_counts <- mutate(capsid_ns_counts, color = ifelse(numIndivs == 8, "alv4", color))
capsid_ns_counts <- mutate(capsid_ns_counts, color = ifelse(numIndivs == 9, "alv5", color))
capsid_ns_counts <- mutate(capsid_ns_counts, color = ifelse(numIndivs == 10, "alv6", color))
capsid_ns_counts <- mutate(capsid_ns_counts, color = ifelse(numIndivs == 11, "alv7", color))
capsid_ns_counts <- mutate(capsid_ns_counts, color = ifelse(numIndivs == 12, "alv8", color))
capsid_ns_counts <- mutate(capsid_ns_counts, color = ifelse(numIndivs > 12, "alv9", color))

capsid_ns_counts <- mutate(capsid_ns_counts, pymol_command = paste0("color ", color, ", chain ", chain, " and resi ", residue_str))
capsid_ns_counts <- mutate(capsid_ns_counts, pymol_command_2 = paste0("set transparency, 0, chain ", chain, " and resi ", residue_str))

capsid_ns_counts <- filter(capsid_ns_counts, numIndivs >= 3)

write.csv(capsid_ns_counts, "data/processed/PymolCapsidColor.csv") # then copy commands column into text file and run.

# Colors:
viridis::viridis(9, option = "plasma") -> palette_plasma

t(data.frame(col2rgb(palette_plasma))) -> palette_plasma_df

# Feature info
site_info <- read.csv("data/reference/OPV2_AntigenicSites_Absolute.csv")

site_info <- mutate(site_info, chain = "0")
site_info <- mutate(site_info, chain = ifelse(gene_antigenic == "VP4", "4", chain))
site_info <- mutate(site_info, chain = ifelse(gene_antigenic == "VP2", "2", chain))
site_info <- mutate(site_info, chain = ifelse(gene_antigenic == "VP3", "3", chain))
site_info <- mutate(site_info, chain = ifelse(gene_antigenic == "VP1", "1", chain))

site_info <- mutate(site_info, color = "0")
site_info <- mutate(site_info, color = ifelse(name == "receptorBinding", "col_receptor", color))
site_info <- mutate(site_info, color = ifelse(name == "NAg1", "col_nag1", color))
site_info <- mutate(site_info, color = ifelse(name == "NAg2", "col_nag2", color))
site_info <- mutate(site_info, color = ifelse(name == "NAg3a", "col_nag3a", color))
site_info <- mutate(site_info, color = ifelse(name == "NAg3b", "col_nag3b", color))

site_info <- mutate(site_info, pymol_command = paste0("color ", color, ", chain ", chain, " and resi ", aa_position))

write.csv(site_info, "data/reference/PDB/Pymol_Commands_CapsidFeatures.csv")
