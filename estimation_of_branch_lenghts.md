# 
For one of my projects, I am trying understand whether annual and perennial species have differences in their rate of molecular evolution. Annual species have complete their life cycle in a single year, whereas perennials have 
Therefore, given a fixed amount of time, annual species should have gone through several rounds of completed life cycles, whereas perennials should not. Given that mutations that get into the next generation
occur during meiosis (the reproductive stage), annuals should have an acceleration in molecular evolution when compared to perennials. Given that mutation rate is difficult to estimate, one way to do so is by comparing the
amount of mutations accumulated since divergence from their most closely related perennial species.

One way to make this comparison is by using a phylogenetic approach where one would have a phylogeny of one (or multiple) gene(s). In theory, the branch length in a phylogenetic tree should be proportional to the amount of mutations accumulated since the Most Recent Common Ancestor (MRCA).
So, one could use the branch lenght at the terminal branch of each species and compare that between annual and perennial species. A problems comes when thinking that different sites in a gene evolve under different
evolutionary forces. In theory, only synonymous sites should be evolve free of selection. 

## Generate summary file
```sh
# Output file
output_file="branch_lengths_output.tsv"
echo -e "Species\tgene\ttype\tbranch_length" > "$output_file"

# Read species list
mapfile -t species_list < ./species_list

# Loop over JSON files
for json_file in ./results_FITTER/*.json; do
  gene_id=$(basename "$json_file" .aligned.trimAl.PAML.fasta.FITTER.json)

  for species in "${species_list[@]}"; do
    # Extract the block for the species (assumes each block is on its own lines)
    species_block=$(awk "/\"$species\"[[:space:]]*:\\{/ , /\\}/" "$json_file")

    # Extract synonymous and nonsynonymous values
    syn=$(echo "$species_block" | grep '"synonymous"' | sed -E 's/.*: *([0-9.eE+-]+).*/\1/')
    nonsyn=$(echo "$species_block" | grep '"nonsynonymous"' | sed -E 's/.*: *([0-9.eE+-]+).*/\1/')

    if [[ -n "$syn" && -n "$nonsyn" ]]; then
      echo -e "$species\t$gene_id\tsynonymous\t$syn" >> "$output_file"
      echo -e "$species\t$gene_id\tnon_synonymous\t$nonsyn" >> "$output_file"
    fi
  done
done
```
## Analyse data in R
### Prepare the data
```R
library(tidyverse)
library(AICcmodavg)
library(phylolm)

setwd("~/r_analysis")

# Read the output of the summary script above.

branch_length <- read.table("branch_lengths_output.tsv", header = T, sep = '\t')

# Read additional file. This file should have at least two columns: 'Species' and 'life_cycle'. The latter, should be 'annual' or 'perennial'.
life_cycle <- read.csv("species_annual_perennial.csv", header = T)

# Merge files by the colum species
branch_length_merged <- merge(branch_length, life_cycle, by = "Species")

branch_length_merged <- branch_length_merged %>%
  mutate(life_cycle_binary = ifelse(life_cycle == "annual", 1, 0)) # Perennials coded as '0' and annuals as '1'

rm(branch_length, life_cycle)

# HyPhy estimates branch length and both synonymous and non-synonymous sites. Both types of sites evolve a different rates and non-synonymous sites may be contrained by natural selection. Therefore, we need to split the file and keep only the synonymous sites. 

branch_lengths_syn <- filter(branch_length_merged, type != c("non_synonymous"))

# Filter outlier branch lengths per 'Species' by using the Interquantile Range (IQR) method. This is following Smith and Donoghue (2008) DOI: 10.1126/science.1163197.

branch_lengths_syn_filtered_per_species <- branch_lengths_syn %>%
  group_by(Species) %>%
  filter({
    Q1 <- quantile(branch_length, 0.25, na.rm = TRUE)
    Q3 <- quantile(branch_length, 0.75, na.rm = TRUE)
    IQR_val <- Q3 - Q1#Estimate Interquantile Range (IQR)
    lower_bound <- Q1 - (1.5 * IQR_val)
    upper_bound <- Q3 + (1.5 * IQR_val)
    branch_length >= lower_bound & branch_length <= upper_bound & branch_length >= 0.001 #I need to establish the lower bound manually to eliminate branches that are near zero
  }) %>%
  ungroup()

rm(branch_lengths_syn,branch_lengths_nonsyn)
```
### Statistical analyses
There are multiuple ways to analyse this data. For example, we might want to just colapse all the synonymous branch lengths of all genes per species, and therefore just compare the the raw differences between annual and perennial species. This might be done like this:
```R
# Read phylogenetic tree. The analysis will require a phylogenetic tree given that we will use the 'phylolm' package that accounts for the phylogenetic history of the species.
 tree_angios <- read.tree("/angiossperms_tree.tree")

 branch_lengths_syn_filtered_per_species$Species <- as.character(branch_lengths_syn_filtered_per_species$Species)

 branch_lengths_syn_filtered_per_species_summary <- branch_lengths_syn_filtered_per_species %>% group_by(life_cycle,Species,life_cycle_binary,SI_SC_binary,genus) %>%
   summarise_at(vars(branch_length),
                list(branch_length = sum))

## Prune tree. Eliminate potential taxa for which the branch length wa snot estimated, or that simply, must not be compared.

 species_in_tree <- intersect(branch_lengths_syn_filtered_per_species_summary$Species, tree_angios$tip.label)
 species_to_remove <- setdiff(tree_angios$tip.label, species_in_tree)
 pruned_tree <- drop.tip(tree_angios, species_to_remove)

 # Set the row names to be the Species column
 rownames(branch_lengths_syn_filtered_per_species_summary) <- branch_lengths_syn_filtered_per_species_summary$Species

 # Now reorder to match the tree tip labels exactly
 branch_lengths_syn_filtered_per_species_summary <- branch_lengths_syn_filtered_per_species_summary[pruned_tree$tip.label, ]

 # Set the row names to be the Species column
 rownames(branch_lengths_syn_filtered_per_species_summary) <- branch_lengths_syn_filtered_per_species_summary$Species

 # Verify the match
 all(rownames(branch_lengths_syn_filtered_per_species_summary) == pruned_tree$tip.label)

 model_full <- phylolm(branch_length ~ life_cycle_binary * SI_SC_binary,
                       data = branch_lengths_syn_filtered_per_species_summary,
                       phy = pruned_tree,
                       model="BM",
                       boot = 0)

 model_additive <- phylolm(branch_length ~ life_cycle_binary + SI_SC_binary,
                           data = branch_lengths_syn_filtered_per_species_summary,
                           phy = pruned_tree,
                           model="BM",
                           boot = 0)

 model_null <- phylolm(branch_length ~ 1,
                       data = branch_lengths_syn_filtered_per_species_summary,
                       phy = pruned_tree,
                       model="BM",
                       boot = 0)

 # Calculate AICc for model selection
 n <- model_additive$n

 # AICc for additive model
 k_add <- length(coef(model_additive))
 aicc_add <- AIC(model_additive) + (2 * k_add * (k_add + 1)) / (n - k_add - 1)

 # AICc for full model
 k_full <- length(coef(model_full))
 aicc_full <- AIC(model_full) + (2 * k_full * (k_full + 1)) / (n - k_full - 1)

 # Select model based on AICc
 # ΔAICc > 2 indicates support for full model

 aicc_full < aicc_add - 2

 ggplot(data = branch_lengths_nonsyn_filtered_per_species, aes(x = type, y = branch_length)) +
   geom_boxplot(aes(fill = life_cycle), width = 0.8) + theme_bw() + ylim (c(0, 0.1)) +
   labs(x = "", y = "Branch lengths")
```
Another way of doing the analysis is just by doing a Wilcoxon Signed-Rank test by comparing the cases in which branch length is indeed bigger in annuals than in perennials. The code below assumes that you have an additional column in the 'life_cycle' dataframe that includes the genus of both species.
```R
 # Separate annual and perennial
 annual <- branch_lengths_syn_filtered_per_species_summary %>% filter(life_cycle == "annual") %>% 
   select(genus, branch_length) %>% rename(annual_branch_length = branch_length)
 
 perennial <- branch_lengths_syn_filtered_per_species_summary %>% filter(life_cycle == "perennial") %>% 
   select(genus, branch_length) %>% rename(perennial_branch_length = branch_length)
 
 # Merge by genus
 paired <- merge(annual, perennial, by = "genus")
 
 # Compute difference (annual - perennial)
 paired$diff <- paired$annual_branch_length - paired$perennial_branch_length
 
 # Wilcoxon signed-rank test (one-tailed)
 
 wilcox_test <- wilcox.test(paired$annual_branch_length, paired$perennial_branch_length, 
                            paired = TRUE, alternative = "greater")
 print(wilcox_test)

 # Paired plot (lines connecting annual and perennial for each genus)
 paired_long <- paired %>%
   tidyr::pivot_longer(cols = c(annual_branch_length, perennial_branch_length), 
                       names_to = "life_stage", 
                       values_to = "measure") %>%
   mutate(life_stage = ifelse(life_stage == "annual_branch_length", "Annual", "Perennial"))
 
 ggplot(paired_long, aes(x = life_stage, y = measure, group = genus, color = genus)) +
   geom_point(size = 3) +
   geom_line() +
   labs(title = "Paired measures: Annual vs Perennial by Genus",
        x = "Life form", y = "Measure") +
   theme_minimal() +
   theme(legend.position = "right")
 
 # Boxplot of differences
 ggplot(paired, aes(y = diff)) +
   geom_boxplot(fill = "lightblue") +
   geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
   labs(title = "Distribution of differences (Annual - Perennial)",
        y = "Difference") +
   theme_minimal()
 
 # Boxplot of differences
 ggplot(paired, aes(y = diff)) +
   geom_boxplot(fill = "lightblue") +
   geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
   labs(title = "Distribution of differences (Annual - Perennial)",
        y = "Difference") +
   theme_minimal()
```
