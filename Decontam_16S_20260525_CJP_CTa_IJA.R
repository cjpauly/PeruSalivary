##### 16S Decontam CJP running CTa's script #####
## Initial 2026 May 25, Current 2026 July 29
## R 4.6.1

### Load Libraries ###
install.packages("installr"); library(installr); packageVersion("installr") ## 0.23.4
if (!require("BiocManager", quietly = TRUE)){install.packages("BiocManager")} ## Bioconductor 3.23, BiocManager 1.30.27
BiocManager::install("decontam"); library(decontam); packageVersion("decontam") ## 1.32.0
BiocManager::install("phyloseq"); library(phyloseq); packageVersion("phyloseq") ## 1.56.0
library(ggplot2); packageVersion("ggplot2") ##‘4.0.3’
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}; library(devtools); packageVersion("devtools") #2.5.2
pak::pak("jbisanz/qiime2R"); library(qiime2R); packageVersion("qiime2R") ## 0.99.6
install.packages("tidyverse"); library(tidyverse); packageVersion("tidyverse") ## 2.0.0
install.packages("ggthemes"); library(ggthemes); packageVersion("ggthemes") #7.0.0
install.packages("extrafont"); library(extrafont); packageVersion("extrafont") #0.20
install.packages("scales"); library(scales); packageVersion("scales") #1.4.0
library(reshape2); packageVersion("reshape2") #1.4.5
library(RColorBrewer); packageVersion("RColorBrewer") #1.1-3
install.packages("ggpubr"); library(ggpubr); packageVersion("ggpubr") #1.0.0
install.packages("melt"); library(melt); packageVersion("melt") #1.11.4
pak::pak("vmikk/metagMisc@621e9cf"); library(metagMisc); packageVersion("metagMisc") ## ‘0.5.0’
install.packages("data.table"); library(data.table); packageVersion("data.table") #1.18.6.1
library(plyr); packageVersion("plyr")  ## ‘1.8.9’
install.packages("metamisc"); library(metamisc); packageVersion("metamisc") ##‘0.4.0’
install.packages("writexl"); library(writexl); packageVersion("writexl") #2.0.1
library(biomformat); packageVersion("biomformat") #1.40.0
#renv::update() to check for new versions for the WHOLE PROJECT
#renv::snapshot() if it still works, renv::restore() if it doesn't

#Workspace management
rm(list = ls(all.names = TRUE))
## Define subdirectory in console: ".", "Peru", or "Global"
## Ex.: subwd <- "Peru"

###### Part 1: To line 155 #####
# load metadata
METADATA <- read.table(file.path(subwd, "sample_metadata.tsv"), header = TRUE, sep = "\t", quote = "", fill = TRUE)

# Ensure metadata has unique sample names and no empty names
METADATA <- METADATA %>%
  filter(!is.na(id) & id != "") %>%
  distinct(id, .keep_all = TRUE)

# Convert to sample_data object
sample_data <- sample_data(METADATA)
rownames(sample_data) <- sample_data$id

# Import species table
Species <- import_biom(file.path(subwd, "json.biom"))
#Species.matrix<- as.matrix(Species) ## not needed?

# Import taxonomy
S_taxonomy <- read_qza(file.path(subwd, "silva_taxonomy.qza"))
G_taxonomy <- read_qza(file.path(subwd, "greengenes_taxonomy.qza"))
taxatable_S <- S_taxonomy$data %>%
  as_tibble() %>%
  separate(Taxon, sep = ";", into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
taxatable_G <- G_taxonomy$data %>%
  as_tibble() %>%
  separate(Taxon, sep = ";", into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

# Ensure the taxonomy table has OTU IDs as rownames
taxatable_S <- taxatable_S %>%
  column_to_rownames("Feature.ID")
taxatable_G <- taxatable_G %>% column_to_rownames("Feature.ID")

# Get sample names from both metadata and species table
metadata_samples <- rownames(sample_data)
species_samples <- sample_names(Species)

# Find common samples
common_samples <- intersect(metadata_samples, species_samples)

# Subset metadata and species table to include only common samples between the metadata and biom table
sample_data <- sample_data[common_samples, ]
Species <- prune_samples(common_samples, Species)

# Convert taxatable to a matrix to avoid coercion warning
taxatable_S <- as.matrix(taxatable_S)
taxatable_G <- as.matrix(taxatable_G)

#Species.OTU=otu_table(Species.matrix, taxa_are_rows = TRUE) ##not needed?

phyloseq_species_S <- phyloseq(Species, sample_data, tax_table(taxatable_S))
phyloseq_species_S
phyloseq_species_G <- phyloseq(Species, sample_data, tax_table(taxatable_G))
phyloseq_species_G

# I proceeded to identify the prevalence of contaminants
sample_data(phyloseq_species_S)$is.neg <- sample_data(phyloseq_species_S)$SampleOrControl == "Control"
sample_data(phyloseq_species_G)$is.neg <- sample_data(phyloseq_species_G)$SampleOrControl == "Control"
## Alternative syntax: sample_data(phyloseq_species_S)$is.neg <- sample_data(phyloseq_species_S)$SampleType != "Sample"


# Identify contaminants with different thresholds
contamdf_prevS_50 <- isContaminant(phyloseq_species_S, method = "prevalence", neg = "is.neg", threshold = 0.5)
contamdf_prevS_55 <- isContaminant(phyloseq_species_S, method = "prevalence", neg = "is.neg", threshold = 0.55)
contamdf_prevS_60 <- isContaminant(phyloseq_species_S, method = "prevalence", neg = "is.neg", threshold = 0.60)
contamdf_prevS_65 <- isContaminant(phyloseq_species_S, method = "prevalence", neg = "is.neg", threshold = 0.65)
contamdf_prevS_75 <- isContaminant(phyloseq_species_S, method = "prevalence", neg = "is.neg", threshold = 0.75)
contamdf_prevS_90 <- isContaminant(phyloseq_species_S, method = "prevalence", neg = "is.neg", threshold = 0.9)
contamdf_prevG_50 <- isContaminant(phyloseq_species_G, method = "prevalence", neg = "is.neg", threshold = 0.5)
contamdf_prevG_55 <- isContaminant(phyloseq_species_G, method = "prevalence", neg = "is.neg", threshold = 0.55)
contamdf_prevG_60 <- isContaminant(phyloseq_species_G, method = "prevalence", neg = "is.neg", threshold = 0.60)
contamdf_prevG_65 <- isContaminant(phyloseq_species_G, method = "prevalence", neg = "is.neg", threshold = 0.65)
contamdf_prevG_75 <- isContaminant(phyloseq_species_G, method = "prevalence", neg = "is.neg", threshold = 0.75)
contamdf_prevG_90 <- isContaminant(phyloseq_species_G, method = "prevalence", neg = "is.neg", threshold = 0.9)

# Ensure taxonomy is in a dataframe format with OTU IDs as row names
taxatable_df_S <- as.data.frame(taxatable_S)
taxatable_df_S <- tibble::rownames_to_column(taxatable_df_S, "OTU ID")
taxatable_df_G <- as.data.frame(taxatable_G)
taxatable_df_G <- tibble::rownames_to_column(taxatable_df_G, "OTU ID")

# Merge contamination data with taxonomy data
contamdf_prevS_50$taxon<-rownames((contamdf_prevS_50))
contamdf_prevS_55$taxon<-rownames((contamdf_prevS_55))
contamdf_prevS_60$taxon<-rownames((contamdf_prevS_60))
contamdf_prevS_65$taxon<-rownames((contamdf_prevS_65))
contamdf_prevS_75$taxon<-rownames((contamdf_prevS_75))
contamdf_prevS_90$taxon<-rownames((contamdf_prevS_90))
contamdf_prevG_50$taxon<-rownames((contamdf_prevG_50))
contamdf_prevG_55$taxon<-rownames((contamdf_prevG_55))
contamdf_prevG_60$taxon<-rownames((contamdf_prevG_60))
contamdf_prevG_65$taxon<-rownames((contamdf_prevG_65))
contamdf_prevG_75$taxon<-rownames((contamdf_prevG_75))
contamdf_prevG_90$taxon<-rownames((contamdf_prevG_90))

contamdf_prev_taxa_S_50 <- merge(contamdf_prevS_50, taxatable_df_S, by.x = "row.names", by.y = "OTU ID")
contamdf_prev_taxa_S_55 <- merge(contamdf_prevS_55, taxatable_df_S, by.x = "row.names", by.y = "OTU ID")
contamdf_prev_taxa_S_60 <- merge(contamdf_prevS_60, taxatable_df_S, by.x = "row.names", by.y = "OTU ID")
contamdf_prev_taxa_S_65 <- merge(contamdf_prevS_65, taxatable_df_S, by.x = "row.names", by.y = "OTU ID")
contamdf_prev_taxa_S_75 <- merge(contamdf_prevS_75, taxatable_df_S, by.x = "row.names", by.y = "OTU ID")
contamdf_prev_taxa_S_90 <- merge(contamdf_prevS_90, taxatable_df_S, by.x = "row.names", by.y = "OTU ID")
contamdf_prev_taxa_G_50 <- merge(contamdf_prevG_50, taxatable_df_G, by.x = "row.names", by.y = "OTU ID")
contamdf_prev_taxa_G_55 <- merge(contamdf_prevG_55, taxatable_df_G, by.x = "row.names", by.y = "OTU ID")
contamdf_prev_taxa_G_60 <- merge(contamdf_prevG_60, taxatable_df_G, by.x = "row.names", by.y = "OTU ID")
contamdf_prev_taxa_G_65 <- merge(contamdf_prevG_65, taxatable_df_G, by.x = "row.names", by.y = "OTU ID")
contamdf_prev_taxa_G_75 <- merge(contamdf_prevG_75, taxatable_df_G, by.x = "row.names", by.y = "OTU ID")
contamdf_prev_taxa_G_90 <- merge(contamdf_prevG_90, taxatable_df_G, by.x = "row.names", by.y = "OTU ID")

table(contamdf_prevS_50$contaminant)
table(contamdf_prevS_55$contaminant)
table(contamdf_prevS_60$contaminant)
table(contamdf_prevS_65$contaminant)
table(contamdf_prevS_75$contaminant)
table(contamdf_prevS_90$contaminant)
table(contamdf_prevG_50$contaminant)
table(contamdf_prevG_55$contaminant)
table(contamdf_prevG_60$contaminant)
table(contamdf_prevG_65$contaminant)
table(contamdf_prevG_75$contaminant)
table(contamdf_prevG_90$contaminant)

### Compiling decontam table results into one excel workbook
contamlist_S<- write_xlsx(list(S.Threshold.50=contamdf_prev_taxa_S_50, S.Threshold.55=contamdf_prev_taxa_S_55, S.Threshold.60=contamdf_prev_taxa_S_60, S.Threshold.65=contamdf_prev_taxa_S_65, S.Threshold.75=contamdf_prev_taxa_S_75, S.Threshold.90=contamdf_prev_taxa_S_90), path = file.path(subwd, "decontam_results_S.xlsx"), col_names = TRUE, format_headers = FALSE, use_zip64 = FALSE)
contamlist_G<- write_xlsx(list(G.Threshold.50=contamdf_prev_taxa_G_50, G.Threshold.55=contamdf_prev_taxa_G_55, G.Threshold.60=contamdf_prev_taxa_G_60, G.Threshold.65=contamdf_prev_taxa_G_65, G.Threshold.75=contamdf_prev_taxa_G_75, G.Threshold.90=contamdf_prev_taxa_G_90), path = file.path(subwd, "decontam_results_G.xlsx"), col_names = TRUE, format_headers = FALSE, use_zip64 = FALSE)

###### Part 2: To end #####
###Define classifier and threshold in console
###Ex. Silva @ 60: classifier <- "S" AND threshold <- "60"
write.table(get(paste0("contamdf_prev_taxa_", classifier, "_", threshold)), file.path(subwd, paste0("contamdf_prev_taxa_", classifier, "_", threshold, ".txt")), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

###Visualizing contaminants
ps_pa <- transform_sample_counts(get(paste0("phyloseq_species_", classifier)), function(abund) 1 * (abund > 0))
ps_pa_neg <- prune_samples(sample_data(ps_pa)$SampleOrControl == "Control", ps_pa)
ps_pa_pos <- prune_samples(sample_data(ps_pa)$SampleOrControl == "Sample", ps_pa)
## Alternate syntax:
## ps_pa_neg <- prune_samples(sample_data(ps_pa)$SampleType != "Sample", ps_pa)
## ps_pa_pos <- prune_samples(sample_data(ps_pa)$SampleType == "Sample", ps_pa)

df_pa <- data.frame(pa.pos = taxa_sums(ps_pa_pos), pa.neg = taxa_sums(ps_pa_neg), contaminant = get(paste0("contamdf_prev", classifier, "_", threshold))$contaminant)
ggplot(data = df_pa, aes(x = pa.neg, y = pa.pos, color = contaminant)) +
  geom_point() +
  xlab("Prevalence (Negative Controls)") +
  ylab("Prevalence (Samples)")

# Setting the data into a dataframe and adding a library size column to the data
proj_df <- as.data.frame(sample_data(get(paste0("phyloseq_species_", classifier))))
proj_df$LibrarySize <- sample_sums(get(paste0("phyloseq_species_", classifier)))

# Creating a visualization of where the controls fall in comparison to the samples of their sequences/library size
proj_df <- proj_df[order(proj_df$LibrarySize), ]
proj_df$index <- seq(nrow(proj_df))
ggplot(data = proj_df, aes(x = index, y = LibrarySize, color = SampleOrControl)) +
  geom_point()
### OR SampleType instead of SampleOrControl, etc.

# Assign taxonomy to the contaminants
proj_taxonomy <- parse_taxonomy(get(paste0(classifier, "_taxonomy"))$data)
proj_taxonomy <- tibble::rownames_to_column(proj_taxonomy, "featureID")
proj_contamdf <- tibble::rownames_to_column(get(paste0("contamdf_prev", classifier, "_", threshold)), "featureID")

# Merging tables to know the taxonomy of the contaminants
proj_contam_taxonomy <- dplyr::semi_join(proj_taxonomy, proj_contamdf, by = "featureID")

# Exporting the taxonomy of contaminants
write.table(proj_contam_taxonomy, file = file.path(subwd, paste0("contamdf_prev_", classifier, "_", threshold, ".csv")), sep = ",", quote = FALSE, col.names = TRUE, row.names = TRUE)

#  Pruning based on set threshold
proj_pruned <- prune_taxa(!df_pa$contaminant, get(paste0("phyloseq_species_", classifier)))
proj_pruned

# Convert to dataframe and export
proj_pruned_df <- phyloseq_to_df(proj_pruned, addtax = FALSE, addtot = FALSE, addmaxrank = FALSE, sorting = "NULL")
colnames(proj_pruned_df) <- gsub(pattern = "\\bX", replacement = "", x = colnames(proj_pruned_df))

# Export the final table for analysis
write.table(proj_pruned_df, file.path(subwd, paste0("species_decontam_", classifier, "_", threshold, ".txt")), sep = "\t")
#Move row 1 one cell to the right (add tab at beginning of file), then delete column A (cut -f 2-) and quotation marks (sed -i '' 's/"//g')
#manually rename ID's as necessary (PCRNeg1.SA to PCRNeg1-SA with sed -i '' 's/\.SA/-SA/g'; Global: add leading X's back to set from 10317.00161094 to 10317.00192706)
