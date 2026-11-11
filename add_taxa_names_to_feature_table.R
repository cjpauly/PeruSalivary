#Add taxon names to feature tables
#Written by Iyunoluwa Ademola-Popoola
#Edited by CJP
#10/29/2026
#R 4.6.1

#Workspace management
rm(list = ls(all.names = TRUE))

# Use this bash code to export feature table and taxonomy from QIIME
# qiime tools export \
# --input-path silva_taxonomy.qza \
# --output-path . &&
# qiime tools export \
# --input-path nocontrol-noEMC.qza \
# --output-path tempbiom &&
# mv tempbiom/feature-table.biom nocontrol-noEMC.biom &&
# rmdir tempbiom &&
# biom convert -i nocontrol-noEMC.biom -o nocontrol-noEMC.tsv --to-tsv


# Import your feature table
# rows = taxa; columns = samples
featuretable <- read.delim("nocontrol-noEMC.tsv",
                           skip = 1,
                           header = TRUE,
                           check.names = FALSE,
                           comment.char = "",
                           quote = ""
)

# Import the taxonomy table (Feature IDs and names)
taxonomy <- read.delim("taxonomy.tsv", header=TRUE)
head(taxonomy) # Taxonomy file should have Feature.ID and Taxon columns

# Merge the feature table with the taxonomy table by Feature ID
colnames(featuretable)[1] <- "Feature.ID"
head(featuretable)
table_with_names <- merge(featuretable, taxonomy, by = "Feature.ID", all.x = TRUE)
head(table_with_names) # Carrying forward Feature.ID and Taxon columns

# Save the merged table to a CSV file for further analysis
write.csv(table_with_names, "nocontrol-noEMC_with_names.csv", row.names = FALSE)
