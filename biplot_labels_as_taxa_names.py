# This Jupyter Notebook script will allow you to create biplots using Emperor in QIIME2 with human-readable feature names
# It was created by Timur Yergaliyev (https://forum.qiime2.org/t/changing-colors-of-arrows-in-rpca-biplot-graphs-show-arrow-legend-on-graph/16152/12) and edited by Cameron Pauly with assistance from Copilot

#import
from pathlib import Path
import pandas as pd

#input
table = Path("nocontrol-noEMC.tsv")
taxa = Path("taxonomy.tsv")

#Verify feature table and taxa table formatted accurately
taxa_df = pd.read_csv(taxa, sep="\t", dtype={"Feature ID": str})
biom_df = pd.read_csv(table, sep="\t", skiprows=1, dtype={"#OTU ID": str})

required_biom_columns = {"#OTU ID"}
required_taxa_columns = {"Feature ID", "Taxon", "Confidence"}

if not required_biom_columns.issubset(biom_df.columns):
    raise ValueError(f"Feature table is missing columns: {required_biom_columns - set(biom_df.columns)}")
if not required_taxa_columns.issubset(taxa_df.columns):
    raise ValueError(f"Taxonomy table is missing columns: {required_taxa_columns - set(taxa_df.columns)}")

if biom_df["#OTU ID"].duplicated().any() or taxa_df["Feature ID"].duplicated().any():
    raise ValueError("Feature IDs must be unique in both input tables.")

#Align taxonomy rows to the feature table using IDs
old_ids = biom_df["#OTU ID"].astype(str)
taxa_df = (
    taxa_df.set_index("Feature ID")
    .reindex(old_ids)
    .rename_axis("Feature ID")
    .reset_index()
)

if taxa_df["Taxon"].isna().any():
    raise ValueError("Some feature table IDs have no matching taxonomy entry.")

#Clean taxonomy text
taxonomy_labels = (
    taxa_df["Taxon"]
    .fillna("")
    .astype(str)
    .str.rsplit("__", n=1)
    .str[-1]
)

taxonomy_labels = taxonomy_labels.mask(taxonomy_labels.eq(""), "unclassified")

#Replace Feature IDs & OTU IDs with unique short names
def make_unique(labels):
    labels = [str(label) for label in labels]
    reserved = set(labels)
    used = set()
    result = []

    for label in labels:
        unique_label = label
        if unique_label in used:
            suffix = 1
            unique_label = f"{label}.{suffix}"
            while unique_label in used or unique_label in reserved:
                suffix += 1
                unique_label = f"{label}.{suffix}"

        result.append(unique_label)
        used.add(unique_label)

    return result

new_ids = make_unique(taxonomy_labels)

if pd.Index(new_ids).has_duplicates:
    raise ValueError("Generated feature IDs are not unique.")

#Apply same IDs to both tables
biom_df["#OTU ID"] = new_ids
taxa_df["Feature ID"] = new_ids

### writing all output files
biom_df.to_csv("nocontrol-noEMC_with_names.tsv", sep="\t", index=False)
taxa_df[["Feature ID", "Taxon", "Confidence"]].to_csv(
    "taxonomy_with_names.tsv", sep="\t", index=False)

### Re-importing feature table and taxonomy with modified hashes
!biom convert -i nocontrol-noEMC_with_names.tsv -o nocontrol-noEMC_with_names.biom --table-type="OTU table" --to-hdf5
!qiime tools import \
    --input-path nocontrol-noEMC_with_names.biom \
    --type 'FeatureTable[Frequency]' \
    --input-format BIOMV210Format \
    --output-path nocontrol-noEMC_with_names.qza
!qiime tools import \
    --input-path taxonomy_with_names.tsv \
    --type 'FeatureData[Taxonomy]' \
    --input-format TSVTaxonomyFormat \
    --output-path taxonomy_with_names.qza

#Re-running qiime gemelli rpca and qiime emperor biplot with the new files
!qiime gemelli rpca \
  --i-table nocontrol-noEMC_with_names.qza \
  --p-min-feature-count 5 \
  --p-min-sample-count 5 \
  --o-biplot rpca/rpca_biplot_with_names.qza \
  --o-distance-matrix rpca/rpca_distance_with_names.qza  
!qiime emperor biplot \
  --i-biplot rpca/rpca_biplot_with_names.qza \
  --m-sample-metadata-file sample_metadata.tsv \
  --m-feature-metadata-file taxonomy_with_names.qza \
  --o-visualization rpca/rpca_biplot_with_names.qzv \
  --p-number-of-features 8
