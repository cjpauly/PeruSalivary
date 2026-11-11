##### NOTE TO SELF: THIS STARTS FROM [75] WHEN I SWITCHED TO DADA2/
##### CHECK BEFORE THIS IN NOTEBOOK AS WELL [also there is a Peru run between importing Ethiopia and Malaysia]
##### PUT ONLY INTRO FILES IN GITHUB (SEE "These data were demultiplexed ..." FOR LIST)
######### Bolivia: Global/demux_bolivia/Bolivia_FASTQ_script.sh and Global/demux_bolivia/Bolivia_FASTQ.tsv
######### Industrial: Global/demux_industrial/FileNames.txt and .qza's from /storage/group/LiberalArts/default/lsw132_collab/09_LabMembersFolders/Collab/AGP_Project/2_ProcessedSequences/PostDemux/
######### Malaysia: Global/demux_malaysia/ENA_report.tsv , Global/demux_malaysia/OrangAsli.csv , and Global/demux_malaysia/samplenames.txt

#!/usr/bin/env bash

# QIIME2 v2026.4 16S Pipeline
# See environment.yml to load Conda environment
# Using DADA2 for denoising

#### Import sequencing data for global analysis ####
cd Global/

# Relevant metadata files:
# sample_metadata.tsv
# dada2/FileNames.txt (includes only the eight projects which were included in the final dataset - no Ecuador or Malaysia)
# XXXX

# These data were demultiplexed and imported to QIIME2 by Iyunoluwa J. Ademola-Popoola (https://github.com/microARCHlab/UgandaGlobalMicrobiome_2025/blob/main/Data_wrangling)
# Nepal data (Ryu et al. 2024): demux_2026/demux_nepal.qza
# Tanzania data (Bisanz et al. 2015): demux_2026/demux_tanzania.qza
# Uganda data (Ademola-Popoola et al., in press): demux_2026/demux_uganda.qza
# Venezuela data (Clemente et al. 2015): demux_2026/demux_venezuela.qza

# These data were wrangled, as below
# Bolivia data (Han et al. 2016): demux_bolivia/
# Ecuador data (Vásquez-Toasa et al. 2025): demux_ecuador/Ecuador_FASTQ.tsv
# Ethiopia data (Mulisa et al. 2026): EOSC_weyrich.tar.gz
# Industrial data (McDonald et al. 2018): demux_industrial/
# Malaysia data (Yeo et al. 2022): demux_malaysia/
# Peru data (Barbieri et al. 2019): Fastq_Raw/

## 1. Industrial ##
cd demux_industrial/

# There is a FileNames.txt in here as well, along with .qza files
# Move .qza files to subdirectory
mkdir filtered-demuxed-qzas/
mv *-filtered.qza filtered-demuxed-qzas/

# Unzip AGP .qza files
mkdir fastq/
while IFS= read -r sample || [ -n "$sample" ];
do
  mkdir -p "fastq/$sample";
  unzip -q "filtered-demuxed-qzas/${sample}-filtered.qza" '*fastq*' -d "fastq/$sample";
done < FileNames.txt

# Make directory for forward-read FASTQ's
mkdir fastq-pooled/

# Find any duplicated file names
find fastq \
  -type f -name '*R1*fastq*' \
  -exec basename {} \; |
sort |
uniq -d

# Move all forward-read FASTQ's with unique names to fastq-pooled/
find fastq/ -type f -name '*R1*fastq*' \
 -exec mv -n {} fastq-pooled/ \;

# Find any file names shared by 2+ files still in fastq/
find fastq \
  -type f -name '*R1*fastq*' \
  -exec basename {} \; |
sort |
uniq -d

# Take output of 'Find any duplicated file names' and put in this file, EXCEPT the output of 'Find any file names shared by 2+ files still in fastq/'
cd fastq/
touch duplicates.txt

# Manually rename and move output of 'Find any file names shared by 2+ files still in fastq/'
find . -type f -name '10317.000001507_0_L001_R1_001.fastq.gz' -print
find . -type f -name '10317.000001619_1_L001_R1_001.fastq.gz' -print
mv -i ./2363-Pace_GI_S1/9696a6c1-e8e3-49c3-8b6a-0dfbf6690aab/data/10317.000001507_0_L001_R1_001.fastq.gz ../fastq-pooled/DUPLICATE1.10317.000001507_0_L001_R1_001.fastq.gz
mv -i ./2201-AG_31_S1_L001/7dc9b5a0-ccb8-4591-8727-c6b23df9392e/data/10317.000001507_0_L001_R1_001.fastq.gz ../fastq-pooled/DUPLICATE2.10317.000001507_0_L001_R1_001.fastq.gz
mv -i ./2253-AG32_S1_L001/4c8bf2c8-60c3-4aab-8cc9-d2d2d86fa01b/data/10317.000001507_0_L001_R1_001.fastq.gz ../fastq-pooled/DUPLICATE3.10317.000001507_0_L001_R1_001.fastq.gz
mv -i ./2201-AG_31_S1_L001/7dc9b5a0-ccb8-4591-8727-c6b23df9392e/data/10317.000001619_1_L001_R1_001.fastq.gz ../fastq-pooled/DUPLICATE1.10317.000001619_1_L001_R1_001.fastq.gz
mv -i ./2253-AG32_S1_L001/4c8bf2c8-60c3-4aab-8cc9-d2d2d86fa01b/data/10317.000001619_1_L001_R1_001.fastq.gz ../fastq-pooled/DUPLICATE2.10317.000001619_1_L001_R1_001.fastq.gz

# Find all remaining duplicate files
while IFS= read -r sample || [ -n "$sample" ];
do
  find . -type f -name "$sample" -print;
done < duplicates.txt

# Put output of 'Find all remaining duplicate files' into this file
touch duplicates-fullpath.txt

# Move duplicates from subdirectories to fastq/
while IFS= read -r fullpath || [ -n "$fullpath" ];
do
  mv -i -- "$fullpath" ./;
done < duplicates-fullpath.txt

# Confirm pulled all duplicates
while IFS= read -r sample || [ -n "$sample" ];
do
  find . -mindepth 2 -type f -name "$sample" -print;
done < duplicates.txt

# Move duplicates to fastq-pooled/
while IFS= read -r sample || [ -n "$sample" ];
do
  mv -i -- "$sample" "../fastq-pooled/DUPLICATE.$sample";
done < duplicates.txt

# Confirm moved all duplicates
find . -type f -name '*R1*fastq*' -print
find . -type f -name '*fastq*' -print

# Verify fastq file are Casava 1.8 demultiplexed format
cd ../fastq-pooled/
find . -type f ! -name '*_*_L001_R1_001.fastq.gz' -print

# Keep only the samples which weren't filtered out in Ademola-Popoola et al.
cd ../
grep '^10317' ../sample_metadata.tsv | cut -f1 > SampleNamesFinal.txt
wc -l SampleNamesFinal.txt
# There are 3289 samples to retain
mkdir fastq-final/
while IFS= read -r sample || [ -n "$sample" ];
do
  mv -i -- "fastq-pooled/${sample}_"* fastq-final/;
done < SampleNamesFinal.txt
sed 's/^/DUPLICATE./' SampleNamesFinal.txt > DuplicatedSampleNames.txt
while IFS= read -r sample || [ -n "$sample" ];
do
  matches=( "fastq-pooled/${sample}_"* )
  if [ -e "${matches[0]}" ]; then
    mv -i -- "${matches[@]}" fastq-final/
  fi
done < DuplicatedSampleNames.txt

# Check if there are any samples from SampleNamesFinal.txt which do not have a matched file in fastq-final/
shopt -s nullglob
while IFS= read -r sample || [ -n "$sample" ];
do
  matches=( fastq-final/"$sample"* )
  if [ "${#matches[@]}" -eq 0 ]; then
    printf '%s\n' "$sample"
  fi
done < SampleNamesFinal.txt

# Check if there are any duplicates in fastq-final/
shopt -s nullglob
while IFS= read -r sample || [ -n "$sample" ];
do
  matches=( "fastq-final/${sample}_"* )
  if [ -e "${matches[0]}" ] && [ "${#matches[@]}" -gt 1 ]; then
    printf '%s\n' "${matches[@]}"
  fi
done < SampleNamesFinal.txt

# Drop duplicate files
# Put filenames (only!) from 'Check if there are any duplicates in fastq-final/' into the below file
touch MultiRun-SampleNames.txt
sed -i '' 's|fastq-final/||g' MultiRun-SampleNames.txt
mkdir fastq-multiruns/
while IFS= read -r sample || [ -n "$sample" ];
do
mv -i -- "fastq-final/$sample" "fastq-multiruns/DUPLICATE.$sample";
done < MultiRun-SampleNames.txt

# Import using Casava 1.8 file format protocol
cd ../
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path demux_industrial/fastq-final/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux_2026/demux_industrial.qza

## 2. Ethiopia ##
# Extract metadata and .fastq.gz files from healthy patients (ID's: EOSC_148 - EOSC_286)
tar -tzf EOSC_weyrich.tar.gz |
awk -F'[_/.]' '
    $0 == "reads_to_weyrich/EOSC_metadata.tsv" ||
    ($1 == "reads" && $2 == "to" && $3 == "weyrich" && $4 == "EOSC" &&
     $5 >= 148 && $5 <= 286 &&
     $6 ~ /^[FR]$/ &&
     $7 == "fastq" && $8 == "gz")
' |
tar -xzf EOSC_weyrich.tar.gz -C ./ -T /dev/stdin &&
mv reads_to_weyrich/ demux_ethiopia/

# Clean up metadata file
cd demux_ethiopia/
grep -v 'Diseased' EOSC_metadata.tsv | cut -f1-2,4-5,7-9 > meta_output.tsv &&
mv meta_output.tsv EOSC_metadata.tsv

# Create file with sample names
grep -v 'SampleID' EOSC_metadata.tsv | cut -f1 > samplenames.txt

# Create QIIME manifest file
printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > qiime_manifest.txt
sed -n 'l' qiime_manifest.txt
while IFS= read -r sample || [[ -n "$sample" ]];
    do
        printf '%s\t$PWD/%s_F.fastq.gz\t$PWD/%s_R.fastq.gz\n' \
            "$sample" "$sample" "$sample";
done < samplenames.txt >> qiime_manifest.txt
head qiime_manifest.txt | sed -n 'l'

# Determine Phred offset
while IFS= read -r sample || [[ -n "$sample" ]];
    do
        vsearch --fastq_chars "${sample}_F.fastq.gz" 2>&1 |
            grep '^Guess:';
        vsearch --fastq_chars "${sample}_R.fastq.gz" 2>&1 |
            grep '^Guess:';
done < samplenames.txt

# Import with Phred+33 offset
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path qiime_manifest.txt \
  --output-path ../demux_2026/demux_ethiopia.qza \
  --input-format PairedEndFastqManifestPhred33V2

## 3. Malaysia ##
#Downloaded MetaOral15032021.csv and QC_Meta09032021.csv from https://doi.org/10.26180/16959583
#Selected relevant QC runs from QC_Meta09032021.csv (BodySite: NA or Oral) and added them to MetaOral15032021.csv, keeping columns SampleID, Code, and SampleType (merging that column with OA.Group from MetaOral15032021.csv and keeping it as SampleType)
#Added "-repeat" to the end of the SampleID for repeated samples; changed "Kit" to "OralKit"; filled out remaining cells with NA; deleted QC_Meta09032021.csv; renamed MetaOral15032021.csv as OrangAsli.csv
#OrangAsli.csv lists all oral samples

#Located .fastq.gz files at https://www.ebi.ac.uk/ena/browser/view/PRJNA778070
#I could not find matching .fastq.gz files for some samples; these are in OrangAsli.csv without a listed run accession
#Sample accessions have sample names (as "OA Group")
#Downloaded ENA_report.tsv with all sample accessions, run accessions, and FASTQ file paths

#Many samples are incompatible with further analysis, including Jehai samples and those QC runs with "LF-"
#Repeated samples were also not carried forward
#Identified relevant samples ONLY, by sample accession, and placed these sample accessions in samplenames.txt

cd ../demux_malaysia/

# Select the relevant rows in ENA_report.tsv from samplenames.txt
touch oralsamples.txt
while IFS= read -r sample || [[ -n "$sample" ]];
do
  grep -F -- "$sample" ENA_report.tsv >> oralsamples.txt;
done < samplenames.txt

# Create download script in ENA format
cut -f3- oralsamples.txt > getoralfastqs.sh
tr ';' '\n' < getoralfastqs.sh > output.txt && mv output.txt getoralfastqs.sh
sed -i '' 's|^|wget -nc ftp://|' getoralfastqs.sh

# Run download script
bash getoralfastqs.sh

# Convert OrangAsli.csv to .tsv and clean up
tr ',' $'\t' < OrangAsli.csv > OA_metadata.tsv
grep -E 'RunAccession|TemiarGM|Temuan|SRR16830035|SRR16830034' OA_metadata.tsv |
    grep -v $'^\t' > meta_output.tsv &&
    mv meta_output.tsv OA_metadata.tsv
cut -f1-2,4- OA_metadata.tsv > output.tsv &&
    mv output.tsv OA_metadata.tsv

# Create QIIME manifest file
printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > qiime_manifest.txt
sed -n 'l' qiime_manifest.txt
while IFS= read -r sample || [[ -n "$sample" ]];
    do
        printf '%s\t$PWD/%s_1.fastq.gz\t$PWD/%s_2.fastq.gz\n' \
            "$sample" "$sample" "$sample";
done < samplenames.txt >> qiime_manifest.txt
head qiime_manifest.txt | sed -n 'l'

# Determine Phred offset
while IFS= read -r sample || [[ -n "$sample" ]];
    do
        vsearch --fastq_chars "${sample}_1.fastq.gz" 2>&1 |
            grep '^Guess:';
        vsearch --fastq_chars "${sample}_2.fastq.gz" 2>&1 |
            grep '^Guess:';
done < samplenames.txt

# Import with Phred+33 offset
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path qiime_manifest.txt \
  --output-path ../demux_2026/demux_malaysia.qza \
  --input-format PairedEndFastqManifestPhred33V2

###### These reads have been pre-trimmed to 125 bp and so I did not carry them forward ######

## 4. Bolivia ##
# Selected only samples MSA12, MSA14, and MSA15
# MSA13 has ample caries and is the sister of MSA15
# Other samples are not expectorated saliva

cd ../demux_bolivia/

# Download FASTQ files
bash Bolivia_FASTQ_script.sh

# Rename files
while IFS=$'\t' read -r accession sample || [[ -n "$accession" || -n "$sample" ]];
do
    for direction in 1 2;
    do
        old="demux_bolivia/${accession}_${direction}.fastq.gz"
        new="demux_bolivia/${sample}_${direction}.fastq.gz"
        mv "$old" "$new"
    done
done < demux_bolivia/Bolivia_FASTQ.tsv

# Create QIIME manifest file
printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > qiime_manifest.txt
sed -n 'l' qiime_manifest.txt
while IFS=$'\t' read -r discard sample || [[ -n "$sample" ]];
    do
        printf '%s\t$PWD/%s_1.fastq.gz\t$PWD/%s_2.fastq.gz\n' \
            "$sample" "$sample" "$sample";
done < Bolivia_FASTQ.tsv >> qiime_manifest.txt
head qiime_manifest.txt | sed -n 'l'

# Determine Phred offset
while IFS=$'\t' read -r discard sample || [[ -n "$sample" ]];
    do
        vsearch --fastq_chars "${sample}_1.fastq.gz" 2>&1 |
            grep '^Guess:';
        vsearch --fastq_chars "${sample}_2.fastq.gz" 2>&1 |
            grep '^Guess:';
done < Bolivia_FASTQ.tsv

# Import with Phred+33 offset
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path qiime_manifest.txt \
  --output-path ../demux_2026/demux_bolivia.qza \
  --input-format PairedEndFastqManifestPhred33V2

## 5. Ecuador ##
# Sample ID's are NEM1-NEM5
# Downloaded from https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1097549 using sra-tools

cd ../demux_ecuador/

# Rename files
while IFS=$'\t' read -r accession sample || [[ -n "$accession" || -n "$sample" ]];
do
    for direction in 1 2;
    do
        old="demux_ecuador/${accession}_${direction}.fastq.gz"
        new="demux_ecuador/${sample}_${direction}.fastq.gz"
        mv "$old" "$new"
    done
done < demux_ecuador/Ecuador_FASTQ.tsv

# Create QIIME manifest file
printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > qiime_manifest.txt
sed -n 'l' qiime_manifest.txt
while IFS=$'\t' read -r discard sample || [[ -n "$sample" ]];
    do
        printf '%s\t$PWD/%s_1.fastq.gz\t$PWD/%s_2.fastq.gz\n' \
            "$sample" "$sample" "$sample";
done < Ecuador_FASTQ.tsv >> qiime_manifest.txt
head qiime_manifest.txt | sed -n 'l'

# Determine Phred offset
while IFS=$'\t' read -r discard sample || [[ -n "$sample" ]];
    do
        vsearch --fastq_chars "${sample}_1.fastq.gz" 2>&1 |
            grep '^Guess:';
        vsearch --fastq_chars "${sample}_2.fastq.gz" 2>&1 |
            grep '^Guess:';
done < Ecuador_FASTQ.tsv

# Import with Phred+33 offset
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path qiime_manifest.txt \
  --output-path ../demux_2026/demux_ecuador.qza \
  --input-format PairedEndFastqManifestPhred33V2

# Most sequences are chimeric; I attempted to save these files with cutadapt
cd ../demux_2026/
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux_ecuador.qza \
  --p-front-f TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG \
  --p-front-r GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC \
  --p-match-adapter-wildcards \
  --p-match-read-wildcards \
  --p-cores 4 \
  --o-trimmed-sequences demux_ecuador_cutadapt.qza

###### While cutadapt successfully removed a partial primer sequence, the taxonomic mapping from these files was very poor and so they were not carried forward ######

## 6. Peru ##
cd ../

# Import files
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path Fastq_Raw/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux_2026/demux_peru.qza

# All projects imported
# Carrying forward eight of them

#### Visualizing demux output and denoising forward reads, trimmed to 150 bp
while IFS= read -r project || [ -n "$project" ];  
do  
qiime demux summarize \
  --i-data demux_2026/demux_$project.qza \
  --o-visualization demux_2026/demux_$project.qzv;  
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux_2026/demux_$project.qza \
  --p-trunc-len 150 \
  --o-table dada2/table_$project.qza \
  --o-representative-sequences dada2/rep-seqs_$project.qza \
  --o-denoising-stats dada2/denoising-stats_$project.qza \
  --o-base-transition-stats dada2/base-transition-stats_$project.qza;  
qiime metadata tabulate \
  --m-input-file dada2/denoising-stats_$project.qza \
  --o-visualization dada2/denoising-stats_$project.qzv;  
qiime feature-table summarize \
  --i-table dada2/table_$project.qza \
  --o-summary dada2/table_$project.qzv \
  --o-sample-frequencies dada2/sample-freqs_$project.qza \
  --o-feature-frequencies dada2/feature-freqs_$project.qza;  
done < dada2/FileNames.txt

#### Merging feature tables and rep seqs

# Merge all raw feature tables
qiime feature-table merge \
  --i-tables dada2/table_bolivia.qza \
  --i-tables dada2/table_ethiopia.qza \
  --i-tables dada2/table_industrial.qza \
  --i-tables dada2/table_nepal.qza \
  --i-tables dada2/table_peru.qza \
  --i-tables dada2/table_tanzania.qza \
  --i-tables dada2/table_uganda.qza \
  --i-tables dada2/table_venezuela.qza \
  --o-merged-table dada2/merged-table_global.qza  

# Filter out samples not retained in final metadata file
qiime feature-table filter-samples \
  --i-table dada2/merged-table_global.qza \
  --m-metadata-file sample_metadata.tsv \
  --o-filtered-table dada2/filtered-table_global.qza  

# Filter out low-frequency samples and summarize final table
qiime feature-table filter-features \
  --i-table dada2/filtered-table_global.qza \
  --p-min-frequency 10 \
  --o-filtered-table final-table_global.qza
qiime feature-table summarize \
  --i-table final-table_global.qza \
  --m-metadata-file sample_metadata.tsv \
  --o-summary final-table_global.qzv \
  --o-sample-frequencies final-sample-freqs_global.qza \
  --o-feature-frequencies final-feature-freqs_global.qza

# Merge representative sequences and filter by final feature table
qiime feature-table merge-seqs \
  --i-data dada2/rep-seqs_bolivia.qza \
  --i-data dada2/rep-seqs_ethiopia.qza \
  --i-data dada2/rep-seqs_industrial.qza \
  --i-data dada2/rep-seqs_nepal.qza \
  --i-data dada2/rep-seqs_peru.qza \
  --i-data dada2/rep-seqs_tanzania.qza \
  --i-data dada2/rep-seqs_uganda.qza \
  --i-data dada2/rep-seqs_venezuela.qza \
  --o-merged-data dada2/merged-rep-seqs_global.qza  
qiime feature-table filter-seqs \
  --i-data dada2/merged-rep-seqs_global.qza \
  --i-table final-table_global.qza \
  --o-filtered-data final-rep-seqs_global.qza
qiime feature-table tabulate-seqs \
  --i-data final-rep-seqs_global.qza \
  --o-visualization final-rep-seqs_global.qzv

# Create phylogeny and taxonomy files (using SILVA: https://www.arb-silva.de/documentation/classification/qiime-2)
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences final-rep-seqs_global.qza \
  --output-dir phylogeny  
qiime feature-classifier classify-sklearn \
  --i-classifier SILVA_144_SSURef_NR99_uniform_classifier_V4-515f-806r.qza \
  --i-reads final-rep-seqs_global.qza \
  --o-classification silva_taxonomy.qza  
qiime metadata tabulate \
  --m-input-file silva_taxonomy.qza \
  --o-visualization silva_taxonomy.qzv

# Export feature table as .biom for decontam
qiime tools export \
  --input-path final-table_global.qza \
  --output-path .
biom convert -i feature-table.biom -o table_from_biom.txt --to-tsv
biom convert -i table_from_biom.txt -o json.biom --table-type="OTU table" --to-json


#### Decontam (see R script)

#### Post-decontam
# Import back to QIIME and filter controls and Eukaryotes
biom convert -i species_decontam_S_60.txt -o decontam_table.biom --table-type="OTU table" --to-hdf5
qiime tools import \
  --input-path decontam_table.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path decontam_table.qza
qiime feature-table filter-samples \
  --i-table decontam_table.qza \
  --m-metadata-file sample_metadata.tsv \
  --p-where "SampleType='Sample'" \
  --o-filtered-table decontam_nocontrol.qza
qiime taxa filter-table \
  --i-table decontam_nocontrol.qza \
  --i-taxonomy silva_taxonomy.qza \
  --p-exclude Eukaryota,mitochondria,chloroplast \
  --o-filtered-table nocontrol-noEMC.qza
qiime feature-table summarize \
  --i-table nocontrol-noEMC.qza \
  --m-metadata-file sample_metadata.tsv \
  --o-summary nocontrol-noEMC.qzv \
  --o-sample-frequencies sample-freqs_nocontrol-noEMC.qza \
  --o-feature-frequencies feature-freqs_nocontrol-noEMC.qza

# Rarefy
# Set max depth to the highest level that retains all projects
qiime diversity alpha-rarefaction \
  --i-table nocontrol-noEMC.qza  \
  --i-phylogeny phylogeny/rooted_tree.qza \
  --p-max-depth 37789 \
  --p-steps 75 \
  --m-metadata-file sample_metadata.tsv \
  --o-visualization alphararefaction.qzv
# Choose to rarefy to 7,379 in order to retain all 3 Bolivian samples

# Analytics
# Taxa barplot
qiime taxa barplot \
  --i-table nocontrol-noEMC.qza \
  --i-taxonomy silva_taxonomy.qza \
  --m-metadata-file sample_metadata.tsv \
  --o-visualization taxonomy_barplot.qzv

# Core alpha- and beta-diversity metrics
# By country, continent, and mode of subsistence
# Additionally by lysis mechanism and sequencer
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny phylogeny/rooted_tree.qza \
  --i-table nocontrol-noEMC.qza \
  --p-sampling-depth 7379 \
  --m-metadata-file sample_metadata.tsv \
  --output-dir coremetrics  
qiime diversity alpha-group-significance \
  --i-alpha-diversity coremetrics/faith_pd_vector.qza \
  --m-metadata-file sample_metadata.tsv \
  --o-visualization coremetrics/faithpd_significance.qzv  
qiime diversity alpha-group-significance \
  --i-alpha-diversity coremetrics/evenness_vector.qza \
  --m-metadata-file sample_metadata.tsv \
  --o-visualization coremetrics/evenness_significance.qzv  
qiime diversity beta-group-significance \
  --i-distance-matrix coremetrics/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample_metadata.tsv \
  --m-metadata-column country \
  --p-pairwise \
  --o-visualization coremetrics/beta_significance-unweighted-unifrac.qzv  
qiime diversity beta-group-significance \
  --i-distance-matrix coremetrics/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample_metadata.tsv \
  --m-metadata-column Continent \
  --p-pairwise \
  --o-visualization coremetrics/beta_significance-unweighted-unifrac_continent.qzv  
qiime diversity beta-group-significance \
  --i-distance-matrix coremetrics/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample_metadata.tsv \
  --m-metadata-column Subsistence \
  --p-pairwise \
  --o-visualization coremetrics/beta_significance-unweighted-unifrac_subsistence.qzv  
qiime diversity beta-group-significance \
  --i-distance-matrix coremetrics/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample_metadata.tsv \
  --m-metadata-column ExtractionMethod \
  --o-visualization coremetrics/beta_significance-unweighted-unifrac_extraction.qzv  
qiime diversity beta-group-significance \
  --i-distance-matrix coremetrics/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample_metadata.tsv \
  --m-metadata-column LysisMechanism \
  --o-visualization coremetrics/beta_significance-unweighted-unifrac_lysis.qzv  
qiime diversity beta-group-significance \
  --i-distance-matrix coremetrics/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample_metadata.tsv \
  --m-metadata-column Sequencing_Method \
  --o-visualization coremetrics/beta_significance-unweighted-unifrac_sequencing_length.qzv  
qiime diversity beta-group-significance \
  --i-distance-matrix coremetrics/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample_metadata.tsv \
  --m-metadata-column Sequencer \
  --p-pairwise \
  --o-visualization coremetrics/beta_significance-unweighted-unifrac_sequencer.qzv  
qiime diversity adonis \
  --i-distance-matrix coremetrics/bray_curtis_distance_matrix.qza \
  --m-metadata-file sample_metadata.tsv \
  --p-formula "country" \
  --o-visualization coremetrics/adonis.qzv  
qiime diversity adonis \
  --i-distance-matrix coremetrics/bray_curtis_distance_matrix.qza \
  --m-metadata-file sample_metadata.tsv \
  --p-formula "Continent" \
  --o-visualization coremetrics/adonis_continent.qzv  
qiime diversity adonis \
  --i-distance-matrix coremetrics/bray_curtis_distance_matrix.qza \
  --m-metadata-file sample_metadata.tsv \
  --p-formula "Subsistence" \
  --o-visualization coremetrics/adonis_subsistence.qzv  
qiime diversity adonis \
  --i-distance-matrix coremetrics/bray_curtis_distance_matrix.qza \
  --m-metadata-file sample_metadata.tsv \
  --p-formula "ExtractionMethod" \
  --o-visualization coremetrics/adonis_extraction.qzv  
qiime diversity adonis \
  --i-distance-matrix coremetrics/bray_curtis_distance_matrix.qza \
  --m-metadata-file sample_metadata.tsv \
  --p-formula "LysisMechanism" \
  --o-visualization coremetrics/adonis_lysis.qzv  
qiime diversity adonis \
  --i-distance-matrix coremetrics/bray_curtis_distance_matrix.qza \
  --m-metadata-file sample_metadata.tsv \
  --p-formula "Sequencing_Method" \
  --o-visualization coremetrics/adonis_sequencing_length.qzv  
qiime diversity adonis \
  --i-distance-matrix coremetrics/bray_curtis_distance_matrix.qza \
  --m-metadata-file sample_metadata.tsv \
  --p-formula "Sequencer" \
  --o-visualization coremetrics/adonis_sequencer.qzv  

# Robust Principal Component Analysis on Aitchison Distances
mkdir rpca/  

# Modify feature table and taxonomy files to generate human-readable biplot labels
# See biplot_labels_as_taxa_names.py

qiime gemelli rpca \
  --i-table nocontrol-noEMC_with_names.qza \
  --p-min-feature-count 5 \
  --p-min-sample-count 5 \
  --o-biplot rpca/rpca_biplot_with_names.qza \
  --o-distance-matrix rpca/rpca_distance_with_names.qza  
qiime emperor biplot \
  --i-biplot rpca/rpca_biplot_with_names.qza \
  --m-sample-metadata-file sample_metadata.tsv \
  --m-feature-metadata-file silva_taxonomy.qza \
  --o-visualization rpca/rpca_biplot_with_names.qzv \
  --p-number-of-features 8  
qiime diversity beta-group-significance \
  --i-distance-matrix rpca/rpca_distance.qza \
  --m-metadata-file sample_metadata.tsv \
  --m-metadata-column country \
  --p-pairwise \
  --o-visualization rpca/beta_significance_gemelli.qzv  
qiime diversity beta-group-significance \
  --i-distance-matrix rpca/rpca_distance.qza \
  --m-metadata-file sample_metadata.tsv \
  --m-metadata-column Continent \
  --o-visualization rpca/beta_significance_gemelli_continent.qzv  
qiime diversity beta-group-significance \
  --i-distance-matrix rpca/rpca_distance.qza \
  --m-metadata-file sample_metadata.tsv \
  --m-metadata-column Subsistence \
  --o-visualization rpca/beta_significance_gemelli_subsistence.qzv  
qiime diversity beta-group-significance \
  --i-distance-matrix rpca/rpca_distance_with_names.qza \
  --m-metadata-file sample_metadata.tsv \
  --m-metadata-column ExtractionMethod \
  --o-visualization rpca/beta_significance_gemelli_extraction.qzv  
qiime diversity beta-group-significance \
  --i-distance-matrix rpca/rpca_distance_with_names.qza \
  --m-metadata-file sample_metadata.tsv \
  --m-metadata-column LysisMechanism \
  --o-visualization rpca/beta_significance_gemelli_lysis.qzv  
qiime diversity beta-group-significance \
  --i-distance-matrix rpca/rpca_distance_with_names.qza \
  --m-metadata-file sample_metadata.tsv \
  --m-metadata-column Sequencing_Method \
  --o-visualization rpca/beta_significance_gemelli_sequencing_length.qzv  
qiime diversity beta-group-significance \
  --i-distance-matrix rpca/rpca_distance_with_names.qza \
  --m-metadata-file sample_metadata.tsv \
  --m-metadata-column Sequencer \
  --p-pairwise \
  --o-visualization rpca/beta_significance_gemelli_sequencer.qzv  
