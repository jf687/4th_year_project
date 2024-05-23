rm(list=ls())
library(dplyr)
# Download TCGA Breast Cancer proteomic data, and data on pathological tumour stage

# Load TCGA BRCA proteomic data --------------------------------------------------------------------------

# We use RPPA (with replicate-base normalization) data
# Data from https://xenabrowser.net/datapages/?dataset=TCGA.BRCA.sampleMap%2FRPPA_RBN&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
options(timeout=200)
utils::download.file("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FRPPA_RBN.gz", "Genomic_application/data/TCGA.BRCA.sampleMap%2FRPPA_RBN.gz")
BRCA_download <- utils::read.table(gzfile("Genomic_application/data/TCGA.BRCA.sampleMap%2FRPPA_RBN.gz"), header = T)
rownames(BRCA_download) = BRCA_download[,1]
BRCA_download = t(BRCA_download[,-1])
dim(BRCA_download) # 747 x 131
# Controlled that all Ids end with '01', i.e., no controls. 

# Load TCGA BRCA clinical data --------------------------------------------------------------------------
BRCA_clindat_download <- utils::read.csv("Genomic_application/data/TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix", sep='\t')
rownames(BRCA_clindat_download) = BRCA_clindat_download[,1] # Make patient ID be rownames
BRCA_clindat_download = BRCA_clindat_download[,-1]
names(BRCA_clindat_download) # Lots of other clinical variables to use instead if we want to
unique(BRCA_clindat_download$pathologic_stage) # The different tumour stages (higher=more severe)

# Get patient ID for tumours of stages I, II, and III
IDs.stageI = rownames(BRCA_clindat_download)[which(BRCA_clindat_download$pathologic_stage %in% c('Stage I', 'Stage IA', 'Stage IB'))]
IDs.stageII = rownames(BRCA_clindat_download)[which(BRCA_clindat_download$pathologic_stage %in% c('Stage II', 'Stage IIA', 'Stage IIB'))]
IDs.stageIII = rownames(BRCA_clindat_download)[which(BRCA_clindat_download$pathologic_stage %in% c('Stage III', 'Stage IIIA', 'Stage IIIB'))]

# If the above does not give anything interesting (too much heterogeneity): can consider only luminal cancers. Then replace the above with this:
#IDs.stageI = rownames(BRCA_clindat_download)[which(BRCA_clindat_download$pathologic_stage %in% c('Stage I', 'Stage IA', 'Stage IB') & 
#                                                     BRCA_clindat_download$PAM50_mRNA_nature2012 %in% c('Luminal A', 'Luminal B'))]
#IDs.stageII = rownames(BRCA_clindat_download)[which(BRCA_clindat_download$pathologic_stage %in% c('Stage II', 'Stage IIA', 'Stage IIB') & 
#                                                      BRCA_clindat_download$PAM50_mRNA_nature2012 %in% c('Luminal A', 'Luminal B'))]
#IDs.stageIII = rownames(BRCA_clindat_download)[which(BRCA_clindat_download$pathologic_stage %in% c('Stage III', 'Stage IIIA', 'Stage IIIB') & 
#                                                       BRCA_clindat_download$PAM50_mRNA_nature2012 %in% c('Luminal A', 'Luminal B'))]


# Get same ID format as proteomic data
IDs.stageI = gsub("-", ".",IDs.stageI)
IDs.stageII = gsub("-", ".",IDs.stageII)
IDs.stageIII = gsub("-", ".",IDs.stageIII)

# Get the corresponding IDs that are present in the proteomic data set (not all of the above patients have proteomic data)
IDs.stageI.rppa = rownames(BRCA_download)[which(rownames(BRCA_download) %in% IDs.stageI)]
IDs.stageII.rppa = rownames(BRCA_download)[which(rownames(BRCA_download) %in% IDs.stageII)]
IDs.stageIII.rppa = rownames(BRCA_download)[which(rownames(BRCA_download) %in% IDs.stageIII)]
length(IDs.stageI.rppa) # 113
length(IDs.stageII.rppa) # 429
length(IDs.stageIII.rppa) # 140

# Get the corresponding data sets
brca_dat_stageI = BRCA_download[which(rownames(BRCA_download) %in% IDs.stageI.rppa),]
brca_dat_stageII = BRCA_download[which(rownames(BRCA_download) %in% IDs.stageII.rppa),]
brca_dat_stageIII = BRCA_download[which(rownames(BRCA_download) %in% IDs.stageIII.rppa),]
dim(brca_dat_stageI) # 113 x 131
dim(brca_dat_stageII) # 429 x 131
dim(brca_dat_stageIII) # 140 x 131

brca_dat = list(brca_dat_stageI,brca_dat_stageII, brca_dat_stageIII)
brca_stages = c(1,2,3)

# Map the names of the antibodies used to identify the proteins to the names of the genes that encode them
# Use file showing which genes the proteins detected by the antibodies are encoded from (can be downloaded from TCGA website)
protein.names = colnames(BRCA_download)
rppa.to.gene = read.table("Genomic_application/data/RPPA_to_gene.txt", sep = "\t", stringsAsFactors = F)
mapping.frame = data.frame(protein = protein.names,gene=rppa.to.gene[match(protein.names,rppa.to.gene[,1]),2])
length(unique(mapping.frame$gene)) # 102 unique proteins

# Save all data
save(brca_dat, brca_stages, mapping.frame, file='Genomic_application/data/brca_dat_formatted.RData')


