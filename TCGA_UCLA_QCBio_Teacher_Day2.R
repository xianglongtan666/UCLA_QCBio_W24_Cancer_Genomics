### Install necessary packages
if(!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

library(BiocManager)

BiocManager::install('BioinformaticsFMRP/TCGAbiolinksGUI.data')
BiocManager::install('BioinformaticsFMRP/TCGAbiolinks')
BiocManager::install('tidyverse')
BiocManager::install('SummarizedExperiment') # for data mining and data structure organizing

BiocManager::install('maftools') # for MAF file analysis

BiocManager::install('DESeq2') # for differential expression analysis
BiocManager::install('EnhancedVolcano') # DEG data visualization

BiocManager::install('umap') # for UMAP projection

BiocManager::install('survival')
BiocManager::install('survminer') # survival analysis

### Load necessary packages for day 2
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(maftools)
library(DESeq2)
library(EnhancedVolcano)

### Set up working directory
path <- '/Users/zhuqian/Work/TCGA_UCLA_QCBio'
setwd(path)
getwd()

### Get table of contents from TCGAbiolinks
gdc_project <- getGDCprojects()
getProjectSummary('TCGA-READ')

### Build a query
snv_query <- GDCquery(project = 'TCGA-READ',
                      data.category = 'Simple Nucleotide Variation')
output_snv_query <- getResults(snv_query)
View(output_snv_query)

### Filter the query for masked somatic mutations
snv_query <- GDCquery(project = 'TCGA-READ',
                      data.category = 'Simple Nucleotide Variation',
                      access = 'open',
                      data.type = 'Masked Somatic Mutation',
                      data.format = 'MAF',
                      barcode = c('TCGA-AF-2689-01A-01D-1989-10,TCGA-AF-2689-10A-01D-1989-10','TCGA-EI-6512-01A-11D-1733-10,TCGA-EI-6512-10A-01D-1733-10','TCGA-AF-2693-01A-02D-1733-10,TCGA-AF-2693-10A-01D-1733-10','TCGA-AF-3913-01A-02W-1073-09,TCGA-AF-3913-11A-01W-1073-09','TCGA-DC-6154-01A-31D-1924-10,TCGA-DC-6154-10A-01D-1924-10'))
output_snv_query <- getResults(snv_query)

### Retrieve those MAF datasets
GDCdownload(snv_query)

### Prepare data into R
maf <- GDCprepare(snv_query)

### Visualize mutational datasets by maftools package
maf <- read.maf(maf = maf)

getSampleSummary(maf)

plotmafSummary(maf = maf,
               rmOutlier = TRUE,
               dashboard = TRUE) # Plot summary plots for somatic mutations

oncoplot(maf = maf,
         top = 20) # plot top mutated genes

titv(maf = maf,
     plot = TRUE,
     useSyn = TRUE) # transition vs. transversion

### Query, retrieve and load txn datasets from TCGA-READ project
getProjectSummary('TCGA-READ')

rna_query <- GDCquery(project = 'TCGA-READ',
                      data.category = 'Transcriptome Profiling',
                      access = 'open',
                      experimental.strategy = 'RNA-Seq',
                      sample.type = c('Primary Tumor','Solid Tissue Normal'))

output_rna_query <- getResults(rna_query)

unique(output_rna_query$sample_type)

GDCdownload(rna_query)

rna <- GDCprepare(rna_query)

rna_matrix <- assay(rna,'unstranded')

### Differential expression analysis by DESeq2
dds_input <- DESeqDataSetFromMatrix(countData = rna_matrix,
                                    colData = colData(rna),
                                    design = ~sample_type)

# Filter low counts
keep <- rowSums(counts(dds_input)) >= 10

dds_input <- dds_input[keep,]

dds <- DESeq(dds_input) # Run DESeq2

res <- results(dds) # Get DESeq2 result

res

summary(res)

res_df <- as.data.frame(res) # Transfer res to a dataframe

write.csv(res_df, file = 'DEA.csv')

### Transfer gene ID
gene_metadata <- as.data.frame(rowData(rna))

res_df$gene_id <- rownames(res_df)

res_df <- left_join(x = res_df,
                    y = gene_metadata,
                    by = 'gene_id')

dim(res_df)

### Data visualization
EnhancedVolcano(toptable = res_df,
                lab = res_df$gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1)
