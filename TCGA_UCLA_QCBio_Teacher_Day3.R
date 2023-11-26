### Load our necessary libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)
library(umap)
library(survival)
library(survminer)

### Set up working directiory
setwd('/Users/zhuqian/Work/TCGA_UCLA_QCBio')
getwd()

### Let's build up query from yesterday's txn datasets
rna_query <- GDCquery(project = 'TCGA-READ',
                      data.category = 'Transcriptome Profiling',
                      access = 'open',
                      experimental.strategy = 'RNA-Seq',
                      sample.type = c('Primary Tumor','Solid Tissue Normal'))

rna <- GDCprepare(rna_query)

rna_matrix <- assay(rna,'unstranded')

### Differential expression analysis by DESeq2
dds_input <- DESeqDataSetFromMatrix(countData = rna_matrix,
                                    colData = colData(rna),
                                    design = ~sample_type)


# Filter low counts
keep <- rowSums(counts(dds_input)) >= 10

dds_input <- dds_input[keep,]

dds <- DESeq(dds_input)

# Variance stable transformation
vsd <- vst(dds, blind = FALSE)

# Dimentionality reduction by PCA
plotPCA(vsd,intgroup = 'sample_type')

# Psedudo single cell analysis by UMAP
vsd_counts <- assay(vsd)

umap <- umap(t(vsd_counts),
             n_neighbors = 15)

umap_res <- as.data.frame(umap$layout)

# Plot umap result
ggplot(data = umap_res,
       aes(x = V1, y = V2)) +
  geom_point() +
  xlab('umap1') +
  ylab('umap2')

# Label umap clusters by different criteria
umap_merged <- cbind(umap_res,colData(vsd))

ggplot(umap_merged, aes(x = V1, y = V2, color = umap_merged$sample_type)) +
  geom_point() +
  xlab('umap1') +
  ylab('umap2')

ggplot(umap_merged, aes(x = V1, y = V2, color = umap_merged$ajcc_pathologic_stage)) +
  geom_point() +
  xlab('umap1') +
  ylab('umap2')

ggplot(umap_merged, aes(x = V1, y = V2, color = umap_merged$prior_malignancy)) +
  geom_point() +
  xlab('umap1') +
  ylab('umap2')

ggplot(umap_merged, aes(x = V1, y = V2, color = umap_merged$gender)) +
  geom_point() +
  xlab('umap1') +
  ylab('umap2')

ggplot(data = umap_merged,
       aes(x = V1, y = V2, color = umap_merged$race)) +
  geom_point() +
  xlab('umap1') +
  ylab('umap2') 

colnames(umap_merged)

ggplot(data = umap_merged,
       aes(x = V1, y = V2, color = umap_merged$ajcc)) +
  geom_point() +
  xlab('umap1') +
  ylab('umap2') 

### Retrieve clinical data
clinical_query <- GDCquery_clinic('TCGA-READ')

# Neat our clinical datasets
all(c('vital_status','days_to_last_follow_up','days_to_death') %in% colnames(clinical_query)) # to check the existance of our columns

which(colnames(clinical_query) %in% c('vital_status','days_to_last_follow_up','days_to_death')) # to check column orders of our data of interest

View(clinical_query[,c(9,107,113)])

table(clinical_query$vital_status)

clinical_query <- clinical_query[clinical_query$vital_status != 'Not Reported',] # filter out outlier datasets

table(clinical_query$vital_status)

clinical_query$deceased <- ifelse(clinical_query$vital_status=='Alive',FALSE,TRUE)

table(clinical_query$deceased)

clinical_query$overall_survival <- ifelse(clinical_query$deceased==FALSE, clinical_query$days_to_last_follow_up, clinical_query$days_to_death)

# Integrate with our vst txn datasets
gene_metadata <- as.data.frame(rowData(rna))

norm_counts <- counts(dds,normalized = T) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene_id') %>%
  gather(key = 'case_id', value = 'counts',-'gene_id') %>%# transform the data structure for downstream analysis
  left_join(.,gene_metadata,by='gene_id')

# Pick our gene of interest
goi <- norm_counts %>%
  filter(gene_name == 'CDKN1A')

# Filter out txn datasets from solid tissue normals
tumor_indices <- grep('-01A.*',goi$case_id)

tumor_goi <- goi[tumor_indices,]

# Separate our tumor group into 2 based on the median txn intensity of goi
median_goi <- median(tumor_goi$counts)
median_goi

tumor_goi$strata <- ifelse(tumor_goi$counts >= median_goi,'High','Low')

table(tumor_goi$strata)

# Merge clinical data
tumor_goi$submitter_id <- gsub('-01A.*','',tumor_goi$case_id) # trim off additional barcodes after patient id

tumor_goi <- merge(tumor_goi,clinical_query, by = 'submitter_id')

# Fitting our survival curve
fit <- survfit(Surv(overall_survival,deceased) ~strata, data = tumor_goi)

ggsurvplot(fit = fit,
           data = tumor_goi,
           pval = T,
           risk.table = T)
