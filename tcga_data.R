setwd("/Users/1Air/Documents/Insight/project")
library(RTCGAToolbox)
#source("https://bioconductor.org/biocLite.R")
#biocLite("TCGAbiolinks")
#source("https://bioconductor.org/biocLite.R")
#biocLite("TCGAbiolinks")
library(TCGAbiolinks)
library(dplyr)
#source("https://bioconductor.org/biocLite.R")
#biocLite("RTCGA")
library(RTCGA)
#install.packages("magrittr")
library(magrittr)

pacman::p_load('vtreat','caret','gbm','doParallel','data.table','dplyr','h2o')



# QUERY -----------------------------------------------

# project A list of valid project (see table below)
# data.category A valid project (see list with getProjectSummary(project))
# data.type A data type to filter the files to download
# sample.type A sample type to filter the files to download (See table below)
# workflow.type GDC workflow type
# barcode A list of barcodes to filter the files to download
# legacy Search in the legacy repository? Default: FALSE
# platform Experimental data platform (HumanMethylation450, AgilentG4502A_07 etc). Used only for legacy repository
# file.type A string to filter files, based on its names. Used only for legacy repository

query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")

query.mirna <- GDCquery(project = "TCGA-BRCA", 
                        data.category = "Transcriptome Profiling", 
                        data.type = "miRNA Expression Quantification",
                        sample.type = c("Primary solid Tumor","Solid Tissue Normal"))

query.exp.hg19 <- GDCquery(project = "TCGA-BRCA",
                           data.category = "Gene expression",
                           data.type = "Gene expression quantification",
                           platform = "Illumina HiSeq", 
                           file.type  = "normalized_results",
                           experimental.strategy = "RNA-Seq",
                           #barcode = c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01"),
                           legacy = TRUE)

getGDCprojects()
#getProjectSummary('TCGA-BRCA')


# CLINICAL -------------------------------------------

clin.query <- GDCquery(project = "TCGA-BLCA", data.category = "Clinical", barcode = "TCGA-FD-A5C0")
json  <- tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical.patient <- GDCprepare_clinic(clin.query, clinical.info = "patient")
clinical.patient.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up")
clinical.index <- GDCquery_clinic("TCGA-BLCA")

clin <- GDCquery_clinic("TCGA-ACC", type = "clinical", save.csv = TRUE)
clin <- GDCquery_clinic("TCGA-ACC", type = "biospecimen", save.csv = TRUE)

query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Clinical", 
                  barcode = c("TCGA-RU-A8FL","TCGA-AA-3972"))
GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info = "patient")
clinical.drug <- GDCprepare_clinic(query, clinical.info = "drug")
clinical.radiation <- GDCprepare_clinic(query, clinical.info = "radiation")
clinical.admin <- GDCprepare_clinic(query, clinical.info = "admin")

query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Biospecimen", 
                  barcode = c("TCGA-RU-A8FL","TCGA-AA-3972"))
GDCdownload(query)
clinical.admin <- GDCprepare_clinic(query, clinical.info = "admin")
clinical.sample <- GDCprepare_clinic(query, clinical.info = "sample")
clinical.slide <- GDCprepare_clinic(query, clinical.info = "slide")
clinical.portion <- GDCprepare_clinic(query, clinical.info = "portion")


# This code will get all clinical indexed data from TCGA
library(TCGAbiolinks)
library(data.table)
clinical <- TCGAbiolinks:::getGDCprojects()$project_id %>% 
  regexPipes::grep("TCGA",value=T) %>% 
  sort %>% 
  plyr::alply(1,GDCquery_clinic, .progress = "text") %>% 
  rbindlist
readr::write_csv(clinical,path = paste0("all_clin_indexed.csv"))

# This code will get all clinical XML data from TCGA
getclinical <- function(proj){
  message(proj)
  while(1){
    result = tryCatch({
      query <- GDCquery(project = proj, data.category = "Clinical")
      GDCdownload(query)
      clinical <- GDCprepare_clinic(query, clinical.info = "patient")
      for(i in c("admin","radiation","follow_up","drug","new_tumor_event")){
        message(i)
        aux <- GDCprepare_clinic(query, clinical.info = i)
        if(is.null(aux)) next
        # add suffix manually if it already exists
        replicated <- which(grep("bcr_patient_barcode",colnames(aux), value = T,invert = T) %in% colnames(clinical))
        colnames(aux)[replicated] <- paste0(colnames(aux)[replicated],".",i)
        if(!is.null(aux)) clinical <- merge(clinical,aux,by = "bcr_patient_barcode", all = TRUE)
      }
      readr::write_csv(clinical,path = paste0(proj,"_clinical_from_XML.csv")) # Save the clinical data into a csv file
      return(clinical)
    }, error = function(e) {
      message(paste0("Error clinical: ", proj))
    })
  }
}
clinical <- TCGAbiolinks:::getGDCprojects()$project_id %>% 
  regexPipes::grep("TCGA",value=T) %>% sort %>% 
  plyr::alply(1,getclinical, .progress = "text") %>% 
  rbindlist(fill = TRUE) %>% setDF %>% subset(!duplicated(clinical))
readr::write_csv(clinical,path = "all_clin_XML.csv")
# result: https://drive.google.com/open?id=0B0-8N2fjttG-WWxSVE5MSGpva1U
# Obs: this table has multiple lines for each patient, as the patient might have several followups, drug treatments,
# new tumor events etc...

?TCGAquery_SampleTypes
?TCGAquery_MatchedCoupledSampleTypes


bar <- c("TCGA-G9-6378-02A-11R-1789-07", "TCGA-CH-5767-04A-11R-1789-07")

S <- TCGAquery_SampleTypes(bar,"TP")
S2 <- TCGAquery_SampleTypes(bar,"NB")

# Retrieve multiple tissue types  NOT FROM THE SAME PATIENTS
SS <- TCGAquery_SampleTypes(bar,c("TP","NB"))

# Retrieve multiple tissue types  FROM THE SAME PATIENTS
SSS <- TCGAquery_MatchedCoupledSampleTypes(bar,c("NT","TP"))


# DOWNLOAD------------------------------------------


query <- GDCquery(project = "TCGA-BRCA", data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment",
                  barcode = c( "TCGA-OR-A5KU-01A-11D-A29H-01", "TCGA-OR-A5JK-01A-11D-A29H-01"))
GDCdownload(query)
data <- GDCprepare(query)

query_brca <- GDCquery(project = "TCGA-BRCA",
                           data.category = "Gene expression",
                           data.type = "Gene expression quantification",
                           platform = "Illumina HiSeq", 
                           file.type  = "normalized_results",
                           experimental.strategy = "RNA-Seq",
                           #barcode = c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01"),
                           legacy = TRUE)
GDCdownload(query_brca)
data <- GDCprepare(query_brca)
#-------------------------------------- 
# Gene expression
#-------------------------------------- 
# brca2
# mRNA pipeline: https://gdc-docs.nci.nih.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
query.brca2 <- GDCquery(project = "TCGA-BRCA", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - Counts"
                           #barcode =  c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01")
                        )
GDCdownload(query.brca2)
data <- GDCprepare(query = query.brca2,
                     save = TRUE, 
                     save.filename = "brca.rda")

# brca3
query.brca3 <- GDCquery(project = "TCGA-BRCA",
                           data.category = "Gene expression",
                           data.type = "Gene expression quantification",
                           platform = "Illumina HiSeq", 
                           file.type  = "normalized_results",
                           experimental.strategy = "RNA-Seq",
                           #barcode = c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01"),
                           legacy = TRUE)
GDCdownload(query.brca3)
data.brca <- GDCprepare(query.brca3,
                        save = TRUE, 
                        directory = GDCdata2,
                        save.filename = "brca2.rda")

# ANALYSIS ---------------------------------------------------

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
BRCARnaseqSE <- GDCprepare(query)

BRCAMatrix <- assay(BRCARnaseqSE,"raw_count") 
# or BRCAMatrix <- assay(BRCARnaseqSE,"raw_count")

# For gene expression if you need to see a boxplot correlation 
# and AAIC plot to define outliers you can run
BRCARnaseq_CorOutliers <- TCGAanalyze_Preprocessing(BRCARnaseqSE)


# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("NT"))

# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), 
                                   typesample = c("TP"))

# Diff.expr.analysis (DEA)
dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                            mat2 = dataFilt[,samplesTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")

# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"Tumor","Normal",
                                          dataFilt[,samplesTP],dataFilt[,samplesNT])

# FROM TXT -------------------------------------------

pacman::p_load('data.table','dplyr','ggplot2')

m_filename = '/Users/1Air/Documents/Insight/project/gdac.broadinstitute.org_BRCA.mRNA_Preprocess_Median.Level_3.2016012800.0.0/BRCA.medianexp.txt'
clin_filename = '/Users/1Air/Documents/Insight/project/gdac.broadinstitute.org_BRCA.Merge_Clinical.Level_1.2016012800.0.0/BRCA.clin.merged.txt'
r_filename = '/Users/1Air/Documents/Insight/project/gdac.broadinstitute.org_BRCA.Merge_rnaseq__illuminahiseq_rnaseq__unc_edu__Level_3__gene_expression__data.Level_3.2016012800.0.0 2/BRCA.rnaseq__illuminahiseq_rnaseq__unc_edu__Level_3__gene_expression__data.data.txt'

# rna gene expression quantification
rna <- fread(
  r_filename,
  header = T,
  sep = '\t',
  strip.white = TRUE,
  quote="",
  stringsAsFactors = FALSE,
  na.strings = c('NA', '')
)

# microarray gene expression quantification
marray <- fread(
  m_filename,
  header = T,
  sep = '\t',
  strip.white = TRUE,
  quote="",
  stringsAsFactors = FALSE,
  na.strings = c('NA', '')
)

# clinical data
clin <- fread(
  clin_filename,
  header = T,
  sep = '\t',
  strip.white = TRUE,
  quote="",
  stringsAsFactors = FALSE,
  na.strings = c('NA', '')
)

dim(rna) # [1] 20533  2635
dim(marray) # [1] 17815   591
dim(clin) # [1] 1097 3719


# CLEAN----------------------------------------------------

clin <- as.data.frame(t(clin))
rna <- as.data.frame(t(rna))
marray <- as.data.frame(t(marray))

# rna
#gene.bcr <- row.names(rna)[-1]
genenames <- as.character(unname(unlist(rna[1,1:ncol(rna)])))
colnames(rna) <- genenames
rna <- rna %>% mutate(pid=row.names(rna)) %>% filter(gene=='median_length_normalized') %>% select(-gene)
rna.bcr <- rna$pid
#gene.bcr <- gene.bcr[-1]
#rna.bcr <- unique(gene.bcr)
rm(gene.bcr.uniq)

# microarray
marray.bcr <- row.names(marray)[-1]
mnames <- as.character(unname(unlist(marray[1,1:ncol(marray)])))
colnames(marray) <- mnames
marray <- marray[-1,-1]
rm(mnames)
marray.bcr <- marray.bcr[-1]

# clinical
varnames <- as.character(unname(unlist(clin[1,1:ncol(clin)])))
colnames(clin) <- varnames
clin <- clin[-1,-1]

save(rna,file='rna.rda' )
save(marray,file='marray.rda' )
save(clin,file='clin.rda' )

# merge rna and marray if possible

# get column for join by deriving patient provenance across both datasets
# extract patient ID/ barcodes
clin.bcr <- aliq.bcra
rm(aliq.bcr)
#pat.bcr <- as.character(unname(unlist(clin[,19])))
clin.bcr <- as.character(unname(unlist(clin[,3040]))) #rna
pat.id <- as.character(unname(unlist(clin[,1347])))
# drg.bcr <- unname(unlist(clin[,546]))
# rad.bcr <- unname(unlist(clin[,1433]))
# samp.bcr <- unname(unlist(clin[,2854]))
# omf.bcr <- unname(unlist(clin[,3564]))

# selection of normal samples "NT"
rna.norm <- TCGAquery_SampleTypes(barcode = rna.bcr,
                                   typesample = c("NT"))
marray.norm <- TCGAquery_SampleTypes(barcode = marray.bcr,
                                  typesample = c("NT"))

# selection of tumor samples "TP"
clin.tumor <- TCGAquery_SampleTypes(barcode = clin.bcr, 
                                   typesample = c("TP"))
rna.tumor <- TCGAquery_SampleTypes(barcode = rna.bcr, 
                                   typesample = c("TP"))
marray.tumor <- TCGAquery_SampleTypes(barcode = marray.bcr, 
                                   typesample = c("TP"))

rna_id <- data.frame(rna.bcr, tolower(substr(rna.bcr,9,12)))
names(rna_id)[2] <- 'pat.id'
rna_id <- arrange(rna_id, pat.id)

mray_id <- data.frame(marray.bcr, tolower(substr(marray.bcr,9,12)))
names(mray_id)[2] <- 'pat.id'
mray_id <- arrange(mray_id, pat.id)

clin_id <- data.frame(pat.bcr, pat.id)
clin_id <- arrange(clin_id, pat.id)

merge_id <- left_join(clin_id,gene_id, by=c('pat.id'='gene.id')) %>% arrange(pat.id)
merge_id_gen <- left_join(gene_id,clin_id, by=c('gene.id'='pat.id')) %>% arrange(gene.id)
# now create join
genexp$patient.patient_id <- tolower(substr(row.names(genexp),9,12))
genedata <- left_join(genexp, clin)

# FEATURE EXTRACTION -----------------------------------------------
genedata$patient.age <- floor(as.numeric(as.character(genedata$patient.days_to_birth))/-365.4)
genedata$TP53 <- as.numeric(as.character(genedata$TP53))
genedata$BRCA1 <- as.numeric(as.character(genedata$BRCA1))
genedata$BRCA2 <- as.numeric(as.character(genedata$BRCA2))
rownames(genedata) <- gene.bcr

# class variable
#genexp$class <- genedata$patient.vital_status
clin$class <- factor( rep("late",nrow(clin)),ordered = F, levels = c('early','late'))
clin$class[clin$patient.stage_event.pathologic_stage=='i'
             | clin$patient.stage_event.pathologic_stage=='ia'
             | clin$patient.stage_event.pathologic_stage=='ib'
             | clin$patient.stage_event.pathologic_stage=='ii'
             | clin$patient.stage_event.pathologic_stage=='iia'
             | clin$patient.stage_event.pathologic_stage=='iib'
             | clin$patient.stage_event.pathologic_stage=='iii'
             | clin$patient.stage_event.pathologic_stage=='iiia'
             ] <- 'early'

table(clin$patient.stage_event.pathologic_stage)
# early/good prognosis i, ia, ib, iia, iib, iiia, iii, ii
# late/concern all else

data <- select(data, -IM_VisitType)

# head(colnames(genedata[,1:17814]))
# head "ELMO2"    "CREB3L1"  "RPS11"    "PNMA1"    "MMP2"     "C10orf90"
# tail "RXFP2"   "PIK3IP1" "SLC39A6" "SNRPD2"  "AQP7"    "CTSC" 
genedata[,1:17814] <- as.numeric(as.character(unlist(genedata[,1:17814])))



# VISUALIZE -----------------------------------------
plot(genedata$patient.vital_status)
plot(genedata$patient.person_neoplasm_cancer_status)
plot(genedata$patient.ethnicity)
plot(genedata$patient.race_list.race)
plot(genedata$patient.vital_status)
plot(genedata$patient.menopause_status)
plot(genedata$patient.histological_type)
plot(genedata$patient.tumor_tissue_site)
plot(genedata$patient.drugs.drug.drug_name)
plot(genedata$patient.radiation_therapy)
plot(genedata$patient.breast_carcinoma_estrogen_receptor_status)
plot(genedata$patient.breast_carcinoma_progesterone_receptor_status)
plot(genedata$patient.breast_carcinoma_surgical_procedure_name)
hist(genedata$patient.age)
plot(genedata$patient.stage_event.pathologic_stage)
# unique or table or levels

# ggplot2

# age histogram
g <- ggplot(genedata, aes(patient.age))
g + geom_histogram()

# boxplot age by vital status
g <- ggplot(genedata, aes(patient.vital_status, patient.age))
g + geom_boxplot(aes(fill=patient.vital_status))

# age by tumor stage

# age and tp53/cytp450 color by vital status
g <- ggplot(genedata, aes(BRCA1, BRCA2))
g + geom_point(aes(color=patient.vital_status))


# QUICK MODEL ---------------------------------

set.seed(123)
inTrain <- createDataPartition(y = data$Readmit, p = .80,list = FALSE)
train <- data[ inTrain,]
test <- data[-inTrain,]
nrow(train) # 57564
nrow(test) # 29654
plot(train$Readmit)
summary(train)

#fit1 <- lm(class~V1410,data=train)
#fit1 <- lm(paste0('class ~ ',paste(preds,collapse = '+')),data=train)
fit1 <- lm(class ~ .,data=train)
summary(fit1)

rf1 <- randomForest(class ~ .,
                    data=train)
print(rf1)
imps <- importance(rf1)


# KEY VARIABLES ----------------------------------

#aliq.bcr <- colnames(data)[2:ncol(data)]
#clin$patient.samples.sample.portions.portion.analytes.analyte-2.aliquots.aliquot-2.bcr_aliquot_barcode
#clin[patient.omfs.omf.omf_staging.pathologic_stage,] # stage/class early or late, preventable/managmt

#patient.drugs.drug-2.drug_name
#patient.drugs.drug-2.regimen_indication
#patient.drugs.drug-2.therapy_types.therapy_type
#patient.drugs.drug-2.days_to_drug_therapy_start
#patient.drugs.drug-2.days_to_drug_therapy_start
# up to 23 drugs per patient

#patient.follow_ups.follow_up.brc_followup_barcode
# patient.race_list.race
# patient.ethnicity
# patient.patientid    --part of barcode
# patient.person_neoplasm_cancer_status
# patient.other_dx
# patient.tissue_source_site  --part of barcode
# patient.vital_status
# patient.menopause_status
# patient.histological_type
# patient.icd_10

# patient.normal_controls

#patient.radiation_therapy
#patient.radiations.radiation.bcr_radiation.barcode
#patient.radiations.radiation.anatomic_treatment_site
# up to 5 rads

# patient.samples.sample.bcr_sample_barcode
# patient.samples.sample.pathology_report_file_name
# patient.samples.sample-2.portions.portion.analytes.analyte-2.aliquots.aliquot-2.bcr_aliquot_barcode
# patient.samples.sample-2.portions.portion.analytes.analyte-2.aliquots.aliquot-2.plate_id   --- part of barcode
# patient.samples.sample.portions.portion.slides.slide.bcr_slide_barcode
# patient.samples.sample.portions.portion.slides.slide.image_file_name

# patient.omfs.omf.bcr_omf_barcode

# patient.breast_carcinoma_estrogen_receptor_status
# patient.breast_carcinoma_primary_surgical_procedure_name
# patient.breast_carcinoma_progesterone_receptor_status
# patient.breast_carcinoma_surgical_procedure_name
# patient.days_to_birth
# patient.days_to_death
# patient.stage_event.pathologic_stage
# patient.stage_event.tnm_categories.pathologic_categories.pathologic_m
# patient.stage_event.tnm_categories.pathologic_categories.pathologic_n
# patient.stage_event.tnm_categories.pathologic_categories.pathologic_t

