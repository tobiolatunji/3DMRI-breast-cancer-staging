setwd("/Users/1Air/Documents/Insight/project")
library(RTCGAToolbox)
#source("https://bioconductor.org/biocLite.R")
#biocLite("TCGAbiolinks")
library(TCGAbiolinks)
library(dplyr)
#source("https://bioconductor.org/biocLite.R")
#biocLite("RTCGA")
library(RTCGA)
#install.packages("magrittr")
library(magrittr)

pacman::p_load('vtreat','caret','gbm','doParallel','data.table','dplyr','h2o',
               'RTCGAToolbox','TCGAbiolinks', 'RTCGA')

# FROM TXT -------------------------------------------

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

save(rna,file='/Users/1Air/Documents/Insight/project/rna.rda' )
save(marray,file='/Users/1Air/Documents/Insight/project/marray.rda' )
save(clin,file='/Users/1Air/Documents/Insight/project/clin.rda' )

load('/Users/1Air/Documents/Insight/project/rna.rda')
load('/Users/1Air/Documents/Insight/project/marray.rda')
load('/Users/1Air/Documents/Insight/project/clin.rda')
rna.bcr  <-  rna$pat.bcr
marray.bcr  <- marray$pat.bcr
clin.bcr  <-  clin$pat.bcr

# CLEAN----------------------------------------------------

clin <- as.data.frame(t(clin))
rna <- as.data.frame(t(rna))
marray <- as.data.frame(t(marray))

# rna
rna.bcr <- row.names(rna)[-1]
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


# merge rna and marray if possible

# get column for join by deriving patient provenance across both datasets

# extract patient ID/ barcodes

#pat.bcr <- as.character(unname(unlist(clin[,19])))
clin.bcr1 <- as.character(unname(unlist(clin[,3040])))
clin.bcr2 <- as.character(unname(unlist(clin[,1662])))
clin.bcr3 <- as.character(unname(unlist(clin[,3125])))
clin.bcr4 <- as.character(unname(unlist(clin[,2129])))
clin.bcr5 <- as.character(unname(unlist(clin[,2161])))
clin.bcr6 <- as.character(unname(unlist(clin[,3272])))
clin.bcr7 <- as.character(unname(unlist(clin[,2577])))
clin.bcr8 <- as.character(unname(unlist(clin[,2278])))
clin.bcr9 <- as.character(unname(unlist(clin[,3058])))
clin.bcr10 <- as.character(unname(unlist(clin[,1559])))
clin.bcr11 <- as.character(unname(unlist(clin[,2593])))
clin.bcr12 <- as.character(unname(unlist(clin[,2960])))

clin$clin.bcr1 <- toupper(clin[,3040])
clin$clin.bcr2 <- toupper(clin[,1662])
clin$clin.bcr3 <- toupper(clin[,3125])
clin$clin.bcr4 <- toupper(clin[,2129])
clin$clin.bcr5 <- toupper(clin[,2161])
#clin$clin.bcr4 <- toupper(clin[,2129])

bcr.master <- c(clin.bcr1,clin.bcr2,clin.bcr3,clin.bcr4,clin.bcr5,clin.bcr6,
                clin.bcr7,clin.bcr8,clin.bcr9,clin.bcr10,clin.bcr11,clin.bcr12)
bcr.master <- toupper(unique(bcr.master))

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


# FEATURE EXTRACTION -----------------------------------------------
clin$patient.age <- floor(as.numeric(as.character(clin$patient.days_to_birth))/-365.4)
marray$TP53 <- as.numeric(as.character(marray$TP53))
marray$BRCA1 <- as.numeric(as.character(marray$BRCA1))
marray$BRCA2 <- as.numeric(as.character(marray$BRCA2))
#rownames(genedata) <- gene.bcr

# class variable
#genexp$class <- genedata$patient.vital_status
clin$class <- factor( rep("late",nrow(clin)),ordered = F, levels = c('early','late'))
clin$class[clin$patient.stage_event.pathologic_stage=='stage i'
             | clin$patient.stage_event.pathologic_stage=='stage ia'
             | clin$patient.stage_event.pathologic_stage=='stage ib'
             | clin$patient.stage_event.pathologic_stage=='stage ii'
             | clin$patient.stage_event.pathologic_stage=='stage iia'
             | clin$patient.stage_event.pathologic_stage=='stage iib'
             | clin$patient.stage_event.pathologic_stage=='stage iii'
             | clin$patient.stage_event.pathologic_stage=='stage iiia'
             ] <- 'early'

table(clin$patient.stage_event.pathologic_stage)
levels(clin$patient.stage_event.pathologic_stage)
# early/good prognosis i, ia, ib, iia, iib, iiia, iii, ii
# late/concern all else
table(clin$class)/nrow(clin)

#clin <- select(clin, -patient.stage_event.pathologic_stage)

# head(colnames(genedata[,1:17814]))
# head "ELMO2"    "CREB3L1"  "RPS11"    "PNMA1"    "MMP2"     "C10orf90"
# tail "RXFP2"   "PIK3IP1" "SLC39A6" "SNRPD2"  "AQP7"    "CTSC" 
genedata[,1:17814] <- as.numeric(as.character(unlist(genedata[,1:17814])))

# now create joins
#row.names(rna) <- rna.bcr
rna$pat.bcr <- rna.bcr
marray$pat.bcr <- marray.bcr
clin$pat.bcr <- toupper(clin.bcr)

master.bcr <- data.frame(bcr.master, 
                      patient.patient_id=tolower(substr(bcr.master,9,12)),
                      patient.bcr_patient_barcode=tolower(substr(bcr.master,1,12)))

head(clin$patient.patient_id)
head(master.bcr)
length(unique(clin$patient.bcr_patient_barcode))
anyDuplicated(clin$patient.bcr_patient_barcode)
anyDuplicated(substr(rna.bcr,1,12)) #100
anyDuplicated(substr(marray.bcr,1,12)) #76

head(clin[,c(19,3040,3722)])
which(colnames(clin)=='class')
testjoin2 <- merge(clin[,c(19,3040,3722)], master.bcr, by='patient.bcr_patient_barcode', all.x=T, sort=F)
sum(is.na(testjoin2$class))
testjoin2 <- left_join(clin[,c(19,3040,3722)], master.bcr)

#clin$pat.bcr <- toupper(clin$pat.bcr)
#row.names(clin) <- clin$pat.bcr
resp <- clin[,c('clin.bcr1', 'class')]
marray[,c('pat.bcr', 'class')]

length(intersect(bcr.master,rna.bcr)) # 874/878
setdiff(rna.bcr,bcr.master)
length(intersect(bcr.master,marray.bcr)) #584/590
setdiff(marray.bcr,bcr.master)

clin$patient.patient_id[unique(row(clin[-1])[clin[-1] == "TCGA-A7-A0CG-01A-12R-A056-07"])]
df[!!rowSums(df==1),]
clin[,which(clin=="TCGA-A7-A0CG-01A-12R-A056-07",arr.ind = T)]
colnames(clin)[which(!is.na(colSums(clin=='TCGA-A7-A0CG-01A-12R-A056-07')))]

l <- sapply(colnames(clin), function(x) grep('TCGA-A7-A0CG-01A-12R-A056-07', clin[,x]))

length(unique(marray$pat.bcr))
resp$pat.bcr
marray$pat.bcr[9] %in% resp$pat.bcr
mray <- marray[,17810:17816]
ray <- rna[,20530:20534]
rna$class <- NULL

testjoin <- merge(ray, testjoin2, by.x='pat.bcr',by.y ='bcr.master' , all.x=T, sort=F)
#head(ray$pat.bcr)
sum(is.na(testjoin$class))

testjoin3 <- merge(mray, testjoin2, by.x='pat.bcr',by.y ='bcr.master' , all.x=T, sort=F)
sum(is.na(testjoin3$class))

marray <- merge(marray, clin[,c('pat.bcr', 'class')], by='pat.bcr', all.x=T, sort=F)
rna <- merge(rna, clin[,c('pat.bcr', 'class')], by='pat.bcr', all.x=T, sort=F)
summary(marray$class)
summary(rna$class)

genomics <- inner_join(rna, marray)
genomics <- full_join(rna, marray)

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
# microarray model
summary(marray$class)
inTrain <- createDataPartition(y = marray$class, p = .80,list = FALSE)
train <- marray[ inTrain,]
test <- marray[-inTrain,]
nrow(train) # 57564
nrow(test) # 29654
plot(train$class)
summary(marray$class)

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

