#Script for running annotation functions as well as defining sRNA loci for chlamy.
# The appropriate thresholds (e.g. FDR) are determined in:
# Segmentation_Analysis.R

##### Load necesary libraries
library(stringr)
library(xtable)
library(rtracklayer)
library(segmentSeq)

library(MASS)
library(RColorBrewer)
library(baySeq)
library(MKmisc) #binomCI, binomial confidence interval
library(readr)
library(dplyr)
library(purrr)

###### setting variables

fdr <- 0.05

baseDir <- "/projects/nick_matthews"
# Specify location of segmentation
segLocation <- file.path(baseDir, "segmentation_2018")
annoDir     <- file.path(baseDir, "resources")
annoFile    <- "chlamy_all_annotations.Rdata"
#set working directory to github repository on cluster
gitdir  <- file.path(baseDir, "chlamy_locus_map_github")

meta    <- read_csv(file.path(gitdir, "Summary_of_Data.csv")) %>%
  filter(InCurrentLociRun == "Yes")
metawt <- meta$Controls %in% c("wt") # only WT libs

source(file.path(gitdir, "Scripts/chlamy_source_code.R"))

# list.files(segLocation, pattern = "segD.*")
inputdata <- "segD_chlamy_segmentation_LociRun2018_multi200_gap100.RData" #14390
aDfile    <- "aD_chlamy_segmentation_LociRun2018_multi200_gap100.RData"

gitfingerprint <- system2("git", args = "rev-parse --short HEAD", stdout = TRUE)

prefix <- str_replace(inputdata, "segD_chlamy_segmentation_(.*).RData", "\\1")
saveLocation <- file.path(segLocation, paste(prefix, gitfingerprint, sep = "_"))
dir.create(saveLocation)


#####Load in and process segmentation, selecting significant loci#####
# loads segD object (lociData class)
load(file.path(segLocation, inputdata))

# Select loci based on some fdr
# perReplicate: If TRUE, selection of loci is done on a replicate by replicate basis. If FALSE, selection will be done on the likelihood that the locus represents a true locus in at least one replicate group.

# which setting are we using? see issue #9
loci <- selectLoci(cD = segD, FDR = fdr, perReplicate = TRUE) # 6164
attr(loci, "parameter") <- prefix
attr(loci, "git") <- gitfingerprint

## creating/exporting coordinates object
gr <- loci@coordinates
# safe attributes in object
attr(gr, "parameter") <- prefix
attr(gr, "git") <- gitfingerprint
# name the loci 'CR' - chlamydomonas reinhardtii, 'SL' - srna locus
names(gr) <- sprintf("CRSL%05.f0", 1:length(gr))

#####These next few functions compute and compile annotation files, this shouldn't need runnding every time#####
#Compute introns - this takes a while, don't run unless necessary
intronCalculate()
#Process transposon file
annotenv <- new.env()
myfile <- file.path(annoDir, "transposon_annotations.Rdata")
if(!file.exists(myfile)) {
  transposonProcess(annoDir, gitdir)
} else {
  message(myfile, " exists, loading..")
  load(myfile, envir = annotenv, verbose = FALSE)
}

myfile <- file.path(annoDir, annoFile)
if(!file.exists(myfile)) {
  #Compute and compiles annotations
  compileAnnotations(annoDir, annoFile = annoFile)
} else {
  message(myfile, " exists, loading..")
  load(myfile, envir = annotenv, verbose = FALSE)
}
# attach annotenv so objects are accesible
attach(annotenv)
ls(annotenv)

#####the next set of functions take the annotated locus object 'gr' and add some more annotation data#####
#The functions are found 'chlamy_source_code.R'

# annotate by size class
gr <- sizeClass(gr, intervals = c(0,100,400,1500,3000,Inf))
pdf(file.path(saveLocation, paste0("sizeLoci_density", prefix, ".pdf")))
plot(density(log(gr$size[gr$size>0])))
abline(v = log(c(0,100,400,1500,3000,Inf)))
dev.off()

# annotate with overlapping features
gr <- expressionClass(locAnn = gr, loci = loci, wt = metawt)
gr <- featureAnn(locAnn = gr, annotations = annotenv)

#methylation function
gr <- methylation(gr, annoDir)

#Extra Current annotations
gr <- strainSpec(gr, loci, meta = meta, gitdir = gitdir)
gr <- lifeCycle(gr, loci, meta = meta, gitdir = gitdir)
gr <- mutantSpec(gr, loci, meta = meta, gitdir = gitdir)

# PhaseTank
gr <- phaseMatch2(gr,annoDir=file.path(gitdir, "Data/PhaseTank_OUTPUT_2018.11.27_18.08"),outputName="Pred_tab_2018.11.27_18.08")
pdf(file.path(saveLocation, paste0("phasing_density", prefix, ".pdf")))
plot(density((gr$phaseScore[gr$phaseScore>0])))
dev.off()

# annotate with counting biases; i.e, is there a higher than average ratio of 21s to 20s, or a higher number of reads starting with As than usual
cl <- makeCluster(24)
gr <- countingBiases(locAnn = gr, cl = cl,
                     segLocation = segLocation,
                     wt = metawt,
                     aDfile = aDfile)
stopCluster(cl)

pdf(file.path(saveLocation, paste0("21v20ratio_density", prefix, ".pdf")))
plot(density(log(gr$ratio21vs20[gr$ratio21vs20>0]), na.rm = T))
abline(v = log(c(0.3, 0.7)))
dev.off()


# what is the dominating size class?
sizedf <- data.frame( "smaller_20bp" = gr$countsSmall,
                      "equal_20bp"   = gr$counts20,
                      "equal_21bp"   = gr$counts21,
                      "larger_21bp"  = gr$countsBig)

# determine index of each row with the max number of sRNAs
idx <- sizedf %>%
  pmap_int(function(...) which.max(c(...)))
gr$predominant_sRNA_sizeClass <- colnames(sizedf)[idx]
# breakdown of how many loci have a specif prevailing sRNA mapping to it:

gr$ratio_strand_class <- classCI(gr$countsplus, gr$countsminus,
                                 probs = c(0.2,0.4,0.6,0.8), comma = 1,
                                 plotname=file.path(saveLocation, paste0("standbias_", prefix, ".pdf")))
levels(gr$ratio_strand_class) <- c("strong_bias", "med_bias", "no_bias", "med_bias", "strong_bias")

# countsallwt: WT reads
# countsnormalwtnorm: WT normalized reads (corrected for multi read count)
# the ratio between the constitutes repetitiveness
gr$repetitivenessClass <- classCI(gr$countsallwt, gr$countsnormalwtnorm,
                                  probs = c(med = 0.6, high = 0.9), comma = 1,
                                  plotname=file.path(saveLocation, paste0("Repetitiveness_", prefix, ".pdf")))

#Phasing
gr$phaseClass <- as.ordered(cut(gr$phaseScore, c(-1, 0, 60, Inf),include.lowest=TRUE,
                                labels= c("none","median","high")))


save(gr, meta, metawt, loci, baseDir, prefix, saveLocation, file = file.path(saveLocation, paste0("gr_fdr", fdr,  ".RData")))
load(file = file.path(saveLocation, paste0("gr_fdr", fdr,  ".RData")))
# export as gff3 file for viewing in browser
export.gff3(gr, con = file.path(saveLocation, paste0("loci_fdr", fdr, prefix, ".gff")))
#Write csv for phasing

write.csv(as.data.frame(gr), file = file.path(saveLocation, paste0("loci_fdr", fdr, prefix, ".csv")))
