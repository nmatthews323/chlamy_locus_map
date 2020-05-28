
library(segmentSeq)
library(readr)
library(dplyr)
library(stringr)
baseDir      <- "/projects/nick_matthews"
segLocation  <- file.path(baseDir, "segmentation_2018")
gitdir       <- file.path(baseDir, "chlamy_locus_map_github")
setwd(file.path(gitdir, "PhaseTank_v1.0"))

# selected libraries as discussed with Nick/Seb/Adrian (july2018)
files <- read_csv(file.path(gitdir, "Summary_of_Data.csv")) %>%
    filter(LociRun2018 == "Yes") %>%
    filter(Controls == "wt")

files$DataCode


libDir <- file.path(baseDir, "sequencing_data/fasta_with_counts")
mylibs <- file.path(libDir, paste0(files$DataCode,"_assembly5_Chlamydomonas_reinhardtii.patman.aligned_reads.fasta"))
head(mylibs)
file.exists(mylibs)

tankcmd <- "perl PhaseTank_v1.0_mod.pl --genome /projects/nick_matthews/resources/Creinhardtii_236.fa --lib"

system(paste(tankcmd, paste(mylibs, collapse = ",")))

#You can check the output files in directory './OUTPUT_2018.11.27_18.08/'.
