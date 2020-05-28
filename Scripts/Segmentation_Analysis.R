#Analyses quality of segmentation producing a variety of graphs

library(segmentSeq)
library(stringr)
#Specify location of the segmentation
baseDir <- "/projects/nick_matthews"
# Specify location of segmentation
segLocation <- file.path(baseDir, "segmentation_2018")
annoDir     <- file.path(baseDir, "resources")
#set working directory to github repository on cluster
gitdir      <- file.path(baseDir, "chlamy_locus_map_github")
annoDir <- "/projects/nick_matthews/resources"
inputdata <- "segD_chlamy_segmentation_LociRun2018_multi200_gap100.RData" #14390
gitfingerprint <- system2("git", args = "rev-parse --short HEAD", stdout = TRUE)
prefix <- str_replace(inputdata, "segD_chlamy_segmentation_(.*).RData", "\\1")
saveLocation <- file.path(segLocation, paste(prefix, gitfingerprint, sep = "_"))
dir.create(saveLocation)

# code for locus summary plots mostly to determine fdr cut-off
# the locus map (not having selected loci yet, but after calculating loci likelihoods) is 'segD'
load(file.path(segLocation, inputdata))

#Choose FDR
myfdr <- 0.05

# check parameter used to generate gr
attr(gr, "parameter")
# [1] "_multi200_gap100"
attr(gr, "git")
# [1] "59a4cbc"

# nct<-segD
                                        
master <- read.csv(file.path(gitdir, "Summary_of_Data.csv"))

#To select very specific control WTs
wt<-subset(master,master[,"Controls"] =="wt")
wtReps<-unique(wt[,"Replicate"]) 

nfdrnum <- NULL

# check how many loci get selected at various FDR levels
for(ii in (1:20 / 4)) {
  fdr <- 10^-ii
  loci <- selectLoci(segD, FDR = fdr, perReplicate = TRUE) #
  nfdrnum <- cbind(nfdrnum, c(fdr, nrow(loci)))
}

pdf(file.path(saveLocation, "fdrNumbers.pdf"))
plot(nfdrnum[1,], nfdrnum[2,], log = "xy", xlab = "FDR", ylab = "# of loci", col = rep(c("red", "blue"), each =20))
abline(v=0.05)
dev.off()

# check how many loci appear in N replicate groups (at various fdr levels)
pdf(file.path(saveLocation, "fdr_hists.pdf"), width = 10, paper = "a4r")
par(mfcol = c(3,5))
for(fdr in c(0.1, 0.05, 10^-(2:5))) {
  loci <- selectLoci(segD, FDR = fdr, perReplicate = TRUE)
  hist(rowSums(exp(loci@locLikelihoods[,wtReps])), breaks = 0:length(wtReps), main = paste("FDR =", fdr), xlab = "Expectation")
  hist(rowSums(exp(loci@locLikelihoods)), breaks = 0:nlevels(segD@replicates), main = paste("FDR =", fdr), xlab = "Expectation")
  # also have a look at locus lenght distribution while we're at it.
  plot(density(log10(width(loci@coordinates))), main = "", xlab = "log width")
}
dev.off()

#########Select loci
loci <- selectLoci(segD, FDR = myfdr, perReplicate = TRUE)


# get an idea of the sequencing depth added by each replicate groups
sumLibsizes <- sapply(levels(loci@replicates), function(rep) sum(libsizes(loci)[loci@replicates == rep]))
meanLibscale <- sapply(levels(loci@replicates), function(rep) mean(libsizes(loci)[loci@replicates == rep]))

# order by increasing sequencing depth added
ordLoc <- order(sumLibsizes, decreasing = FALSE)

# number of additional loci added by each library (with increasing sequencing depth)
cumloc <- sapply(1:length(ordLoc), function(ii) {
  message(ii)
  selLoc <- segD[,segD@replicates %in% ordLoc[1:ii]]
  selLoc@locLikelihoods <- as.matrix(selLoc@locLikelihoods[,ordLoc[1:ii],drop = FALSE])
  z <- try(selectLoci(selLoc, FDR = 0.05 , perReplicate = TRUE))
    if(class(z) == "try-error") return(0) else return(nrow(z))
})

collibs <- rep("black", nlevels(loci@replicates))
collibs[wtReps] <- "red"

# plot number of loci discovered as we add deeper libraries. Ideally, we want to see the number of additional loci tailing off, indicating we've achieved enough sequencing depth/variety to get most loci
pdf(file.path(saveLocation, "CumSeqVolume.pdf"))
plot(x = cumsum(sumLibsizes[ordLoc]), y = cumloc, xlab = "Cumulative sequencing volume", ylab = "Total loci discovered", col = collibs[ordLoc], pch = 19)
dev.off()

# individual numbers of loci per replicate group. Should correlate roughly with library size for lower library sizes, hopefully become more or less constant for higher library sizes as returns diminish.
pdf(file.path(saveLocation, "LibScalingFactor.pdf"))
plot(y = summariseLoci(loci, perReplicate = TRUE), x = meanLibscale, log = "x", ylab = "# of loci", xlab = "Library scaling factor")
dev.off()
