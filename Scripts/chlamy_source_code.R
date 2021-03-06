#Source code file containing all functions applied in other scripts

#####Function set 1 - utility functions#####
#These functions compile annotations or are used in later functions

#Compute annotations loads in annotation files from a provided directory and after processing them outputs and saves them as a Rdata file for further use
#path to annotation files has a default but can also be user specified
compileAnnotations <- function(annoDir = "resources", annoFile) {
  #annoDir = "C:/Users/Nick/Documents/Uni Work/Third Year/Project/segmentMap_II/Old_Annotation_Files"
  #Specify what libraries are required
  require(GenomicRanges)
  require(rtracklayer)

  #Load large transposon file
  load(file.path(annoDir,"transposon_annotations.Rdata"))
  #Load in irs and trs data
  irs <- import.gff3(file.path(annoDir,"irf_output_Creinhardtii_236.gff"))
  irs$type<-irs$source
  trs <- import.gff3(file.path(annoDir,"trf_output_Creinhardtii_236.gff"))
  trs$type<-trs$source

  #Load in miRNA data
  miRNA <- import.gff3(file.path(annoDir,"Vallietal2016_miRNAs.gff3"))
  miRNA$type <- "miRNA"
  #Load in rRNA data
  rRNA <- import.gff3(file.path(annoDir,"Creinhardtii_rRNA.gff3"))
  #Load in MSAT data
  MSAT <- import.gff3(file.path(annoDir,"Creinhardtii_MSAT.gff3"))

  #Load in annotation data from v5.5 chlamy annotation
  anno <- import.gff3(file.path(annoDir,"Creinhardtii_281_v5.5.gene_exons.gff3"))
  #Extract from genes data the relevant subsets
  mRNA          <- anno[anno$type=="mRNA",]
  genes         <- anno[anno$type=="gene",]
  exons         <- anno[anno$type=="exon",]
  CDS           <- anno[anno$type=="CDS",]
  threeprimeUTR <- anno[anno$type=="three_prime_UTR",]
  fiveprimeUTR  <- anno[anno$type=="five_prime_UTR",]

  #Calculate promoter region the 500bp flanking region around gene
  promoter <- promoters(anno[anno$type=="mRNA"],upstream=500,downstream=0)
  promoter$type <- droplevels(promoter$type)
  levels(promoter$type) <- "promoter"

  #Load in intron file, calculated seperately using intronCalculate file
  introns <- import.gff3(file.path(annoDir,"Creinhardtii_281_v5.5.CalculatedIntrons.gff3"))
  levels(introns$type) <- "intron"

  #Save all as one big RData file
  save(genes,anno,
       transposons, repetativeSeq, TE_Class_DNA, TE_Class_RET,
       TE_L1,TE_Gypsy,TE_Copia,TE_hAT,TE_RTE,TE_Novosib,TE_DualenRandI,
       TE_P,TE_Mariner,TE_REM1,TE_EnSpm,TE_DIRS,TE_TOC2,TE_TOC1,TE_Gulliver,TE_TCR1,TE_Harbinger,
       rRNA,miRNA,MSAT,mRNA,irs,trs,promoter,threeprimeUTR,fiveprimeUTR,exons,CDS,
       introns,file=file.path(annoDir, annoFile))
}

#Function which shortens the standard annotation process
#Takes locAnn object and list of 1 or more annotation files, computes overlap and assigns either TRUE/FALSE or named annotation
computeOverlaps <- function(locAnn,annotations,assignNames = FALSE,  namesColumn = "type") {
  if(!(is.environment(annotations) | is.list(annotations)) | is.null(names(annotations))) {
    stop("Requires list or environment input with name of annotation")
  }
  for(annotName in names(annotations)) {
    annotation <- annotations[[annotName]]
    #Ensure sequence levels are equivalent - stops warnings for scaffolds
    seqlevels(annotation) <- seqlevels(locAnn)
    #Calculate overlap
    overlap <- findOverlaps(locAnn,annotation)
    #extract only unique overlaps
    overlapunique <- !duplicated(queryHits(overlap))
    
    #Assign TRUE for those that overlap unless want to specify names of overlap
    if(assignNames == FALSE) {
      #Set up column
      locAnn$temp <- rep(FALSE, length(locAnn))
      locAnn$temp[queryHits(overlap)[overlapunique]] <- TRUE
    } else if(assignNames == TRUE) {
      #If assigning names takes the
      #Set up column
      locAnn$temp <- rep("none", length(locAnn))
      locAnn$temp[queryHits(overlap)[overlapunique]] <- as.character(mcols(annotation)[,namesColumn])[subjectHits(overlap)[overlapunique]]
    } else {
      warning("assignNames must be either TRUE or FALSE")
    }
    #Assign temp column to correct name and remove temporary column
    mcols(locAnn)[annotName] <- locAnn$temp
    locAnn$temp <- NULL
  }
  #Return locAnn
  locAnn
}

#Function calculates introns from the gene annotation v5.5 from p
intronCalculate <- function(annoDir = "resources") {
  #Load in necessary libraries
  require(tidyverse)
  require(GenomicRanges)
  require(rtracklayer)

  #####Load in dataset#####
  allGenes <- import.gff3(file.path(annoDir,"Creinhardtii_281_v5.5.gene_exons.gff3"))
  widths <- data.frame(type = c(), width = c())
  widths <- data.frame(type=c(as.character(unique(allGenes$type))),
                       width=sapply(unique(allGenes$type),function(x) mean(width(allGenes[allGenes$type==x]))))

  #####Extract exonic sequences#####
  genes <- allGenes[allGenes$type=="gene",]
  exons <- allGenes[allGenes$type=="exon",]
  mRNA  <- allGenes[allGenes$type=="mRNA",]

  #####Derive intronic sequences#####
  #remove mRNA which don't have at least two exons overlapping them - therefore wouldn't have an intron by definition
  exons <- unique(exons)
  exonOverlaps <- data.frame(table(queryHits(findOverlaps(mRNA,exons))))
  exonOverlaps$Keep <- rep(FALSE,nrow(exonOverlaps))
  exonOverlaps$Keep[exonOverlaps$Freq > 1] <- TRUE
  mRNAList <- split(mRNA[exonOverlaps$Keep],1:length(mRNA[exonOverlaps$Keep]))
  #Use endoapply to derive introns from areas of mRNA not covered by exons
  intronList <- endoapply(mRNAList, function(x){
    introns <- NULL
    introns <- GenomicRanges::setdiff(x,exons)
    if(length(introns) > 0) {
      #If intron found, add in relevant information
      introns$source <- paste0("CalculatedFrom:",x$source)
      introns$parent <- x$ID
      introns$type   <- "intron"
      introns$id     <- paste0(x$ID,".intron.",1:length(introns))
      introns$pacid  <- x$pacid
    } else {
      #If no introns, just return empty Granges object
      introns = GRanges()
    }
    introns
  })

  #Unlist into one GRanges object
  introns <- unlist(intronList)

  #Save introns as gff3 file
  export.gff3(introns,file.path(annoDir,"Creinhardtii_281_v5.5.CalculatedIntrons.gff3"))
}

#Function classified transposon dataset into standardised classification scheme
#Have to specify whether to keep unknowns, default FALSE
transposonProcess <- function(annoDir = "resources",gitdir,keepUnknowns=FALSE) {
  #Require some libraries
  require(rtracklayer)
  require(GenomicRanges)
  #Read in repeatMasked file and reference file
  transposons <- import.gff3(file.path(annoDir, "Creinhardtii_281_v5.5.repeatmasked_assembly_v5.0.gff3"))
  referenceFile <- read.csv(file.path(gitdir,"transposon_classification_scheme.csv"),header=TRUE, stringsAsFactors = FALSE)
  #Set up empty columns
  transposons$name        <- rep("Unknown", length(transposons))
  transposons$class       <- rep("Unknown", length(transposons))
  transposons$order       <- rep("Unknown", length(transposons))
  transposons$superfamily <- rep("Unknown", length(transposons))
  #Run loop through each row and assign order
  for(ii in 1:nrow(referenceFile)) {
    #Find matches for the search terms
    tempSearch <- grep(referenceFile$SearchTerm[ii],transposons$Name, ignore.case = TRUE)
    #Assign corresponding values
    transposons$name[tempSearch]        <- referenceFile$Name[ii]
    transposons$class[tempSearch]       <- referenceFile$Class[ii]
    transposons$order[tempSearch]       <- referenceFile$Order[ii]
    transposons$superfamily[tempSearch] <- referenceFile$Superfamily[ii]
  }
  #Extract all sequences as repetative sequences
  repetativeSeq <- transposons
  #Extract just identified transposons
  transposons <- transposons[transposons$name != "Unknown"]
  #Extract the two major classes
  TE_Class_DNA <- transposons[transposons$class =="DNA"]
  TE_Class_RET <- transposons[transposons$class =="RET"]
  #Export files
  export.gff3(transposons,file.path(annoDir,"Creinhardtii_281_v5.5.transposonsfromrepeatmasked_assembly_v5.0.gff3"))
  export.gff3(repetativeSeq, file.path(annoDir,"Creinhardtii_281_v5.5.repeatmaskedannotated_assembly_v5.0.gff3"))
  export.gff3(TE_Class_DNA, file.path(annoDir,"Creinhardtii_281_v5.5.DNAtransposonsfromrepeatmasked_assembly_v5.0.gff3"))
  export.gff3(TE_Class_RET, file.path(annoDir,"Creinhardtii_281_v5.5.RETtransposonsfromrepeatmasked_assembly_v5.0.gff3"))

  #extract major Orders
  orders <- unique(transposons$order)
  for(ii in 1:length(orders)) {
    #Extract superfamily
    temp <- transposons[transposons$order == orders[ii]]
    #Export as gff3
    export.gff3(temp,file.path(annoDir,paste0("Creinhardtii_281_v5.5.Order",orders[ii],"transposonsfromrepeatmasked_assembly_v5.0.gff3")))
    #Assign to object in environment for saving
    assign(paste0("TE_Order_",orders[ii]),temp)
  }

  #Now do the same for all the superfamilies
  superfamilies <- unique(transposons$superfamily)
  for(ii in 1:length(superfamilies)) {
    #Extract superfamily
    temp <- transposons[transposons$superfamily == superfamilies[ii]]
    #Export as gff3
    export.gff3(temp,file.path(annoDir,paste0("Creinhardtii_281_v5.5.",superfamilies[ii],"transposonsfromrepeatmasked_assembly_v5.0.gff3")))
    #Assign to object in environment for saving
    assign(paste0("TE_",superfamilies[ii]),temp)
  }

  #Export as consolidated Rdata file
  save(transposons, repetativeSeq, TE_Class_DNA, TE_Class_RET,
       TE_Order_LTR,TE_Order_LINE,TE_Order_TIR,TE_Order_SINE,TE_Order_DIRS,
       TE_L1,TE_Gypsy,TE_Copia,TE_hAT,TE_RTE,TE_Novosib,TE_DualenRandI,
       TE_P,TE_Mariner,TE_REM1,TE_EnSpm,TE_DIRS,TE_TOC2,TE_TOC1,TE_Gulliver,TE_TCR1,TE_Harbinger,
  file=file.path(annoDir,"transposon_annotations.Rdata"))
}

#Function used by classCI function
assignCI <- function(cis, windows, comma) {
  intwin <- cbind(findInterval(cis[1,], windows), findInterval(cis[2,], windows))
  z <- intwin[,2] - intwin[,1]
  asswin <- rep(NA, nrow(intwin))
  if(any(z == 0, na.rm = TRUE))
    asswin[which(z == 0)] <- intwin[which(z == 0), 1]
  if(any(z == 1, na.rm = TRUE))
    asswin[which(z == 1)] <- intwin[cbind(which(z == 1), apply(abs(intwin[which(z == 1),] - (comma + 1)), 1, which.min))]
  classIDs <- names(windows)[asswin]

  classIDs <- ordered(as.factor(classIDs), levels = names(windows)[-length(windows)])
}

#Needed to compute confidence intervals in countingBiases function
classCI <- function(x1, x2, probs, divisions, comma, plot = TRUE,plotname) {
  require(MKmisc) # binomCI function
  if(missing(probs))
    probs <- 1 / (1 + 1/ divisions)
  windows <- c(low = -Inf, probs, infinite = Inf)

  if(is.null(names(probs))) names(windows) <- apply(cbind(c(-Inf, round(probs,3)), round(c(probs, Inf), 3)), 1, paste, collapse = "-")

  ratioCI <- apply(cbind(x1, x2), 1, function(x) binomCI(x[1], sum(x), method = "modified jeffreys")$conf.int)

  classIDs <- assignCI(ratioCI, windows, comma)
  pdf(plotname)
  plot(density(log10(colMeans(ratioCI[,!is.na(classIDs)])), na.rm = TRUE)); abline(v = log10(probs))

  plot(density(colMeans(ratioCI[,!is.na(classIDs)]), na.rm = TRUE, from = 0, to = 1)); abline(v = (probs))
  p2r <- function(p) 1 / ((1 / p) - 1)
  plot(density(log2(p2r(colMeans(ratioCI[,!is.na(classIDs)]))))); abline(v = log2(p2r(probs)))
  plot(density(log10((ratioCI[1,!is.na(classIDs)])), na.rm = TRUE)); abline(v = log10(probs))
  plot(density((ratioCI[1,!is.na(classIDs)]), na.rm = TRUE, from = 0, to = 1)); abline(v = (probs))
  plot(density(log10((ratioCI[2,!is.na(classIDs)])), na.rm = TRUE)); abline(v = log10(probs))
  plot(density((ratioCI[2,!is.na(classIDs)]), na.rm = TRUE, from = 0, to = 1)); abline(v = (probs))

  dev.off()

  classIDs
}

#####Scripts for generating summary tables#####
#This function generates loci identified per chromosome and loci density table
chromosomeLociTable <- function(locusMap, mapLocation = "segmentation2018") {
  #libraries requires
  require(segmentSeq)
  require(rtracklayer)
  #Load in loci
  loci <- locusMap
  #Specify Chlamy assembly v5 characteristics
  summaryTable <-
    data.frame(Chromosome=c("chromosome_1","chromosome_2","chromosome_3","chromosome_4","chromosome_5","chromosome_6","chromosome_7","chromosome_8","chromosome_9","chromosome_10",
                            "chromosome_11","chromosome_12","chromosome_13","chromosome_14","chromosome_15","chromosome_16","chromosome_17","scaffold_18","scaffold_19","scaffold_20",
                            "scaffold_21","scaffold_22","scaffold_23","scaffold_24","scaffold_25","scaffold_26","scaffold_27","scaffold_28","scaffold_29","scaffold_30",
                            "scaffold_31","scaffold_32","scaffold_33","scaffold_34","scaffold_35","scaffold_36","scaffold_37","scaffold_38","scaffold_39","scaffold_40",
                            "scaffold_41","scaffold_42","scaffold_43","scaffold_44","scaffold_45","scaffold_46","scaffold_47","scaffold_48","scaffold_49","scaffold_50",
                            "scaffold_51","scaffold_52","scaffold_53","scaffold_54"),
               Length=c(8033585,9223677,9219486,4091191,3500558,9023763,6421821,5033832,7956127,6576019,
                        3826814,9730733,5206065,4157777,1922860,7783580,7188315,271631,219038,200793,
                        189560,163774,127913,127161,102191,80213,55320,55278,52813,52376,
                        48183,42264,39192,33576,32450,25399,24537,24437,22408,22082,
                        21325,21000,20974,17736,16939,16627,14746,14165,13462,
                        12727,11225,6241,2479,2277))
  #Now calculate loci
  summaryTable$LociIdentified <- apply(summaryTable,1,function(x) sum(seqnames(loci) == x[1]))
  summaryTable$LociDensity <- summaryTable$LociIdentified/summaryTable$Length*1000000
  summaryTable
}

#####Function set 2 - locus feature annotations#####
#This set of functions annotate loci based on their intrinsic features or occurances

#sizeClasses divides the loci into intervals defined through visual analysis of sizes
#intervals which have a default in the function but can also be user specified when the function is called
sizeClass <- function(locAnn, intervals = c(0,30,75,150,1000,Inf)) {
  locAnn$size <- width(locAnn)
  locAnn$sizeclass <- as.ordered(cut(width(locAnn),intervals))
  locAnn
}

#Expression class calculates whether it occurs accross all WTs or not
expressionClass <- function(locAnn, loci, wt) {
  #Require packages
  require(GenomicRanges)
  require(rtracklayer)

  #expression only for wt data
  wtrepgroups <- as.integer(unique(loci@replicates[wt]))
  locAnn$expression <- rowSums(loci@locLikelihoods[,wtrepgroups]>log(0.5))
  locAnn$expressionClass <- as.ordered(cut(locAnn$expression,breaks=c(-Inf,1,5,Inf),include.lowest=TRUE,labels=c("specific","inbetween","common"))) #changed common threshold to more than 5
  locAnn
}

# counting biases - uses the alignment data object used to run the segmentation
countingBiases <- function(locAnn, cl, samplefile, segLocation = "segmentation_2018", aDfile,wt) {
  #Firstly load in raw alignment data
  load(file.path(segLocation, aDfile))
  colnames(values(aD@alignments))[2] <- "multireads"

  #Extract raw widths from the aD object for key wt libraries
  aDwidths <- width(aD@alignments)
  aDnormal <- aD[aDwidths>19 & aDwidths < 22,]
  #Remove any multi-reads
  aDnormalnomulti <- aDnormal[!duplicated(as.character(aDnormal@alignments$tag)),]
  aDwidths <- width(aDnormal@alignments)

  #Looking at first nucleotide biases
  #Find first nucleotide
  firstNuc <- substr(aDnormal@alignments$tag,1,1)
  firstNucnomulti <- substr(aDnormalnomulti@alignments$tag,1,1)
  #proportion of sRNAs with a given 5'nuc
  # find normal ratio of first base nucleotides
  expectedRatio <-  tapply(rowSums(aDnormalnomulti@data),firstNucnomulti,sum)/sum(aDnormalnomulti@data)
  # get counts in each locus
  aDA <- aDnormal[firstNuc=="A",]
  countsA <- getCounts(segments=locAnn,aD=aDA,cl=cl); rm(aDA)
  locAnn$countsA <- rowSums(countsA[,wt])
  aDT <- aDnormal[firstNuc=="T",]
  countsT <- getCounts(segments=locAnn,aD=aDT,cl=cl); rm(aDT)
  locAnn$countsT <- rowSums(countsT[,wt])
  aDC <- aDnormal[firstNuc=="C",]
  countsC <- getCounts(segments=locAnn,aD=aDC,cl=cl); rm(aDC)
  locAnn$countsC <- rowSums(countsC[,wt])
  aDG <- aDnormal[firstNuc=="G",]
  countsG <- getCounts(segments=locAnn,aD=aDG,cl=cl); rm(aDG)
  locAnn$countsG <- rowSums(countsG[,wt])
  #Combine results
  firstBase <- cbind(locAnn$countsA, locAnn$countsC, locAnn$countsG, locAnn$countsT)
  # assumes binomial distribution and looks for significant variation from the overall ratio for each locus
  z <- pbinom(firstBase,
              matrix(rowSums(firstBase), ncol = ncol(firstBase), nrow = nrow(firstBase)),
              prob = matrix(expectedRatio, ncol = ncol(firstBase), nrow = nrow(firstBase), byrow = TRUE), lower.tail = FALSE)
  z[rowSums(firstBase) == 0,] <- NA
  zq <- matrix(p.adjust(z, method = "BH"), ncol = ncol(firstBase))
  zq[rowSums(firstBase) == 0,] <- 1
  pred <- apply(do.call("cbind", lapply(1:4, function(ii) c("", c("A", "C", "G", "T")[ii])[as.integer(zq[,ii] < 0.01) + 1])), 1, paste, collapse = "")
  pred[pred == ""] <- NA
  locAnn$predominant_5prime_letter <- as.factor(pred)

  #Looking at biases in sides of sRNAs
  # counting numbers of 20s and 21s
  aDs <- aD[width(aD@alignments)<20,]
  countsSmall <- getCounts(segments=locAnn,aD=aDs,cl=cl); rm(aDs)
  countsSmallwt <- countsSmall[,wt]
  aD20 <- aD[width(aD@alignments)==20,]
  counts20 <- getCounts(segments=locAnn,aD=aD20,cl=cl); rm(aD20)
  counts20wt <- counts20[,wt]
  aD21 <- aD[width(aD@alignments)==21]
  counts21 <- getCounts(segments=locAnn,aD=aD21,cl=cl); rm(aD21)
  counts21wt <- counts21[,wt]
  aDnorm <- aD[width(aD@alignments)<22 & width(aD@alignments)>19]
  countsnorm <- getCounts(segments=locAnn,aD=aDnorm,cl=cl); rm(aDnorm)
  countsnormwt <- countsnorm[,wt]
  aDb <- aD[width(aD@alignments)>21,]
  countsBig <- getCounts(segments=locAnn,aD=aDb,cl=cl); rm(aDb)
  countsBigwt <- countsBig[,wt]

  #sRNA size ratio
  #this part computes confidence intervals on the ratio of 20s and 21s for each locus, and then uses that to put each locus in a window.
  #This part will produce a plot of the density of the mean (and the log of the mean) of the confidence intervals; choose the 'probs' values in such a way to split the modes of the density plots.

  #Assign counts to locAnn
  locAnn$counts20     <- rowSums(counts20wt)
  locAnn$counts21     <- rowSums(counts21wt)
  locAnn$countsSmall  <- rowSums(countsSmallwt)
  locAnn$countsBig    <- rowSums(countsBigwt)
  locAnn$countsNormal <- rowSums(countsnormwt)

  #For 21 vs 20 ratio
  locAnn$ratio21vs20 <- log2((rowSums(counts21wt))/(rowSums(counts20wt)))
  locAnn$ratio21vs20[(locAnn$counts20+locAnn$counts21)<=5] <- NaN

  #For 20 vs 21 ratio
  locAnn$ratio20vs21 <- log2((rowSums(counts20wt))/(rowSums(counts21wt)))
  locAnn$ratio20vs21[(locAnn$counts20+locAnn$counts21)<=5] <- NaN

  #For small RNAs vs Normal RNAs
  locAnn$ratioSmallvsNormal <- log2((rowSums(countsSmallwt))/(rowSums(countsnormwt)))
  locAnn$ratioSmallvsNormal[(locAnn$countsSmall+locAnn$countsNormal)<=5] <- NaN

  #For big RNAs vs Normal RNAs
  locAnn$ratioBigvsNormal <- log2((rowSums(countsBigwt))/(rowSums(countsnormwt)))
  locAnn$ratioBigvsNormal[(locAnn$countsBig+locAnn$countsNormal)<=5] <- NaN

  #Calculate strand biases - phased loci should be strongly biased
  # as above, but now looking at strand biases (i.e., whether most of the sRNAs are on the same strand, or evenly split across strands)
  countsminus   <- getCounts(segments=locAnn,aD=aDnormal[strand(aDnormal@alignments)=="-",],cl=cl)
  countsplus    <- getCounts(segments=locAnn,aD=aDnormal[strand(aDnormal@alignments)=="+",],cl=cl)
  countsall     <- getCounts(segments=locAnn,aD=aDnormal, cl=cl)
  countsallwt   <- countsall[,wt]
  countsminuswt <- countsminus[,wt]
  countspluswt  <- countsplus[,wt]

  locAnn$countsall    <- rowSums(countsall)
  #strand bias
  locAnn$countsallwt  <- rowSums(countsallwt)
  locAnn$countsminus  <- rowSums(countsminuswt)
  locAnn$countsplus   <- rowSums(countspluswt)
  locAnn$ratio_strand <- log2((rowSums(countsminuswt))/(rowSums(countspluswt)))
  locAnn$ratio_strand[(locAnn$countsminus+locAnn$countsplus)<=5] <- NaN
  locAnn$ratio_strand_abs <- abs(locAnn$ratio_strand)


  ##computing repetetivness for each loci (total reads div by multi match corrected) - use full setof aD not aDnormal
  matches               <- aDnormal@alignments$multireads
  aDnormal@data         <- aDnormal@data/matches
  countsnormalwtnorm    <- getCounts(segments=locAnn,aD=aDnormal[,wt],cl=cl)
  locAnn$countsnormalwtnorm <- rowSums(countsnormalwtnorm)
  locAnn$repetitiveness <- 1-(rowSums(countsnormalwtnorm)/rowSums(countsallwt))

  locAnn
}

##Compares loci from different life cycles
lifeCycle <- function(locAnn, loci, FDR=0.05, meta, gitdir) {
  source(file.path(gitdir, "Scripts/selectLoci.R"))
  selLoc <- selectLoci(loci, FDR = FDR, perReplicate = TRUE, returnBool = TRUE)
  locAnn$vegetative <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(meta$LifeCycle == "Vegetative")]))]) > 0
  locAnn$zygote <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(meta$LifeCycle == "Zygote")]))]) > 0
  locAnn$vegetativespecific <- locAnn$vegetative & !locAnn$zygote
  locAnn$zygotespecific <- locAnn$zygote & !locAnn$vegetative
  locAnn
}

##Compares loci from different strains
strainSpec <- function(locAnn, loci, FDR=0.05, meta, gitdir) {
  source(file.path(gitdir, "Scripts/selectLoci.R"))
  selLoc <- selectLoci(loci, FDR = FDR, perReplicate = TRUE, returnBool = TRUE)
  #Work out which strains the loci are present in
  locAnn$CC1883 <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(meta$Genotype == "CC1883")]))]) > 0
  locAnn$CC125 <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(meta$Genotype == "CC125")])),drop=FALSE]) > 0
  locAnn$CC4350 <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(meta$Genotype == "CC4350")])),drop=FALSE]) > 0
  locAnn$J <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(meta$Genotype == "J")])),drop=FALSE]) > 0

  #Assign any that are specific
  locAnn$CC1883specific <- locAnn$CC1883 & !(locAnn$CC125 | locAnn$J | locAnn$CC4350)
  locAnn$CC4350specific <- locAnn$CC4350 & !(locAnn$CC1883 | locAnn$CC125 | locAnn$J)
  locAnn$CC125specific <- locAnn$CC125 & !(locAnn$CC1883 | locAnn$J | locAnn$CC4350)
  locAnn$Jspecific <- locAnn$J & !(locAnn$CC1883 | locAnn$CC125 | locAnn$CC4350)
  locAnn
}

##Compares loci from mutant experiments - NOTE: selects only the wts from Adrian's mutant experiments for comparison!
mutantSpec <- function(locAnn, loci,FDR=0.05, meta, gitdir) {
  #Uses development version of segmentseq code
  source(file.path(gitdir, "Scripts/selectLoci.R"))
  selLoc <- selectLoci(loci, FDR = FDR, perReplicate = TRUE, returnBool = TRUE)

  locAnn$wtAdrian <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(meta$AdrianExp == "wt")]))]) > 0
  locAnn$dcl3Adrian <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(meta$AdrianExp == "dcl3")]))]) > 0
  locAnn$ago3Adrian <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(meta$AdrianExp == "ago3")]))]) > 0
  locAnn$mutantAdrian <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(!meta$AdrianExp == "wt" & !meta$AdrianExp == "FALSE")]))]) > 0

  #Specific to wild type from Adrian's experiments
  locAnn$wtAdrianspecific <- locAnn$wtAdrian & !locAnn$mutantAdrian
  #Specific to Adrian's mutants
  locAnn$mutantAdrianspecific <- locAnn$mutantAdrian & !locAnn$wtAdrian

  #Specific to dcl3 related mutants
  locAnn$notDCL3dependent <- locAnn$dcl3Adrian & !locAnn$wtAdrian
  #Specifically not present in dcl3 mutants
  locAnn$DCL3dependent <- locAnn$wtAdrian & !locAnn$dcl3Adrian
  #Specific to AGO3 related mutants
  locAnn$AGO3dependent <- locAnn$wtAdrian & !locAnn$ago3Adrian
  locAnn$notAGO3dependent <- !locAnn$wtAdrian & locAnn$ago3Adrian
  locAnn
}

#####Function set 3 - annotation overlap functions#####
#This set of function computes overlaps with given annotation sets

#This function loads in the annotation data and computes overlaps with given locAnn object
#There is a default location for the annotation files but can also be user specified
featureAnn <- function(locAnn, annotations) {
  #Require packages
  require(rtracklayer)
  require(GenomicRanges)
  attach(annotations)

  #Compute all simple annotations
  #Then compute annotations
  locAnn <- computeOverlaps(locAnn,annotations)

  #Also compute intergenic regions
  geneoverlap <- findOverlaps(locAnn,genes)
  geneoverlapunique <- unique(queryHits(geneoverlap))
  locAnn$intergenic <- rep(FALSE,length(locAnn))
  locAnn$intergenic[-geneoverlapunique] <- TRUE

  #Assign miRNAs information as seperate columns
  seqlevels(miRNA) <- seqlevels(locAnn)
  miRNAoverlap <- findOverlaps(locAnn,miRNA)
  #miRNA overlap should be unique anyway, but just in case, take first overlap
  miRNAoverlapunique <- !rev(duplicated(rev(queryHits(miRNAoverlap))))
  #Assign geneID
  locAnn$miRNAGeneID <- rep("none", length(locAnn))
  locAnn$miRNAGeneID[queryHits(miRNAoverlap)[miRNAoverlapunique]] <- as.character(miRNA$GeneID[subjectHits(miRNAoverlap)[miRNAoverlapunique]])
  #AssignLocation
  locAnn$miRNALocation <- rep("none", length(locAnn))
  locAnn$miRNALocation[queryHits(miRNAoverlap)[miRNAoverlapunique]] <- as.character(miRNA$location[subjectHits(miRNAoverlap)[miRNAoverlapunique]])
  #Assign previously described name
  locAnn$miRNAPrevDes <- rep("none", length(locAnn))
  locAnn$miRNAPrevDes[queryHits(miRNAoverlap)[miRNAoverlapunique]] <- as.character(miRNA$previouslyDescribed[subjectHits(miRNAoverlap)[miRNAoverlapunique]])

  #Do overlapType and transposon superfamily for cluster analysis
  TE_Class_DNA$type <- "TE_Class_DNA"
  TE_Class_RET$type <- "TE_Class_RET"
  elementsOrdered <- c(miRNA[,"type"],rRNA[,"type"],fiveprimeUTR[,"type"],threeprimeUTR[,"type"],
                       exons[,"type"],introns[,"type"],TE_Class_DNA[,"type"],TE_Class_RET[,"type"],
                       promoter[,"type"],MSAT[,"type"],irs[,"type"],trs[,"type"])
  transposonsOrdered <- c(TE_TCR1,TE_TOC2,TE_hAT,TE_Novosib,TE_Gulliver,TE_Harbinger,TE_Mariner,TE_P,TE_EnSpm,
                          TE_RTE,TE_L1,TE_DualenRandI,      
                          TE_Gypsy,TE_Copia,TE_TOC1,TE_REM1,TE_DIRS)
  transposonsOrdered$type <- transposonsOrdered$superfamily
  auxAnnotations <- list(overlapType=elementsOrdered,transposonType=transposonsOrdered)
  locAnn <- computeOverlaps(locAnn,auxAnnotations,assignNames = TRUE)

  locAnn
}

#Methylation overlaps using new methylation loci from Tom
methylation <- function(locAnn,annoDir = "resources") {
  #Require certain packages
  require(GenomicRanges)
  require(rtracklayer)
  #Load in files
  methCG=import.gff3(file.path(annoDir,"meth_data/chlamy_CGmeth.gff3"))
  methCHH=import.gff3(file.path(annoDir,"meth_data/chlamy_CHHmeth.gff3"))
  methCHG=import.gff3(file.path(annoDir,"meth_data/chlamy_CHGmeth.gff3"))
  #Ensure sequence levels are matched to the reference
  seqlevels(methCG) <- seqlevels(methCHG) <- seqlevels(methCHH) <- seqlevels(locAnn)
  #Create additional object that merged the three methylation files
  methAll <- c(methCG,methCHH,methCHG)
  #Make a list of methylation files for overlap computation
  methylation = list(methCG = methCG, methCHH = methCHH, methCHG = methCHG, methAll = methAll)
  #Compute overlaps
  locAnn <- computeOverlaps(locAnn,methylation)
  locAnn
}

#Phasing process and annotate
phaseMatch2 <- function(locAnn, annoDir = "resources",outputName="Pred_tab_2018.11.27_18.08") {
  #Load in packages
  require(GenomicRanges)
  require(rtracklayer)

  #Read in output
  phaseOutput <- read.table(file.path(annoDir,outputName),header=TRUE,stringsAsFactors = FALSE)
  #Build granges file out of table
  #find seqnames
  seqnamessplit <- strsplit(phaseOutput$ID,"_",fixed=TRUE)
  seqnamestemp <- paste0(lapply(seqnamessplit,"[[",1),"_",lapply(seqnamessplit,"[[",2))
  #Extract beginning and end coordinates
  begin <- as.numeric(lapply(strsplit(phaseOutput$Beg.End,":",fixed=TRUE),"[[",1))
  end <- as.numeric(lapply(strsplit(phaseOutput$Beg.End,":",fixed=TRUE),"[[",2))
  #generate GRanges object of phased loci for overlap computation
  phased <- GRanges(seqnames = seqnamestemp,IRanges(start=begin,end=end),
                    length=phaseOutput$Length,phasedRatio=phaseOutput$Phased_Ratio,phasedAbundance=phaseOutput$Phased_Abundance,
                    phasedNumber=phaseOutput$Phased_Number,phasedScore=phaseOutput$Phased_Score)
  
  # exporting to inspect in genome browser
  export.gff3(phased, con = file.path(annoDir, paste0(outputName, ".gff")))
  #Compute matches
  locAnn <- computeOverlaps(locAnn,list(phased=phased))
  overlaps <- findOverlaps(locAnn,phased)
  qHits <- queryHits(overlaps)
  sHits <- subjectHits(overlaps)
  #Default phaseScore is 0 (i.e. no phasing)
  locAnn$phaseScore <- rep(0,length(locAnn))
  #Assign phase scored where only one phased locus is present
  locAnn$phaseScore[qHits[table(qHits)==1]] <- phased$phasedScore[sHits[table(qHits)==1]]
  #For areas with more than one phased locus, assign the average phasing score 
  for(index in as.numeric(names(which(table(qHits)>1)))) {
    sHits == index
    locAnn$phaseScore[index] <- mean(phased$phasedScore[sHits[qHits == index]])
  }

  #Return locAnn
  locAnn
}
