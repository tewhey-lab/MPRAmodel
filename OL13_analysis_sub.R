
attr_FADS <- read.delim("/projects/tewhey-lab/projects/UKBB/OL13/setup/FADS_tile_snp.20190214.attributes", stringsAsFactors=F)
count_FADS <- read.delim("/projects/tewhey-lab/projects/UKBB/OL13/ReplicateCount_output/OL13_20210804.count", stringsAsFactors=F)

total_reps <- 4

cond_FADS <- as.data.frame(c(rep("DNA",total_reps), rep("K562",total_reps)), stringsAsFactors=F)
colnames(cond_FADS) <- "condition"
rownames(cond_FADS) <- colnames(count_FADS)[7:ncol(count_FADS)]

source("/projects/tewhey-lab/projects/UKBB/MPRA_count_analysis/tag_analysis_package.R")

#setwd("/projects/tewhey-lab/projects/UKBB/OL13/anlaysis_output/")

tagWrapper(count_FADS, attr_FADS, cond_FADS, filePrefix="OL13_20210811")

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
