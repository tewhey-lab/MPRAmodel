
#Read in Count and Attributes tables
attr_PROJ <- read.delim("full/path/to/attributes/file", stringsAsFactors=F)
count_PROJ <- read.delim("/full/path/to/count/table", stringsAsFactors=F)

#Create Condition Table if not made in the MPRAcount pipeline
#This is an example of what it would look like with 

#total_reps <- 4
#cond_PROJ <- as.data.frame(c(rep("DNA",total_reps), rep("K562",total_reps)), stringsAsFactors=F)
#colnames(cond_PROJ) <- "condition"
#rownames(cond_PROJ) <- colnames(count_PROJ)[7:ncol(count_PROJ)]

source("/full/path/to/repository/tag_analysis_package.R")

#Run the analysis
tagWrapper(count_PROJ, attr_PROJ, cond_PROJ, filePrefix="example_prefix", projectName="your_base_project_name")

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
