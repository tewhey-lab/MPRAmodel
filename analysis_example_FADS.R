### This script is an example of how to use the Tag Analysis package

## The package requires attributes data, condition data, and count data for a project
## Using the FADS data located on the box as an example

attr_FADS <- read.delim("/Users/deweyh/Box/FADS/FADS_tile_snp.20190214.attributes")
count_FADS <- read.delim("/Users/deweyh/Box/FADS/mpra_data/OL13_FADS_K562_48h-Counts_20190828.out")
cond_FADS <- read.delim("/Users/deweyh/Box/FADS/mpra_data/OL13_FADS_K562_48h-Counts_20190828.cond", header = F, row.names = 1)

## Set the working directory for where you want your plots and results files to be written
setwd("</your/desired/working/directory/here>")

### Example input
## input includes the attribute, condition, and count data for the project, as well as the normalization method to be used.
## filePrefix should be the name of the project with some description
tagWrapper(count_FADS, attr_FADS, cond_FADS, filePrefix = 'FADS_example_remove_outliers', plotSave = T, method = 'ro')

# Please note that rendering the QC plots is very processor intensive with large data sets #
