# MPRA tag sequencing analysis pipeline

tag_analysis_package.R:
      A series of functions to assist in analysis of MPRA count tables

   * `oligoIsolate` - aggregates the count table based on barcodes resulting in one count for each replicate per oligo, if error, cigar, MD and alignment columns are present these are dropped for the purposes of counts
       * INPUTS:  counts table
       * OUTPUTS: oligo count data with the oligo as row names and the aggregate count data
   * `conditionStandard` - Standardize condition data for later use
       * INPUTS:  table of conditions (2 columns w/no header column 1 is replicates as found in the count table, column 2 indicates cell type)
       * OUTPUTS: condition table standardized for use later in the pipeline
   * `processAnalysis` - Performs an initial DESeq analysis
       * INPUTS:  
