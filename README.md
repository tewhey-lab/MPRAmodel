# MPRA count analysis pipeline

tag_analysis_package.R:
      A series of functions to assist in analysis of MPRA count tables

There are 3 files that this pipeline needs for input. <br>
   * Counts Table: Table of counts with header. Should be either the output of the [ReplicateCount](https://github.com/tewhey-lab/tag_analysis_WDL) pipeline, or in the same format. <br>
   * Attributes Table: Attributes for each oligo including oligo name, SNP, chromosome, position, reference allele, alt allele, Allele (ref/alt), window, strand, project, haplotype. For Oligos named in the form chr:pos:ref:alt:allele:window(:haplotype) the script [here](https://github.com/tewhey-lab/tag_analysis_WDL/blob/master/scripts/make_attributes_oligo.pl) can be used to generate the attributes table. <br>
   * Condition Table: 2 columns w/no header column 1 is replicates as found in the count table, column 2 indicates cell type. This can be done in the program of your choice.

Description of Functions:
   * `addHaplo` - Make sure haplotype column is in the attribute data, and check/resolve samples from multiple projects
       * INPUTS:  Attributes Table, strings indicating negative controls, positive controls, and the overall project name
   * `oligoIsolate` - aggregates the count table based on barcodes resulting in one count for each replicate per oligo, if error, cigar, MD and alignment columns are present these are dropped for the purposes of counts
       * INPUTS:  Counts Table (see above description for details)
       * OUTPUTS: oligo count data with the oligo as row names and the aggregate count data
   * `conditionStandard` - Standardize condition data for later use
       * INPUTS:  Condition Table (see above description for details)
       * OUTPUTS: condition table standardized for use later in the pipeline
   * `processAnalysis` - Performs an initial DESeq analysis
       * INPUTS:  aggregated counts table, standardized condition table
       * OUTPUTS: initial dds_results object
   * `tagNorm` - Normalize DESeq results and plot normalized densities for each cell type
       * INPUTS:  aggregated counts table, standardized condition table, attributes table (see above description for details), normalization method ('ss' - summit shift, 'ro' - remove outliers, 'nc' - negative control only, 'mn' - median of ratios i.e. use DESeq normalization method), string indicating what you call negative controls in the attributes table
       * OUTPUTS: normalized dds_results object
   * `tagSig` - Called by `tagNorm`. Replaces dispersions in dds_results with cell type specific dispersions.
   * `dataOut` - Organize output data for future functions and writes results files. **NB** This function calls `oligoIsolate`, `conditionStandard`, `processAnalysis` and `tagNorm` it is not necessary to call any of those functions before this one.
       * INPUTS:  Counts Table, Attributes Table, Condition Table, altRef (logical), file prefix, normalization method, string indicating negative controls in attributes table
       * OUTPUTS: Standard results file for each cell type, cell type specific t-test results files
   * `expandDups` - expand Ids that denote duplicate oligos in count/DESeq results. Called by `dataOut`
   * `cellSpecificTtest` - performs TTests on inidividual cell types comparing alt/ref haplotypes. Called by `dataOut`
   * `panel.cor`, `panel.lm`, `panel.nlm` - functions which set up for correlation scatter plots. Called by `tagWrapper`.
   * `mpraScatter` - function to produce scatter plot of counts data for visualization of correlation between two samples. Called by `tagWrapper`
       * INPUTS:  standardized condition table, normalized count data, strings of sample names, maximum x and y values to include in plots
       * OUTPUTS: Scatter plots comparing two samples of either the same or different cell types.
   * `plot_logFC` - Produces plots showing expression fold change vs. normalized tag counts. Called by `tagWrapper`
   * `tagWrapper` - Function runs the whole analysis
       * INPUTS:  Counts Table, Attributes Table, Condition Table, file prefix, overall project name, negative control name, positive control name, normalization method, plot save (logical), altRef (logical)
       * OUTPUTS: plots normalization curves, writes standard results file for each cell type and cell type specific t-test results files, if plot save is TRUE the correlation tables, scatter plots and log fold change plots are saved.
