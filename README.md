# MPRA count analysis pipeline

MPRAmodel.R:
      A series of functions to assist in analysis of MPRA count tables

## Quick Start:
```
	source("MPRAmodel.R")
	mpra_out <- MPRAmodel(<countsTable>, <attributesData>, <conditionData>, <filePrefix>, <negCtrlName>, <posCtrlName>, <projectName>)
```

## Arguments

### There are 3 files that this pipeline needs for input. <br>
   * `countsTable` : Table of counts with header. Should be either the output of the [MPRAcount](https://github.com/tewhey-lab/tag_analysis_WDL) pipeline, or in the same format. <br>
   * `attributesData` : Attributes for each oligo including oligo name, SNP, chromosome, position, reference allele, alt allele, Allele (ref/alt), window, strand, project, haplotype. For Oligos named in the form chr:pos:ref:alt:allele:window(:haplotype) the scripts [here](https://github.com/tewhey-lab/tag_analysis_WDL/blob/master/scripts/make_infile.py) and [here](https://github.com/tewhey-lab/tag_analysis_WDL/blob/master/scripts/make_attributes_oligo.pl) can be used to generate the attributes table. <br>
   * `conditionData` : 2 columns w/no header column 1 is replicates as found in the count table, column 2 indicates cell type. This can be done in the program of your choice, or taken from the output of [MPRAcount](https://github.com/tewhey-lab/tag_analysis_WDL).

### Other arguments needed <br>
  * `exclList` : ARRAY; List of cell types, in string format, that should be excluded from the analysis. This array is empty by default.
  * `filePrefix` : STRING; All written files will have this string included in the file name, if running the pipeline multiple times with  different settings this is a good place to differentiate them
  * `plotSave` : LOGICAL; Default `TRUE`, indicator of whether or not to save non-normalization QC plots
  * `altRef` :  LOGICAL; Default `TRUE`, indicator of how to sort alleles for `cellSpecificTtest`
  * `method` : STRING; Default "ss", indicator of method to be used to normalize the data
      * `'ss'` : Summit Shift - shifts the l2fc density peak to line up with 0
      * `'ssn'` : Summit Shift (Negative Controls Only) - shifts the peak of negative controls to 0
      * `'ro'` : Remove outliers - remove oligos that don't have a p-value or have a p-value > 0.001
      *`'nc'` : Negative Controls - normalize only the negative controls
  * `negCtrlName` : STRING; Default "negCtrl", how negative controls are indicated in the project column of the `attributesData` file
  * `posCtrlName` : STRING; Default "expCtrl", how positive controls are indicated in the project column of the `attributesData` file
  * `projectName` : STRING; Default "MPRA_PROJ", a generalized name for the overall project
  * `tTest` : LOGICAL; Default `TRUE`, perform cell type specific tTest to determine emVARs
  * `DEase` : LOGICAL; Default `TRUE`, use Differential Expression methods to detect Allele Specific Expression (this method can be found [here](http://rstudio-pubs-static.s3.amazonaws.com/275642_e9d578fe1f7a404aad0553f52236c0a4.html))
  * `correction` : STRING; Default "BH", indicator of whether to use Benjamini Hochberg ("BH") or Bonferroni ("BF") for p-value correction
  * `cutoff` : INT; Default 0.01, significance cutoff for including alleles for skew calculation (tTest only)
  * `upDisp` : LOGICAL; Default `TRUE`, update dispersions with cell type specific calculations
  * `prior` : LOGICAL; Default `TRUE`, use `betaPrior=T` when calculating the cell type specific dispersions


## Functions (in alphabetical order) <br>
  * `addHaplo` -  Add haplotype column to attribute data and resolve oligos with multiple projects, if oligos with multiple projects are listed separately, they will be collapsed into a single row, if oligos are identical but the SNP is different for them, an error will be thrown <br>
     _use_ :
      ```
      addHaplo(attributesData, negCtrlName, posCtrlName, projectName)
      ```
  * `cellSpecificTtest` - Function to perform TTest on individual cell types, determine the allelic skew, and identify emVARs <br>
     _use_ :
      ```
      cellSpecificTtest(attributesData, counts_norm, dups_output, ctrl_mean, exp_mean, ctrl_cols, exp_cols, altRef, correction, cutoff)
      ```
  * `conditionStandard` - Standardize condition data <br>
    _use_ :
    ```
    conditionStandard(conditionData)
    ```
  * `dataOut` - Retrieve output data for future functions - if only looking for results and not the plots this is the only function that needs to be called. Most functions should only need an output from here <br>
    _use_ :
    ```
    dataOut(countsData, attributesData, conditionData, exclList, altRef, file_prefix, method, negCtrlName, tTest, Dease, correction, cutoff, upDisp, prior)
    ```
  * `DESkew` - Function to determine the allelic skew of oligos and identify emVARs using the DESeq method to determine allelic skew <br>
    _use_ :
    ```
    DESkew(conditionData, counts_norm, attributesData, celltype, dups_output, prior)
    ```
  * `expandDups` - Expand IDs that denote duplicate oligos in count/DESeq results, or any dataframe with rownames that are separated by a ";" <br>
    _use_ :
    ```
    expandDups(output)
    ```
  * `fileDate` - Generate the date in YYYYMMDD format for separation of runs in the project <br>
    _use_ :
    ```
    fileDate()
    ```
  * `MPRAmodel` - Overall function which will run everything once called <br>
    _use_ :
    ```
    MPRAmodel(countsData, attributesData, conditionData, exclList, filePrefix, plotSave, altRef, method, negCtrlName, posCtrlName, projectName, tTest, Dease, correction, cutoff, upDisp, prior, …)
    ```
  * `mpraScatter` - Function to produce scatter plot of counts data for visualization of correlation between two samples <br>
      _use_ :
      ```
      mpraScatter(conditionData, countsOut, sampleX, sampleY, xmax, ymax, plotSave)
      ```
  * `oligoIsolate` - Remove Error, CIGAR, MD and position columns if necessary; aggregate count data with relation to the oligo, writes a copy of the collapsed raw counts  <br>
     _use_ :
    ```
    oligoIsolate(countsData, file_prefix)
    ```
  * `panel.cor` - Internal function to produce correlation scatter plots <br>
      _use_ :
      ```
      panel.cor(x, y, digits, prefix, cex.cor, …)
      ```
  * `panel.lm` - Internal function to produce correlation scatter plots <br>
      _use_ :
      ```
      panel.lm(x, y, pch, col.lm, …)
      ```
  * `panel.nlm` - Internal function to produce correlation scatter plots <br>
      _use_ :
      ```
      panel.nlm(x, y, pch, col.lm, …)
      ```
  * `plor_logFC` - Function to produce plots showing the expression fold change vs. normalized tag counts <br>
      _use_ :
      ```
      plot_logFC(full_output, sample, negCtrlName, posCtrlName)
      ```
  * `processAnalysis` - Initial processing of files via DESeq analysis <br>
      _use_ :
      ```
      processAnalysis(countsData, conditionData, exclList) <br>
      ```
  * `tagNorm` - Normalize DESeq results and plot normalized densities for each cell type  <br>
      _use_ :
      ```
      tagNorm(countsData, conditionData, attributesData, exclList, method, negCtrlName, upDisp, prior) <br>
      ```
  * `tagSig` - Replace dispersions of normalized dds with cell type specific dispersions, called within tagNorm does not need to be called by itself <br>
      _use_ :
      ```
      tagSig(dds_results, dds_rna, cond_data, exclList, prior)
      ```
