### A package to perform MPRA Tag Analysis

# Required packages
# install.packages("ggplot2")
# if (!requireNamespace("BiocManager", quietly = TRUE)){
#   install.packages("BiocManager")}
# BiocManager::install("DESeq2")
# install.packages("reshape2")
library(DESeq2)
library(ggplot2)
library(reshape2)

### Generate the date in YYYYMMDD format for separation of runs in the project
fileDate <- function(){
  a <- Sys.Date()
  year <- format(a, "%Y")
  month <- format(a, "%m")
  day <- format(a, "%d")
  today <- paste0(year, month, day)
  return(today)
}

### Add Haplotype column to attribute data and resolve Oligos with multiple projects
### If oligos with multiple projects are listed separately, they will be collapsed into a single row
### If oligos are identical but the SNP is different for them, an error will be thrown
## INPUTS:
# attributesData  : table of full attributes, columns should include: ID, SNP, Project, Window, Strand, Allele,
  # Haplotype, Bash
# negCtrlName     : String which represents negative controls in your attributes table
## Returns: attributes table with ctrl_exp column which resolves oligos which serve as negative and positive controls
addHaplo <- function(attributesData,negCtrlName="negCtrl", posCtrlName="expCtrl", projectName="MPRA_PROJ"){
  if("Haplotype" %in% colnames(attributesData)){
    names(attributesData)[names(attributesData) == "Haplotype"] <- "haplotype"
  }
  if("Project" %in% colnames(attributesData)){
    names(attributesData)[names(attributesData) == "Project"] <- "project"
  }
  if("Allele" %in% colnames(attributesData)){
    names(attributesData)[names(attributesData) == "Allele"] <- "allele"
  }
  if(!("haplotype" %in% colnames(attributesData))){
    attributesData$haplotype <- "NA"
  }
  if(!("ref_allele" %in% colnames(attributesData))){
    attributesData$ref_allele <- "NA"
  }
  if(!("alt_allele" %in% colnames(attributesData))){
    attributesData$alt_allele <- "NA"
  }
  if("snp_pos" %in% colnames(attributesData)){
    names(attributesData)[names(attributesData) == "snp_pos"] <- "pos"
  }
  if(!("pos" %in% colnames(attributesData))){
    attributesData$pos <- "NA"
  }
  #Check for duplicate IDs, and combine projects if they exist

  id_freq <- as.data.frame(table(attributesData$ID))
  dup_ids <- id_freq[which(id_freq$Freq > 1),]
  if(nrow(dup_ids) > 0){
    for(id in dup_ids$Var1){
      tmp <- attributesData[which(attributesData$ID==id),]
      proj_new <- paste0(tmp$proj,collapse=",")
      attributesData$project[which(attributesData$ID==id)] <- proj_new
    }

    attributesData <- unique(attributesData)
    id_freq <- as.data.frame(table(attributesData$ID))
    dup_ids <- id_freq[which(id_freq$Freq > 1),]
    if(nrow(dup_ids) > 1){
      stop("Duplicate IDs found in Attributes file. Please check for mismatching information.")
    }
  }


  multi_project_check <- colsplit(attributesData$project, ",", names = c("proj1","proj2"))
  if(any(is.na(multi_project_check$proj2))){
    attributesData$ctrl_exp <- attributesData$project
  }
  else{
    attributesData$ctrl_exp <- projectName
    attributesData[grep(negCtrlName, attributesData$project),]$ctrl_exp <- negCtrlName
    attributesData[grep(posCtrlName, attributesData$project),]$ctrl_exp <- posCtrlName

  }
  return(attributesData)
}


### Remove Error, CIGAR, MD and position columns if necessary; aggregrate cound data with relation to the oligo
## INPUT:
  # countsData      : table of tag counts, columns should include: Barcode, Oligo, Sample names
## OUTPUT:
  # oligo count data with the oligo as row names and the aggregate count data
oligoIsolate <- function(countsData, file_prefix){
  if("Error" %in% colnames(countsData)){
    countsData <- countsData[,c(1,2,7:dim(countsData)[2])]
  }
  tag_counts <- aggregate(. ~Oligo, data=countsData[,-1], FUN = sum)
  counts_oligo <- tag_counts[,-1]
  rownames(counts_oligo) <- tag_counts[,1]
  write.table(counts_oligo, paste0("results/", file_prefix, "_", fileDate(), "_counts.out"), quote = F, sep="\t")
  return(counts_oligo)
}

### Standardize condition data
## INPUTS:
# conditionData  : table of conditions, 2 columns no header align to the variable column headers of countsData
  # and celltype
# data should initially be read in such that samples are row names, if reading in the condition table produced by MPRAcount
  # simply set row.names=1 when reading in the file.
## OUTPUT: standardized condition data, factorizing the cell types and categorizing DNA as 1 and RNA as 0
conditionStandard <- function(conditionData){
  cond_data <- as.data.frame(conditionData)
  colnames(cond_data)[1] <- "condition"
  cond_data[,1] <- factor(cond_data[,1])
  cond_data$condition <- relevel(cond_data$condition, "DNA")

  ## Not entirely sure this is necessary to keep in
  cond_data$DNA <- 0
  cond_data[cond_data$condition=="DNA",]$DNA <- 1
  ##

  return(cond_data)
}

### Initial processing of files via DESeq analysis
## INPUT:
# countsData      : table of tag counts, columns should include: Barcode, Oligo, Sample names
# conditionData   : table of conditions, 2 columns no header align to the variable column headers of countsData
  # and celltype
# exclList        : list of celltypes that should be excluded from the analysis. Empty by default
## OUTPUT: dds_results (initial)
processAnalysis <- function(countsData, conditionData, exclList=c()){
  # bring in count data and condition data from previously defined functions
  count_data <- countsData #oligoIsolate(countsData)
  cond_data <- conditionStandard(conditionData)

  # make sure that the column names and row names of the count and condition data match
  colnames(count_data) <- row.names(cond_data)

  # perform DESeq analysis
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = cond_data, design = ~condition)
  dds$condition <- relevel(dds$condition, "DNA")
  dds_results <- DESeq(dds, fitType = 'local', minReplicatesForReplace=Inf)

  return(list(dds,dds_results))
}

### Normalize DESeq results and plot normalized densities for each cell type
## INPUT:
# countsData      : table of tag counts, columns should include: Barcode, Oligo, Sample names
# attributesData  : table of full attributes, columns should include: ID, SNP, Project, Window, Strand, Allele,
  # Haplotype, Bash
# conditionData   : table of conditions, 2 columns no header align to the variable column headers of countsData
  # and celltype
# exclList        : list of celltypes that should be excluded from the anlalysis. Empty by default
# method          : normalization method. Default to 'ss'
  # 'ss' : Summit Shift - shifts the l2fc density peak to line up with 0
  # 'ssn': Summit Shift (Negative Controls Only) - shifts the peak of negative controls to 0
  # 'ro' : Remove outliers - remove oligos that don't have a p-value or have a p-value > 0.001
  # 'nc' : Negative Controls - normalize only the negative controls
## OUTPUT: dds_results (normalized), plots normalization curves
tagNorm <- function(countsData, conditionData, attributesData, exclList = c(), method = 'ss', negCtrlName="negCtrl", upDisp=T, prior=F){
  process <- processAnalysis(countsData, conditionData, exclList)
  dds_results <- process[[2]]
  dds <- process[[1]]
  dds_results_orig <- dds_results
  attribute_ids <- (attributesData[attributesData$ctrl_exp==negCtrlName,])$ID
  full_attribute_ids <- attributesData$ID
  #count_data <- oligoIsolate(countsData)
  count_data <- countsData
  cond_data <- conditionStandard(conditionData)
  colnames(count_data) <- row.names(cond_data)

  for(celltype in levels(cond_data$condition)){
    if(celltype=="DNA" | celltype %in% exclList) next
    message(celltype)
    # Create temporary results table for each cell type
    temp_outputA <- results(dds_results, contrast = c("condition", celltype, "DNA"), cooksCutoff=F, independentFiltering=F)
    # Summit shift normalization
    if(method == "ss"){
      message("summit shift")
      summit <- which.max(density(temp_outputA$log2FoldChange, na.rm=T)$y)
      log_offset <- 2^(density(temp_outputA$log2FoldChange, na.rm=T)$x[summit])
      sizeFactors(dds_results)[which(cond_data$condition == celltype)] <- sizeFactors(dds_results)[which(cond_data$condition == celltype)]*(log_offset)
    }
    # Summit shift normalization - negative controls only
    if(method == "ssn"){
      temp_outputA_neg <- temp_outputA[attribute_ids,]
      summit <- which.max(density(temp_outputA_neg$log2FoldChange, na.rm=T)$y)
      log_offset <- 2^(density(temp_outputA_neg$log2FoldChange, na.rm=T)$x[summit])
      sizeFactors(dds_results)[which(cond_data$condition == celltype)] <- sizeFactors(dds_results)[which(cond_data$condition == celltype)]*(log_offset)
    }
    # Remove outliers for normalization
    if(method == "ro"){
      temp_outputA_neg <- temp_outputA[full_attribute_ids,]
      attribute_ids <- row.names(temp_outputA_neg[!is.na(temp_outputA_neg$pvalue) & temp_outputA_neg$pvalue>0.001,])
    }
  }
  # Remove outliers or negative controls only
  if(method == "ro" | method == "nc"){
    dds_results_tmp<-estimateSizeFactors(dds[attribute_ids])
    sizeFactors(dds_results)<-sizeFactors(dds_results_tmp)
  }

  # dds_rna <- list()
  # # Celltype based DESeq analysis
  # for(celltype in levels(cond_data$condition)){
  #   if(celltype=="DNA" | celltype %in% exclList) next
  #   rna_cols <- cond_data[which(cond_data$condition==celltype),]
  #   rna_count <- count_data[,rownames(rna_cols)]
  #   dds_rna_temp <- DESeqDataSetFromMatrix(countData = rna_count, colData = rna_cols, design = ~1)
  #   sizeFactors(dds_rna_temp) <- sizeFactors(dds_results)[rownames(rna_cols)]
  #   dds_rna_temp <- estimateDispersions(dds_rna_temp)
  #   dds_rna[[celltype]] <- dds_rna_temp
  # }

  # # Replace dispersions in normalized dds with the celltype specific dispersions
  # if(upDisp==T){
  #   cellsp_dds <- tagSig(dds_results, dds_rna, cond_data, exclList, prior)
  # }
  # 
  # # Plot normalized density for each cell type -
  # for (celltype in levels(cond_data$condition)) {
  #   
  #   if(upDisp=T){
  #     dds_results <- cellsp_dds[[celltype]]
  #   }
  #   
  #   if(celltype == "DNA" | celltype %in% exclList) next
  # 
  #   message(celltype)
  #   temp_outputB <- results(dds_results_orig, contrast=c("condition",celltype,"DNA"), cooksCutoff=F, independentFiltering=F)
  # 
  #   outputA <- results(dds_results, contrast=c("condition",celltype,"DNA"), cooksCutoff=F, independentFiltering=F)
  # 
  #   message("Plotting Normalization Curves")
  #   pdf(paste0("plots/Normalized_FC_Density_",celltype,".pdf"),width=10,height=10)
  #   plot(density(temp_outputB[attribute_ids,"log2FoldChange"],na.rm=TRUE),xlim=c(-3,3),ylim=c(0,1.5),col="grey",main=paste0("Normalization - ",celltype))
  #   lines(density(temp_outputB$log2FoldChange,na.rm=TRUE),xlim=c(-3,3),col="black")
  #   lines(density(outputA$log2FoldChange,na.rm=TRUE),xlim=c(-3,3),col="red")
  #   lines(density(outputA[attribute_ids,"log2FoldChange"],na.rm=TRUE),xlim=c(-3,3),col="salmon")
  #   text(1.5,0.4,adj=c(0,0),labels="All - baseline",col="black")
  #   text(1.5,0.35,adj=c(0,0),labels="All - corrected",col="red")
  #   text(1.5,0.3,adj=c(0,0),labels=paste0(negCtrlName," - baseline"),col="grey")
  #   text(1.5,0.25,adj=c(0,0),labels=paste0(negCtrlName," - corrected"),col="salmon")
  #   abline(v=0)
  #   dev.off()
  # }
  return(list(dds_results,dds_results_orig))
}

### Replace dispersions of normalized dds with celltype specific dispersions
## INPUT:
# dds_results     : Normalized dds_results
# dds_rna         : List of cell type specific dds results
# cond_data       : Standardized condition data
# exclList        : list of celltypes that should be excluded from the analysis. Empty by default
# prior           : LOGICAL default T. Use betaPrior=T when running the nbinomWaldTest from DESeq2. Applies shrinkage to outlers
## OUTPUT: dds_results (normalized and celltype specific)
tagSig <- function(dds_results, dds_rna, cond_data, exclList=c(), prior=F){
  cellsp_dds <- list()
  for(celltype in levels(cond_data$condition)){
    if(celltype == "DNA" | celltype %in% exclList) next
    message("updating dispersions for:")
    message(celltype)
    dispersions(dds_rna[[celltype]])[which(is.na(dispersions(dds_rna[[celltype]])))] <- 10 #max(dispersions(dds_results))
    mcols(dds_results)$dispersion <- dispersions(dds_rna[[celltype]])
    dds_results <- nbinomWaldTest(dds_results, betaPrior = prior)
    cellsp_dds[[celltype]] <- dds_results
  }
  message(paste(dim(dds_results), collapse = "\t"))
  return(cellsp_dds)
}

### Retrieve output data for future functions - if only looking for results and not the plots this is the only function that needs to be called
  # Any subsequent functions should only need an output from here
## INPUT
# countsData      : table of tag counts, columns should include: Barcode, Oligo, Sample names
# attributesData  : table of full attributes, columns should include: ID, SNP, Project, Window, Strand, Allele,
  # Haplotype, Bash. **NB** If you are just running this function make sure to pass your attributes table through the addHaplo function
# conditionData   : table of conditions, 2 columns no header align to the variable column headers of countsData
  # and celltype
# exclList        : list of celltypes that should be excluded from the analysis. Empty by default
# altRef          : LOGICAL default T, indicating the order to sort alleles for allelic skew. alt/ref default
# file_prefix     : String to indicate what written file names include
# method          : Method to use for normalization
# negCtrlName     : String indicating what negative controls are called in the attributes table
# tTest           : LOGICAL default T, identify emVARs using the tTest method
# DEase           : LOGICAL default T, identify emVARS using the DESeq method of determining allelic skew
# correction      : String indicating whether to use Benjamini Hochberg ("BH", default) or Bonferroni ("BF") for p-value correction
# cutoff          : significance cutoff for including alleles for skew calculation (tTest only)
# upDisp          : LOGICAL default T, update dispersions with celltype specific calculations
# prior           : LOGICAL default T, use betaPrior=T when calculating the celltype specific dispersions.
## OUTPUT: writes duplicate output and ttest files for each celltype
dataOut <- function(countsData, attributesData, conditionData, exclList = c(), altRef = T, file_prefix, method = 'ss',negCtrlName="negCtrl",
                    tTest=T, DEase=T, cSkew=T, correction="BH", cutoff=0.01, upDisp=T, prior=F, paired=F){
  counts_data <- countsData[,c("Barcode","Oligo",rownames(conditionData))]
  message(paste0(colnames(counts_data), collapse = "\t"))
  count_data <- oligoIsolate(counts_data, file_prefix)
  message("Oligos isolated")
  dds_results_all <- tagNorm(count_data, conditionData, attributesData, exclList, method, negCtrlName, upDisp, prior)
  dds_results <- dds_results_all[[1]]
  dds_results_orig <- dds_results_all[[2]]
  attribute_ids <- (attributesData[attributesData$ctrl_exp==negCtrlName,])$ID
  full_attribute_ids <- attributesData$ID
  message("Tags Normalized")
  cond_data <- conditionStandard(conditionData)
  message("condition data standardized")
  colnames(count_data) <- row.names(cond_data)

  full_output<-list()
  full_output_var<-list()

  condition_table <- as.data.frame(cond_data)
  
  dds_rna <- list()
  # Celltype based DESeq analysis
  for(celltype in levels(cond_data$condition)){
    if(celltype=="DNA" | celltype %in% exclList) next
    rna_cols <- cond_data[which(cond_data$condition==celltype),]
    rna_count <- count_data[,rownames(rna_cols)]
    dds_rna_temp <- DESeqDataSetFromMatrix(countData = rna_count, colData = rna_cols, design = ~1)
    sizeFactors(dds_rna_temp) <- sizeFactors(dds_results)[rownames(rna_cols)]
    dds_rna_temp <- estimateDispersions(dds_rna_temp, fitType='local')
    dds_rna[[celltype]] <- dds_rna_temp
  }
  
  # Replace dispersions in normalized dds with the celltype specific dispersions
  if(upDisp==T){
    cellsp_dds <- tagSig(dds_results, dds_rna, cond_data, exclList, prior)
  }
  
  counts_norm_all <- counts(dds_results, normalized = T)
  
  write.table(counts_norm_all,paste0("results/", file_prefix, "_", fileDate(),"_normalized_counts.out"), quote = F, sep = "\t")
  

  for (celltype in levels(cond_data$condition)) {
    if(celltype == "DNA" | celltype %in% exclList) next
    message(celltype)
    
    if(upDisp==T){
      dds_results <- cellsp_dds[[celltype]]
    }
    
    counts_norm <- counts(dds_results, normalized = T)
    counts_norm <- counts_norm[,colnames(counts_norm) %in% rownames(cond_data)[which(condition_table$condition %in% c("DNA",celltype))]]
    
    write.table(counts_norm,paste0("results/", file_prefix, "_", fileDate(),"_",celltype, "_normalized_counts.out"), quote = F, sep = "\t")
    
    if(DEase==T){
      message("Removing count duplicates")
      counts_DE <- counts(dds_results, normalized=F)
      counts_norm_DE <- expandDups(counts_DE)
    }

    outputA <- results(dds_results, contrast=c("condition",celltype,"DNA"), cooksCutoff=F, independentFiltering=F)

    temp_outputB <- results(dds_results_orig, contrast=c("condition",celltype,"DNA"), cooksCutoff=F, independentFiltering=F)
    
    message("Plotting Normalization Curves")
    pdf(paste0("plots/",file_prefix,"_",fileDate(),"_Normalized_FC_Density_",celltype,".pdf"),width=10,height=10)
    plot(density(temp_outputB[attribute_ids,"log2FoldChange"],na.rm=TRUE),xlim=c(-3,3),ylim=c(0,1.5),col="grey",main=paste0("Normalization - ",celltype))
    lines(density(temp_outputB$log2FoldChange,na.rm=TRUE),xlim=c(-3,3),col="black")
    lines(density(outputA$log2FoldChange,na.rm=TRUE),xlim=c(-3,3),col="red")
    lines(density(outputA[attribute_ids,"log2FoldChange"],na.rm=TRUE),xlim=c(-3,3),col="salmon")
    text(1.5,0.4,adj=c(0,0),labels="All - baseline",col="black")
    text(1.5,0.35,adj=c(0,0),labels="All - corrected",col="red")
    text(1.5,0.3,adj=c(0,0),labels=paste0(negCtrlName," - baseline"),col="grey")
    text(1.5,0.25,adj=c(0,0),labels=paste0(negCtrlName," - corrected"),col="salmon")
    abline(v=0)
    dev.off()
    
    message("Results of dds_results recieved")

    ctrl_cols <- row.names(condition_table[condition_table$condition=="DNA",])
    exp_cols <- row.names(condition_table[condition_table$condition==celltype,])

    message("Control and Experiment Columns Set")

    DNA_mean <- rowMeans(count_data[, colnames(count_data) %in% ctrl_cols], na.rm=T)
    ctrl_mean <- rowMeans(counts_norm[, colnames(counts_norm) %in% ctrl_cols], na.rm = T)
    exp_mean <- rowMeans(counts_norm[, colnames(counts_norm) %in% exp_cols], na.rm = T)
    output_2 <- cbind(DNA_mean,ctrl_mean,exp_mean,outputA[,-1])

    counts_norm_sp <- expandDups(counts_norm)

    message("Removing Duplicates")
    dups_output<-expandDups(output_2)
    # message(paste0(colnames(dups_output),collapse = "\t"))

    message("Writing Standard Results File")
    full_outputA<-merge(attributesData, as.data.frame(dups_output), by.x="ID", by.y="row.names", all.x=TRUE, no.dups=F)
    full_output[[celltype]]<-full_outputA
    write.table(full_outputA, paste0("results/", file_prefix, "_", celltype, "_", fileDate(), ".out"), row.names=F, col.names=T, sep="\t", quote=F)

    if(tTest==T){
      message("Writing T-Test Results File")
      outA<-cellSpecificTtest(attributesData, counts_norm_sp, dups_output, ctrl_mean, exp_mean, ctrl_cols, exp_cols, altRef, correction, cutoff, prior)
      full_output_var[[celltype]]<-outA
      write.table(outA,paste0("results/", file_prefix, "_", celltype, "_emVAR_ttest_", fileDate(),".out"), row.names=F, col.names=T, sep="\t", quote=F)
      # write.table(attributesData, "results/ttest_attributes.txt", row.names = F, quote = F, sep = "\t")
      # write.table(counts_norm_sp, paste0("results/ttest_",celltype,"_norm_counts.txt"), quote = F, sep = "\t")
      # write.table(dups_output, "results/ttest_dups_out.txt", quote = F, sep = "\t")
    }

    if(DEase==T){
      message("Writing DESeq Allelic Skew Results File")
      if(paired==F){
        cond_pass <- condition_table
      }
      if(paired==T){
        plas_row <- length(conditionData$condition[which(conditionData=="DNA")])
        message(plas_row)
        cell_row <- length(conditionData$condition[which(conditionData$condition==celltype)])
        message(cell_row)
        if(plas_row==cell_row){
          cond_pass <- condition_table
        }
        else if(plas_row!=cell_row){
          if(plas_row > cell_row){
            drop_num <- plas_row - cell_row
            message(paste0("dropping ", drop_num," DNA replicates for paired GLM analysis"))
            plas_ids <- data.frame(matrix(ncol = 2, nrow=plas_row))
            colnames(plas_ids) <- c("rep","bc_count")
            plas_ids$rep <- rownames(conditionData)[which(conditionData$condition=="DNA")]
            
            for(id in plas_ids$rep){
              # message(id)
              # message(sum(counts_data[,id] > 0))
              
              plas_ids[plas_ids$rep==id,"bc_count"] <- sum(counts_data[,id] > 0)
            }
            plas_ids <- plas_ids[order(plas_ids$bc_count),]
            
            plas_drop <- plas_ids$rep[1:drop_num]
            
            message(paste0("Dropping: ",plas_drop, collapse = "\t"))
            cond_pass <- condition_table[which(!rownames(condition_table)%in%plas_drop),]
          }
          else if(cell_row > plas_row){
            drop_num <- cell_row - plas_row
            message(paste0("dropping ", drop_num," ",celltype, " replicates for paired GLM analysis"))
            cell_ids <- data.frame(matrix(ncol = 2, nrow=cell_row))
            colnames(cell_ids) <- c("rep","bc_count")
            cell_ids$rep <- rownames(conditionData)[which(conditionData$condition==celltype)]
            message(paste0(cell_ids$rep, collapse = "\t"))
          
            message(class(cell_ids))
            
            for(id in cell_ids$rep){
             
              cell_ids[cell_ids$rep==id,"bc_count"] <- sum(counts_data[,id] > 0)
            }
            cell_ids <- cell_ids[order(cell_ids$bc_count),]
            
            cell_drop <- cell_ids$rep[1:drop_num]
            
            message(paste0("Dropping: ",cell_drop, collapse = "\t"))
            cond_pass <- condition_table[which(!rownames(condition_table)%in%cell_drop),]
          }
        }
        counts_norm_DE <- counts_norm_DE[,which(colnames(counts_norm_DE) %in% rownames(cond_pass))]
      }
      outB <- DESkew(cond_pass, counts_norm_DE, attributesData, celltype, dups_output,prior, cutoff, paired)
      if(paired==F){
        write.table(outB,paste0("results/", file_prefix, "_", celltype, "_emVAR_glm_", fileDate(),".out"), row.names=F, col.names=T, sep="\t", quote=F)
      }
      if(paired==T){
        write.table(outB,paste0("results/", file_prefix, "_", celltype, "_emVAR_glm_paired_", fileDate(),".out"), row.names=F, col.names=T, sep="\t", quote=F)
      }
    }

    message("Writing bed File")
    full_bed_outputA<-merge(attributesData, as.matrix(dups_output),by.x="ID",by.y="row.names",all.x=TRUE,no.dups=FALSE)
    # message(paste0(colnames(full_bed_outputA), collapse = "\t"))
    #printbed<-full_bed_outputA[,c("chr","start","stop","ID","strand","log2FoldChange","Ctrl.Mean","Exp.Mean","pvalue","padj","lfcSE","cigar","md-tag","project")]
    if(!(c("start","stop") %in% colnames(full_bed_outputA))){
      full_bed_outputA$start <- ""
      full_bed_outputA$stop <- ""
    }
    printbed<-full_bed_outputA[,c("chr","start","stop","ID","strand","log2FoldChange","ctrl_mean","exp_mean","pvalue","padj","lfcSE","project")]
    printbed$score<-"."
    #printbed<-printbed[,c("chr","start","stop","ID","score","strand","log2FoldChange","Ctrl.Mean","Exp.Mean","pvalue","padj","lfcSE","cigar","md-tag","project")]
    #colnames(printbed)<-c("chr","start","stop","id","score","strand","log2fc","input-count","output-count","log10pval","log10fdr","lfc-se","cigar","md-tag","project")
    printbed<-printbed[,c("chr","start","stop","ID","score","strand","log2FoldChange","ctrl_mean","exp_mean","pvalue","padj","lfcSE","project")]
    colnames(printbed)<-c("chr","start","stop","id","score","strand","log2fc","input-count","output-count","log10pval","log10fdr","lfc-se","project")
    printbed$strand[printbed$strand=="fwd"]="+"
    printbed$strand[printbed$strand=="rev"]="-"
    printbed$log10pval=-log10(printbed$log10pval)
    printbed$log10fdr=-log10(printbed$log10fdr)

    write.table(printbed,paste0("results/",file_prefix,"_",celltype,"_",fileDate(),".bed"),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
  }
  # if(cSkew ==T){
  #   cellSkew(cond_data, counts_norm_DE, attributesData,file_prefix)
  # }
  return(c(full_output, dds_results))
}

### Expand IDs that denote duplicate oligos in count/DESeq results
# Requires the output_2 of dataOut function - This is a df containing the control (DNA) means and experimental (RNA) means and
  # normalized-celltype specific DESeq results
## Returns: Final output table with no duplicates - Table in the same format as above, except with IDs that were previously several oligos
  ## separated by semi-colons have been separated into separate rows.
expandDups <- function(output){
  output_orig <- output
  if(class(output_orig)=="matrix"){
    output_orig <- as.data.frame(output_orig)
  }
  output_new <- cbind(rownames(output_orig), output_orig)
  colnames(output_new)[1] <- "Row.names"
  # Identify duplicates, if any exist
  message("identifying duplicates")
  dups <- output_new[grep("\\(.*\\)$",output_new$Row.names),]
  dups$Row.names <- gsub("^\\((.*)\\)$","\\1",dups$Row.names)
  # Add everything but the duplicates to the final output
  message("resolving duplicates")
  output_final <- output_new[-(grep("\\(.*\\)$",output_new$Row.names)),]
  output_final <- output_final[!(is.na(output_final$Row.names)),]
  if(nrow(dups) > 0) {
    for(i in 1:nrow(dups)){
      dup_id <- unlist(strsplit(dups$Row.names[i],";"))
      dup_exp <- dups[rep(i, length(dup_id)), ]
      dup_exp$Row.names <- dup_id
      output_final <- rbind(dup_exp,output_final)
    }
    rownames(output_final) <- output_final[,1]
    output_final <- output_final[,-1]
  }else {
    output_final <- output_orig
  }
  return(output_final)
}

### Function to perform TTest on individual cell types
# attributesData  : table of full attributes, columns should include: ID, SNP, Project, Window, Strand, Allele,
  # Haplotype, Bash
# counts_norm     : passed from the dataOut function
# dups_output     : passed from the expandDups function
# ctrl_mean       : passed from the dataOut function
# exp_mean        : passed from the dataOut function
# ctrl_cols       : passed from the dataOut function
# exp_cols        : passed from the dataOut function
# altRef          : Logical, default T indicating sorting by alt/ref, if sorting ref/alt set to F
# correction      : String indicating what correction method to use for significance cutoff, "BH" for Benjamini Hochburg and "BF" for Bonferroni
# cutoff          : Integer for significance cutoff to use.
cellSpecificTtest<-function(attributesData, counts_norm, dups_output, ctrl_mean, exp_mean, ctrl_cols, exp_cols, altRef = T, correction="BH", cutoff=0.01, prior=F){
  snp_data <- subset(attributesData,allele=="ref" | allele=="alt")
  snp_data$comb <- paste(snp_data$SNP,"_",snp_data$window,"_",snp_data$strand,"_",snp_data$haplotype,sep="")
  tmp_ct <- as.data.frame(table(snp_data$comb))

  snp_data_pairs <- snp_data[snp_data$comb %in% tmp_ct[tmp_ct$Freq==2,]$Var1,]

  snp_data_rejected <- snp_data[snp_data$comb %in% tmp_ct[tmp_ct$Freq!=2,]$Var1,]

  snp_data_ctdata_pairs <- merge(snp_data_pairs, counts_norm, by.x="ID", by.y="row.names", all.x=T, no.dups=F)
  snp_data_ctdata_pairs <- snp_data_ctdata_pairs[order(snp_data_ctdata_pairs$SNP, snp_data_ctdata_pairs$window, snp_data_ctdata_pairs$strand, snp_data_ctdata_pairs$haplotype, snp_data_ctdata_pairs$allele),]
  snp_data_expdata_pairs <- merge(snp_data_pairs,dups_output, by.x="ID", by.y="row.names", all.x=T, no.dups=F)
  snp_data_expdata_pairs <- snp_data_expdata_pairs[order(snp_data_expdata_pairs$SNP, snp_data_expdata_pairs$window, snp_data_expdata_pairs$strand, snp_data_expdata_pairs$haplotype, snp_data_expdata_pairs$allele),]

  if(altRef==T){
    evens <- seq(2, nrow(snp_data_pairs), by = 2)
    odds <- seq(1, nrow(snp_data_pairs), by = 2)
  }else{
    evens <- seq(1, nrow(snp_data_pairs), by=2)
    odds <- seq(2, nrow(snp_data_pairs), by=2)
  }

  out <- snp_data_pairs[which(snp_data_pairs$allele=="ref"),c("ID","SNP","chr","pos","ref_allele","alt_allele","allele","window","strand","haplotype")]
  out$comb <- paste0(out$SNP,"_",out$window,"_",out$strand,"_",out$haplotype)
  
  out <- merge(out, snp_data_expdata_pairs[which(snp_data_expdata_pairs$allele=="ref"),c("ID","comb","ctrl_mean","exp_mean","log2FoldChange","lfcSE","pvalue","padj")], by = c("comb","ID"))
  colnames(out)[12:17] <- c("A_Ctrl_Mean","A_Exp_Mean","A_log2FC","A_log2FC_SE","A_logP","A_logPadj_BH")
  
  out$A_logPadj_BH <- -log10(out$A_logPadj_BH)
  out$A_logPadj_BH[out$A_logPadj_BH < 0] <- 0
  out$A_logPadj_BH[out$A_logPadj_BH == Inf] <- max(out$A_logPadj_BH[is.finite(out$A_logPadj_BH)])
  out$A_logPadj_BF <- -log10(out$A_logP*(nrow(snp_data_expdata_pairs)/2))
  out$A_logP <- -log10(out$A_logP)
  out$A_logP[is.na(out$A_logP)] <- 0
  out$A_logP[out$A_logP == Inf] <- max(out$A_logP[is.finite(out$A_logP)])
  out$A_logPadj_BF[out$A_logPadj_BF < 0] <- 0
  out$A_logPadj_BF[out$A_logPadj_BF == Inf] <- max(out$A_logPadj_BF[is.finite(out$A_logPadj_BF)])
  out <- merge(out, snp_data_expdata_pairs[which(snp_data_expdata_pairs$allele=="alt"),c("comb","ctrl_mean","exp_mean","log2FoldChange","lfcSE","pvalue","padj")], by = "comb")
  colnames(out)[19:24] <- c("B_Ctrl_Mean","B_Exp_Mean","B_log2FC","B_log2FC_SE","B_logP","B_logPadj_BH")
  
  out$B_logPadj_BH <- -log10(out$B_logPadj_BH)
  out$B_logPadj_BH[out$B_logPadj_BH < 0] <- 0
  out$B_logPadj_BH[out$B_logPadj_BH == Inf] <- max(out$B_logPadj_BH[is.finite(out$B_logPadj_BH)])
  out$B_logPadj_BF <- -log10(out$B_logP*(nrow(snp_data_expdata_pairs)/2))
  out$B_logP <- -log10(out$B_logP)
  out$B_logP[is.na(out$B_logP)] <- 0
  out$B_logP[out$B_logP == Inf] <- max(out$B_logP[is.finite(out$B_logP)])
  out$B_logPadj_BF[out$B_logPadj_BF < 0] <- 0
  out$B_logPadj_BF[out$B_logPadj_BF == Inf] <- max(out$B_logPadj_BF[is.finite(out$B_logPadj_BF)])

  out2 <- out#[,c(1:12, 16:19, 26:28, 23:25)]

   # Don't try to do the t test for ones with all zeros.
  ignore_idx <- which(rowMeans(snp_data_ctdata_pairs[odds,ctrl_cols]) < 10 | rowMeans(snp_data_ctdata_pairs[odds, exp_cols]) < 10 |
                        is.na(rowMeans(snp_data_ctdata_pairs[odds,ctrl_cols]))  | is.na(rowMeans(snp_data_ctdata_pairs[odds, exp_cols])) |
                        rowMeans(snp_data_ctdata_pairs[evens,ctrl_cols]) < 10 | rowMeans(snp_data_ctdata_pairs[evens, exp_cols]) < 10 |
                        is.na(rowMeans(snp_data_ctdata_pairs[evens,ctrl_cols]))  | is.na(rowMeans(snp_data_ctdata_pairs[evens, exp_cols])) )

  # For the numerator, set zero values to 1 so that the log-ratio is defined.
  counts1 <- snp_data_ctdata_pairs
  counts1[counts1 == 0] <- 1

  # t test
  ratios_A <- log2((counts1[evens, exp_cols]) / rowMeans(snp_data_ctdata_pairs[evens, ctrl_cols]))
  ratios_B <- log2((counts1[odds, exp_cols]) / rowMeans(snp_data_ctdata_pairs[odds, ctrl_cols]))

  ratios_list <- list(ratios_A, ratios_B)

  t_pvalue <- sapply(1:nrow(ratios_A), function(i) if(i %in% ignore_idx){NA} else{
        t.test(as.numeric(ratios_A[i,]), as.numeric(ratios_B[i,]), var.equal=F, paired=T)$p.value})
  t_stat <- sapply(1:nrow(ratios_A), function(i) if(i %in% ignore_idx){NA} else{
        t.test(as.numeric(ratios_A[i,]), as.numeric(ratios_B[i,]), var.equal=F, paired=T)$statistic})

  if(prior==T){
    mean_ratios_A <- log2(rowMeans(counts1[evens,exp_cols], na.rm = T) / rowMeans(snp_data_ctdata_pairs[evens, ctrl_cols], na.rm = T))
    mean_ratios_B <- log2(rowMeans(counts1[odds,exp_cols], na.rm = T) / rowMeans(snp_data_ctdata_pairs[odds, ctrl_cols], na.rm = T))

    out2$Log2Skew <- mean_ratios_B - mean_ratios_A
  }

  if(prior==F){
    out2$Log2Skew <- out2$B_log2FC - out2$A_log2FC
  }
  out2$LogSkew_SE <- sqrt(out2$A_log2FC_SE^2+out2$B_log2FC_SE^2)
  out2$Skew_logP <- ifelse(is.na(t_pvalue), 0, -log10(t_pvalue))

  OE_threshold <- -log10(cutoff)
  if(correction=="BF"){
    is_OE <- out2$A_logPadj_BF >= OE_threshold | out2$B_logPadj_BF >= OE_threshold
  }
  if(correction=="BH"){
    is_OE <- out2$A_logPadj_BH >= OE_threshold | out2$B_logPadj_BH >= OE_threshold
  }
  out2$Skew_logFDR <- -log10(p.adjust(t_pvalue, method = "BH"))
  out2$Skew_logFDR_act <- rep(NA, nrow(out))
  q_idx <- intersect(which(is_OE), which(!is.na(t_pvalue)))
  out2$Skew_logFDR_act[q_idx] <- -log10(p.adjust(t_pvalue[q_idx],method="BH"))

  return(out2)
}

### Function to perform DESeq version of Allelic Skew

DESkew <- function(conditionData, counts_norm, attributesData, celltype, dups_output,prior, cutoff, paired=F){
  
  ds_cond_data <- as.data.frame(conditionData[which(conditionData$condition=="DNA" | conditionData$condition==celltype),,drop=F])
  
  # Prepare the sample table
  message("Preparing Sample Table")
  dna_reps <- nrow(as.data.frame(ds_cond_data[which(ds_cond_data$condition=="DNA"),]))
  rna_reps <- nrow(as.data.frame(ds_cond_data[which(ds_cond_data$condition==celltype),]))
  avg_reps <- (dna_reps+rna_reps)/2
  total_cond <- length(unique(ds_cond_data$condition))
  samps <- data.frame(material=factor(c(rep("DNA",dna_reps*total_cond),rep("RNA",rna_reps*total_cond))),
                      allele=factor(rep(c("ref","alt"),((dna_reps+rna_reps))), levels = c("ref","alt")),
                      sample=factor(rep(c(rownames(ds_cond_data)[which(ds_cond_data$condition=="DNA")], rownames(ds_cond_data)[which(ds_cond_data$condition==celltype)]),each=total_cond)))
  
  
  # Reorganize the count data
  message("reorganizing count data")
  snp_data <- subset(attributesData,allele=="ref" | allele=="alt")

  snp_data$comb <- paste(snp_data$SNP,"_",snp_data$window,"_",snp_data$strand,"_",snp_data$haplotype,sep="")

  tmp_ct <- as.data.frame(table(snp_data$comb))
  
  snp_data_pairs <- snp_data[snp_data$comb %in% tmp_ct[tmp_ct$Freq==2,]$Var1,]

  pairs_ids <- snp_data_pairs$ID
  dups_ids <- rownames(dups_output)
  
  snp_data_pairs <- merge(snp_data_pairs,dups_output, by.x="ID", by.y="row.names", all=T, no.dups=F)
  
  message("COUNT OUT")
  out <- snp_data_pairs[which(snp_data_pairs$allele=="ref"),c("ID","SNP","chr","pos","ref_allele","alt_allele","allele","window","strand","haplotype")]
  out$comb <- paste0(out$SNP,"_",out$window,"_",out$strand,"_",out$haplotype)
  
  out <- merge(out, snp_data_pairs[which(snp_data_pairs$allele=="ref"),c("ID","comb","ctrl_mean","exp_mean","log2FoldChange","lfcSE","pvalue","padj")], by = c("comb","ID"))
  colnames(out)[12:17] <- c("A_Ctrl_Mean","A_Exp_Mean","A_log2FC","A_log2FC_SE","A_logP","A_logPadj_BH")

  out$A_logPadj_BH <- -log10(out$A_logPadj_BH)
  out$A_logPadj_BH[out$A_logPadj_BH < 0] <- 0
  out$A_logPadj_BH[out$A_logPadj_BH == Inf] <- max(out$A_logPadj_BH[is.finite(out$A_logPadj_BH)])
  out$A_logPadj_BF <- -log10(out$A_logP*(nrow(snp_data_pairs)/2))
  out$A_logP <- -log10(out$A_logP)
  out$A_logP[is.na(out$A_logP)] <- 0
  out$A_logP[out$A_logP == Inf] <- max(out$A_logP[is.finite(out$A_logP)])
  out$A_logPadj_BF[out$A_logPadj_BF < 0] <- 0
  out$A_logPadj_BF[out$A_logPadj_BF == Inf] <- max(out$A_logPadj_BF[is.finite(out$A_logPadj_BF)])
  out <- merge(out, snp_data_pairs[which(snp_data_pairs$allele=="alt"),c("comb","ctrl_mean","exp_mean","log2FoldChange","lfcSE","pvalue","padj")], by = "comb")
  colnames(out)[19:24] <- c("B_Ctrl_Mean","B_Exp_Mean","B_log2FC","B_log2FC_SE","B_logP","B_logPadj_BH")
  
  out$B_logPadj_BH <- -log10(out$B_logPadj_BH)
  out$B_logPadj_BH[out$B_logPadj_BH < 0] <- 0
  out$B_logPadj_BH[out$B_logPadj_BH == Inf] <- max(out$B_logPadj_BH[is.finite(out$B_logPadj_BH)])
  out$B_logPadj_BF <- -log10(out$B_logP*(nrow(snp_data_pairs)/2))
  out$B_logP <- -log10(out$B_logP)
  out$B_logP[is.na(out$B_logP)] <- 0
  out$B_logP[out$B_logP == Inf] <- max(out$B_logP[is.finite(out$B_logP)])
  out$B_logPadj_BF[out$B_logPadj_BF < 0] <- 0
  out$B_logPadj_BF[out$B_logPadj_BF == Inf] <- max(out$B_logPadj_BF[is.finite(out$B_logPadj_BF)])
  
  counts_allele_id <- merge(snp_data_pairs[,c("ID","comb","allele","SNP","window","strand","haplotype")], counts_norm, by.x="ID", by.y="row.names", all.x=T)

  rownames(counts_allele_id) <- counts_allele_id$ID
  
  message("done")
  
  message("COUNTS REF")
  
  counts_ref_all <- counts_allele_id[which(counts_allele_id$allele=="ref"),,drop=F]
  counts_ref <- counts_ref_all[,colnames(counts_ref_all) %in% rownames(ds_cond_data),drop=F]
 
  colnames(counts_ref) <- paste0(colnames(counts_ref),"_ref")
  counts_ref$comb <- paste0(counts_ref_all$SNP,"_",counts_ref_all$window,"_",counts_ref_all$strand,"_",counts_ref_all$haplotype)
  counts_ref$ID <- rownames(counts_ref)

  counts_ref <- merge(counts_ref, snp_data_pairs[which(snp_data_pairs$ID %in% rownames(counts_ref)),c("ID","SNP","chr","pos","ref_allele","alt_allele","allele","window","strand","haplotype","comb")],by=c("ID","comb"),all.x=T)
  message("done")  
  
  message("COUNTS ALT")

  counts_alt_all <- counts_allele_id[which(counts_allele_id$allele=="alt"),,drop=F]
  counts_alt <- counts_alt_all[,colnames(counts_alt_all) %in% rownames(ds_cond_data),drop=F]
  
  colnames(counts_alt) <- paste0(colnames(counts_alt),"_alt")
  counts_alt$comb <- paste0(counts_alt_all$SNP,"_",counts_alt_all$window,"_",counts_alt_all$strand,"_",counts_alt_all$haplotype)
  counts_alt$ID <- rownames(counts_alt)
  counts_alt <- merge(counts_alt, snp_data_pairs[which(snp_data_pairs$ID %in% rownames(counts_alt)),c("ID","SNP","chr","pos","ref_allele","alt_allele","allele","window","strand","haplotype","comb")],by=c("ID","comb"),all.x=T)
  
  message("done")
  # message(paste0(colnames(counts_alt), collapse = "\t"))
  # message(paste0(colnames(counts_ref), collapse = "\t"))
 
  counts_ref_alt <- merge(counts_ref, counts_alt, by=c("SNP","chr","pos","ref_allele","alt_allele","window","strand","haplotype","comb"), all=T)
  
  message("counts merged") 

  column_order <- data.frame(allele=factor(rep(c("ref","alt"),((dna_reps+rna_reps))), levels = c("ref","alt")),
                             sample=factor(rep(c(rownames(ds_cond_data)[which(ds_cond_data$condition=="DNA")], rownames(ds_cond_data)[which(ds_cond_data$condition==celltype)]),each=total_cond)))
  column_order$order <- paste0(column_order$sample,"_",column_order$allele)
  
  message("order set")
  message(paste0(colnames(counts_ref_alt), collapse = "\t"))
  counts_ref_alt <- counts_ref_alt[,c("ID.x","SNP","chr","pos","ref_allele","alt_allele","allele.x","window","strand", "haplotype",column_order$order)]
  colnames(counts_ref_alt) <- c("ID","SNP","chr","pos","ref_allele","alt_allele","allele","window","strand","haplotype",column_order$order)
  message(paste0(dim(counts_ref_alt), collapse = "\t"))

  counts_comp <- counts_ref_alt[complete.cases(counts_ref_alt[,column_order$order]),]
  counts_mat <- as.matrix(counts_comp[,column_order$order])
  ids_comp <- counts_comp$ID
  
  message("ids pulled")
 
  # Set Design definition
  design <- ~sample + allele
  
  # Run DESeq analysis
  message("Running DESeq")
  dds <- DESeqDataSetFromMatrix(counts_mat, samps, design)
  
  if(paired==F){
    design(dds) <- ~material + allele + material:allele

    dds <- DESeq(dds, fitType = "local", minReplicatesForReplace = Inf)
    
    res.diff <- results(dds, name="materialRNA.allelealt",cooksCutoff=F,independentFiltering=F)

  }
  
  if(paired==T){
    message("running paired samples")
    
    sample_lets <- c(rep(LETTERS[1:dna_reps], each=total_cond), rep(LETTERS[1:rna_reps], each=total_cond))
    dds$sample.n <- as.factor(sample_lets)
    
    design(dds) <- ~material + material:sample.n + material:allele
    sizeFactors(dds) <- rep(1, (dna_reps+rna_reps)*total_cond)
    
    # if(dna_reps != rna_reps){
    #   warning("Number of DNA replicates and RNA replicates unequal. Using the lower number of replicates to run paired samples. If you want to avoid this run MPRAmodel with `paired=F`")
    #   mm <- model.matrix(~material + material:sample.n + material:allele, colData(dds))
    #   col_mm <- ncol(mm)
    #   mm <- mm[,c(1:((min(dna_reps,rna_reps))*2),(col_mm-1),col_mm)]
    #   dds <- DESeq(dds, full = mm, fitType = "local", minReplicatesForReplace=Inf)
    # }
    # else{
    dds <- DESeq(dds, fitType = "local", minReplicatesForReplace = Inf)
    # }

    #Get the skew results
    # cell_res <- paste0("condition",celltype,".countalt")
    message(paste0(resultsNames(dds), collapse = "\t"))
    
    res.diff <- results(dds, contrast=list("materialRNA.allelealt","materialDNA.allelealt"), cooksCutoff=F, independentFiltering=F)
  }
  
  # res.diff <- results(dds, name="materialRNA.allelealt",cooksCutoff=F,independentFiltering=F)
  
  res.diff <- as.data.frame(res.diff)[,-1]
  colnames(res.diff) <- c("Log2Skew","Skew_SE","skewStat","Skew_logP","Skew_logFDR")
  
  res.diff$Skew_logP <- -log10(as.data.frame(res.diff)$Skew_logP)
  res.diff$Skew_logFDR <- -log10(as.data.frame(res.diff)$Skew_logFDR)
  message(paste0("res_diff samples: ", nrow(res.diff)))
  res.diff$ID <- ids_comp
  
  message("combining data")
  res_comp <- merge(out, res.diff, by="ID", all.y=T)
  
  OE_threshold <- -log10(cutoff)

  glm_stats <- res_comp
  glm_stats$p_orig <- 10^(-glm_stats$Skew_logP)
  glm_stats$OE <- ifelse(glm_stats$A_logPadj_BH > OE_threshold | glm_stats$B_logPadj_BH > OE_threshold, T, F)
  
  glm_update <- glm_stats[which(glm_stats$OE),c("ID","p_orig")]
  glm_update$Skew_logFDR_act <- -log10(p.adjust(glm_update$p_orig, method = "BH"))
  res_comp <- merge(res_comp, glm_update[,c("ID","Skew_logFDR_act")], by="ID", all=T)
  
  return(res_comp)
}


### Set up for correlation scatter plot functions.
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
###
panel.lm <- function (x, y,  pch = par("pch"), col.lm = "red",  ...) {
  ymin <- min(y)
  ymax <- max(y)
  xmin <- min(x)
  xmax <- max(x)
  ylim <- c(min(ymin,xmin),max(ymax,xmax))
  xlim <- ylim
  points(x, y, pch = pch,ylim = ylim, xlim= xlim,col=rgb(144,144,144,75,maxColorValue=255),...)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok))
    abline(lm(y[ok]~ x[ok]),
           col = col.lm, ...)
}
###
panel.nlm <- function (x, y,  pch = par("pch"), col.lm = "red",  ...) {
  ymin <- min(y)
  ymax <- max(y)
  xmin <- min(x)
  xmax <- max(x)
  ylim <- c(min(ymin,xmin),max(ymax,xmax))
  xlim <- ylim
  points(x, y, pch = pch,ylim = ylim, xlim= xlim,col=rgb(144,144,144,75,maxColorValue=255),...)
}

### Function to produce scatter plot of counts data for visualization of correlation between two samples
# countsOut       : normalized count data
  ## Can be acquired from tagNorm function (ex: countsOut <- counts(tagNorm(countsData, conditionData, attributesData)))
# sampleX         : string of sample name (column name from counts data)
# sampleY         : string of different sample name (column name from counts data)
# xmax            : maximum x value to include in plot
# ymax            : maximum y value to include in plot
mpraScatter<-function(conditionData, countsOut, sampleX, sampleY,xmax,ymax, plotSave=T, file_prefix) {
  cond_data <- conditionStandard(conditionData)
  count_df<-as.data.frame(countsOut)
  ggplot_output<-ggplot(count_df, aes_string(sampleX,sampleY)) +
    theme_minimal() +
    geom_point(alpha = .3,size=1) +
    xlab(sampleX) + ylab(sampleY) +
    coord_fixed(ratio = 1,xlim=c(0,xmax), ylim=c(0,ymax)) +
    geom_abline(intercept = 0, slope = 1,linetype = 2, size=.75, color=rgb(255,140,0,150,maxColorValue=255))

  if(plotSave==T){
    for(name in row.names(cond_data)){
      if(name == sampleX){
        sample_typeX <- cond_data[name,1]
      }
    }
    for(name in row.names(cond_data)){
      if(name == sampleY){
        sample_typeY <- cond_data[name,1]
      }
    }
    if(sample_typeX == sample_typeY){
      ggsave(paste0("plots/",file_prefix,"_",fileDate(),"_",sample_typeX, "_cor.png"),ggplot_output,units="in",width=4,height=4,device="png")
    }
    if(sample_typeX != sample_typeY){
      ggsave(paste0("plots/",file_prefix,"_",fileDate(),"_",sample_typeX, "_", sample_typeY, "_cor.png"),ggplot_output,units="in",width=4,height=4,device="png")
    }
  }

  return(ggplot_output)
}

### Function to produce plots showing the expression fold change vs. normalized tag counts
# full_output     : Output from dataOut function
# sample          : Cell Type as string
plot_logFC <- function(full_output, sample, negCtrlName="negCtrl", posCtrlName="expCtrl", color_table) {
  message(paste(dim(full_output), collapse = "\t"))
  message(paste(summary(full_output$ctrl_mean), collapse = "\t"))
  message(class(full_output))
  exp_values<-(full_output[full_output$ctrl_mean > 10 & !is.na(full_output$ctrl_mean),])
  message(paste(dim(exp_values), collapse = "\t"))
  exp_values$exp_mean[is.na(exp_values$exp_mean)]<-1
  message(paste(dim(exp_values), collapse = "\t"))
  exp_values$log2FoldChange[is.na(exp_values$log2FoldChange)]<-0
  message(paste(dim(exp_values), collapse = "\t"))
  exp_values$padj[is.na(exp_values$padj)]<-1
  exp_values$sig<-"Not Significant"
  exp_values$sig[exp_values$padj <= 0.00001]<-"Significant (BH adjusted p <= 0.00001)"
  levels(exp_values$sig)<-c("Not Significant", "Significant (BH adjusted p <= 0.00001)")
  exp_values$sig<-factor(exp_values$sig,levels=c("Not Significant", "Significant (BH adjusted p <= 0.00001)"))

  tmp_plotA<-ggplot(exp_values,aes(x=ctrl_mean,y=log2FoldChange,color=sig)) +
    theme_bw() + theme(panel.grid.major = element_line(size = .25,colour = rgb(0,0,0,75,maxColorValue=255)), panel.grid.minor = element_blank()) +
    scale_colour_manual(values=c("Not Significant"=rgb(0,0,0,200,maxColorValue=255),"Significant (BH adjusted p <= 0.00001)"=rgb(55,126,184,255,maxColorValue=255))) +
    geom_point(alpha = .3,size=1) +
    scale_x_log10() +
    #coord_cartesian(xlim = c(10, 1000),ylim = c(-1.5,7.5)) +
    xlab("Normalized Tag Count - Plasmids") + ylab(paste0(sample," Expression Fold Change log2(RNA/Plasmid)")) +
    theme(legend.position = c(.15, .90),
          legend.key = element_blank(),
          legend.background = element_rect(color=rgb(0,0,0,150,maxColorValue=255), fill = "white", size = .5, linetype = "solid")) +
    guides(colour = guide_legend(override.aes = list(size=3,alpha=.7), title=NULL)) +
    geom_abline(intercept = 0, slope = 0,linetype = 1, size=.75, color=rgb(255,140,0,150,maxColorValue=255))

  # message("colors set")
  # cbPalette <- c("#56B4E9","#F84763","#009E73", "#CAA674", "#0072B2", "#D55E00", "#CC79A7","#8057BB","#FBAD12","#999999")
  # cbPalette <- c("#3F47C9","#4274CE","#4F97BB","#64AC99","#7EB976","#9EBE5A","#BEBB48","#D9AE3E","#E69036","#E35F2D","#DB2823")
  # cbPalette <- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F")
  cbPalette <- c("#E6194B","#F58321","#FFE119","#BFEF45","#3CB44B","#42D4F4","#4363D8","#911EB4","#F032E6","#000075","#A9A9A9","#000000")
  
  tmp_plotB<-1
  tmp_plotB<-ggplot(exp_values,aes(x=ctrl_mean,y=log2FoldChange,color=ctrl_exp)) +
    theme_bw() + theme(panel.grid.major = element_line(size = .25,colour = rgb(0,0,0,75,maxColorValue=255)), panel.grid.minor = element_blank()) +
    scale_colour_manual(values=cbPalette) +
    geom_point(alpha = .2,size=1) +
    geom_point(data = subset(exp_values, ctrl_exp == negCtrlName), aes(x=ctrl_mean,y=log2FoldChange),alpha = .8,size=2) +
    geom_point(data = subset(exp_values, ctrl_exp == posCtrlName), aes(x=ctrl_mean,y=log2FoldChange),alpha = .8,size=2) +
    # geom_point(data = subset(exp_values, project == emvarCtrlName), aes(x=ctrl_mean,y=log2FoldChange),alpha = .8,size=2) +
    scale_x_log10() +
    #coord_cartesian(xlim = c(10, 1000),ylim = c(-1.5,7.5)) +
    xlab("Normalized Tag Count - Plasmids") + ylab(paste0(sample," Expression Fold Change log2(RNA/Plasmid)")) +
    theme(legend.position = c(.15, .80),
          legend.key = element_blank(),
          legend.background = element_rect(color=rgb(0,0,0,150,maxColorValue=255), fill = "white", size = .5, linetype = "solid")) +
    guides(colour = guide_legend(override.aes = list(size=3,alpha=.7), title=NULL)) +
    geom_abline(intercept = 0, slope = 0,linetype = 1, size=.75, color=rgb(0,0,0,0,maxColorValue=255))

  return(list(tmp_plotA,tmp_plotB))
}

### Function which runs EVERYTHING
# countsData      : table of tag counts, columns should include: Barcode, Oligo, Sample names
# attributesData  : table of full attributes, columns should include: ID, SNP, Project, Window, Strand, Allele,
  # Haplotype, Bash
# conditionData   : table of conditions, 2 columns no header align to the variable column headers of countsData
  # and celltype
# exclList        : List of cell types to be excluded, default empty list
# filePrefix      : Name of project, it is suggested if normalization methods are being compared to include the normalization method in filePrefix
# plotSave        : Logical, default T indicating that plots will be saved automatically
# altRef          : Logical, default T indicating sorting by alt/ref, if sorting ref/alt set to F
# method          : Method to be used to normalize the data. 4 options - summit shift normalization 'ss', remove the outliers before DESeq normalization 'ro'
  # perform normalization for negative controls only 'nc', median of ratios method used by DESeq 'mn'
MPRAmodel <- function(countsData, attributesData, conditionData, filePrefix, negCtrlName="negCtrl", posCtrlName="expCtrl", projectName="MPRA_PROJ", exclList=c(), plotSave=T, altRef=T, method = 'ss', tTest=T, DEase=T, cSkew=T, correction="BH", cutoff=0.01, upDisp=T, prior=F, raw=T, paired=F, color_table, ...){
  file_prefix <- filePrefix
  # Make sure that the plots and results directories are present in the current directory
  mainDir <- getwd()
  dir.create(file.path(mainDir, "plots"), showWarnings = FALSE)
  dir.create(file.path(mainDir, "results"), showWarnings = FALSE)
  # Resolve any multi-project conflicts, run normalization, and write celltype specific results files
  attributesData <- addHaplo(attributesData, negCtrlName, posCtrlName, projectName)
  message("running DESeq")
  analysis_out <- dataOut(countsData, attributesData, conditionData, altRef=altRef, exclList, file_prefix, method, negCtrlName, tTest, DEase, cSkew, correction, cutoff, upDisp, prior, paired)
  cond_data <- conditionStandard(conditionData)
  n <- length(levels(cond_data$condition))
  full_output <- analysis_out[1:(n-1)]
  dds_results <- analysis_out[[n]]

  # Plot correlation tables using the functions initialized above.
  message("Plotting correlation tables")
  counts_out <- counts(dds_results, normalized=T)
  if(plotSave==F){
    cor_mat_log <- pairs(counts_out,upper.panel=panel.cor,lower.panel=panel.lm,log="xy",pch=16)
    cor_mat1 <- pairs(counts_out,upper.panel=panel.cor,lower.panel=panel.lm,pch=16)
    cor_mat2 <- pairs(counts_out,upper.panel=panel.cor,lower.panel=panel.lm,xlim=c(0,2000),ylim=c(0,2000),pch=16)
  }

  if(plotSave==T){
    png(file=paste0("plots/",file_prefix,"_",fileDate(),"_Cor_mat_log.png"),width=3000,height=3000)
    pairs(counts_out,upper.panel=panel.cor,lower.panel=panel.lm,log="xy",pch=16)
    dev.off()
    png(file=paste0("plots/",file_prefix,"_",fileDate(),"_Cor_mat.png"),width=3000,height=3000)
    pairs(counts_out,upper.panel=panel.cor,lower.panel=panel.lm,pch=16)
    dev.off()
    png(file=paste0("plots/",file_prefix,"_",fileDate(),"_Cor_mat_2.png"),width=3000,height=3000)
    pairs(counts_out,upper.panel=panel.cor,lower.panel=panel.lm,xlim=c(0,2000),ylim=c(0,2000),pch=16)
    dev.off()
  }

  #Prepare for mpraScatter
  message("Plotting Correlation Scatter Plots")
  rep1_loc <- c()
  level_checked <- 0
  level_count <- 0
  for(celltype in levels(cond_data$condition)){
    level_checked <- level_checked + 1
    for(i in 1:dim(cond_data)[1]){
      if(celltype == cond_data[i,1]){
        rep1_loc <- append(rep1_loc, i)
        level_count <- level_count + 1
      }
      if(level_count == level_checked){
        break
      }
    }
  }

  rep1_loc <- grep("_r1", colnames(counts_out))
  replicate_list <- colnames(counts_out)[rep1_loc]
  cell_combinations <- combn(replicate_list,m=2)
  xmax = quantile(counts_out,0.99)
  ymax = quantile(counts_out,0.99)
  for(i in rep1_loc){
    sampleX <- colnames(counts_out)[i]
    sampleY <- colnames(counts_out)[i+1]
    mpraScatter(conditionData = cond_data, countsOut = counts_out, sampleX, sampleY, xmax = xmax, ymax = ymax, plotSave, file_prefix)
  }
  for(combo in dim(cell_combinations)[2]){
    sampleX <- cell_combinations[1,combo]
    sampleY <- cell_combinations[2,combo]
    mpraScatter(conditionData = cond_data, countsOut = counts_out, sampleX, sampleY, xmax = xmax, ymax = ymax, plotSave, file_prefix)
  }

  #Prepare for plot_logFC
  message("Plotting log Fold Change plots")
  for (celltype in levels(cond_data$condition)) {
    if(celltype=="DNA" | celltype %in% exclList ) next
    message(celltype)
    output_tmp<-full_output[[celltype]]
    # message(paste(dim(output_tmp), collapse="\t"))
    plot_list<-plot_logFC(output_tmp, celltype, negCtrlName, posCtrlName, color_table)
    if(plotSave==F){
      plot_list[[1]]
      plot_list[[2]]
    }
    if(plotSave==T){
      ggsave(paste0("plots/",file_prefix,"_",fileDate(),"_logFC_",celltype,".pdf"),plot_list[[1]],units="in",width=8,height=6,device="pdf")
      ggsave(paste0("plots/",file_prefix,"_",fileDate(),"_logFC_",celltype,"_controls.pdf"),plot_list[[2]],units="in",width=8,height=6,device="pdf")
    }
  }
  return(counts_out)
}

