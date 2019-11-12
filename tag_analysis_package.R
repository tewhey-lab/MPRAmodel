### A package to perform MPRA Tag Analysis

# Required packages
install.packages("ggplot2")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)
library(ggplot2)

###Generate information needed for file outputs later on
fileDate <- function(){
  a <- Sys.Date()
  year <- format(a, "%Y")
  month <- format(a, "%m")
  day <- format(a, "%d")
  today <- paste0(year, month, day)
  return(today)
}

### Add Haplotype column to attribute data if needed
# attributesData  : table of full attributes, columns should include: ID, SNP, Project, Window, Strand, Allele,
  # Haplotype, Bash
addHaplo <- function(attributesData){
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
  return(attributesData)
}

### Isolate oligo data
# countsData      : table of tag counts, columns should include: Barcode, Oligo, Sample names
# returns oligo count data with the oligo as row names and the aggregate count data
oligoIsolate <- function(countsData){
  if("Error" %in% colnames(countsData)){
    countsData <- countsData[,c(1,2,7:dim(countsData)[2])]
  }
  tag_counts <- aggregate(. ~Oligo, data=countsData[,-1], FUN = sum)
  counts_oligo <- tag_counts[,-1]
  rownames(counts_oligo) <- tag_counts[,1]
  return(counts_oligo)
}

### Standardize condition data
# conditionData  : table of conditions, 2 columns no header align to the variable column headers of countsData
  # and celltype
# data should initially be read in such that samples are row names
## Returns: cond_data
conditionStandard <- function(conditionData){
  cond_data <- conditionData
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
# countsData      : table of tag counts, columns should include: Barcode, Oligo, Sample names
# conditionData   : table of conditions, 2 columns no header align to the variable column headers of countsData
  # and celltype
## Returns: dds_results (initial)
processAnalysis <- function(countsData, conditionData){
  # bring in count data and condition data from previously defined functions
  count_data <- oligoIsolate(countsData)
  cond_data <- conditionStandard(conditionData)
  
  # make sure that the column names and row names of the count and condition data match
  colnames(count_data) <- row.names(cond_data)
  
  # perform DESeq analysis
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = cond_data, design = ~condition)
  dds$condition <- relevel(dds$condition, "DNA")
  dds_results <- DESeq(dds, fitType = 'local')
  
  return(dds_results)
}

### Normalize DESeq results and plot normalized densities for each cell type
# countsData      : table of tag counts, columns should include: Barcode, Oligo, Sample names
# attributesData  : table of full attributes, columns should include: ID, SNP, Project, Window, Strand, Allele,
  # Haplotype, Bash
# conditionData  : table of conditions, 2 columns no header align to the variable column headers of countsData
  # and celltype
## Returns: dds_results (normalized)
tagNorm <- function(countsData, conditionData, attributesData, exclList = c(), method = 'ss'){
  dds_results <- processAnalysis(countsData, conditionData)
  dds_results_orig <- processAnalysis(countsData, conditionData)
  attr_data <- addHaplo(attributesData)
  attribute_ids <- (attr_data[attr_data$project == "negCtrl",])$ID
  attribute_names <- colnames(attribute_ids)
  message(attribute_names)
  count_data <- oligoIsolate(countsData)
  cond_data <- conditionStandard(conditionData)
  colnames(count_data) <- row.names(cond_data)
  
  for(celltype in levels(cond_data$condition)){
    if(celltype=="DNA" | celltype %in% exclList) next
    message(celltype)
    # Create temporary results table for each cell type
    temp_outputA <- results(dds_results, contrast = c("condition", celltype, "DNA"))
    # Summit shift normalization
    if(method == "ss"){
      summit <- which.max(density(temp_outputA$log2FoldChange)$y)
      log_offset <- 2^(density(temp_outputA$log2FoldChange)$x[summit])
      sizeFactors(dds_results)[which(cond_data$condition == celltype)] <- sizeFactors(dds_results)[which(cond_data$condition == celltype)]*(log_offset)
    }
    
    if(method == "ro"){
      attribute_ids<-row.names(temp_outputA[!is.na(temp_outputA$pvalue) & temp_outputA$pvalue>0.001,])
    }
  }
  if(method == "ss"){
    dds_results <- estimateDispersions(dds_results,fitType='local')
    dds_results <- nbinomWaldTest(dds_results)
  }
  if(method == "ro" | method == "nc"){
    dds_results_tmp<-estimateSizeFactors(dds_results[attribute_ids])
    sizeFactors(dds_results)<-sizeFactors(dds_results_tmp)
  }
  if(method == "mn"){
    dds_results <- dds_results_orig
  }

  # Plot normalized density for each cell type - 
  for (celltype in levels(cond_data$condition)) {
    if(celltype == "DNA" | celltype %in% exclList) next
    
    temp_outputB <- results(dds_results_orig, contrast=c("condition",celltype,"DNA"))
    outputA <- results(dds_results, contrast=c("condition",celltype,"DNA"))
    
    message("Plotting Normalization Curves")
    pdf(paste0("plots/Normalized_FC_Density_",celltype,".pdf"),width=10,height=10)
    plot(density(temp_outputB[attribute_ids,]$log2FoldChange,na.rm=TRUE),xlim=c(-3,3),col="grey",main=paste0("Normalization - ",celltype))
    lines(density(temp_outputB$log2FoldChange,na.rm=TRUE),xlim=c(-3,3),col="black")
    lines(density(outputA$log2FoldChange,na.rm=TRUE),xlim=c(-3,3),col="red")
    lines(density(outputA[attribute_ids,]$log2FoldChange,na.rm=TRUE),xlim=c(-3,3),col="salmon")
    text(1.5,0.4,adj=c(0,0),labels="All - baseline",col="black")
    text(1.5,0.35,adj=c(0,0),labels="All - corrected",col="red")
    text(1.5,0.3,adj=c(0,0),labels="NegCtrl - baseline",col="grey")
    text(1.5,0.25,adj=c(0,0),labels="NegCtrl - corrected",col="salmon")
    abline(v=0)
    dev.off()
  }
  return(dds_results)
}

### Retrieve output data for future functions - ideally this is the function that will be called.
  # Any subsequent functions should only need an output from here
# countsData      : table of tag counts, columns should include: Barcode, Oligo, Sample names
# attributesData  : table of full attributes, columns should include: ID, SNP, Project, Window, Strand, Allele,
  # Haplotype, Bash
# conditionData  : table of conditions, 2 columns no header align to the variable column headers of countsData
  # and celltype
## writes duplicate output and ttest files for each celltype
dataOut <- function(countsData, attributesData, conditionData, exclList = c(), altRef = T, file_prefix, method = 'ss'){
  dds_results <- tagNorm(countsData, conditionData, attributesData, exclList, method)
  attributesData <- addHaplo(attributesData)
  message("Tags Normalized")
  count_data <- oligoIsolate(countsData)
  message("Oligos isolated")
  cond_data <- conditionStandard(conditionData)
  message("condition data standardized")
  colnames(count_data) <- row.names(cond_data)
  counts_norm <- counts(dds_results, normalized = T)
  
  full_output<-list()
  full_output_var<-list()
  
  for (celltype in levels(cond_data$condition)) {
    if(celltype == "DNA" | celltype %in% exclList) next
    message(celltype)
    
    outputA <- results(dds_results, contrast=c("condition",celltype,"DNA"))
    
    message("Results of dds_results recieved")
    
    condition_table <- as.data.frame(cond_data)
    
    ctrl_cols <- row.names(condition_table)[condition_table$condition=="DNA"]
    exp_cols <- row.names(condition_table)[condition_table$condition==celltype]
    
    message("Control and Experiment Columns Set")
    
    ctrl_mean <- rowMeans(counts_norm[, colnames(counts_norm) %in% ctrl_cols], na.rm = T)
    exp_mean <- rowMeans(counts_norm[, colnames(counts_norm) %in% exp_cols], na.rm = T)
    output_2 <- cbind(ctrl_mean,exp_mean,outputA[,-1])
   
    message("Removing Duplicates")
    dups_output<-expandDups(output_2)
    
    message("Writing Standard Results File")
    full_outputA<-merge(attributesData, as.matrix(dups_output), by.x="ID", by.y="row.names", all.x=TRUE, no.dups=F)
    full_output[[celltype]]<-full_outputA
    write.table(full_outputA, paste0("results/", file_prefix, "_", celltype, "_", fileDate(), ".out"), row.names=F, col.names=T, sep="\t", quote=F)

    message("Writing T-Test Results File")
    outA<-cellSpecificTtest(attributesData, counts_norm, dups_output, ctrl_mean, exp_mean, ctrl_cols, exp_cols, altRef)
    full_output_var[[celltype]]<-outA
    write.table(outA,paste0("results/", file_prefix, "_", celltype, "_emVAR_", fileDate(),".out"), row.names=F, col.names=T, sep="\t", quote=F)    
  }
  return(c(full_output, dds_results))
}

### Expand IDs that denote duplicate oligos in count/DESeq results
# Requires the output_2 of dataOut function. 
  # i.e: dataOut(countsData, attributesData, conditionData)$output_2
## Returns: Final output table with no duplicates
expandDups <- function(output){
  output_orig <- output
  output_new <- cbind(rownames(output), output)
  colnames(output_new)[1] <- "Row.names"
  # Identify duplicates, if any exist
  dups <- output_new[grep("\\(.*\\)$",output_new$Row.names),]
  dups$Row.names <- gsub("^\\((.*)\\)$","\\1",dups$Row.names)
  # Add everything but the duplicates to the final output
  output_final <- output_new[-(grep("\\(.*\\)$",output_new$Row.names)),] 
  output_final <- output_final[!(is.na(output_final$Row.names)),]
  # If there are 1 or more duplicates.... Not clear what this step is doing????
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
cellSpecificTtest<-function(attributesData, counts_norm, dups_output, ctrl_mean, exp_mean, ctrl_cols, exp_cols, altRef = T){
  snp_data <- subset(addHaplo(attributesData),allele=="ref" | allele=="alt")
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
  
  out <- cbind(
    snp_data_expdata_pairs[evens,c(1,2,3,4,5,6,7,9)],
    within(data.frame(
      A_Ctrl_Mean <- snp_data_expdata_pairs[evens, "ctrl_mean"],
      A_Exp_Mean <- snp_data_expdata_pairs[evens, "exp_mean"],
      A_log2FC <- snp_data_expdata_pairs[evens, "log2FoldChange"],
      A_log2FC_SE <- -log10(snp_data_expdata_pairs[evens, "lfcSE"]),
      A_logP <- -log10(snp_data_expdata_pairs[evens, "pvalue"]),
      A_logPadj_BH <- -log10(snp_data_expdata_pairs[evens, "padj"]),    #BF Correction
      A_logPadj_BF <- -log10(snp_data_expdata_pairs[evens, "pvalue"]*(nrow(snp_data_expdata_pairs)/2)),    #BF Correction
      B_Ctrl_Mean <- snp_data_expdata_pairs[odds, "ctrl_mean"],
      B_Exp_Mean <- snp_data_expdata_pairs[odds, "exp_mean"],
      B_log2FC <- snp_data_expdata_pairs[odds, "log2FoldChange"],
      B_log2FC_SE <- -log10(snp_data_expdata_pairs[odds, "lfcSE"]),
      B_logP <- -log10(snp_data_expdata_pairs[odds, "pvalue"]),
      B_logPadj_BH <- -log10(snp_data_expdata_pairs[odds, "padj"]),    #BF Correction
      B_logPadj_BF <- -log10(snp_data_expdata_pairs[odds, "pvalue"]*(nrow(snp_data_expdata_pairs)/2))),{   #BF Correction
        A_logP[is.na(A_logP)] <- 0
        A_logP[A_logP == Inf] <- max(A_logP[is.finite(A_logP)])
        A_logPadj_BH[A_logPadj_BH < 0]<-0
        A_logPadj_BH[A_logPadj_BH == Inf] <- max(A_logPadj_BH[is.finite(A_logPadj_BH)])
        A_logPadj_BF[A_logPadj_BF < 0]<-0
        A_logPadj_BF[A_logPadj_BF == Inf] <- max(A_logPadj_BF[is.finite(A_logPadj_BF)])
        B_logP[is.na(B_logP)] <- 0
        B_logP[B_logP == Inf] <- max(B_logP[is.finite(B_logP)])
        B_logPadj_BH[B_logPadj_BH < 0]<-0
        B_logPadj_BH[B_logPadj_BH == Inf] <- max(B_logPadj_BH[is.finite(B_logPadj_BH)])
        B_logPadj_BF[B_logPadj_BF < 0]<-0
        B_logPadj_BF[B_logPadj_BF == Inf] <- max(B_logPadj_BF[is.finite(B_logPadj_BF)])
      }))
  
  out2 <- out[,c(1:12, 16:19, 23:28)]
  colnames(out2) <- c("ID", "SNP", "chr", "snp_pos", "ref_allele", "alt_allele", "allele", "strand", "A_Ctrl_Mean", "A_Exp_Mean", "A_log2FC", "A_log2FC_SE", "B_Ctrl_Mean", "B_Exp_Mean", "B_log2FC", "B_log2FC_SE",
                      "B_logPadj_BF", "B_logPadj_BH", "B_logP", "A_logPadj_BF", "A_logPadj_BH", "A_logP")
  
  # Don't try to do the t test for ones with all zeros.
  ignore_idx <- which(rowMeans(snp_data_ctdata_pairs[odds,ctrl_cols]) < 10 | rowMeans(snp_data_ctdata_pairs[odds, exp_cols]) < 10 | 
                        is.na(rowMeans(snp_data_ctdata_pairs[odds,ctrl_cols]))  | is.na(rowMeans(snp_data_ctdata_pairs[odds, exp_cols])) | 
                        rowMeans(snp_data_ctdata_pairs[evens,ctrl_cols]) < 10 | rowMeans(snp_data_ctdata_pairs[evens, exp_cols]) < 10 | 
                        is.na(rowMeans(snp_data_ctdata_pairs[evens,ctrl_cols]))  | is.na(rowMeans(snp_data_ctdata_pairs[evens, exp_cols])) )
  
  # For the numerator, set zero values to 1 so that the log-ratio is defined.
  counts1 <- snp_data_ctdata_pairs
  counts1[counts1 == 0] <- 1
  
  # t test
  ratios_A <- log((counts1[evens, exp_cols]) / rowMeans(snp_data_ctdata_pairs[evens, ctrl_cols]))
  ratios_B <- log((counts1[odds, exp_cols]) / rowMeans(snp_data_ctdata_pairs[odds, ctrl_cols]))
  
  ratios_list <- list(ratios_A, ratios_B)

  t_pvalue <- sapply(1:nrow(ratios_A), function(i) if(i %in% ignore_idx){NA} else{
        t.test(as.numeric(ratios_A[i,]), as.numeric(ratios_B[i,]), var.equal=F, paired=T)$p.value})
  t_stat <- sapply(1:nrow(ratios_A), function(i) if(i %in% ignore_idx){NA} else{
        t.test(as.numeric(ratios_A[i,]), as.numeric(ratios_B[i,]), var.equal=F, paired=T)$statistic})
  
  out2$LogSkew <- out2$B_log2FC - out2$A_log2FC
  out2$Skew_logP <- ifelse(is.na(t_pvalue), 0, -log10(t_pvalue))

  OE_threshold <- -log10(.01)
  is_OE <- out2$A_logPadj_BF >= OE_threshold | out2$B_logPadj_BF >= OE_threshold
  out2$Skew_logFDR <- rep(NA, dim(out)[1])
  q_idx <- intersect(which(is_OE), which(!is.na(t_pvalue)))
  out2$Skew_logFDR[q_idx] <- -log10(p.adjust(t_pvalue[q_idx],method="BH"))
  
  return(out2)
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
mpraScatter<-function(conditionData, countsOut, sampleX, sampleY,xmax=quantile(countsOut,.99),ymax=quantile(countsOut,.99), plotSave=T) {
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
      ggsave(paste0("plots/", sample_typeX, "_cor.png"),ggplot_output,units="in",width=4,height=4,device="png")
    }
    if(sample_typeX != sample_typeY){
      ggsave(paste0("plots/", sample_typeX, "_", sample_typeY, "_cor.png"),ggplot_output,units="in",width=4,height=4,device="png")
    }
  }
  
  return(ggplot_output)
}

### Function to produce plots showing the expression fold change vs. normalized tag counts
# full_output     : Output from dataOut function
# sample          : Cell Type as string
plot_logFC <- function(full_output, sample) {
  exp_values<-(full_output[full_output[["ctrl_mean"]] > 10 & !is.na(full_output[["ctrl_mean"]]),])
  exp_values$exp_mean[is.na(exp_values$exp_mean)]<-1
  exp_values$log2FoldChange[is.na(exp_values$log2FoldChange)]<-0
  exp_values$padj[is.na(exp_values$padj)]<-1
  exp_values$sig<-"Not Significant"
  exp_values$sig[exp_values$padj <= 0.00001]<-"Active"
  levels(exp_values$sig)<-c("Not Significant", "Active")
  exp_values$sig<-factor(exp_values$sig,levels=c("Not Significant", "Active"))
  
  tmp_plotA<-ggplot(exp_values,aes(x=ctrl_mean,y=log2FoldChange,color=sig)) +
    theme_bw() + theme(panel.grid.major = element_line(size = .25,colour = rgb(0,0,0,75,maxColorValue=255)), panel.grid.minor = element_blank()) +
    scale_colour_manual(values=c("Not Significant"=rgb(0,0,0,200,maxColorValue=255),"Active"=rgb(55,126,184,255,maxColorValue=255))) +
    geom_point(alpha = .3,size=1) +
    scale_x_log10() +
    #coord_cartesian(xlim = c(10, 1000),ylim = c(-1.5,7.5)) +
    xlab("Normalized Tag Count - Plasmids") + ylab(paste0(sample," Expression Fold Change log2(RNA/Plasmid)")) +
    theme(legend.position = c(.15, .90), 
          legend.key = element_blank(),
          legend.background = element_rect(color=rgb(0,0,0,150,maxColorValue=255), fill = "white", size = .5, linetype = "solid")) +
    guides(colour = guide_legend(override.aes = list(size=3,alpha=.7), title=NULL)) + 
    geom_abline(intercept = 0, slope = 0,linetype = 1, size=.75, color=rgb(255,140,0,150,maxColorValue=255))
  
  cbPalette <- c("#56B4E9", "#999999", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  tmp_plotB<-1
  tmp_plotB<-ggplot(exp_values,aes(x=ctrl_mean,y=log2FoldChange,color=project)) +
    theme_bw() + theme(panel.grid.major = element_line(size = .25,colour = rgb(0,0,0,75,maxColorValue=255)), panel.grid.minor = element_blank()) +
    scale_colour_manual(values=cbPalette) +
    geom_point(alpha = .2,size=1) +
    geom_point(data = subset(exp_values, project == 'negCtrl'), aes(x=ctrl_mean,y=log2FoldChange),alpha = .8,size=2) +
    geom_point(data = subset(exp_values, project == 'expCtrl'), aes(x=ctrl_mean,y=log2FoldChange),alpha = .8,size=2) +
    geom_point(data = subset(exp_values, project == 'emVarCtrl'), aes(x=ctrl_mean,y=log2FoldChange),alpha = .8,size=2) +
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
tagWrapper <- function(countsData, attributesData, conditionData, exclList=c(), filePrefix, plotSave=T, altRef=T, method = 'ss', ...){
  file_prefix <- filePrefix
  analysis_out <- dataOut(countsData, attributesData, conditionData, altRef=altRef, exclList, file_prefix, method)
  cond_data <- conditionStandard(conditionData)
  full_output <- analysis_out[1]
  dds_results <- analysis_out[[2]]
  
  message("Plotting correlation tables")
  counts_out <- counts(dds_results, normalized=T)
  if(plotSave==F){
    cor_mat_log <- pairs(counts_out,upper.panel=panel.cor,lower.panel=panel.lm,log="xy",pch=16)
    cor_mat1 <- pairs(counts_out,upper.panel=panel.cor,lower.panel=panel.lm,pch=16)
    cor_mat2 <- pairs(counts_out,upper.panel=panel.cor,lower.panel=panel.lm,xlim=c(0,2000),ylim=c(0,2000),pch=16)
  }
    
  if(plotSave==T){
    png(file="plots/Cor_mat_log.png",width=3000,height=3000) 
    pairs(counts_out,upper.panel=panel.cor,lower.panel=panel.lm,log="xy",pch=16)
    dev.off()
    png(file="plots/Cor_mat.png",width=3000,height=3000) 
    pairs(counts_out,upper.panel=panel.cor,lower.panel=panel.lm,pch=16)
    dev.off()
    png(file="plots/Cor_mat_2.png",width=3000,height=3000) 
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
  for(i in rep1_loc){
    sampleX <- colnames(counts_out)[i]
    sampleY <- colnames(counts_out)[i+1]
    mpraScatter(conditionData = cond_data, countsOut = counts_out, sampleX, sampleY, plotSave)
  }
  for(combo in dim(cell_combinations)[2]){
    sampleX <- cell_combinations[1,combo]
    sampleY <- cell_combinations[2,combo]
    mpraScatter(conditionData = cond_data, countsOut = counts_out, sampleX, sampleY, plotSave)
  }
  
  #Prepare for plot_logFC
  message("Plotting log Fold Change plots")
  for (celltype in levels(cond_data$condition)) {
    if(celltype=="DNA" | celltype %in% exclList ) next
    message(celltype)
    output_tmp<-full_output[[celltype]]
    message(class(output_tmp))
    plot_list<-plot_logFC(output_tmp, celltype)
    if(plotSave==F){
      plot_list[[1]]
      plot_list[[2]]
    }
    if(plotSave==T){
      ggsave(paste0("plots/logFC_",celltype,".pdf"),plot_list[[1]],units="in",width=8,height=6,device="pdf")
      ggsave(paste0("plots/logFC_",celltype,"_controls.pdf"),plot_list[[2]],units="in",width=8,height=6,device="pdf")
    }
  }
}



