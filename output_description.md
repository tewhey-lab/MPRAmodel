# MPRAmodel Output File Descriptions:
 - `<proj>_<date>_counts.out`
      - Oligo level counts table
      - Sum of all individual barcode counts by oligo
 - `<proj>_<date>_<cell>_normalized_counts.out`
      - Normalized counts data (default summit shift)
      - Produced when running `dataOut` function
 - `<proj>_<cell>_<date>.out`
      + <details>
          <summary>Standard output including the following columns (* are from the provided attributes table): </summary>
          <ul>
               <li> ID*  </li>
               <li> SNP* </li>
               <li> chr* </li>
               <li> pos* </li>
               <li> ref_allele* </li>
               <li> alt_allele* </li>
               <li> allele* </li>
               <li> window* </li>
               <li> strand* </li>
               <li> project* </li>
               <li> haplotype* </li>
               <li> ctrl_exp -- minimized project list (negative controls, positive controls, other; used for plotting) </li>
               <li> DNA_mean -- mean of oligo level plasmid counts (pre-normalization) </li>
               <li> ctrl_mean -- mean of normalized plasmid counts </li>
               <li> exp_mean -- mean of normalized RNA counts </li>
               <li> log2FoldChange -- Oligo activity as log2(RNA/DNA) </li>
               <li> lfcSE -- Standard error of log2FoldChange </li>
               <li> stat -- Wald statistic of activity </li>
               <li> pvalue -- unadjusted pvalue of oligo activity </li>
               <li> padj -- BH adjusted pvalue of oligo activity  </li>
        </ul>
      </details>
 - `<proj>_<cell>_emVAR_ttest_<date>.out`
      - Expression-Modulating Variants (emVAR) calculated by tTest
      + <details>
          <summary>Columns include: </summary>
          <ul>
            <li> comb -- <code> SNP_window_strand_haplotype </code> </li>
                <li> ID </li>
                <li> SNP </li>
                <li> chr </li>
                <li> pos </li>
                <li> ref_allele </li>
                <li> alt_allele </li>
                <li> allele </li>
                <li> window </li>
                <li> strand </li>
                <li> haplotype </li>
                <li> A_Ctrl_Mean -- mean of normalized plasmid counts for the reference allele </li>
                <li> A_Exp_Mean -- mean of normalized RNA counts for the reference allele </li>
                <li> A_log2FC -- Reference allele activity as log2(RNA/DNA) </li>
                <li> A_log2FC_SE -- standard error of the log2FC for the reference allele </li>
                <li> A_logP -- -log10(pvalue) for the reference allele activity </li>
                <li> A_logPadj_BH -- -log10 Benjamini Hochberg adjusted pvalue for the reference allele activity </li>
                <li> A_logPadj_BF -- -log10 Bonferroni adjusted pvalue for the reference allele activity </li>
                <li> B_Ctrl_Mean -- mean of normalized plasmid counts for the alternate allele </li>
                <li> B_Exp_Mean -- mean of normalized RNA counts for the alternate allele </li>
                <li> B_log2FC -- Alternate allele activity as log2(RNA/DNA) </li>
                <li> B_log2FC_SE -- standard error of the log2FC for the alternate allele </li>
                <li> B_logP -- -log10(pvalue) for the alternate allele activity </li>
                <li> B_logPadj_BH -- -log10 Benjamini Hochberg adjusted pvalue for the alternate allele activity </li>
                <li> B_logPadj_BF -- -log10 Bonferroni adjusted pvalue for the alternate allele activity </li>
                <li> Log2Skew -- Allelic Skew `B_log2FC - A_log2FC` </li>
                <li> LogSkew_SE -- `sqrt(A_log2FC_SE^2 + B_log2FC_SE^2)` </li>
                <li> Skew_logP -- -log10(pvalue) from the t-test </li>
                <li> Skew_logFDR -- -log10 Benjamini Hochberg adjusted pvalue from the t-test </li>
                <li> Skew_logFDR_act -- -log10 Benjamini Hochberg adjusted pvalue of the active sequences (A or B must have significance) </li>
        </ul>
      </details>
          <ul> <li>Produced when running <code>dataOut</code> function if <code>tTest=T</code> is used. </li> </ul>
 - `<proj>_<cell>_emVAR_glm_[paired]_<date>.out`
      - emVAR calculated by DESeq2 glm
      + <details>
            <summary>Columns include: </summary>
            <ul>
                <li> ID </li>
                <li> comb -- <code> SNP_window_strand_haplotype </code> </li>
                <li> SNP </li>
                <li> chr </li>
                <li> pos </li>
                <li> ref_allele </li>
                <li> alt_allele </li>
                <li> allele </li>
                <li> window </li>
                <li> strand </li>
                <li> haplotype </li>
                <li> A_Ctrl_Mean -- mean of normalized plasmid counts for the reference allele </li>
                <li> A_Exp_Mean -- mean of normalized RNA counts for the reference allele </li>
                <li> A_log2FC -- Reference allele activity as log2(RNA/DNA)  </li>
                <li> A_log2FC_SE -- standard error of the log2FC for the reference allele </li>
                <li> A_logP -- -log10(pvalue) for the reference allele activity </li>
                <li> A_logPadj_BH -- -log10 Benjamini Hochberg adjusted pvalue for the reference allele activity </li>
                <li> A_logPadj_BF -- -log10 Bonferroni adjusted pvalue for the reference allele activity </li>
                <li> B_Ctrl_Mean -- mean of normalized plasmid counts for the alternate allele </li>
                <li> B_Exp_Mean -- mean of normalized RNA counts for the alternate allele </li>
                <li> B_log2FC -- Alternate allele activity as log2(RNA/DNA)  </li>
                <li> B_log2FC_SE -- standard error of the log2FC for the alternate allele </li>
                <li> B_logP -- -log10(pvalue) for the alternate allele activity </li>
                <li> B_logPadj_BH -- -log10 Benjamini Hochberg adjusted pvalue for the alternate allele activity </li>
                <li> B_logPadj_BF -- -log10 Bonferroni adjusted pvalue for the alternate allele activity </li>
                <li> Log2Skew -- Allelic Skew calculated by the DESeq2 GLM </li>
                <li> LogSkew_SE -- Standard Error of allelic skew, calculated by the GLM </li>
                <li> skewStat -- Wald statistic of allelic skew </li>
                <li> Skew_logP -- -log10(pvalue) from the Wald test </li>
                <li> Skew_logFDR -- -log10 Benjamini Hochberg adjusted pvalue from the Wald test </li>
                <li> Skew_logFDR_act -- -log10 Benjamini Hochberg adjusted pvalue of the active sequences (A or B must have significance) </li>
        </ul>
      </details>
          <ul> <li>Produced when running <code>dataOut</code> function if <code>DEase=T</code> is used. </li> </ul>
