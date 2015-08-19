library(DESeq);

setwd('C:/Users/slowsmile/Desktop/HT29-RNAseq/Simulation');
folds <- c(2, 0.5, log2(2), log2(0.5));
###################### read in data ############################################;
counts <- read.csv('seq.csv', header = TRUE, row.names = 1, stringsAsFactors = FALSE);
group <- c(rep(0,5), rep(1,5));

########################################################################
################### Set up design for Comparison #######################
########################################################################
avsd.design <- data.frame(
  row.names = colnames(counts),
  condition = rep(c('control','case'),each = 5),
  libType = rep('pair-end',5));

avsd.conds <- factor(avsd.design$condition);
avsd.cds <- newCountDataSet(counts.avsd, avsd.conds);
avsd.cds <- estimateSizeFactors(avsd.cds);
sizeFactors(avsd.cds);
avsd.cds <- estimateDispersions(avsd.cds);
avsd.res <- nbinomTest(avsd.cds, "control", "case");

ad.up <- avsd.res[avsd.res$foldChange >= folds[1] & avsd.res$padj <= 0.05, ];
ad.dn <- avsd.res[avsd.res$foldChange <= folds[2] & avsd.res$padj <= 0.05, ];

ad.up.out <- data.frame(GeneNames = ad.up$id, CTRL = ad.up$baseMeanB, CASE = ad.up$baseMeanA, Fold_Change = ad.up$foldChange, Pvalue = ad.up$padj);
ad.dn.out <- data.frame(GeneNames = ad.dn$id, CTRL = ad.dn$baseMeanB, CASE = ad.dn$baseMeanA, Fold_Change = ad.dn$foldChange, Pvalue = ad.dn$padj);

######################### OUtput result files ################################
write.table(ad.up.out, 'C:/Users/slowsmile/Desktop/HT29-RNAseq/Simulation/Results/rseq_Deseq_up.txt', 
            row.names = FALSE, sep = '\t');
write.table(ad.dn.out, 'C:/Users/slowsmile/Desktop/HT29-RNAseq/Simulation/Results/rseq_Deseq_dn.txt', 
            row.names = FALSE, sep = '\t');
