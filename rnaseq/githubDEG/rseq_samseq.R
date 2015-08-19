library(samr);

setwd('C:/Users/slowsmile/Desktop/HT29-RNAseq/Simulation');
folds <- c(2, 0.5, log2(2), log2(0.5));
###################### read in data ############################################;
counts <- read.csv('seq.csv', header = TRUE, row.names = 1, stringsAsFactors = FALSE);
group <- c(rep(1,5), rep(2,5));

# Calculate the normalized gene counts, one can use getLibsizes from BaySeq to calculate size factor
# for each sample
size.factor <- read.table(file = 'C:/Users/slowsmile/Desktop/HT29-RNAseq/Size_factor.txt', 
                          row.names = 1, header = FALSE);
counts.norm <- counts;
for (i in 1:dim(counts)[2]){
  counts.norm[,i] <- counts[,i]/size.factor[i,1];
}

replicates <- group;
grp <- list(rep(1,10), replicates);

####################### Set up Comparison ###########################################
samseq.ac <- SAMseq(counts, group, resp.type = "Two class unpaired", random.seed = 820, 
                    genenames = rownames(counts), nperms = 1000, fdr.output = 0.05);

avsc.up <- samseq.ac$siggenes.table$genes.up[as.numeric(samseq.ac$siggenes.table$genes.up[,4]) > folds[1],];
avsc.dn <- samseq.ac$siggenes.table$genes.lo[as.numeric(samseq.ac$siggenes.table$genes.lo[,4]) < folds[2],];

avsc.up.out <- data.frame(GeneName = avsc.up[,1], CTRL = rowMeans(counts.norm[as.numeric(avsc.up[,2]), 1:5]),
                          CASE = rowMeans(counts.norm[as.numeric(avsc.up[,2]), 6:10]),
                          Fold_Change = as.numeric(avsc.up[,4]), Pvalue = as.numeric(avsc.up[,5])/100);
avsc.dn.out <- data.frame(GeneName = avsc.dn[,1], CTRL = rowMeans(counts.norm[as.numeric(avsc.dn[,2]), 1:5]),
                          CASE = rowMeans(counts.norm[as.numeric(avsc.dn[,2]), 6:10]),
                          Fold_Change = as.numeric(avsc.dn[,4]), Pvalue = as.numeric(avsc.dn[,5])/100);
# OUtput result files 
write.table(avsc.up.out, 'C:/Users/slowsmile/Desktop/HT29-RNAseq/Simulation/Results/rseq_SAMseq_up.txt', 
            row.names = FALSE, sep = '\t');
write.table(avsc.dn.out, 'C:/Users/slowsmile/Desktop/HT29-RNAseq/Simulation/Results/rseq_SAMseq_dn.txt', 
            row.names = FALSE, sep = '\t');
