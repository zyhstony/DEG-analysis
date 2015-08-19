library(baySeq);

setwd('C:/Users/slowsmile/Desktop/HT29-RNAseq/Simulation');
folds <- c(2, 0.5, log2(2), log2(0.5));
###################### read in data ############################################;
counts <- read.csv('seq.csv', header = TRUE, row.names = 1, stringsAsFactors = FALSE);

replicates <- c(1,1,1,1,1,2,2,2,2,2);
grp <- list(rep(1,10), replicates);

############################################################################
#################### Create DGE object for Comparison ######################
############################################################################
#avsc.lib <- lib.size[c(1:3,7:9),];
avsc.CD <- new("countData", data = counts, replicates = replicates, groups = grp);
avsc.CD@libsizes <- getLibsizes(avsc.CD);
avsc.CD@annotation <- data.frame(name = rownames(counts.avsc));

cl <- NULL;
avsc.CDP.NB <- getPriors.NB(avsc.CD, samplesize = 1000, takemean = TRUE, iterations = 1000, cl = cl);
avsc.CDP.NB@priors;

avsc.CDPost.NB <- getLikelihoods.NB(avsc.CDP.NB, pET = "BIC", bootStraps = 3, cl = cl);
avsc.CDPost.NB@estProps;
#avsc.CDPost.NB@posteriors[1:10, ];

selTags <- order(avsc.CDPost.NB@posteriors[, 2], decreasing = TRUE);
avsc.sig.genes <- topCounts(avsc.CDPost.NB, group = 2, FDR = 0.05, normaliseData = TRUE);

CTRL <- rowMeans(avsc.sig.genes[,2:6]);
CASE <- rowMeans(avsc.sig.genes[,7:11]);
fc <- CASE/CTRL;

avsc.up <- avsc.sig.genes[fc >= folds[1],];
avsc.dn <- avsc.sig.genes[fc <= folds[2],];

avsc.up.out <- data.frame(GeneName = avsc.up$name, CTRL = CTRL[fc >= folds[1]], 
                          CASE = CASE[fc >= folds[1]], Fold_Change = fc[fc >= folds[1]], Pvalue = avsc.up$FDR);
avsc.dn.out <- data.frame(GeneName = avsc.dn$name, CTRL = CTRL[fc <= folds[2]], 
                          CASE = CASE[fc <= folds[2]], Fold_Change = fc[fc <= folds[2]], Pvalue = avsc.dn$FDR);

############################## OUtput result files  ############################################
write.table(avsc.up.out, 'C:/Users/slowsmile/Desktop/HT29-RNAseq/Simulation/Results/rseq_Bayseq_up.txt', 
            row.names = FALSE, sep = '\t');
write.table(avsc.dn.out, 'C:/Users/slowsmile/Desktop/HT29-RNAseq/Simulation/Results/rseq_Bayseq_dn.txt', 
            row.names = FALSE, sep = '\t');
