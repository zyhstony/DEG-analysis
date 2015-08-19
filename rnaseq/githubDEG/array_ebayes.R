library(limma)
library(hgu133plus2.db)

folds <- c(2, 0.5, log2(2), log2(0.5));
############################ Read in microarray data ############################
raw.dat <- read.csv('C:/Users/slowsmile/Desktop/HT29-RNAseq/Simulation/array.csv', 
                    header = TRUE, stringsAsFactors = FALSE);
dat.log <- log2(data + 0.01);

cl<- c(rep(1,5), rep(0,5));
########### quantile normalization ###########
avsc.norm <- normalizeBetweenArrays(as.matrix(dat.log), method = "quantile");

##################################################################
############### Comparing between CTRL and CASE ##################
##################################################################
design <- as.matrix(cbind(rep(1,10), cl));
rownames(design) <- colnames(avsc.norm);
colnames(design) <- c('control', 'treat-control');

fit <- lmFit(avsc.norm, design);
avsc.mean <- cbind(rowMeans(avsc.norm[,1:5]), rowMeans(avsc.norm[,6:10]));
fitted <- eBayes(fit);
sig.genes <- topTable(fitted, coef = 'treat-control', number = dim(fitted), adjust="BH", p.value = 0.05);

up.reg <- sig.genes[sig.genes$logFC > folds[3],];
dn.reg <- sig.genes[sig.genes$logFC < folds[4],];

up.reg.out <- data.frame(GeneNames = up.reg$ID, CTRL = 2^avsc.mean[as.numeric(rownames(up.reg)),1],
                         CASE = 2^avsc.mean[as.numeric(rownames(up.reg)),2],
                         Fold_Change = 2^up.reg$logFC, Pvalue = up.reg$adj.P.Val);
dn.reg.out <- data.frame(GeneNames = dn.reg$ID, CTRL = 2^avsc.mean[as.numeric(rownames(dn.reg)),1],
                         CASE = 2^avsc.mean[as.numeric(rownames(dn.reg)),2],
                         Fold_Change = 2^dn.reg$logFC, Pvalue = dn.reg$adj.P.Val);

######################### OUtput result files ################################
write.table(up.reg.out, 'C:/Users/slowsmile/Desktop/HT29-RNAseq/Simulation/Results/array_ebayes_up.txt', 
            row.names = FALSE, sep = '\t');
write.table(dn.reg.out, 'C:/Users/slowsmile/Desktop/HT29-RNAseq/Simulation/Results/array_ebayes_dn.txt', 
            row.names = FALSE, sep = '\t');