# Array data t-test test pipeline
library(limma);

folds <- c(2, 0.5, 1, -1);
############################ Read in microarray data ############################
Path <- 'C:/Users/slowsmile/Desktop/HT29-RNAseq/Simulation/';
Out <- 'C:/Users/slowsmile/Desktop/HT29-RNAseq/Simulation/Results/';
data <- read.csv(paste(Path, 'array.csv', sep = ''), header = TRUE, stringsAsFactors = FALSE, row.name = 1);
dat.log <- log2(data + 0.01);

########### quantile normalization ###########
avsc.norm <- normalizeBetweenArrays(as.matrix(dat.log), method = "quantile");

################################################################
########### t test analysis for A vs C #########################
################################################################
p.cutoff <- 0.05;
input <- avsc.norm;  
t.rec <- as.numeric();
for (i in 1:dim(input)[1]){         
  buffer <- input[i,];
  t.out <- t.test(buffer[1:5], buffer[6:10]);
  t.buf <- c(t.out$p.value, mean(buffer[1:5]), mean(buffer[6:10]), mean(buffer[1:5]) - mean(buffer[6:10]));
  t.rec <- rbind(t.rec, t.buf);
}

colnames(t.rec) <- c('p-value', 'CTRL', 'CASE', 'log2(FC)');
rownames(t.rec) <- rownames(input);

#Multiple testing correction
t.rec[,1] <- p.adjust(t.rec[,1], method = 'BH');
t.rec <- data.frame(t.rec);

# Select significant genes
up.reg <- t.rec[t.rec[,1] <= p.cutoff & t.rec[,4] > folds[3],];
dn.reg <- t.rec[t.rec[,1] <= p.cutoff & t.rec[,4] < folds[4],];

if (dim(up.reg)[1] > 0){
  up.reg.out <- data.frame(GeneNames = rownames(up.reg), CTRL = 2^up.reg[,2], CASE = 2^up.reg[,3], Fold_Change = 2^up.reg[,4], Pvalue = up.reg[,1]);
} else {
  up.reg.out <- data.frame(GeneNames = NULL, CTRL = NULL, CASE = NULL, Fold_Change = NULL, Pvalue = NULL);
}

if (dim(dn.reg)[1] > 0){
  dn.reg.out <- data.frame(GeneNames = rownames(dn.reg), CTRL = 2^dn.reg[,2], CASE = 2^dn.reg[,3], Fold_Change = 2^dn.reg[,4], Pvalue = dn.reg[,1]);
} else {
  dn.reg.out <- data.frame(GeneNames = NULL, CTRL = NULL, CASE = NULL, Fold_Change = NULL, Pvalue = NULL);
}
######################### OUtput result files ################################
write.table(up.reg.out, paste(Out,'array_ttest_up.txt', sep = ''), 
            row.names = FALSE, sep = '\t');
write.table(dn.reg.out, paste(Out, 'array_ttest_dn.txt', sep = ''), 
            row.names = FALSE, sep = '\t');