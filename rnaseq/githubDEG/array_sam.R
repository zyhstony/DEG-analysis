library(hgu133plus2.db);
library(limma);
library(affy);
library(siggenes);

folds <- c(2, 0.5, 1, -1);

############################ Read in microarray data ############################
raw.dat <- read.csv('C:/Users/slowsmile/Desktop/HT29-RNAseq/Simulation/array.csv', 
                    header = TRUE, stringsAsFactors = FALSE);
dat.log <- log2(data + 0.01);
group <- c(rep(1,5), rep(0,5));
########### quantile normalization ###########
avsc.norm <- normalizeBetweenArrays(as.matrix(dat.log), method = "quantile");

########### SAM analysis ######################
best_delta <- function(sam.out, step = 0.01, alpha = 0.05) {
  warning("The input must be a object produced by sam function");
  deltas <- sam.out@mat.fdr[,1];
  sequence <- rev(seq(min(deltas), max(deltas), by = step));
  for (j in 1:length(sequence)){
    buffer <- summary(sam.out, sequence[j]);
    buf.fdr <- buffer@mat.fdr[5];
    if (buf.fdr < alpha) {
      delta <- sequence[j];
    } else {break};
  }  
  return(delta);	
}


sam.avsc <- sam(avsc.norm, group, rand = 820);
glist.avsc <- summary(sam.avsc, best_delta(sam.avsc));

sam.avsc.down <- glist.avsc@mat.sig[glist.avsc@mat.sig$R.fold < folds[2],];
sam.avsc.up <- glist.avsc@mat.sig[glist.avsc@mat.sig$R.fold > folds[1],];

avsc.down <- data.frame(GeneName=rownames(sam.avsc.down), CTRL = 2^rowMeans(avsc.input[sam.avsc.down$Row, 1:5]), 
                        CASE = 2^rowMeans(avsc.input[sam.avsc.down$Row, 6:10]), 
                        Fold_change=sam.avsc.down$R.fold, P_value=sam.avsc.down$q.value);

avsc.up <- data.frame(GeneName=rownames(sam.avsc.up), CTRL = 2^rowMeans(avsc.input[sam.avsc.up$Row, 1:5]), 
                      CASE = 2^rowMeans(avsc.input[sam.avsc.up$Row, 6:10]),
                      Fold_change=sam.avsc.up$R.fold, P_value=sam.avsc.up$q.value);

########## write output to files #############
write.table(avsc.down, 'C:/Users/slowsmile/Desktop/HT29-RNAseq/Simulation/Results/array_sam_dn.txt', row.names = FALSE, sep = '\t');
write.table(avsc.up, 'C:/Users/slowsmile/Desktop/HT29-RNAseq/Simulation/Results/array_sam_up.txt', row.names = FALSE, sep = '\t');