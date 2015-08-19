########simulate data and generate DEG table for the algorithms as described in Xu, X., Zhang, Y., Williams, J., Antoniou, E., McCombie, W. R., Wu, S., ... & Li, E. (2013). Parallel comparison of Illumina RNA-Seq and Affymetrix microarray platforms on transcriptomic profiles generated from 5-aza-deoxy-cytidine treated HT-29 colon cancer cells and simulated datasets. BMC bioinformatics, 14(Suppl 9), S1.
library(limma)
library(hgu133plus2.db)
library(affy);
library(siggenes);
library(baySeq);
library(DESeq);
library(samr);
source('generate.r')
names<-c("table.ttest","table.eb","table.sam","table.bayseq","table.DESeq","table.SAMseq")
table.all<-array(0,c(6,3,10))
for(eff in 1:6){
for(run in 1:10){
de=0.5*eff  ####### mean fold change
all<-gen(i=2,k=5,j=10000,percent=0.1,de=de) #generate data detail see generate.r


###################### Microarray ############################################;

array.raw<-cbind(all[[2]][1,,],all[[2]][2,,])
seq.raw<-cbind(all[[1]][1,,],all[[1]][2,,])
DE<-all[[4]]
dv<-rep(0,10)


write.table(array.raw,paste("c:/rnaseq/simulation/array.raw.csv",eff,run,sep=""),sep=",")
write.table(seq.raw,paste("c:/rnaseq/simulation/seq.raw.csv",eff,run,sep=""),sep=",")
array<-array.raw
log2.diff<-log2(rowMeans(array[,6:10])/rowMeans(array[,1:5]))
#ttest
tpv<-rep(1,10000)
for(i in 1:10000)
{
if(sum(array.raw[i,])!=0){
tpv[i]<-t.test(array[i,1:5],array[i,6:10])$p.value }
}

table.ttest<-cbind(log2.diff,tpv,p.adjust(tpv,method="BH"),DE)
colnames(table.ttest)<-c("log2.diff","pv","FDR","DE")
write.table(table.ttest,paste("c:/rnaseq/simulation/ttest.csv",eff,run,sep=""),sep=",")

#ebayes
cl<- c(rep(0,5), rep(1,5));
design <- as.matrix(cbind(rep(1,10), cl));
rownames(design) <- colnames(array);
colnames(design) <- c('control', 'treat-control');
fit <- lmFit(array, design);
fitted <- eBayes(fit);
table.eb<-cbind(log2.diff,fitted$p.value[,2],p.adjust(fitted$p.value[,2],method="BH"),DE)
colnames(table.eb)<-c("log2.diff","pv","FDR","DE")
write.table(table.eb,paste("c:/rnaseq/simulation/eBayes.csv",eff,run,sep=""),sep=",")

#SAM
group <- c(rep(0,5), rep(1,5));
samfit <- sam(array, group, rand = 820);
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
delta<-best_delta(samfit)
index<-summary(samfit,delta)@row.sig.genes
flag<-rep(0,10000)
flag[index]<-1
table.sam<-cbind(log2.diff,samfit@q.value,flag,DE)
colnames(table.sam)<-c("log2.diff","Local.FDR","significant","DE")
write.table(table.sam,paste("c:/rnaseq/simulation/sam.csv",eff,run,sep=""),sep=",")



###################### RNA_SEQ ############################################;
counts <-seq.raw

replicates <- c(1,1,1,1,1,2,2,2,2,2);
grp <- list(rep(1,10), replicates);


#bayseq
avsc.CD <- new("countData", data = counts, replicates = replicates, groups = grp);
avsc.CD@libsizes <- getLibsizes(avsc.CD);
avsc.CD@annotation <- data.frame(name = c(1:10000));
cl <- NULL;
avsc.CDP.NB <- getPriors.NB(avsc.CD, samplesize = 500, takemean = TRUE, iterations = 500, cl = cl);
avsc.CDPost.NB <- getLikelihoods.NB(avsc.CDP.NB, pET = "BIC", bootStraps = 3, cl = cl);
avsc.sig.genes <- topCounts(avsc.CDPost.NB, group = 2, number = 10000, normaliseData = TRUE);
fdr<-avsc.sig.genes[,ncol(avsc.sig.genes)]
flag<-rep(0,10000)
flag[fdr<0.05]<-1
order=order(avsc.sig.genes$name)
flag<-flag[order]
log2.diff<-log2(rowMeans(avsc.sig.genes[,7:11])/rowMeans(avsc.sig.genes[,2:6]))
log2.diff<-log2.diff[order]
table.bayseq<-cbind(log2.diff,avsc.sig.genes$Likelihood[order],flag,DE)
colnames(table.bayseq)<-c("log2.diff","posterior","significant","DE")
write.table(table.bayseq,paste("c:/rnaseq/simulation/bayseq.csv",eff,run,sep=""),sep=",")

#DESeq
group <- c(rep(0,5), rep(1,5))
avsd.design <- data.frame(
  row.names = c(1:10),
  condition = rep(c('control','case'),each = 5),
  libType = rep('pair-end',5));
avsd.conds <- factor(avsd.design$condition);
avsd.cds <- newCountDataSet(counts, avsd.conds);
avsd.cds <- estimateSizeFactors(avsd.cds);
sizeFactors(avsd.cds);
avsd.cds <- estimateDispersions(avsd.cds);
avsd.res <- nbinomTest(avsd.cds, "control", "case");
table.DESeq<-cbind(avsd.res$log2FoldChange,avsd.res$pval,avsd.res$padj,DE)
colnames(table.DESeq)<-c("log2.diff","pvalue","FDR","DE")
write.table(table.DESeq,paste("c:/rnaseq/simulation/DESeq.csv",eff,run,sep=""),sep=",")

#SAMseq
group <- c(rep(1,5), rep(2,5));
samseq.ac <- SAMseq(counts, group, resp.type = "Two class unpaired", random.seed = 820, geneid=c(1:10000),nperms = 500,fdr.output =0.05);
table<-rbind(samseq.ac$siggenes.table$genes.up,samseq.ac$siggenes.table$genes.lo)
table<-table[,c(2,4)]
name<-as.numeric(table[,1])
flag<-rep(0,10000)
flag[name]<-1
log2.diff<-log2(samseq.ac$samr.obj$foldchange)
table.SAMseq<-cbind(log2.diff,flag,DE)
colnames(table.SAMseq)<-c("log2.diff","significant","DE")
write.table(table.SAMseq,paste("c:/rnaseq/simulation/SAMseq.csv",eff,run,sep=""),sep=",")

#result
temp<-getresult(get(names[1]))
for(i in 2:6)
{
temp<-rbind(temp,getresult(get(names[i])))
}

table.all[,,run]<-temp

}

###generate output for summarizing algorithm comparison given different mean fold change
name<-paste("res",de,sep="")
assign(name,array(0,c(6,3)))
temp<-get(name)

for(i in 1:6)
for(j in 1:3)

{ {
temp[i,j]<-mean(table.all[i,j,])
}  }

assign(name,temp)
} 
