#############defined function#####################
info<-function(index1,index2){ 
#pairwise comparison of DEG result with DEG table provided
#index1 and index2 are the code of simulation dataset
method=c("ttest","ebayes","sam","bayseq","DESeq","SAMseq")
name<-paste("c:/rnaseq/simulation/",method,".csv",index1,index2,sep="")

for(i in 1:6){assign(paste("r",i,sep=""),do.call(read.table,list(name[i],sep=",",header=T)))}
true.gene<-getnames(r1)[[1]]
false.gene<-getnames(r1)[[2]]
all.gene<-getnames(r1)[[3]]
for(i in 2:6){temp<-get(paste("r",i,sep=""));
true.gene<-cbind(true.gene,getnames(temp)[[1]])
false.gene<-cbind(false.gene,getnames(temp)[[2]])
all.gene<-cbind(all.gene,getnames(temp)[[3]])}
rownames(true.gene)<-c(1:10000)
colnames(true.gene)<-method
rownames(false.gene)<-c(1:10000)
colnames(false.gene)<-method
rownames(all.gene)<-c(1:10000)
colnames(all.gene)<-method

name1<-paste("c:/rnaseq/simulation/","array.raw",".csv",index1,index2,sep="")
name2<-paste("c:/rnaseq/simulation/","seq.raw",".csv",index1,index2,sep="")
assign(paste("array",index1,index2,sep=""),do.call(read.table,list(name1,sep=",",header=T)))
assign(paste("seq",index1,index2,sep=""),do.call(read.table,list(name2,sep=",",header=T)))
temp.array<-get(paste("array",index1,index2,sep=""))
temp.seq<-get(paste("seq",index1,index2,sep=""))
temp.all<-cbind(temp.array[,6:10],temp.seq[,6:10],all.gene)
return(temp.all)
}
}

getnames<-function(table){
###format DEG table to compute sensitivity 
flag<-rep(0,10000)
flag.t<-rep(0,10000)
flag.f<-rep(0,10000)
if(colnames(table)[ncol(table)-1]=="significant"){
flag[table[,ncol(table)-1]==1 & abs(table[,1])>=1]<-1
flag.t<-flag
flag.t[table[,ncol(table)]==0]<-0
flag.f<-flag
flag.f[table[,ncol(table)]!=0]<-0
}
if(colnames(table)[ncol(table)-1]=="FDR"){
flag[table[,ncol(table)-1]<0.05 & abs(table[,1])>=1]<-1
flag.t<-flag
flag.t[table[,ncol(table)]==0]<-0
flag.f<-flag
flag.f[table[,ncol(table)]!=0]<-0
}

return(list(flag.t,flag.f,flag))
}

######################################################################################################################

###############Compare DEG Table on different simulation runs###########################
for(i in 1:8 )
{
for(j in 1:10)
{
assign(paste("info",i,j,sep=""),info(i,j))
}
}

###################characterize genes that methods agree on and disagree on in pariwise############################################
for(level in 1:8)
{
assign(paste("matrix",level,sep=""),matrix(0,9,18))
now<-get(paste("matrix",level,sep=""))
for(row in 1:3)
{
for(col in 1:3)
{
assign(paste("a.array",row,col,sep=""),NULL)
name.a<-paste("a.array",row,col,sep="")
assign(paste("b.array",row,col,sep=""),NULL)
name.b<-paste("b.array",row,col,sep="")
assign(paste("c.array",row,col,sep=""),NULL)
name.c<-paste("c.array",row,col,sep="")
assign(paste("a.seq",row,col,sep=""),NULL)
name.as<-paste("a.seq",row,col,sep="")
assign(paste("b.seq",row,col,sep=""),NULL)
name.bs<-paste("b.seq",row,col,sep="")
assign(paste("c.seq",row,col,sep=""),NULL)
name.cs<-paste("c.seq",row,col,sep="")
for(rep in 1:10)
{
name.info<-paste("info",level,rep,sep="")
temp.array<-as.matrix(get(name.info)[get(name.info)[,10+row]==1 & get(name.info)[,13+col]==0,1:10] )
temp.seq<-as.matrix(get(name.info)[get(name.info)[,13+col]==1 & get(name.info)[,10+row]==0,1:10] )
temp.all<-as.matrix(get(name.info)[get(name.info)[,13+col]*get(name.info)[,10+row]==1,1:10])
assign(name.a,c(get(name.a),as.vector(temp.array[,1:5])))
assign(name.b,c(get(name.b),as.vector(temp.seq[,1:5])))
assign(name.c,c(get(name.c),as.vector(temp.all[,1:5])))
assign(name.as,c(get(name.as),as.vector(temp.array[,6:10])))
assign(name.bs,c(get(name.bs),as.vector(temp.seq[,6:10])))
assign(name.cs,c(get(name.cs),as.vector(temp.all[,6:10])))
}
now[(row-1)*3+col,]<-c(mean(get(name.a)),median(get(name.a)),sqrt(var(get(name.a))),mean(get(name.as)),median(get(name.as)),sqrt(var(get(name.as))),mean(get(name.b)),median(get(name.b)),sqrt(var(get(name.b))),mean(get(name.bs)),median(get(name.bs)),sqrt(var(get(name.bs))),mean(get(name.c)),median(get(name.c)),sqrt(var(get(name.c))),mean(get(name.cs)),median(get(name.cs)),sqrt(var(get(name.cs))))
}}
assign(paste("matrix",level,sep=""),now)
} #gather the information

method1=c("ttest","ebayes","sam")
method2=c("bayseq","DESeq","SAMseq")
name.method<-c(do.call(paste,list(method1[1],"*",method2,sep="")),do.call(paste,list(method1[2],"*",method2,sep="")),do.call(paste,list(method1[3],"*",method2,sep="")))
method3<-c("array-only","seq-only","both")
method1<-c("mean","median","std")
method2<-c("intensity","count")
colna<-NULL
for(i in 1:3){for(j in 1:2){for(k in 1:3){colna<-c(colna,paste(method3[i],method2[j],method1[k],sep="*"))}}} # generate column names for out

for(level in 1:8)
{
assign(paste("matrix",level,sep=""),as.data.frame(get(paste("matrix",level,sep=""))))
temp<-get(paste("matrix",level,sep=""))
rownames(temp)<-name.method
colnames(temp)<- colna
assign(paste("matrix",level,sep=""),temp)
write.table(temp,paste("c:/rnaseq/simulation/",paste("matrix.case",level,sep=""),".csv",sep=""),sep=",")
} #output














