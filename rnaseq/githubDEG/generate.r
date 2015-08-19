seqdata<-function(i,j,k,lambda,delta,phi)
{
  #generate RNA-Seq data as described in method section, arg detail see gen
  Y<-array(0,c(i,j,k))
  for(a in 1:i)
  {for(b in 1:j)
  {for(c in 1:k)
  {
    Y[a,b,c]=rnbinom(1,size=1/phi[b],mu=lambda[b]*exp((a-1)*delta[b]))
  }}}
  Y
}

arraydata<-function(i,j,k,sigma1,u,delta)
{
  #generate Microarray data as described in method section, arg detail see gen
  Y<-array(0,c(i,j,k))
  for(a in 1:i)
  {for(b in 1:j)
  {for(c in 1:k)
  {
    level=(u[b])
    temp=8.9^2/(exp(sigma1^2)*(exp(sigma1^2)-1))
    yma<-log(max((exp(2*level)-temp),0)/(2*exp(level)))+rnorm(1,0,sigma1)+(a-1)*delta[b]
    Y[a,b,c]=max(exp(yma)+rnorm(1,23.2,8.9),0)
  }}}
  Y
}

getresult<-function(table){
  #show sensitivity and specificity of each method from the generated DEG table
  flag<-rep(0,10000)
  if(colnames(table)[ncol(table)-1]=="significant"){
    flag[table[,ncol(table)-1]==1 & abs(table[,1])>=1]<-1
    true.positive<-nrow(matrix(table[flag==1 & table[,ncol(table)]!=0,],nc=ncol(table)))
    false.positive<-nrow(matrix(table[flag==1  & table[,ncol(table)]==0,],nc=ncol(table)))
    undetected<-nrow(matrix(table[flag==0 & abs(table[,ncol(table)])>1,],nc=ncol(table)))
    
  }
  if(colnames(table)[ncol(table)-1]=="FDR"){
    flag[table[,ncol(table)-1]<0.05 & abs(table[,1])>=1]<-1
    true.positive<-nrow(matrix(table[flag==1 & table[,ncol(table)]!=0,],nc=ncol(table)))
    false.positive<-nrow(matrix(table[flag==1  & table[,ncol(table)]==0,],nc=ncol(table)))
    undetected<-nrow(matrix(table[flag==0 & abs(table[,ncol(table)])>1,],nc=ncol(table)))
  }
  cbind(true.positive,false.positive,undetected)
}
gen<-function(alpha1=0.462,beta1=1/0.00467,alpha2=0.85,beta2=0.5,mu0=5.4,sigma0=1.08,de=0.5,sigma=0.7,sigma1=0.467,i=2,j=10000,k=2,percent=0.2)
{
# genrate parallele Microarray and RNA-Seq as decribed in method section. 
#alpha1 and beta1 are gamma parameters for mean distribution of RNA-Seq 
#alpha2 and beta2 are gamma parameters for dispersion distribution of RNA-Seq 
#mu0 and sigma0 are normal parameters for mean distribution of Microarray
#de and sigma are normal parameters for fold change distribution of DEGs.
#sigma1 is the multiplicative white noise variance as described in Durbin, Blythe P., et al. "A variance-stabilizing transformation for gene-expression microarray data." Bioinformatics 18.suppl 1 (2002): S105-S110.
# i is number of groups
# j is number of genes
# k is number of subjects in each group
# percent is the expected percentage of DEG
lambda<-rgamma(j,alpha1,scale=beta1)
temp1<-lambda[lambda>400]
temp2<-rgamma(length(temp1),7.724,5.716)
temp3<-temp1*sort(temp2)[rank(temp1)]
lambda[lambda>400]<-temp3
delta<-rep(0,j)

index<-sample(j,percent*j)
prob<-rbinom(index,1,0.5)
delta[index]=prob*rnorm(length(index),de*log(2),sigma)+(1-prob)*rnorm(length(index),-de*log(2),sigma)
table<-cbind(index,delta[index])
colnames(table)<-c("DE gene","log fold change")

phi<-rgamma(j,alpha2,scale=beta2)
temp<-rnorm(j,mu0,sigma0)
mu<-abs(temp-min(temp))+min(temp)
mean<-rep(0,j)
temp<-sort(mu)
tempp<-floor(rank(lambda))
for(l in 1:j)
{mean[l]=temp[tempp[l]]}

seq<-seqdata(i,j,k,lambda,delta,phi)
array<-arraydata(i,j,k,sigma1,u=mean,delta)
return(list(seq,array,table,delta))

}
