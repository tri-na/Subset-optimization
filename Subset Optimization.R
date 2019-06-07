
# preparation of data and features
rm(list=ls())
set.seed(5)
mydata=read.csv('Liion_comp_528.csv')
# define the grouping from hierachical clustering
mydata$hgp=NA
mydata[1:55,'hgp']=1
mydata[56:104,'hgp']=2
mydata[105:108,'hgp']=3
mydata[109:167,'hgp']=4
mydata[168:235,'hgp']=5
mydata[236:269,'hgp']=6
mydata[270:528,'hgp']=7


# we use all groups and then groups 4-7 to simultaneously reduce the within cluster sum to improve clustering
alldata=mydata
subdata=mydata[109:528,]
Xall_ini=as.matrix(mydata[,5:4505])
Xsub_ini=as.matrix(subdata[,5:4505])
yall=mydata$log10_cond
ysub=subdata$log10_cond

# function returning the variance for all groups
myobjall=function(X){
  ncut_vec=seq(2,10)
  hc_df=mydata[,c('formula_id','index','log10_cond')]
  hc_df[,as.character(ncut_vec)]=NA
  hc_run <- hclust(dist(X),method="ward.D2")
  for (cut_dex in 1:length(ncut_vec)){
    ncuts=ncut_vec[cut_dex]
    mycut=cutree(hc_run,ncuts)
    hc_df[,as.character(ncuts)]=mycut
  }
  hc_sum=data.frame(matrix(NA,nrow=length(ncut_vec),ncol=1))
  colnames(hc_sum)=c('label_var')
  rownames(hc_sum)=as.character(ncut_vec)
  # use the measure of summing both labeled and unlabeled
  for (group_meth in as.character(ncut_vec)){
    aggcount=aggregate(hc_df[,'log10_cond'],list(hc_df[,group_meth]),function(x)sum(!is.na(x)))
    out_stat=matrix(NA,nrow=nrow(aggcount),ncol=3)
    colnames(out_stat)=c('gp','count','var')
    aggvar=aggregate(na.omit(hc_df)[,'log10_cond'],list(na.omit(hc_df)[,group_meth]),var)
    out_stat[,'gp']=aggcount[,1]
    out_stat[,'count']=aggcount[,2]
    out_stat[match(aggvar[,1],out_stat[,'gp']),'var']=aggvar[,2]
    myout3=sum(out_stat[,'count']*out_stat[,'var'],na.rm=T)
    hc_sum[group_meth,'label_var']=myout3
  }
  return((hc_sum))
}

# function returning the variance for groups 4-7
myobjsub=function(X){
  ncut_vec=seq(2,10)
  hc_df=subdata[,c('formula_id','index','log10_cond')]
  hc_df[,as.character(ncut_vec)]=NA
  hc_run <- hclust(dist(X),method="ward.D2")
  for (cut_dex in 1:length(ncut_vec)){
    ncuts=ncut_vec[cut_dex]
    mycut=cutree(hc_run,ncuts)
    hc_df[,as.character(ncuts)]=mycut
  }
  hc_sum=data.frame(matrix(NA,nrow=length(ncut_vec),ncol=1))
  colnames(hc_sum)=c('label_var')
  rownames(hc_sum)=as.character(ncut_vec)
  # use the measure of summing both labeled and unlabeled
  for (group_meth in as.character(ncut_vec)){
    aggcount=aggregate(hc_df[,'log10_cond'],list(hc_df[,group_meth]),function(x)sum(!is.na(x)))
    out_stat=matrix(NA,nrow=nrow(aggcount),ncol=3)
    colnames(out_stat)=c('gp','count','var')
    aggvar=aggregate(na.omit(hc_df)[,'log10_cond'],list(na.omit(hc_df)[,group_meth]),var)
    out_stat[,'gp']=aggcount[,1]
    out_stat[,'count']=aggcount[,2]
    out_stat[match(aggvar[,1],out_stat[,'gp']),'var']=aggvar[,2]
    myout3=sum(out_stat[,'count']*out_stat[,'var'],na.rm=T)
    hc_sum[group_meth,'label_var']=myout3
  }
  return((hc_sum))
}


# initial run with all features weighted to be 1.
wvec=rep(1,4501)
Xwall=t(t(Xall_ini)*wvec)
Xwsub=t(t(Xsub_ini)*wvec)
benchobjall=myobjall(Xwall)
benchobjsub=myobjsub(Xwsub)
for (featdex in seq(1,4501,by=1)){
  if (wvec[featdex]==1){
    wvec[featdex]=0
    Xwsubn=t(t(Xsub_ini)*wvec)
    newobjsub=myobjsub(Xwsubn)
    
    if (all(newobjsub<=benchobjsub)){ 
      Xwalln=t(t(Xall_ini)*wvec)
      newobjall=myobjall(Xwalln)
      if (all(newobjall<=benchobjall)){
        wvec[featdex]=0
        benchobjall=newobjall
        benchobjsub=newobjsub
      }
      else {
        wvec[featdex]=1
      }
      
    }
    else {
      wvec[featdex]=1
      
    }
    print (benchobjall)
    print(wvec[featdex])
    print(featdex)
  }
}
write.table(wvec,'simultaneous_decreased_1.txt')


# run iteratively and save the weight vectors for each dimension of the mxrd
for (dex in seq(1,22)){
inputname=paste0('simultaneous_decreased_',dex,'.txt')
outputname=paste0('simultaneous_decreased_',dex+1,'.txt')

wvec=read.table(inputname,header = T)$x
Xwall=t(t(Xall_ini)*wvec)
Xwsub=t(t(Xsub_ini)*wvec)
benchobjall=myobjall(Xwall)
benchobjsub=myobjsub(Xwsub)
for (featdex in seq(1,4501,by=1)){
  if (wvec[featdex]==1){
    wvec[featdex]=0
    Xwsubn=t(t(Xsub_ini)*wvec)
    newobjsub=myobjsub(Xwsubn)
    
    if (all(newobjsub<=benchobjsub)){ 
      Xwalln=t(t(Xall_ini)*wvec)
      newobjall=myobjall(Xwalln)
      if (all(newobjall<=benchobjall)){
        wvec[featdex]=0
        benchobjall=newobjall
        benchobjsub=newobjsub
      }
      else {
        wvec[featdex]=1
      }
      
    }
    else {
      wvec[featdex]=1
      
    }
    print (benchobjall)
    print(wvec[featdex])
    print(featdex)
  }
}

write.table(wvec,outputname)

}