utils_loadObject<-function
### loads an R object when you don't know the name
(fname
 ### file name
){
  x<-load(fname);
  get(x);
}

utils_stripwhite<-function
### strip whitespace from a string
(string
 #### string
 ){
  gsub("^\\s+|\\s+$", "", string)
}

utils_myDate<-function
### print date
()
{
  format(Sys.time(), "%b_%d_%Y");
}

GEP_makeMean<-function
### return a dataframe of mean or median-ed data based on given groupings
(exp,
 ### exp df
 groupings,
 ### vector of groupings
 type='mean'
 ### mean or median
){
  
  ##<<note colnames become the column name of the first sample in each group from the original data
  ans<-data.frame();
  grps<-unique(groupings);
  if(type=='mean'){
    for(grp in grps){
      gi<-which(groupings==grp);
      if(length(gi)==1){
        cat("\nhere")
        if(nrow(ans)==0){
          ans<-data.frame(exp[,gi]);
        }else{
          ans<-cbind(ans, exp[,gi]);
        }
      }
      else{
        xxx<-apply(exp[,gi],1,mean);
        if(nrow(ans)==0){
          ans<-data.frame(xxx);
        }
        else{
          ans<-cbind(ans, xxx);
        }
      }
    }
  }
  else{
    for(grp in grps){
      gi<-which(groupings==grp);
      xxx<-apply(exp[,gi],1,median);
      if(nrow(ans)==0){
        ans<-data.frame(xxx);
      }
      else{
        ans<-cbind(ans, xxx);
      }
    }
  }
  
  colnames(ans)<-grps;
  ans;
  ### data.frame of mean or median-ed data based on given groupings
}


zscore<-function
### compute zscore
(x,
 ### numeric vector
 meanVal, 
 ### mean of distribution to compute zscore of x against 
 sdVal
 ### standard deviation of distribution to compute zscore of x agains
 ){ 
  (x-meanVal)/sdVal;
  ### zscore
}

cn_zscoreVect<-function
### Compute the mean zscore of given genes in each sample
(genes,
 ### genes
 xvals,
 ### named vector
 tVals,
 ### tvals
 ctt
 ### ctt
 ){
  ans<-vector();
  for(gene in genes){
    ans<-append(ans, zscore(xvals[gene], tVals[[ctt]][['mean']][[gene]], tVals[[ctt]][['sd']][[gene]]));
  }
  ans;
  ### zscore vector
}

cn_reduceMatLarge<-function
### reduce a data.frame
(datFrame,
 ### df
 valCol="score",
 ### value to merge
 cName="description",
 ### column to reduce on 
 iterOver="subNet"
 ### iterate over
 ){
  
  iterOvers<-unique(as.vector(datFrame[,iterOver]));
  ans<-data.frame();
  for(io in iterOvers){
    #  cat(io,"\n");
    xi<-which(datFrame[,iterOver]==io);
    dfTmp<-datFrame[xi,];
    x<- utils_reduceMat(dfTmp,valCol=valCol,cName=cName);
    x<-cbind(x, subNet=rep(io, nrow(x)));
    ans<-rbind(ans, x);
  }
  ans;
  ### ans
}

cn_extract_SN_DF<-function
### returns a DF of: sample_id, description, ctt, subnet_name, score
(scores,
 ### a matrix of subNet scores
 sampTab,
 ### sample table
 dLevel,
 ### column name of sampTab to group values by
 rnames=NULL
 ### rownames to extract
){
  
  if(is.null(rnames)){
    rnames<-rownames(scores);
    #cat("GOT NULL\n");
  }
  tss<-scores[rnames,];
  if(length(rnames)==1){
    tss<-t(as.matrix(scores[rnames,]));
    rownames(tss)<-rnames;
  #  cat(dim(tss),"\n")
  }
  nSamples<-ncol(tss);
  stTmp<-sampTab[colnames(tss),]; ####
  snNames<-rownames(tss);
  num_subnets<-length(snNames);    
  snNames<-unlist(lapply(snNames, rep, times=nSamples));
  sample_ids<-rep(as.vector(stTmp$sample_id), num_subnets);
  descriptions<-rep(as.vector(stTmp[,dLevel]), num_subnets);
    # myCtts<-rep(ctt, length(snNames));
  scores<-as.vector(t(tss));
  data.frame(sample_id=sample_ids,
             description=descriptions,
             #         ctt=myCtts,
             subNet = snNames, 
             score=scores);
  ### data.frame
}

utils_reduceMat<-function
### reduce the data.matrix values by averaging and getting st dvs
(datFrame,
 ### data frame to reduce
 valCol,
 ### value column
 cName='ann'
 ### cname to reduce on e.g. ann
){
    
  mids<-unique(as.vector(datFrame[,cName]));
  means<-rep(0, length(mids));
  sds<-rep(1, length(mids));
  indexI<-rep(0, length(mids)); # index to extract other columns from original data.frame
  
  for(i in seq(length(mids))){
    mid<-mids[i]
    xi<-which(datFrame[,cName]==mid);
    tmpDat<-datFrame[xi,];
    means[i]<-mean(tmpDat[,valCol]);
    sds[i]<-sd(tmpDat[,valCol]);
  }  
  data.frame(grp_name=mids, mean=means, stdev=sds);
  ### df of grp_name, mean, sd
}


utils_stderr<-function
### calculate standard error
(x){
  sqrt(var(x)/length(x));
  ### stderr
}

utils_concord<-function 
### concordance is a useful measure of set overalp 
(vect1,
 ### vector 1
 vect2
 ### vector 2
 ){
  length(intersect(vect1, vect2))/length(union(vect1, vect2));
  ### length(intersection) / length(union)
}

cellnetr_log<-function
#### make a CellNet logging string
(str
 ### str to log
 ){
  lstring1<-paste("CellNet [",format(Sys.time(), "%Y_%b_%d_%H:%M:%S"),"]: ", sep='');
  paste(lstring1, str, "\n", sep='');
  ### output str
}

