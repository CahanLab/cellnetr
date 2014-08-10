#
# Copyright (C) 2014 Patrick Cahan
#

# to run this file from command line:
# Rscript ./cellnetr_website.R config_mouse4302.txt

logFname<-paste("log_", format(Sys.time(), "%Y_%b_%d_%H_%M_%S"),".txt", sep='');

args <- commandArgs(trailingOnly = TRUE)
configFname<-args[1];

# path to CellNet objects. Unsure right now where these files will be stored on Orchestra.
path_CN_obj<-"~/CellNet/CN_Objects/";
path_CN_code<-"~/git_repos/cellnetr/"
#.libPaths(c("~/installed_R_libs","/home/pc87/R/x86_64-unknown-linux-gnu-library/2.15","/opt/R-2.15.3/lib64/R/library"));
### require(cellnetr);
source( paste(path_CN_code,"cellnetr_apply.R", sep=''));
source( paste(path_CN_code,"cellnetr_cnres_ops.R", sep=''));
source( paste(path_CN_code,"cellnetr_graphics.R", sep=''));
source( paste(path_CN_code,"cellnetr_preprocess.R", sep=''));
source( paste(path_CN_code,"cellnetr_utils.R", sep=''));

# path to custom R packages. On Orchestra, this should include CellNet, randomForest, and any organism/platform annotations that differ from mine
### path_R_packages<-""
###  .libPaths(path_R_packages);

logString<-"";
lstring<-"Start web-based analysis";
cat(cellnetr_log(lstring), file=logFname, append=TRUE);

lstring<-"Load libraries";
cat(cellnetr_log(lstring), file=logFname, append=TRUE);

require(affy);
require(ggplot2);
require(gplots);
require(randomForest);


# USER supplied arguments, extracted from configFame
# - analysis_name
# - ct_target
# - fname_sampTab
# - platform
# - dLevel (optional)

lstring<-"Read user-supplied parameters";
cat(cellnetr_log(lstring), file=logFname, append=TRUE);

raw <- scan(configFname, what="", sep="\n");
aList <- strsplit(raw, "=");
names(aList) <- sapply(aList, `[[`, 1);
aList <- lapply(aList, `[`, -1);
#aList[['sOrder']]<-strsplit(aList[['sOrder']], ",")[[1]];

mydate<-utils_myDate();

# load sample table and correct CEL file names
# assumes that the CEL files have been uploaded to the current directory
# if not, then prepend the correct path to stQuery$file_name 
lstring<-"Read and correct sample table";
cat(cellnetr_log(lstring), file=logFname, append=TRUE);

stQuery<-expr_readSampTab(aList[['fname_sampTab']]);
stQuery<-geo_fixNames(stQuery);

lstring<-paste("Found ", nrow(stQuery), " expression profiles", sep='');
cat(cellnetr_log(lstring), file=logFname, append=TRUE);

# load and normalize the query expression data
lstring<-"Load and normalize expression data";
cat(cellnetr_log(lstring), file=logFname, append=TRUE);

expQuery<-Norm_cleanPropRaw(stQuery, aList[['platform']]);  
lstring<-paste("Data load complete. ", nrow(expQuery), " genes",sep='');
cat(cellnetr_log(lstring), file=logFname, append=TRUE);

# load CellNet object for this platform. 
cnObjName<-switch(aList[['platform']],
                  mouse4302 = paste(path_CN_obj_obj, "cnProc_mouse4302_062414.R",sep=''),
                  mogene10stv1 = paste(path_CN_obj, "cnProc_mogene_062414.R",sep=''),
                  ###                    mouseIllumina8v2 = paste(path_CN_obj, "",sep=''), ### Add this platform later
                  hgu133plus2 = paste(path_CN_obj, "cnProc_Hg1332_062414.R",sep=''),
                  hugene10stv1 = paste(path_CN_obj, "cnProc_Hugene_062414.R",sep=''));

lstring<-paste("Load CellNet ",cnObjName,sep='');
cat(cellnetr_log(lstring), file=logFname, append=TRUE);

cnProc<-utils_loadObject(cnObjName);

# Run CellNet
lstring<-"Run CellNet";
cat(cellnetr_log(lstring), file=logFname, append=TRUE);
tmpAns<-cn_apply(expQuery, stQuery, cnProc, dLevelQuery=aList[['dLevel']]);

# Score transcription factors
lstring<-paste("Score transcription factors in GRN ", aList[['ct_target']],sep='');
cat(cellnetr_log(lstring), file=logFname, append=TRUE);
tfScores<-cn_nis_all(tmpAns, cnProc, aList[['ct_target']])

# Run standard output script
lstring<-"Make standard CellNet plots";
cat(cellnetr_log(lstring), file=logFname, append=TRUE);
plotFile<-cn_stdout(tmpAns, cnProc, tfScores, aList[['analysis_name']], aList[['ct_target']]);  
  
# Write tables
lstring<-"Write data: classification values, GRN scores, normalized expression, and NIS values";
cat(cellnetr_log(lstring), file=logFname, append=TRUE);
datFiles<-cn_outputRes(tmpAns, tfScores, aList[['analysis_name']]);

# close log file, tar and compress everything
lstring<-"Done";
cat(cellnetr_log(lstring), file=logFname, append=TRUE);

fnameOut<-paste(aList[['analysis_name']],"_analysis.tar", sep='');
cmd<-paste("tar -c ", logFname, " ",plotFile, " ", datFiles[1], " ",datFiles[2]," ", datFiles[3]," ", datFiles[4]," > ", fnameOut, sep='');
system(cmd);
cmd<-paste("gzip -f ", fnameOut, sep='');
system(cmd);
