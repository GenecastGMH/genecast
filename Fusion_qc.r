
#================读取参数=========================#
arg <- commandArgs(T)
if ((arg[1]=="-h")|(length(arg)==0)) {
  cat("Usage : Rscript Fusion_qc.r input_patient_id sample_list (output_dir)\n")
  quit("no")
}
# sample_list : 
#   patient time2 tissueType tissue
#   P27 2 tumor2 tumor2_1a

#================定义常量=========================#
patient_data <- "/work/shared/GeneticTest/"
cnv_dir <- "/work/user/zhanggh/Work/CNVresult/CNVreport/"
sample_dir <- "/work/shared/Analysis/"
sample_raw_file <- "/work/data/projects/sequencing_data_list.txt"
# arg <- c("P331","/work/user/gemh/Tools/mytools/QC_jikou/sample_list_P331.txt","/work/user/gemh/Tools/mytools/QC_jikou/P331/")
#================加载包===========================#
tryCatch({
  library(dplyr,quietly = T,verbose = F,warn.conflicts = F)
  library(stringr,quietly = T,verbose = F,warn.conflicts = F)
  library(data.table,quietly = T,verbose = F,warn.conflicts = F)
},error = function(e){
  cat("Need R packages: dplyr,stringr,data.table\n")
  quit("no")
})
if (length(arg)==2) {
  arg[3] <- paste(patient_data,arg[1],sep = "")
}

result_dir <- paste(arg[3],"/qc_data/",sep = "")
if (!file.exists(result_dir)){
  dir.create(result_dir,recursive = TRUE)
}
sample_list <- read.table(arg[2],sep = "\t",header = F,stringsAsFactors = F)
#==============
tryCatch({
  fusion_file_dir <-  paste(patient_data,arg[1],"/database/",arg[1],
                         "_",sample_list$V3,"_",sample_list$V2,".fusion.txt",sep = "")
  fusion_data <- c()
  for (i in 1:length(fusion_file_dir)) {
    if (!file.exists(fusion_file_dir[i])) {
      cat(paste("Missing ",arg[1]," file : ",fusion_file_dir[i],sep = ""),"\n")
      quit("no")
    }
    if(file.info(fusion_file_dir[i])$size!=0){
      temp_data <- read.table(fusion_file_dir[i],sep = "\t",fill = T,stringsAsFactors = F,header = F)
      fusion_data <- rbind(fusion_data,
                           select(temp_data,patient_id=V2,submission_time=V3,
                                  tissue=V26,chrom=V4,chrom2=V7,pos=V5,pos2=V8,strand=V6,
                                  strand2=V9,gene=V10,gene2=V15,transcript=V11,
                                  transcript2=V16,exon=V12,exon2=V17,
                                  anno=V14,anno2=V19,reads=V20,
                                  depth1=V21,depth2=V22,freq=V23,freq1=V24,
                                  freq2=V25))
         }
  }
  if (is.null(fusion_data)) {
    cat(paste(arg[1]," fusion is empty! " ,sep = ""),"\n")
    quit("no") 
  }
  fusion_data[,18:23][fusion_data[,18:23]=="-"] <- ""
  write.table(fusion_data,paste(result_dir,"fusion_qc.txt",sep = ""),sep = "\t",quote = F,row.names = F)
  cat("successful!\n")
  
},error = function(e){
  cat("Bug!\n")
  quit("no")})












