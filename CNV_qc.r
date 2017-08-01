
#================读取参数=========================#
arg <- commandArgs(T)
if ((arg[1]=="-h")|(length(arg)==0)) {
  cat("Usage : Rscript CNV_qc.r input_patient_id sample_list (output_dir)\n")
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
# arg <- c("P4008","/work/user/gemh/Tools/mytools/QC_jikou/sample_list_P4008.txt","/work/user/gemh/Tools/mytools/QC_jikou/")
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

sample_list <- read.table(arg[2],sep = "\t",header = F,stringsAsFactors = F) %>%
  filter(V3!="blood")
#===========
tryCatch({
  cnv_file_dir <-  paste(patient_data,arg[1],"/database/",arg[1],
                            "_",sample_list$V3,"_",sample_list$V2,".cnv.txt",sep = "")
  cnv_data <- c()
  for (i in 1:length(cnv_file_dir)) {
    if (!file.exists(cnv_file_dir[i])) {
      cat(paste("Missing ",arg[1]," file : ",cnv_file_dir[i],sep = ""),"\n")
      quit("no")
    }
      temp_data <- read.table(cnv_file_dir[i],sep = "\t",fill = T,stringsAsFactors = F,header = F)
      temp_data[,11] <- (2^(temp_data$V8)) *2
      cnv_data <- rbind(cnv_data,
                           select(temp_data,patient_id=V2,
                                  submission_time=V3,chrom=V4,
                                  start=V5,end=V6,gene=V7,log2_ratio=V8,
                                  tissue=V9,copy_number=V11))

  }

  cnv_data[cnv_data=="\\N"] <- ""
  write.table(cnv_data,paste(result_dir,"cnv_qc.txt",sep = ""),sep = "\t",quote = F,row.names = F)
  cat("successful!\n")
  
},error = function(e){
  cat("Bug!\n")
  quit("no")})










