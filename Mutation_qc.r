
#================读取参数=========================#
arg <- commandArgs(T)
if ((arg[1]=="-h")|(length(arg)==0)) {
  cat("Usage : Rscript Mutation_qc.r input_patient_id sample_list (output_dir)\n")
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

sample_list <- read.table(arg[2],sep = "\t",header = F,stringsAsFactors = F)
#=================================================#
tryCatch({
  snp_file_dir <-  paste(patient_data,arg[1],"/database/",arg[1],
                         "_",unique(sample_list[,c(1:2)])$V2,".snv.txt",sep = "")
  snp_data <- c()
  for (i in 1:length(snp_file_dir)) {
    if (!file.exists(snp_file_dir[i])) {
      cat(paste("Missing ",arg[1]," file : ",snp_file_dir[i],sep = ""),"\n")
      quit("no")
    }
    varscan_file_dir <- paste(patient_data,arg[1],"/varscan/",
                             sample_list[sample_list$V2 == (unique(sample_list[,c(1:2)])$V2[i]),]$V4,
                              ".cns.vcf.gz",sep = "")
    varscan_data <- c()
    for (j in 1:length(varscan_file_dir)) {
      if (!file.exists(varscan_file_dir[j])) {
        cat(paste("Missing ",arg[1]," file : ",varscan_file_dir[j],sep = ""),"\n")
        quit("no")
      }
      varscan_data <- rbind(varscan_data,fread(paste("zcat ",varscan_file_dir[j],sep = ""),sep = "\t",header = T)[,c(1,2,4,5,10)] %>%
        mutate(id = sample_list[sample_list$V2 == (unique(sample_list[,c(1:2)])$V2[i]),]$V4[j]))
    }
    varscan_data <-  cbind(varscan_data[,-6], data.frame(str_split(varscan_data$id,"_",simplify = T)))
    varscan_data[,7] <- as.numeric(as.character(str_extract(varscan_data[,7], "\\d+")))
    varscan_data <- cbind(varscan_data[,-5],data.frame(str_split(varscan_data$Sample1,":",simplify = T)[,c(6,4)],stringsAsFactors = F))
    colnames(varscan_data)[c(1,5:8)] <- c("chr","tissuetype","time2","reads","depth")
    temp_data <- fread(snp_file_dir[i],sep = "\t",header = F) %>% data.frame()
    temp_data <- left_join(temp_data,varscan_data,by=c("V18"="tissuetype",
                                                       "V3"="time2",
                                                       "V4"="chr","V5"="POS",
                                                       "V6"="REF","V7"="ALT"))
    snp_data <- rbind(snp_data,
                      select(temp_data,patient_id=V2,submission_time=V3,tissue=V18,
                             chrom=V4,pos=V5,ref=V6,alt=V7,gene=V8,types=V9,
                             transcript=V10,exon=V12,mutation=V14,protein_change=V13,
                             mutation_c=V23,dbsnp=V15,pop_mut_freq_eas=V16,
                             pop_mut_freq_eur=V17,frequency=V19,reads,depth))
  }
  snp_data$frequency <- as.numeric(as.character(str_split(snp_data$frequency,"%",simplify = T)[,1]))/100
  snp_data[is.na(snp_data)] <- ""
  snp_data[snp_data=="\\N"] <- ""
  
  write.table(snp_data,paste(result_dir,"mutation_qc.txt",sep = ""),sep = "\t",quote = F,row.names = F)
  cat("successful!\n")
  
},error = function(e){
  cat("Bug!\n")
  quit("no")})