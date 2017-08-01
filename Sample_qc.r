
#================读取参数=========================#
arg <- commandArgs(T)
if ((arg[1]=="-h")|(length(arg)==0)) {
  cat("Usage : Rscript Sample_qc.r input_patient_id sample_list (output_dir)\n")
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
#================创建输出文件夹===============#
if (length(arg)==2) {
  arg[3] <- paste(patient_data,arg[1],sep = "")
}

result_dir <- paste(arg[3],"/qc_data/",sep = "")
if (!file.exists(result_dir)){
  dir.create(result_dir,recursive = TRUE)
}
#=================================================#
tryCatch({
  #====Genecast/list.txt
  if(!file.exists(paste(patient_data,"list.txt",sep = ""))){
    cat(paste("Missing "," list.txt !",sep = ""),"\n")
    quit("no")            
  }
  if(!file.exists(arg[2])){
    cat(paste("Missing ",arg[2],sep = ""),"\n")
    quit("no")            
  }
  sample_list <- read.table(arg[2],sep = "\t",header = F,stringsAsFactors = F)
  patient_list <- left_join( read.table(paste(patient_data,"list.txt",sep = ""),sep = "\t",
                                        stringsAsFactors = F,header = T,fill =T) %>% filter(patient==arg[1]), 
                             read.table(sample_raw_file,sep = "\t",fill=T,header = F,stringsAsFactors = F,
                                        col.names=c("sequencingId","sample_dir")),
                             by=c("sequencingId"))
  patient_list <- semi_join(patient_list,sample_list,by=c("patient"="V1","time2"="V2","tissueType"="V3","tissue"="V4"))
  #==========查看文件完整性========#
  if(!file.exists(paste(patient_data,arg[1],"/mpileup/summary.txt.gz",sep = ""))){
    cat(paste("Missing ",arg[1]," file : ",
              paste(patient_data,arg[1],"/mpileup/summary.txt.gz",sep = ""),sep = ""),"\n")
    quit("no")
  }
  if(!file.exists(paste(patient_data,arg[1],"/mapped/summary.txt.gz",sep = ""))){
    cat(paste("Missing ",arg[1]," file : ",
              paste(patient_data,arg[1],"/mapped/summary.txt.gz",sep = ""),sep = ""),"\n")
    quit("no")
  }
  for (i in 1:length(unique(patient_list$sample_dir))) {
    fqstat_dir <- paste("/work/shared/Analysis/",
                        unique(patient_list$sample_dir)[i],
                        "/fqstat/summary.txt.gz",sep = "")
    if(!file.exists(fqstat_dir)){
      cat(paste("Missing ",arg[1]," file : ",fqstat_dir,sep = ""),"\n")
      quit("no")
    }
  }
  for (i in 1:nrow(patient_list)) {
    duplica_dir <- paste("/work/shared/Analysis/",
                         patient_list$sample_dir[i],
                         "/mapped/",patient_list$sampleId[i],
                         "_",patient_list$indexId[i],
                         ".sorted.filtered.mkdup.metrics",sep = "")
    if(!file.exists(duplica_dir)){
      cat(paste("Missing ",arg[1]," file : ",duplica_dir,sep = ""),"\n")
      quit("no")
    }
  }
  #==========fqstat/summary========#
  fqstat_summary <- c()
  for (i in 1:length(unique(patient_list$sample_dir))) {
    fqstat_dir <- paste("/work/shared/Analysis/",
                        unique(patient_list$sample_dir)[i],
                        "/fqstat/summary.txt.gz",sep = "")
    fqstat_sample <-  filter(patient_list,sample_dir==unique(patient_list$sample_dir)[i]) %>%
      select(patient_id= patient,id=tissue,sampleId,indexId)
    
    fqstat_sample$sampleId <- paste(fqstat_sample$sampleId,"_",fqstat_sample$indexId,sep = "")
    
    fqstat_summary <- rbind(fqstat_summary,left_join(fqstat_sample[,-4] ,
                                                     fread(paste("zcat ",fqstat_dir,sep = ""),
                                                           sep = "\t",header = T,integer64 = "numeric"),by=c("sampleId"="id")))
  }
  #=修改成小数
  fqstat_summary[,6] <- as.numeric(as.character(str_split(fqstat_summary[,6],"%",simplify = T)[,1]))/100
  fqstat_summary[,7] <- as.numeric(as.character(str_split(fqstat_summary[,7],"%",simplify = T)[,1]))/100
  fqstat_summary[,8] <- as.numeric(as.character(str_split(fqstat_summary[,8],"%",simplify = T)[,1]))/100
  fqstat_summary[,9] <- as.numeric(as.character(str_split(fqstat_summary[,9],"%",simplify = T)[,1]))/100
  fqstat_summary[,10] <- as.numeric(as.character(str_split(fqstat_summary[,10],"%",simplify = T)[,1]))/100

  #==============duplicate==============================#
  duplica_summary <- patient_list[,c(1,4)]
  for (i in 1:nrow(patient_list)) {
    duplica_dir <- paste("/work/shared/Analysis/",
                         patient_list$sample_dir[i],
                         "/mapped/",patient_list$sampleId[i],
                         "_",patient_list$indexId[i],
                         ".sorted.filtered.mkdup.metrics",sep = "")
    if(!file.exists(duplica_dir)){
      cat(paste("Missing ",arg[1]," ",duplica_dir," duplica-mkdup.metrics !",sep = ""),"\n")
      quit("no")
    }
    duplica_summary[i,3] <- str_split(readLines(duplica_dir)[8],"\t",simplify = T)[,8]
  }
  #================mpileup/summary============================#
  mpileup_summary <- fread(paste("zcat ",patient_data,arg[1],"/mpileup/summary.txt.gz",sep = ""),
                           sep = "\t",header = T,integer64 = "numeric") %>% data.frame()
  #=修改成小数
  mpileup_summary[,6] <-  as.numeric(as.character(str_split(mpileup_summary[,6],"%",simplify = T)[,1]))/100
  mpileup_summary[,9] <-  as.numeric(as.character(str_split(mpileup_summary[,9],"%",simplify = T)[,1]))/100
  mpileup_summary[,11] <-  as.numeric(as.character(str_split(mpileup_summary[,11],"%",simplify = T)[,1]))/100
  mpileup_summary[,13] <-  as.numeric(as.character(str_split(mpileup_summary[,13],"%",simplify = T)[,1]))/100
  mpileup_summary[,15] <-  as.numeric(as.character(str_split(mpileup_summary[,15],"%",simplify = T)[,1]))/100
  mpileup_summary[,17] <-  as.numeric(as.character(str_split(mpileup_summary[,17],"%",simplify = T)[,1]))/100
  mpileup_summary[,19] <-  as.numeric(as.character(str_split(mpileup_summary[,19],"%",simplify = T)[,1]))/100
  
  #================mapped/summary============================#
  
  mapped_summary <- fread(paste("zcat ",patient_data,arg[1],"/mapped/summary.txt.gz",sep = ""),
                          sep = "\t",header = T,integer64 = "numeric") %>% data.frame()
  #=修改成小数
  mapped_summary[,4] <-  as.numeric(as.character(str_split( mapped_summary[,4],"%",simplify = T)[,1]))/100
  mapped_summary[,6] <-  as.numeric(as.character(str_split( mapped_summary[,6],"%",simplify = T)[,1]))/100
  mapped_summary[,8] <-  as.numeric(as.character(str_split( mapped_summary[,8],"%",simplify = T)[,1]))/100
  mapped_summary[,10] <-  as.numeric(as.character(str_split( mapped_summary[,10],"%",simplify = T)[,1]))/100
  
  #==============combind=======================================#
  result <- full_join(mapped_summary,mpileup_summary,by=c("id")) %>% mutate(patient_id=arg[1])
  result <- full_join(result,fqstat_summary[,-3],by=c("patient_id","id"))
  result <- full_join(result,duplica_summary,by=c("patient_id"="patient","id"="tissue"))
  
  ID <- str_split(result$id,"_",simplify = T) %>% data.frame()
  result <- mutate(result,submission_time=str_extract(ID$X2, "\\d+"),
                   tissue = ID$X1)
  
  #===============sex===================================#
  cnv_file <- paste(cnv_dir,arg[1],".",filter(result,tissue!="blood")$id,".cnvreport.txt",sep = "")
  sex <- c()
  for (i in 1:length(cnv_file)) {
    if(!file.exists(cnv_file[i])){
      sex <- c(sex,"-")
      cat(paste("Missing ",arg[1],"-",filter(result,tissue!="blood")$id[i],"-",cnv_file," cnv-file !",sep = ""),"\n")
    } 
    
    sex <- c(sex , scan(cnv_file[i],what = "character",sep="\t")[7])
  }
  sex <-  substr(sex, 1, 1)
  result <- left_join(result,data.frame(id= filter(result,tissue!="blood")$id,
                                        sex,stringsAsFactors = F),
                      by="id")
  #======= V44:"germline_mutation_concordance",V45:"avg_reads"======#
  result[,44] <- ""
  result[,45] <- ""
  #==============整理=====================#
  colnames(result) <- gsub("%","",colnames(result))
  result <- select(result,patient_id,submission_time,tissue,panel,panel_size=panelSize,
                   gender=sex,total_reads=total,avg_reads=V44,gc_ratio=GC,
                   Q20_ratio=Q20,Q30_ratio=Q30,mapped_reads=mapped,mapping_ratio=mappedRatio,
                   paired_mapping_reads=paired,paired_mapping_ratio=pairedRatio,
                   duplication_ratio=V3,targeted_ratio=targetedRatio,
                   average_sequencing_depth_on_target=depth,panel_coverage_ratio_1X=coverage,
                   panel_coverage_ratio_100X=coverage100x,panel_coverage_ratio_1000X=coverage1000x,
                   panel_coverage_ratio_5000X=coverage5000x,n_ratio=N,germline_mutation_concordance=V45
  )
  result[is.na(result)] <- ""
  #============写文件==============#
  write.table(result,paste(result_dir,"sample_qc.txt",sep = ""),sep = "\t",quote = F,row.names = F)
  cat("successful!\n")
},error = function(e){
  cat("Bug!\n")
  quit("no")})