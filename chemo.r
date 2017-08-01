
#=====================================#
arg <- commandArgs(T)
if ((arg[1]=="-h")|(length(arg)==0)) {
  cat("Usage : Rscript chemo.r input_patient_id sample_list (output_dir)\n")
  quit("no")
}
#=======================================#
options(warn =-1)
patient_data <- "/work/shared/GeneticTest/"
cnv_dir <- "/work/user/zhanggh/Work/CNVresult/CNVreport/"
sample_dir <- "/work/shared/Analysis/"
sample_raw_file <- "/work/data/projects/sequencing_data_list.txt"
# arg <- c("P4008","/work/user/gemh/Tools/mytools/QC_jikou/sample_list_P4008.txt","/work/user/gemh/Tools/mytools/QC_jikou/")
#================
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
  filter(V3=="blood")

chemical_rs_data <- fread("/work/user/gemh/Tools/mytools/QC_tools/chemical_rs_site.txt",sep = "\t",header = F,fill = T)

#=================================
tryCatch({
 chemo_file_dir <-  paste(patient_data,arg[1],"/varcheck.mkdup/",
                          sample_list$V4,".varcheck.gz",sep = "")
 chemo_data <- c()
 change_format <- function(varcheck_line){
   if(nchar(varcheck_line[5])<=1 & nchar(varcheck_line[8])==0){
     return(data.table(t(varcheck_line))[,1:7])
   }
    data.table(V1=varcheck_line[1],
               V2=varcheck_line[2],
               V3=varcheck_line[3],
               V4=varcheck_line[4],
    V5=c(str_split(varcheck_line[5],"[|]",simplify = T),
         str_split(varcheck_line[8],"[|]",simplify = T)),
    V6= c(str_split(varcheck_line[6],"[|]",simplify = T),
          str_split(varcheck_line[9],"[|]",simplify = T)),
    V7= c(str_split(varcheck_line[7],"[|]",simplify = T),
    str_split(varcheck_line[10],"[|]",simplify = T))) %>%
      filter(V5!="") %>% return()
 }
 arrange_format <- function(loc_data) {
   temp_result <- filter(temp_data,V1==loc_data[1] & V2==loc_data[2])
   if(nrow(temp_result)==1){
     temp_result[,8] <- paste(temp_result[,3],temp_result[,4],sep = "")
     return(temp_result[,c(1,2,8)])
   }
   
   temp_result <- arrange(temp_result,desc(V6))
   temp_result[,8] <- paste(temp_result$V5[1],temp_result$V5[2],sep = "")
   return(temp_result[1,c(1,2,8)])
 }
for (i in 1:length(chemo_file_dir)) {
  if (!file.exists(chemo_file_dir[i])) {
    cat(paste("Missing ",arg[1]," file : ",chemo_file_dir[i],sep = ""),"\n")
    quit("no")
  }
temp_data  <- fread(paste("zcat ",chemo_file_dir[i],sep = ""),sep = "\t",header = F) 
temp_data <- semi_join(temp_data,chemical_rs_data,by=c("V1"="V4","V2"="V6"))
temp_data  <-  apply(as.matrix(temp_data), 1, change_format) 
temp_data[is.na(temp_data)] <- NULL
temp_data <- rbindlist(temp_data) 
temp_data <- cbind(temp_data[,-c(4,6,7)],
                   data.frame(apply(as.matrix(temp_data)[,c(4,6,7)], 2, as.numeric),stringsAsFactors = F))
temp_data <- filter(temp_data,V6 >=10 & V4>=10)
chemo_data_tmp <- apply(as.matrix(unique(temp_data[,c(1,2)])),1,arrange_format) %>% rbindlist()
chemo_data_tmp$V2 <- as.numeric(chemo_data_tmp$V2)
chemo_data_tmp <- left_join(chemo_data_tmp,unique(chemical_rs_data[,-9]),by=c('V1'="V4","V2"="V6")) %>%
  mutate(pid=sample_list[i,1],submission_times=sample_list[i,2],no = "") %>%
  select(no,pid,gene=V2.y,locus=V3,detection_result=V8.x,tumor_type=V1.y,submission_times,chrom=V1,pos=V2,ref=V7,alt=V8.y)
chemo_data <- rbind(chemo_data,chemo_data_tmp)
}
 
  write.table(chemo_data,paste(result_dir,"chemo_qc.txt",sep = ""),sep = "\t",quote = F,row.names = F)
  
},error = function(e){
  cat("Bug!\n")
  quit("no")})





