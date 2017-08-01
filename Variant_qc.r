
#================读取参数=========================#
arg <- commandArgs(T)
if ((arg[1]=="-h")|(length(arg)==0)) {
  cat("Usage : Rscript Variant_qc.r input_patient_id sample_list (output_dir)\n")
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
#=====================
tryCatch({
  recheck_file_dir <-  paste(patient_data,arg[1],"/recheck/",
                         sample_list$V4,".filtered.sites.gz",sep = "")
  recheck_data <- c()
 for (i in 1:length(recheck_file_dir)) {
   if (!file.exists(recheck_file_dir[i])) {
     cat(paste("Missing ",arg[1]," file : ",recheck_file_dir[i],sep = ""),"\n")
     quit("no")
   }
   temp_data <- fread(paste("zcat ",recheck_file_dir[i],sep=""),
                      sep="\t",header=T) %>% mutate(patient_id=sample_list[i,1],
                                                    submission_time=sample_list[i,2],
                                                    tissue=sample_list[i,3],
                                                    variant_bases_sd_quality="",
                                                    mkdup_variant_bases_sd_quality="",
                                                    mismatch_count_in_support_reads="" )
   
   recheck_data <- rbind(recheck_data,
                         select(temp_data,patient_id,submission_time,tissue,
                                chrom,pos,ref,alt,gene,exon,mutation=mutation.c,
                                protein_change=mutation.p,variant_vaf=freq,
                                supports_reads=reads,total_depth=depth,variant_positive=reads_plus,
                                variant_negative=reads_minus,variant_head=reads_begin,variant_tail=reads_end,
                                variant_bases_mean_quality=avg_qual,variant_bases_sd_quality,
                                mkdup_variant_vaf=freq.mkdup ,mkdup_supports_reads=reads.mkdup,
                                mkdup_total_depth=depth.mkdup,mkdup_variant_positive=reads_plus.mkdup,
                                mkdup_variant_negative=reads_minus.mkdup,mkdup_variant_head=reads_begin.mkdup,
                                mkdup_variant_tail=reads_end.mkdup,mkdup_variant_bases_mean_quality=avg_qual.mkdup,
                                mkdup_variant_bases_sd_quality,mismatch_count_in_support_reads,
                                recheck=confident
                                ))
 }
  write.table(recheck_data,paste(result_dir,"variant_qc.txt",sep = ""),sep = "\t",quote = F,row.names = F)
  cat("successful!\n")
},error = function(e){
  cat("Bug!\n")
  quit("no")})






