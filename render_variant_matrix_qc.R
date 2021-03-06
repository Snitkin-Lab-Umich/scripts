# SNP and/or Indel variant matrix QC
# Argument 1: SNP_mat_code.csv or Indel_mat_code.csv
# Argument 2: prefix for output files (optional)
# Output: html document with QC statistics and plots

# Get command line arguments
args = commandArgs(trailingOnly = TRUE)

if(length(args) == 2){
  pref=paste0(getwd(), "/", format(Sys.time(), "%Y-%m-%d"),'_',args[[2]])
}else{
  pref=paste0(getwd(), "/", format(Sys.time(), "%Y-%m-%d"))
}

if(grepl('SNP',args[[1]],ignore.case=T)){
  pref=paste0(pref,'_SNP')
  tmpdir = paste0('tmp_',args[[2]],'_SNP')
  title='SNP Matrix QC'
}else if(grepl('Indel',args[[1]],ignore.case=T)){
  pref=paste0(pref,'_Indel')
  tmpdir = paste0('tmp_',args[[2]],'_Indel')
  title='Indel Matrix QC'
}

mat = normalizePath(args[[1]])

rmarkdown::render('/nfs/esnitkin/bin_group/pipeline/Github/scripts/variant_matrix_qc.R', 
                  params=list(mat=mat,title=title),
                  output_format='html_document',
                  output_file=paste0(pref,"_variant_matrix_qc.html"),
                  intermediates_dir=tmpdir)

unlink(tmpdir,recursive=T)
