suppressPackageStartupMessages(library(haven))
root <- file.path(Sys.getenv("DXA_DATA_ROOT", file.path(getwd(), "data_raw")), "charls_official/raw_extracted/2011")
rows <- list()
for(fn in c("health_status_and_functioning.dta","family_information.dta","health_care_and_insurance.dta")) {
  d <- read_dta(file.path(root,fn), n_max=0)
  labels <- vapply(d,function(x) {
    z<-attr(x,"label");if(is.null(z)) "" else as.character(z)
  },character(1))
  keep <- grepl("fracture|steroid|rheuma|arthritis|alcohol|liquor|wine|beer|parent",labels,ignore.case=TRUE)
  rows[[fn]] <- data.frame(file=rep(fn,sum(keep)),variable=names(d)[keep],label=labels[keep])
}
x <- do.call(rbind,rows)
print(x,row.names=FALSE)
write.csv(x,file.path(Sys.getenv("DXA_PROJECT_DIR", file.path(getwd(), "DXA_transportability_results")), "reviewer1_revision/analysis_results/R1_FRAX_raw_variable_audit.csv"),row.names=FALSE)
