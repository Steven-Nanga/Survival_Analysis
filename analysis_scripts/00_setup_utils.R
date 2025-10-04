#!/usr/bin/env Rscript
# Utility functions shared by the dataset analysis scripts
load_or_install <- function(pkgs){
  to_install <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  if(length(to_install)) {
    message("Installing missing packages: ", paste(to_install, collapse=", "))
    install.packages(to_install, repos = "https://cloud.r-project.org")
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

safe_mkdir <- function(dir){
  if(!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

# export a gtsummary table to a Word (.docx) file using flextable/officer
export_gtsummary_docx <- function(tbl, file){
  if(!('flextable' %in% installed.packages()[,'Package']) || !('officer' %in% installed.packages()[,'Package'])){
    install.packages(c('flextable','officer'), repos = 'https://cloud.r-project.org')
  }
  ft <- tryCatch(as_flex_table(tbl), error = function(e) as_flex_table(tbl %>% as_flextable()))
  doc <- officer::read_docx()
  doc <- officer::body_add_flextable(doc, value = ft)
  print(doc, target = file)
}
