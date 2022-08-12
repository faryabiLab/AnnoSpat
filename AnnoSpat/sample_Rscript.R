#filename='/mnt/data1/aanchal/data/IMC_T1D/labels_ssc/trte_labels_ssc_ELM_withImmuneCelltypes_withNegMarkers_d-adaptive2.csv'

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  for (filename in args) {
    cat("---R script running-------", sep="\n")
    dat <- read.csv(file = filename, header = TRUE)
    cat(nrow(dat) , ncol(dat) , sep="\n")
  }
}

main()