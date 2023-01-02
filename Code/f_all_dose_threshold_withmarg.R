# Function that recodes allele info (from allele names to number of doses of effect allele)
# The effect allele is the first element of each allele vector
f_all_dose_threshold_withmarg <- function(all){
  .allEff <- all[1]
  .all <- all[-1]
  
  .allRec <- rep(0, length(.all))
  .allRec[which(.all == 0)] <- NA
  .allRec[which(.all == .allEff)] <- 1

  return(.allRec)
}


