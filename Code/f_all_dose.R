# Function that recodes allele info (from allele names to number of doses of effect allele)
# The effect allele is the first element of each allele vector
f_all_dose <- function(all1, all2){
  .allEff <- all1[1]
  .all1 <- all1[-1]
  .all2 <- all2[-1]
  
  .allRec1 <- rep(0, length(.all1))
  .allRec1[which(.all1 == 0)] <- NA
  .allRec1[which(.all1 == .allEff)] <- 1
  .allRec2 <- rep(0, length(.all2))
  .allRec2[which(.all2 == 0)] <- NA
  .allRec2[which(.all2 == .allEff)] <- 1
  
  .dose <- .allRec1 + .allRec2
  
  return(.dose)
}


