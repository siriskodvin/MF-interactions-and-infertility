# Function that identifies effect alleles (= major alleles)
f_eff_all_major <- function(pAll1, pAll2){
  .pAll1 <- pAll1
  .pAll2 <- pAll2
  .allCount <- as.data.frame(table(c(.pAll1, .pAll2)))
  
  if(0 %in% .allCount$Var1){
    .allCount <- .allCount[-which(.allCount$Var1 == 0),]
  }
  .allCount <- .allCount[order(.allCount$Freq, decreasing = TRUE, na.last = TRUE),]
  .allEff <- as.character(.allCount$Var1[1])
  
  return(.allEff)
}


