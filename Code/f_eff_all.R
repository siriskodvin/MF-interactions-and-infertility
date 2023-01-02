# Function that identifies effect alleles (= minor alleles)
f_eff_all <- function(pAll1, pAll2){
  .pAll1 <- pAll1
  .pAll2 <- pAll2
  .allCount <- as.data.frame(table(c(.pAll1, .pAll2)))
  
  if(0 %in% .allCount$Var1){
    .allCount <- .allCount[-which(.allCount$Var1 == 0),]
  }
  .allCount <- .allCount[order(.allCount$Freq, decreasing = TRUE, na.last = TRUE),]
  .allEff <- as.character(.allCount$Var1[2])
  
  return(.allEff)
}


