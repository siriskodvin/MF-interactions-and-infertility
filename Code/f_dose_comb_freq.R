# Function that summarizes the frequencies of ART and nonART for all parental dose combinations
f_dose_comb_freq <- function(mDose, fDose, pheno){
  .mDose <- mDose
  .fDose <- fDose
  .ART <- pheno

  # Initiate dataframe to store results
  .modResSing <- data.frame(tmp = 0)

  # Calculate and save allele combination frequencies
  .mDoseART <- .mDose[which(.ART == 1)]
  .fDoseART <- .fDose[which(.ART == 1)]
  .ARTfreq <- table(.mDoseART, .fDoseART)
  .ARTfreq <- as.data.frame(.ARTfreq)
  if(0 %in% .ARTfreq$.mDoseART){
    if(0 %in% .ARTfreq$.fDoseART){
      .modResSing$ART_0_0 <- .ARTfreq$Freq[which(.ARTfreq$.mDoseART == 0 & .ARTfreq$.fDoseART == 0)]
    }else{
      .modResSing$ART_0_0 <- NA
    }
    if(1 %in% .ARTfreq$.fDoseART){
      .modResSing$ART_0_1 <- .ARTfreq$Freq[which(.ARTfreq$.mDoseART == 0 & .ARTfreq$.fDoseART == 1)]
    }else{
      .modResSing$ART_0_1 <- NA
    }
    if(2 %in% .ARTfreq$.fDoseART){
      .modResSing$ART_0_2 <- .ARTfreq$Freq[which(.ARTfreq$.mDoseART == 0 & .ARTfreq$.fDoseART == 2)]
    }else{
      .modResSing$ART_0_2 <- NA
    }
  }else{
    .modResSing$ART_0_0 <- NA
    .modResSing$ART_0_1 <- NA
    .modResSing$ART_0_2 <- NA
  }
  if(1 %in% .ARTfreq$.mDoseART){
    if(0 %in% .ARTfreq$.fDoseART){
      .modResSing$ART_1_0 <- .ARTfreq$Freq[which(.ARTfreq$.mDoseART == 1 & .ARTfreq$.fDoseART == 0)]
    }else{
      .modResSing$ART_1_0 <- NA
    }
    if(1 %in% .ARTfreq$.fDoseART){
      .modResSing$ART_1_1 <- .ARTfreq$Freq[which(.ARTfreq$.mDoseART == 1 & .ARTfreq$.fDoseART == 1)]
    }else{
      .modResSing$ART_1_1 <- NA
    }
    if(2 %in% .ARTfreq$.fDoseART){
      .modResSing$ART_1_2 <- .ARTfreq$Freq[which(.ARTfreq$.mDoseART == 1 & .ARTfreq$.fDoseART == 2)]
    }else{
      .modResSing$ART_1_2 <- NA
    }
  }else{
    .modResSing$ART_1_0 <- NA
    .modResSing$ART_1_1 <- NA
    .modResSing$ART_1_2 <- NA
  }
  if(2 %in% .ARTfreq$.mDoseART){
    if(0 %in% .ARTfreq$.fDoseART){
      .modResSing$ART_2_0 <- .ARTfreq$Freq[which(.ARTfreq$.mDoseART == 2 & .ARTfreq$.fDoseART == 0)]
    }else{
      .modResSing$ART_2_0 <- NA
    }
    if(1 %in% .ARTfreq$.fDoseART){
      .modResSing$ART_2_1 <- .ARTfreq$Freq[which(.ARTfreq$.mDoseART == 2 & .ARTfreq$.fDoseART == 1)]
    }else{
      .modResSing$ART_2_1 <- NA
    }
    if(2 %in% .ARTfreq$.fDoseART){
      .modResSing$ART_2_2 <- .ARTfreq$Freq[which(.ARTfreq$.mDoseART == 2 & .ARTfreq$.fDoseART == 2)]
    }else{
      .modResSing$ART_2_2 <- NA
    }
  }else{
    .modResSing$ART_2_0 <- NA
    .modResSing$ART_2_1 <- NA
    .modResSing$ART_2_2 <- NA
  }
  .mDosenonART <- .mDose[which(.ART == 0)]
  .fDosenonART <- .fDose[which(.ART == 0)]
  .nonARTfreq <- table(.mDosenonART, .fDosenonART)
  .nonARTfreq <- as.data.frame(.nonARTfreq)
  if(0 %in% .nonARTfreq$.mDosenonART){
    if(0 %in% .nonARTfreq$.fDosenonART){
      .modResSing$nonART_0_0 <- .nonARTfreq$Freq[which(.nonARTfreq$.mDosenonART == 0 & .nonARTfreq$.fDosenonART == 0)]
    }else{
      .modResSing$nonART_0_0 <- NA
    }
    if(1 %in% .nonARTfreq$.fDosenonART){
      .modResSing$nonART_0_1 <- .nonARTfreq$Freq[which(.nonARTfreq$.mDosenonART == 0 & .nonARTfreq$.fDosenonART == 1)]
    }else{
      .modResSing$nonART_0_1 <- NA
    }
    if(2 %in% .nonARTfreq$.fDosenonART){
      .modResSing$nonART_0_2 <- .nonARTfreq$Freq[which(.nonARTfreq$.mDosenonART == 0 & .nonARTfreq$.fDosenonART == 2)]
    }else{
      .modResSing$nonART_0_2 <- NA
    }
  }else{
    .modResSing$nonART_0_0 <- NA
    .modResSing$nonART_0_1 <- NA
    .modResSing$nonART_0_2 <- NA
  }
  if(1 %in% .nonARTfreq$.mDosenonART){
    if(0 %in% .nonARTfreq$.fDosenonART){
      .modResSing$nonART_1_0 <- .nonARTfreq$Freq[which(.nonARTfreq$.mDosenonART == 1 & .nonARTfreq$.fDosenonART == 0)]
    }else{
      .modResSing$nonART_1_0 <- NA
    }
    if(1 %in% .nonARTfreq$.fDosenonART){
      .modResSing$nonART_1_1 <- .nonARTfreq$Freq[which(.nonARTfreq$.mDosenonART == 1 & .nonARTfreq$.fDosenonART == 1)]
    }else{
      .modResSing$nonART_1_1 <- NA
    }
    if(2 %in% .nonARTfreq$.fDosenonART){
      .modResSing$nonART_1_2 <- .nonARTfreq$Freq[which(.nonARTfreq$.mDosenonART == 1 & .nonARTfreq$.fDosenonART == 2)]
    }else{
      .modResSing$nonART_1_2 <- NA
    }
  }else{
    .modResSing$nonART_1_0 <- NA
    .modResSing$nonART_1_1 <- NA
    .modResSing$nonART_1_2 <- NA
  }
  if(2 %in% .nonARTfreq$.mDosenonART){
    if(0 %in% .nonARTfreq$.fDosenonART){
      .modResSing$nonART_2_0 <- .nonARTfreq$Freq[which(.nonARTfreq$.mDosenonART == 2 & .nonARTfreq$.fDosenonART == 0)]
    }else{
      .modResSing$nonART_2_0 <- NA
    }
    if(1 %in% .nonARTfreq$.fDosenonART){
      .modResSing$nonART_2_1 <- .nonARTfreq$Freq[which(.nonARTfreq$.mDosenonART == 2 & .nonARTfreq$.fDosenonART == 1)]
    }else{
      .modResSing$nonART_2_1 <- NA
    }
    if(2 %in% .nonARTfreq$.fDosenonART){
      .modResSing$nonART_2_2 <- .nonARTfreq$Freq[which(.nonARTfreq$.mDosenonART == 2 & .nonARTfreq$.fDosenonART == 2)]
    }else{
      .modResSing$nonART_2_2 <- NA
    }
  }else{
    .modResSing$nonART_2_0 <- NA
    .modResSing$nonART_2_1 <- NA
    .modResSing$nonART_2_2 <- NA
  }
  .modResSing <- subset(.modResSing, select = -tmp)
  return(.modResSing)
}


