# Recode and logistic regression of the complementary model
f_log_reg_m68 <- function(mDose, fDose, pheno, fam, cov){
  require(lme4)
  .mfDose <- rep(0, length(mDose))
  .mfDose[which((mDose == 2 & fDose == 0) | (mDose == 0 & fDose == 2))] <- 1
  .ART <- pheno
  .fam_big <- fam$fam_big
  .fam_unique <- fam$fam_unique
  .cov <- cov

  # General linear mixed effects model
  if(is.null(.cov)){
    .mod <- glmer(.ART ~ .mfDose + (1|.fam_unique), family = binomial(link = 'logit'), nAGQ = 0)
  }else{
    .form_main <- ".ART ~ .mfDose + (1|.fam_unique)"
    .form_cov <- paste0(paste0(".cov$", names(.cov)), collapse = " + ")
    .mod <- glmer(paste0(.form_main, " + ", .form_cov), family = binomial(link = 'logit'), nAGQ = 0)
  }
  .modSum <- summary(.mod)

  # Adding missing estimates when there are no interaction estimates
  if(!(".mfDose" %in% rownames(.modSum$coefficients))){
    rowno <- dim(.modSum$coefficients)[1]
    rownam <- c(rownames(.modSum$coefficients)[1], ".mfDose", rownames(.modSum$coefficients)[2:rowno])
    .modSum$coefficients <- rbind(.modSum$coefficients[1,], rep(NA, 4), .modSum$coefficients[2:rowno,])
    rownames(.modSum$coefficients) <- rownam
  }

  .modResSing <- data.frame(N = length(.modSum$residuals))

  .modSumCoef <- .modSum$coefficients
  .modSumCoefVec <- as.vector(t(.modSumCoef))
  names(.modSumCoefVec) <- paste0(rep(rownames(.modSumCoef), each = 4), colnames(.modSumCoef))  
  .modSumCoefVec <- as.data.frame(t(.modSumCoefVec))

  .modResSing <- cbind(.modResSing, .modSumCoefVec)
  return(.modResSing)
}


