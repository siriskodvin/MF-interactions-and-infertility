# Logistic dose-response regression function
f_log_reg <- function(mDose, fDose, pheno, fam, cov){
  require(lme4)

  .mDose <- mDose
  .fDose <- fDose
  .ART <- pheno
  .fam_big <- fam$fam_big
  .fam_unique <- fam$fam_unique
  .cov <- cov

  # General linear mixed effects model
  if(is.null(.cov)){
    .mod <- glmer(.ART ~ .mDose + .fDose + .mDose:.fDose + (1|.fam_unique), family = binomial(link = 'logit'), nAGQ = 0)
  }else{
    .form_main <- ".ART ~ .mDose + .fDose + .mDose:.fDose + (1|.fam_unique)"
    .form_cov <- paste0(paste0(".cov$", names(.cov)), collapse = " + ")
    .mod <- glmer(paste0(.form_main, " + ", .form_cov), family = binomial(link = 'logit'), nAGQ = 0)
  }
  .modSum <- summary(.mod)

  # Adding missing estimates when there are no interaction estimates
  if(!(".mDose:.fDose" %in% rownames(.modSum$coefficients))){
    .modSum$coefficients <- rbind(.modSum$coefficients, rep(NA, 4))
    rownames(.modSum$coefficients)[dim(.modSum$coefficients)[1]] <- ".mDose:.fDose"
  }

  .modResSing <- data.frame(N = length(.modSum$residuals))

  # Check deficiency warning
  if(length(.modSum$fitMsgs) > 0){
    .modSumCoefVec <- NA
  }else{
    .modSumCoef <- .modSum$coefficients
    .modSumCoefVec <- as.vector(t(.modSumCoef))
    names(.modSumCoefVec) <- paste0(rep(rownames(.modSumCoef), each = 4), colnames(.modSumCoef))  
    .modSumCoefVec <- as.data.frame(t(.modSumCoefVec))
  }
  .modResSing <- cbind(.modResSing, .modSumCoefVec)
  return(.modResSing)
}


