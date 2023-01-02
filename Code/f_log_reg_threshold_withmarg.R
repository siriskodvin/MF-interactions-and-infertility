# Logistic regression function for the threshold model with marginal effects
f_log_reg_threshold_withmarg <- function(mAll1, mAll2, fAll1, fAll2, mfThresh, pheno, fam, cov){
  require(lme4)
  .mAll1 <- mAll1
  .mAll2 <- mAll2
  .mAll <- .mAll1 + .mAll2
  .fAll1 <- fAll1
  .fAll2 <- fAll2
  .fAll <- .fAll1 + .fAll2
  .mfThresh <- mfThresh
  .ART <- pheno
  .fam_big <- fam$fam_big
  .fam_unique <- fam$fam_unique
  .cov <- cov

  # General linear mixed effects model (with threshold w/ marginal effects)
  if(is.null(.cov)){
    .mod <- glmer(.ART ~ .mAll + .mAll1:.mAll2 + .fAll + .fAll1:.fAll2 + .mfThresh + (1|.fam_unique), family = binomial(link = 'logit'), nAGQ = 0)
  }else{
    .form_main <- ".ART ~ .mAll + .mAll1:.mAll2 + .fAll + .fAll1:.fAll2 + .mfThresh + (1|.fam_unique)"
    .form_cov <- paste0(paste0(".cov$", names(.cov)), collapse = " + ")
    .mod <- glmer(paste0(.form_main, " + ", .form_cov), family = binomial(link = 'logit'), nAGQ = 0)
  }
  .modSum <- summary(.mod)

  # Adding missing estimates when there are no effect estimates
  if(!(".mfThresh" %in% rownames(.modSum$coefficients))){
    if(!(".mAll1:.mAll2" %in% rownames(.modSum$coefficients))){
      if(!(".fAll1:.fAll2" %in% rownames(.modSum$coefficients))){
        rowno <- dim(.modSum$coefficients)[1]
        rownam <- c(rownames(.modSum$coefficients)[1:3], ".mfThresh", rownames(.modSum$coefficients)[4:rowno])
        .modSum$coefficients <- rbind(.modSum$coefficients[1:3,], rep(NA, 4), .modSum$coefficients[4:rowno,])
        rownames(.modSum$coefficients) <- rownam
        .modSum$coefficients <- rbind(.modSum$coefficients, rep(NA, 4))
        rownames(.modSum$coefficients)[dim(.modSum$coefficients)[1]] <- ".mAll1:.mAll2"
        .modSum$coefficients <- rbind(.modSum$coefficients, rep(NA, 4))
        rownames(.modSum$coefficients)[dim(.modSum$coefficients)[1]] <- ".fAll1:.fAll2"
      }else{
        rowno <- dim(.modSum$coefficients)[1]
        rownam <- c(rownames(.modSum$coefficients)[1:3], ".mfThresh", rownames(.modSum$coefficients)[4:rowno])
        .modSum$coefficients <- rbind(.modSum$coefficients[1:3,], rep(NA, 4), .modSum$coefficients[4:rowno,])
        rownames(.modSum$coefficients) <- rownam
        rowno <- dim(.modSum$coefficients)[1]
        rownam <- c(rownames(.modSum$coefficients)[1:14], ".mAll1:.mAll2", rownames(.modSum$coefficients)[rowno])
        .modSum$coefficients <- rbind(.modSum$coefficients[1:14,], rep(NA, 4), .modSum$coefficients[rowno,])
        rownames(.modSum$coefficients) <- rownam
      }
    }else if(!(".fAll1:.fAll2" %in% rownames(.modSum$coefficients))){
      rowno <- dim(.modSum$coefficients)[1]
      rownam <- c(rownames(.modSum$coefficients)[1:3], ".mfThresh", rownames(.modSum$coefficients)[4:rowno])
      .modSum$coefficients <- rbind(.modSum$coefficients[1:3,], rep(NA, 4), .modSum$coefficients[4:rowno,])
      rownames(.modSum$coefficients) <- rownam
      .modSum$coefficients <- rbind(.modSum$coefficients, rep(NA, 4))
      rownames(.modSum$coefficients)[dim(.modSum$coefficients)[1]] <- ".fAll1:.fAll2"
    }else{
      rowno <- dim(.modSum$coefficients)[1]
      rownam <- c(rownames(.modSum$coefficients)[1:3], ".mfThresh", rownames(.modSum$coefficients)[4:rowno])
      .modSum$coefficients <- rbind(.modSum$coefficients[1:3,], rep(NA, 4), .modSum$coefficients[4:rowno,])
      rownames(.modSum$coefficients) <- rownam
    }
  }else if(!(".mAll1:.mAll2" %in% rownames(.modSum$coefficients))){
    if(!(".fAll1:.fAll2" %in% rownames(.modSum$coefficients))){
      .modSum$coefficients <- rbind(.modSum$coefficients, rep(NA, 4))
      rownames(.modSum$coefficients)[dim(.modSum$coefficients)[1]] <- ".mAll1:.mAll2"
      .modSum$coefficients <- rbind(.modSum$coefficients, rep(NA, 4))
      rownames(.modSum$coefficients)[dim(.modSum$coefficients)[1]] <- ".fAll1:.fAll2"
    }else{
      rowno <- dim(.modSum$coefficients)[1]
      rownam <- c(rownames(.modSum$coefficients)[1:14], ".mAll1:.mAll2", rownames(.modSum$coefficients)[rowno])
      .modSum$coefficients <- rbind(.modSum$coefficients[1:14,], rep(NA, 4), .modSum$coefficients[rowno,])
      rownames(.modSum$coefficients) <- rownam
    }
  }else if(!(".fAll1:.fAll2" %in% rownames(.modSum$coefficients))){
    .modSum$coefficients <- rbind(.modSum$coefficients, rep(NA, 4))
    rownames(.modSum$coefficients)[dim(.modSum$coefficients)[1]] <- ".fAll1:.fAll2"
  }

  .modResSing <- data.frame(N = length(.modSum$residuals))
    
  .modSumCoef <- .modSum$coefficients
  .modSumCoefVec <- as.vector(t(.modSumCoef))
  names(.modSumCoefVec) <- paste0(rep(rownames(.modSumCoef), each = 4), colnames(.modSumCoef))  
  .modSumCoefVec <- as.data.frame(t(.modSumCoefVec))

  .modResSing <- cbind(.modResSing, .modSumCoefVec)
  return(.modResSing)
}


