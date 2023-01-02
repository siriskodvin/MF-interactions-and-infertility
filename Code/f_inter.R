# Function that reads a ped file (data.frame) and loads the mfs structure,
# recodes and performs the mf-interaction analysis.
# Use comment/uncomment to select the desired type of interaction model.
f_inter <- function(filenames, mfs, fam, cov = NULL){
  require(data.table)
  require(dplyr)
  source("/mnt/work/workbench/sirinaerland.skodvin/MFinter/Code/funk/f_eff_all.R")
  source("/mnt/work/workbench/sirinaerland.skodvin/MFinter/Code/funk/f_eff_all_major.R")
  source("/mnt/work/workbench/sirinaerland.skodvin/MFinter/Code/funk/f_all_dose.R")
  source("/mnt/work/workbench/sirinaerland.skodvin/MFinter/Code/funk/f_dose_comb_freq.R")
  source("/mnt/work/workbench/sirinaerland.skodvin/MFinter/Code/funk/f_all_dose_dominant.R")
  source("/mnt/work/workbench/sirinaerland.skodvin/MFinter/Code/funk/f_all_dose_recessive.R")
  source("/mnt/work/workbench/sirinaerland.skodvin/MFinter/Code/funk/f_all_dose_threshold_withmarg.R")
  source("/mnt/work/workbench/sirinaerland.skodvin/MFinter/Code/funk/f_log_reg.R")
  source("/mnt/work/workbench/sirinaerland.skodvin/MFinter/Code/funk/f_log_reg_threshold_withmarg.R")
  source("/mnt/work/workbench/sirinaerland.skodvin/MFinter/Code/funk/f_log_reg_m68.R")

  # Read ped file
  .ped <- fread(file = filenames[1])
  .ped <- as.data.frame(.ped)
  .ped <- .ped[, -c(3:6)]
  
  # Checking if the data contains only one SNP
  if(dim(.ped)[2] == 4){.sing <- TRUE}else{.sing <- FALSE}
  
  # Read map file and adjust to store results (res)
  .res <- fread(file = filenames[2])
  .res <- as.data.frame(.res)
  .res <- .res[, -3]
  colnames(.res) <- c("CHR", "SNP", "BP")

  # Merging mfs structure with mat and pat genetics
  .mGen <- left_join(mfs, .ped, by = c("fam" = "V1", "m_id" = "V2"))
  .mGen <- .mGen[,-c(1:4)]
  .fGen <- left_join(mfs, .ped, by = c("fam" = "V1", "f_id" = "V2"))
  .fGen <- .fGen[,-c(1:4)]

  ## OPTION 1 - EFFECT ALLELE = MINOR ALLELE
  # Identifying the minor alleles
  if(.sing){
    .effAllAll <- f_eff_all(.ped[,3], .ped[,4])
  }else{
    .col <- seq(3, ncol(.ped), by = 2)
    .effAllAll <- mapply(f_eff_all, .ped[, .col], .ped[, -c(1, 2, .col)])
  }

#  ## OPTION 2 - EFFECT ALLELE = MAJOR ALLELE
#  # Identifying the major alleles
#  if(.sing){
#    .effAllAll <- f_eff_all_major(.ped[,3], .ped[,4])
#  }else{
#    .col <- seq(3, ncol(.ped), by = 2)
#    .effAllAll <- mapply(f_eff_all_major, .ped[, .col], .ped[, -c(1, 2, .col)])
#  }

  # Saving the effect alleles in .res
  .res$A1 <- .effAllAll

  # Calculate effect allele frequencies
  .res$A1_freq <- NA
  for(i in 1:length(.effAllAll)){
    .alls <- c(.ped[,2*i + 1], .ped[,2*i + 2])
    .res$A1_freq[i] <- length(.alls[which(.alls == .effAllAll[i])])/
      (length(.alls[which(.alls == .effAllAll[i])]) +
         length(.alls[which(.alls != .effAllAll[i] & .alls != 0)]))
  }

  # Remove abundant objects
  rm(.ped)

  # Saving the effect alleles in the first row of the m/f ped structures (for use in analyses)
  .effAllAll <- rep(.effAllAll, each = 2)
  .mGen <- rbind(.effAllAll, .mGen)
  .fGen <- rbind(.effAllAll, .fGen)

  ## RECODE TO DOSE = 0,1,2
  # Call for allele recode function
  if(.sing){
    .mDoseAll <- f_all_dose(.mGen[,1], .mGen[,2])
    .fDoseAll <- f_all_dose(.fGen[,1], .fGen[,2])
  }else{
    .col <- seq(1, ncol(.mGen), by = 2)
    .mDoseAll <- as.data.frame(mapply(f_all_dose, .mGen[, .col], .mGen[, -.col]))
    .fDoseAll <- as.data.frame(mapply(f_all_dose, .fGen[, .col], .fGen[, -.col]))
  }

  # Retrieve phenotype
  .ART <- mfs$pheno

  # Add dose combination frequencies to .res
  if(.sing){
    .combFreqs <- f_dose_comb_freq(.mDoseAll, .fDoseAll, .ART)
  }else{
    .combFreqs <- mapply(f_dose_comb_freq, .mDoseAll, .fDoseAll, MoreArgs = list(.ART), SIMPLIFY = FALSE)
    .combFreqs <- do.call("rbind", .combFreqs)
  }
  .res <- cbind(.res, .combFreqs)

#  ## OPTIONAL - RECODE TO DOSE = 0,1 (DOMINANT)
#  # Call for allele recode function
#  if(.sing){
#    .mDoseAll <- f_all_dose_dominant(.mGen[,1], .mGen[,2])
#    .fDoseAll <- f_all_dose_dominant(.fGen[,1], .fGen[,2])
#  }else{
#    .col <- seq(1, ncol(.mGen), by = 2)
#    .mDoseAll <- as.data.frame(mapply(f_all_dose_dominant, .mGen[, .col], .mGen[, -.col]))
#    .fDoseAll <- as.data.frame(mapply(f_all_dose_dominant, .fGen[, .col], .fGen[, -.col]))
#  }

#  ## OPTIONAL - RECODE TO DOSE = 0,1 (RECESSIVE)
#  # Call for allele recode function
#  if(.sing){
#    .mDoseAll <- f_all_dose_recessive(.mGen[,1], .mGen[,2])
#    .fDoseAll <- f_all_dose_recessive(.fGen[,1], .fGen[,2])
#  }else{
#    .col <- seq(1, ncol(.mGen), by = 2)
#    .mDoseAll <- as.data.frame(mapply(f_all_dose_recessive, .mGen[, .col], .mGen[, -.col]))
#    .fDoseAll <- as.data.frame(mapply(f_all_dose_recessive, .fGen[, .col], .fGen[, -.col]))
#  }

#  ## OPTIONAL - RECODE TO SEPARATE ALLELIC DOSES (THRESHOLD MODEL WITH THRESHOLD = k)
#  k <- 3
#  # Split alleles
#  if(.sing){
#    .mAll1raw <- .mGen[,1]
#    .mAll2raw <- .mGen[,2]
#    .fAll1raw <- .fGen[,1]
#    .fAll2raw <- .fGen[,2]
#  }else{
#    .mAll1raw <- .mGen[,seq(1, dim(.mGen)[2], by = 2)]
#    .mAll2raw <- .mGen[,seq(2, dim(.mGen)[2], by = 2)]
#    .fAll1raw <- .fGen[,seq(1, dim(.fGen)[2], by = 2)]
#    .fAll2raw <- .fGen[,seq(2, dim(.fGen)[2], by = 2)]
#  }
#  # Call for allele recode function
#    if(.sing){
#    .mAll1all <- f_all_dose_threshold_withmarg(.mAll1raw)
#    .mAll2all <- f_all_dose_threshold_withmarg(.mAll2raw)
#    .fAll1all <- f_all_dose_threshold_withmarg(.fAll1raw)
#    .fAll2all <- f_all_dose_threshold_withmarg(.fAll2raw)
#  }else{
#    .mAll1all <- as.data.frame(apply(.mAll1raw, 2, f_all_dose_threshold_withmarg))
#    .mAll2all <- as.data.frame(apply(.mAll2raw, 2, f_all_dose_threshold_withmarg))
#    .fAll1all <- as.data.frame(apply(.fAll1raw, 2, f_all_dose_threshold_withmarg))
#    .fAll2all <- as.data.frame(apply(.fAll2raw, 2, f_all_dose_threshold_withmarg))
#  }
#  # Generate combination variables and recode according to threshold
#  .mfThresh <- .mAll1all + .mAll2all + .fAll1all + .fAll2all
#  .mfThresh[.mfThresh < k] <- 0
#  .mfThresh[.mfThresh >= k] <- 1

  # Remove abundant objects
  rm(.mGen, .fGen)

  ## OPTION 1 - STANDARD INTERACTION MODEL (use this model with the standard 0,1,2 recoding, or dominant/recessive recoding)
  # Call for logistic regression function, returning the interaction coefficient and p value
  if(.sing){
    .modRes <- f_log_reg(.mDoseAll, .fDoseAll, .ART, fam, cov)
  }else{
    .modRes <- mapply(f_log_reg, .mDoseAll, .fDoseAll, MoreArgs = list(.ART, fam, cov), SIMPLIFY = FALSE)

    # Set equal dimensions for .modRes elements
    varno <- NULL
    for(i in 1:length(.modRes)){
      varno <- c(varno, dim(.modRes[[i]])[2])
    }
    varmax <- max(varno)
    for(j in 1:length(.modRes)){
      if(varno[j] < varmax){
        .modRes[[j]] <- cbind(.modRes[[j]], matrix(NA, dim(.modRes[[j]])[1], varmax - dim(.modRes[[j]])[2]))
        colnames(.modRes[[j]]) <- colnames(.modRes[[which.max(varno)]])
      }
    }
    .modRes <- do.call("rbind", .modRes)
  }

#  ## OPTION 2 - THRESHOLD MODEL WITH MARGINAL EFFECTS (use this model with threshold recoding)
#  # Call for logistic regression for threshold model
#  if(.sing){
#    .modRes <- f_log_reg_threshold_withmarg(.mAll1all, .mAll2all, .fAll1all, .fAll2all, .mfThresh, .ART, fam, cov)
#  }else{
#    .modRes <- mapply(f_log_reg_threshold_withmarg, .mAll1all, .mAll2all, .fAll1all, .fAll2all, .mfThresh, MoreArgs = list(.ART, fam, cov), SIMPLIFY = FALSE)
#    .modRes <- do.call("rbind", .modRes)
#  }

#  ## OPTION 3 - COMPLEMENTARY MODEL (use this model with the standard 0,1,2 recoding)
#  # Call for logistic regression for the complementary model
#  if(.sing){
#    .modRes <- f_log_reg_m68(.mDoseAll, .fDoseAll, .ART, fam, cov)
#  }else{
#    .modRes <- mapply(f_log_reg_m68, .mDoseAll, .fDoseAll, MoreArgs = list(.ART, fam, cov), SIMPLIFY = FALSE)
#    .modRes <- do.call("rbind", .modRes)
#  }

  .res <- cbind(.res, .modRes)

  return(.res)
}