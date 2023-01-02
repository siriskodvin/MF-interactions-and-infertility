#### -------------------------------- Script header -------------------------------- ####
# Date:         09.12.2022                                                              #
# Author:       Siri Naerland Skodvin                                                   #
# Filename:     MFinter.R                                                               #
# Description:  Reads genetic data (ped- and map-files), a key file to identify case    #
#               and control mother-father dyads, and a covariate data frame. The        #
#               interaction analyses are done using generalized linear mixed-effects    #
#               models with logistic regression.                                        #
#### --------------------------------------------------------------------------------####



#### INPUT, DATA AND PACKAGES/FUNCTIONS ####

# Paths
cpath <- /PATH_TO_CODE_AND_RESULTS/
dpath <- /PATH_TO_TEMPORARY_STORAGE_AND_MAIN_R_SCRIPTS/

library(parallel)     # To run analyses in parallel
source(paste0(cpath, "/f_inter.R"))

# Set chromosome number (two characters) - CHANGE HERE
chr <- "01"

# Load mf structure file and extract family and couple ids
mfs <- readRDS(file = /COUPLE_KEY_FILE.rds/)
fam <- mfs[,c("fam_big", "fam_unique")]
mfs <- mfs[, 1:4]

# Load covariate matrix (optional - alternatively use cov <- NULL)
cov <- readRDS(file = /COV_FILE.rds/)



#### ANALYSES ####

# List of ped and map file paths
fileAll <- list.files(dpath)
fileInt <- cbind(paste0(dpath, "/", fileAll[endsWith(fileAll, ".ped")]),
                 paste0(dpath, "/", fileAll[endsWith(fileAll, ".map")]))
fileInt <- split(fileInt, seq(nrow(fileInt)))

# Initializing cluster for parallel processing
coresNo <- length(fileInt)
cl <- makeCluster(coresNo)

# Call for interaction regression function
res <- parLapply(cl = cl, X = fileInt, fun = f_inter, mfs = mfs, fam = fam, cov = cov)

# Set equal dimensions for .modRes elements
varno <- NULL
for(i in 1:length(res)){
  varno <- c(varno, dim(res[[i]])[2])
}
varmax <- max(varno)
for(j in 1:length(res)){
  if(varno[j] < varmax){
    res[[j]] <- cbind(res[[j]], matrix(NA, dim(res[[j]])[1], varmax - dim(res[[j]])[2]))
    colnames(res[[j]]) <- colnames(res[[which.max(varno)]])
  }
}
res <- do.call("rbind", res)
rownames(res) <- NULL

# Save results and log with time stamp
t <- Sys.time()
t <- gsub(" ", "_", t)
resFile <- paste0(t, "_res_chr", chr, ".rds")
saveRDS(res, file = paste0(cpath, "/", resFile))