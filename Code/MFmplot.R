#### MERGE AND PLOT RESULTS ####

library(qqman)
library(Haplin)

# Input
chr <- "02"

# Read res-files
path <- /PATH_TO_CODE_AND_RESULTS/
resFiles <- list.files(path)
resFiles <- resFiles[grepl("res", resFiles, fixed = TRUE)]
resFiles <- resFiles[endsWith(resFiles, "rds")]
resList <- list()
if(length(resFiles) != 0){
  for(j in 1:length(resFiles)){
    pathFil <- paste0(path, "/", resFiles[j])
    res <- readRDS(file = pathFil)
    resList[[length(resList) + 1]] <- res
  }
}

# Check dimensions of res files
varno <- NULL
for(i in 1:length(resList)){
  varno <- c(varno, dim(resList[[i]])[2])
}
varmax <- max(varno)
for(j in 1:length(resList)){
  if(varno[j] < varmax){
    resList[[j]] <- cbind(resList[[j]], matrix(NA, dim(resList[[j]])[1], varmax - dim(resList[[j]])[2]))
    colnames(resList[[j]]) <- colnames(resList[[which.max(varno)]])
  }
}

res <- do.call("rbind", resList)
dinds <- which(duplicated(res))
if(length(dinds) != 0){
  res <- res[-dinds,]
}
saveRDS(res, file = paste0(path, "/res", chr, ".rds"))
writeLines(as.character(dim(res)[1]), paste0(path, "/numvar", chr, ".txt"))

# Manhattan plot
res <- res[!is.na(res$`.mDose:.fDosePr(>|z|)`),]
ymax <- max((-log10(min(res$`.mDose:.fDosePr(>|z|)`)) + 1), 8)
jpeg(file = paste(path, "/manh", chr, ".jpeg", sep = ""),
     width = 1200, height = 600, quality = 100)
par(mar = c(4, 18, 11, 20))
manhattan(res, # xlim = c(xmin, xmax),
          p = ".mDose:.fDosePr(>|z|)",					# Change in accordance with type of interaction model
          ylim = c(0, ymax),
          col = c("blue", "red", "green"), suggestiveline = FALSE,
          annotatePval = 5e-08, annotateTop = FALSE)
dev.off()

# QQ plot Haplin
pdf(file = paste0(path, "/qq", chr, ".pdf"))
pQQ(res$`.mDose:.fDosePr(>|z|)`)						# Change in accordance with type of interaction model
dev.off()