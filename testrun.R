# LiblineaR version 2.10-12 ist das package was wir zum vorhersagen brauchen

# navigate to diretory where you downloaded the github files
setwd("~/Desktop/RandomForest/3rd/Classifier_paper/general/github_files/")
outFun <- function(Counts, ID_class, sep) {
# 1. preprocessing ############################################################
# load ML models ###############################################################
# load ID conversion table 
  ID_conv <- read.csv("ID_conversion.tsv", sep = sep, stringsAsFactors = F)
# load count data, where the first column should be gene identifiers
  Counts <- read.csv(Counts, sep = "\t", stringsAsFactors = F, row.names = 1)
  cat("counts loaded...")
  
if (length(rownames(Counts)) == length(which(rownames(Counts) == as.character(1:nrow(Counts))))) {
  cat("Error: symbol, ensemble or entrez are not provided in the first column")
} else {
  
# select the genes used for classifier trainig
  ma <- match(ID_conv[,match(ID_class, colnames(ID_conv))], rownames(Counts))
  Counts <- Counts[ma[!is.na(ma)],]

# convert to symbol (classifier was trained on symbols)
  ma <- match(rownames(Counts), ID_conv[,match(ID_class, colnames(ID_conv))])
  rownames(Counts) <- ID_conv$symbol[ma]
  
# normalize data and scale between 0 and 1
  Counts.norm <- Counts+1
  Counts.norm <- apply(Counts.norm, 2, log10)
  Counts.norm <- apply(Counts.norm, 2, scale)
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  for (i in 1:ncol(Counts.norm)) {
    Counts.norm[,i] <- range01(Counts.norm[,i])
  }

# transpose data
  Counts.norm <- as.data.frame(t(Counts.norm))
  colnames(Counts.norm) <- rownames(Counts)
  colnames(Counts.norm) <- make.names(colnames(Counts.norm))

# find genes not provided by user
  ma <- match(ID_conv$symbol, rownames(Counts))
  GenesNoFound <- ID_conv$symbol[is.na(ma)]

# Print number of missing genes
  cat(paste0(length(GenesNoFound)), " of ", nrow(ID_conv), " genes not found\n")
  if ((length(GenesNoFound)) > 0) {
  cat("impute missing genes...\n")
  }
  
# impute missing genes
  GenesNoFound_df <- matrix(ID_conv$norm_exp[match(GenesNoFound, ID_conv$symbol)], nrow = nrow(Counts.norm), ncol = length(GenesNoFound))
  colnames(GenesNoFound_df) <- GenesNoFound
  Counts.norm <- cbind(Counts.norm, GenesNoFound_df)

# 2. classification ############################################################
# load ML models ###############################################################
  load("svm0.3_models.RData")
# predict 
  cat("predict...\n")
  preds <- list()
  for (i in 1:length(svm0.3)) {
    preds[[i]] <- as.character(predict(svm0.3[[i]], Counts.norm))
  }
# combine preds   
  preds <- do.call("cbind",preds)
  rownames(preds) <- rownames(Counts.norm)
  mat <- matrix(ncol =  length(unique(c(preds[,1],preds[,2],preds[,3],preds[,4],preds[,5],
                                      preds[,6],preds[,7],preds[,8],preds[,9],preds[,10]))))
  colnames(mat) <- unique(c(preds[,1],preds[,2],preds[,3],preds[,4],preds[,5],
                          preds[,6],preds[,7],preds[,8],preds[,9],preds[,10]))
  mat <- as.data.frame(mat)
  # count number of machines predicting a subtype 
  types <- ncol(mat)
  for (i in 1:nrow(preds)) {
    TAB <- table(as.character(preds[i,]))
    mat[i,] <- TAB[match(colnames(mat),names(TAB))]
  }
  mat[is.na(mat)] <- 0
  # count maximal number of machines predicting a subtype
  Max <- apply(mat, 1, max)
  Prediction <- colnames(mat)[apply(mat, 1, function(x) which.max(x))]
  mat$max_prob <- Max
  # make vales "probalistic"
  mat <- mat/10
  # assign predictions
  mat$prediction <- Prediction
  mat$prediction[which(mat$max<=0.5)] <- "Unclassified"
  mat$sample <- rownames(Counts.norm)
  cat("predictions saved in:", getwd())
  # save predictions
  write.table(mat,"predictions.tsv", sep = "\t", row.names = F)
}
}
outFun(Counts = "StJude_data_missing_genes.tsv", ID_class = "ensemble_ID", sep = "\t")
# Counts: /path/to/the/counts/file.tsv
# possible ID_class: "symbol", "ensemble_ID" or "entrez_ID"
# possible sep: "\t", ",", ";", " ")

