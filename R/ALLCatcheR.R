#' @title Classifiction
#'
#' @description Prediction of B-ALL leukemia subtypes based on expression data
#' @export
#'
#' @param Counts.file count data
#' @param ID_class gene ids
#' @param sep file seperator
#' @return data.frame containing class predictions
#' @examples
#' outFun()
#'
outFun <- function(Counts.file=NA, ID_class="symbol", sep="\t") {
  # 1. preprocessing ############################################################
  # load count data, where the first column should be gene identifiers
  if(is.na(Counts.file)){
    Counts <- test_data
    cat("test counts loaded...\n")
  }else{
    Counts <- read.csv(Counts.file, sep = "\t", stringsAsFactors = F, row.names = 1)
    cat("counts loaded...\n")
  }
 
  if (length(rownames(Counts)) == length(which(rownames(Counts) == as.character(1:nrow(Counts))))) {
    stop("Error: symbol, ensemble or entrez are not provided in the first column")
  }
        
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
  
  # 2. classification using ML############################################################
  # load ML models ###############################################################
  #load("/media/tbeder/platte/models_L_Europe_StJude_comb_7030_26.RData")
  # predict 
  cat("ML prediction...\n")
  predsFin <- list()
  for (j in 1:6) {
    preds <- list()
    for (i in 1:length(models_L[[j]][[1]])) {
      preds[[i]] <- as.character(predict(models_L[[j]][[1]][[i]], Counts.norm))
    }
    predsFin[[j]] <- preds <- do.call("cbind",preds)
  }
  predsFin <- do.call("cbind",predsFin)
  preds <- predsFin
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
  Max <- as.numeric(apply(mat, 1, max)/60)
  Prediction <- colnames(mat)[apply(mat, 1, function(x) which.max(x))]
  
  cat("Classification using ssGSEA...\n")
  # 3. classification using ssGSEA ###########################################
  # calculate tpmps
  # load cds lengths
  cds_length <-  cds_length$cds_length[match(rownames(Counts), cds_length$gene)]
  # calculate mean cds length for missing values
  cds_length[is.na(cds_length)] <- mean(cds_length[is.na(cds_length)], na.rm = T)
  # calculate tpmps itself # https://support.bioconductor.org/p/91218/
  tpm3 <- function(counts,len) {
    x <- counts/len
    return(t(t(x)*1e6/colSums(x)))
  }
  tpms <- tpm3(Counts[,1:ncol(Counts)],cds_length)
  tpms <- as.data.frame(tpms)
  
  # rank data
  rankData <- singscore::rankGenes(tpms)
  
  # load gene sets used for ssGSEA and transform them to a list
  Genes_up <- list()
  for (i in 1:length(unique(GeneSets_up$subtype))) {
    Genes_up[[i]] <- GeneSets_up$gene[which(GeneSets_up$subtype == unique(GeneSets_up$subtype)[i])]
  }
  names(Genes_up) <- unique(GeneSets_up$subtype)
  
  Genes_dn <- list()
  for (i in 1:length(unique(GeneSets_dn$subtype))) {
    Genes_dn[[i]] <- GeneSets_dn$gene[which(GeneSets_dn$subtype == unique(GeneSets_dn$subtype)[i])]
  }
  names(Genes_dn) <- unique(GeneSets_dn$subtype)
  
  scoredf <- singscore::simpleScore(rankData, upSet = Genes_up[[1]], downSet = Genes_dn[[1]])
  head(scoredf)
  TotalScore <- scoredf[,1, drop = FALSE]
  
  for (i in 1:length(Genes_up)) {
    scoredf <- singscore::simpleScore(rankData, upSet = Genes_up[[i]], downSet = Genes_dn[[i]])
    TotalScore[,i] <- scoredf[,1]
  }
  
  # assign pathways names to data
  colnames(TotalScore) <- names(Genes_up)
  head(TotalScore)
  
  # scale total enrichment scores
  TotalScore_scaled <- as.data.frame(t(apply(TotalScore, 1, scale)))
  for (i in 1:nrow(TotalScore_scaled)) {
    TotalScore_scaled[i,] <- range01(TotalScore_scaled[i,])
  }
  colnames(TotalScore_scaled) <- colnames(TotalScore)
  
  # load training data with known subtype assignments
  colnames(trainData)
  colnames(TotalScore_scaled)
  TrainingSubtypes <- trainData$subtype
  rownames(trainData) <- paste0("train_", rownames(trainData))
  names(TrainingSubtypes) <- rownames(trainData)
  
  trainData <- trainData[,match(colnames(TotalScore_scaled), colnames(trainData))]
  #################
  TotalScore_scaled_ALL <- rbind(TotalScore_scaled,trainData)
  DIST_ALL <- dist(TotalScore_scaled_ALL)
  DIST_ALL <- as.data.frame(as.matrix(DIST_ALL))
  
  trainPos <- (nrow(TotalScore_scaled)+1):ncol(DIST_ALL)
  testPos <- 1:nrow(TotalScore_scaled)

  # calculate ssGSEA NN results and the final prediction score 
  cat("combine ML and ssGSEA...\n")
  ML_KNN <- c()
  geneSetPreds_df <- as.data.frame(matrix(NA,ncol = length(sort(unique(TrainingSubtypes)))))
  colnames(geneSetPreds_df) <- sort(unique(TrainingSubtypes))
  i <- 4
  for (i in 1:length(testPos)) {
    NN10 <- colnames(DIST_ALL)[trainPos][order(as.numeric(DIST_ALL[testPos[i],trainPos]))][2:11]
    NN10 <- as.character(TrainingSubtypes[match(NN10, names(TrainingSubtypes))])
    for (j in 1:ncol(geneSetPreds_df)) {
      geneSetPreds_df[i,j] <- length(which(NN10 == colnames(geneSetPreds_df)[j]))
    }
    
    if (Prediction[i] == "Ph.pos") {
      ML_KNN[i] <- mean(c(Max[i], length(which(NN10 == "Ph.group"))/10
      ))
    }
    else if (Prediction[i] == "Ph.like") {
      ML_KNN[i] <- mean(c(Max[i], length(which(NN10 == "Ph.group"))/10
      ))
      
    }
    else if (Prediction[i] == "ETV6.RUNX1") {
      ML_KNN[i] <- mean(c(Max[i], length(which(NN10 == "ETV6.RUNX1.group"))/10
      ))
      
    }
    else if (Prediction[i] == "ETV6.RUNX1.like") {
      ML_KNN[i] <- mean(c(Max[i], length(which(NN10 == "ETV6.RUNX1.group"))/10
      ))
      
    } else {
      ML_KNN[i] <- mean(c(Max[i], length(which(NN10 == Prediction[i]))/10
   ))
    }
  }
  geneSet_preds <- colnames(geneSetPreds_df)[apply(geneSetPreds_df, 1, which.max)]
  # make final preditions 
  final_pred <- rep("Unclassified",nrow(Counts.norm))
  # load ML_KNN cutoffs for each subtype
  for (i in 1:nrow(Counts.norm)) {
    if(ML_KNN[i] >= ML_KNN_cutoffs$cutoffs[match(Prediction[i], ML_KNN_cutoffs$subtype)])
      final_pred[i] <- Prediction[i]
  }
  table(final_pred)
  
  # 4. generate output 
  colnames(mat) <- paste0(colnames(mat), "_ML")
  colnames(geneSetPreds_df) <- paste0(colnames(geneSetPreds_df), "_geneSets")
  
  output <- cbind("sample" = rownames(Counts.norm),
        mat/60,
        "ML_prediction" = Prediction,
        geneSetPreds_df,
        "GeneSet_prediction" = geneSet_preds,
        score = ML_KNN,
        final_prediction = final_pred)
  
  cat("predictions saved in:", getwd())
  # save predictions
  cat("Writing output file:",paste0(getwd(), "/predictions.tsv"),"...\n")
  write.table(output,"predictions.tsv", sep = "\t", row.names = F)
  return(invisible(output))
}
