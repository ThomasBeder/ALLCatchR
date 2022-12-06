globalVariables(c("test_data","models_20","NH","BC_model_GMALL","BC_model_MLL","models_L_sex","models_L_Immuno"))

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
#' allcatch()
#'

allcatch <- function(Counts.file=NA, ID_class="symbol", sep="\t") {
  # 1. preprocessing ############################################################
  # load count data, where the first column should be gene identifiers
  if(is.na(Counts.file)){
    Counts <- test_data
    cat("test counts loaded...\n")
  }else{
    Counts <- utils::read.csv(Counts.file, sep = sep, stringsAsFactors = F, row.names = 1)
    cat("counts loaded...\n")
  }
  
  if (length(rownames(Counts)) == length(which(rownames(Counts) == as.character(1:nrow(Counts))))) {
    stop("Error: symbol, ensemble or entrez are not provided in the first column")
  }
    ID_conv <- ID_conversion
  # select the genes used for classifier trainig
  ma <- match(ID_conv[,match(ID_class, colnames(ID_conv))], rownames(Counts))
  Counts <- Counts[ma[!is.na(ma)],]
  
  # convert to symbol (classifier was trained on symbols)
  ma <- match(rownames(Counts), ID_conv[,match(ID_class, colnames(ID_conv))])
  Counts <- Counts[!is.na(ma),]
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
  # predict 
  i <- 1
  cat("ML prediction...\n")
  preds <- list()
  for (i in 1:10) {
    preds[[i]] <- as.character(stats::predict(models_20[[i]], Counts.norm))
  }
  predsFin20 <- do.call("cbind",preds)
  
  preds <- list()
  for (i in 1:10) {
    preds[[i]] <- as.character(stats::predict(NH[[i]], Counts.norm))
  }
  predsFinNH <- do.call("cbind",preds)
  
  rownames(predsFin20) <- rownames(Counts.norm)
  mat20 <- matrix(ncol =  21)
  colnames(mat20) <- cutoffs$class[1:21]
  mat20 <- as.data.frame(mat20)
  # count number of machines predicting a subtype 
  types <- ncol(mat20)
  for (i in 1:nrow(predsFin20)) {
    TAB <- table(as.character(predsFin20[i,]))
    mat20[i,] <- TAB[match(colnames(mat20),names(TAB))]
  }
  mat20[is.na(mat20)] <- 0
  # count maximal number of machines predicting a subtype
  mat20 <- mat20/10
  rownames(mat20) <- rownames(Counts.norm)
  Max20 <- apply(mat20, 1, max)
  Prediction20 <- colnames(mat20)[apply(mat20, 1, function(x) which.max(x))]
  #############################################################################
  rownames(predsFinNH) <- rownames(Counts.norm)
  matNH <- matrix(ncol =  2)
  colnames(matNH) <- c("No_NH","Near.haploid")
  matNH <- as.data.frame(matNH)
  # count number of machines predicting a subtype 
  types <- ncol(matNH)
  for (i in 1:nrow(predsFinNH)) {
    TAB <- table(as.character(predsFinNH[i,]))
    matNH[i,] <- TAB[match(colnames(matNH),names(TAB))]
  }
  matNH[is.na(matNH)] <- 0
  # count maximal number of machines predicting a subtype
  matNH <- matNH/10
  rownames(matNH) <- rownames(Counts.norm)
  MaxNH <- apply(matNH, 1, max)
  PredictionNH <- colnames(matNH)[apply(matNH, 1, function(x) which.max(x))]

  mat20$Near.haploid <- 0
  
 if (length(which(PredictionNH == "Near.haploid" & Prediction20 == "NH.HeH")) > 0 ) {
    ma <- which(PredictionNH == "Near.haploid" & Prediction20 == "NH.HeH")
    mat20[ma,which(colnames(mat20) == "NH.HeH")] <- 0
    mat20[ma,which(colnames(mat20) == "Near.haploid")] <- MaxNH[ma]
    Prediction20[ma] <- "Near.haploid"
    }
  
  cat("Classification using ssGSEA...\n")
  # 3. classification using ssGSEA ###########################################
  # calculate tpmps
  # load cds lengths
  cds_lengthNum <-  cds_length$cds_length[match(rownames(Counts), cds_length$gene)]
  # calculate mean cds length for missing values
  cds_lengthNum[is.na(cds_lengthNum)] <- mean(cds_lengthNum[is.na(cds_lengthNum)], na.rm = T)
  # calculate tpmps itself # https://support.bioconductor.org/p/91218/
  tpm3 <- function(counts,len) {
    x <- counts/len
    return(t(t(x)*1e6/colSums(x)))
  }
  tpms <- tpm3(Counts[,1:ncol(Counts)],cds_lengthNum)
  tpms <- as.data.frame(tpms)
  
  # rank data
  rankData <- singscore::rankGenes(tpms)
  
  # load gene sets used for ssGSEA and transform them to a list
  Genes_up <- Genes_up20
  Genes_dn <- Genes_dn20
  scoredf <-  suppressWarnings({singscore::simpleScore(rankData, upSet = Genes_up[[1]], downSet = Genes_dn[[1]])})
  utils::head(scoredf)
  TotalScore <- scoredf[,1, drop = FALSE]
  
  for (i in 1:length(Genes_up)) {
    scoredf <- suppressWarnings({singscore::simpleScore(rankData, upSet = Genes_up[[i]], downSet = Genes_dn[[i]])})
    TotalScore[,i] <- scoredf[,1]
  }
  
  # assign pathways names to data
  colnames(TotalScore) <- names(Genes_up)
  utils::head(TotalScore)
  
  # scale total enrichment scores
  TotalScore_scaled <- as.data.frame(t(apply(TotalScore, 1, scale)))
  for (i in 1:nrow(TotalScore_scaled)) {
    TotalScore_scaled[i,] <- range01(TotalScore_scaled[i,])
  }
  colnames(TotalScore_scaled) <- colnames(TotalScore)
  
  # load training data with known subtype assignments
  trainData <- TotalScore_scaled_training
  colnames(trainData)
  table(trainData$subtype)
  TrainingSubtypes <- trainData$subtype
  rownames(trainData) <- paste0("train_", rownames(trainData))
  names(TrainingSubtypes) <- rownames(trainData)
  trainData <- trainData[,match(colnames(TotalScore_scaled), colnames(trainData))]
  #################
  TotalScore_scaled_ALL <- rbind(TotalScore_scaled,trainData)
  DIST_ALL <- stats::dist(TotalScore_scaled_ALL)
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
    ML_KNN[i] <- mean(c(Max20[i], length(which(NN10 == Prediction20[i]))/10
      ))
    }
  geneSetPreds_df <- geneSetPreds_df/10
  geneSet_preds <- colnames(geneSetPreds_df)[apply(geneSetPreds_df, 1, which.max)]
  colnames(geneSetPreds_df) <- paste0("NN_",colnames(geneSetPreds_df))
  colnames(mat20) <- paste0("ML_",colnames(mat20))
 
################################################################################   
##### NN classifer NH ##########################################################
################################################################################
  # load gene sets used for ssGSEA and transform them to a list
  Genes_up <- Genes_upNH
  Genes_dn <- Genes_dnNH
  length(Genes_up)
  length(Genes_dn)
  scoredf <- suppressWarnings({singscore::simpleScore(rankData, upSet = Genes_up[[1]], downSet = Genes_dn[[1]])})
  utils::head(scoredf)
  TotalScore <- scoredf[,1, drop = FALSE]
  
  for (i in 1:length(Genes_up)) {
    scoredf <- suppressWarnings({singscore::simpleScore(rankData, upSet = Genes_up[[i]], downSet = Genes_dn[[i]])})
    TotalScore[,i] <- scoredf[,1]
  }
  
  # assign pathways names to data
  colnames(TotalScore) <- names(Genes_up)
  utils::head(TotalScore)
  
  # scale total enrichment scores
  TotalScore_scaled <- as.data.frame(t(apply(TotalScore, 1, scale)))
  for (i in 1:nrow(TotalScore_scaled)) {
    TotalScore_scaled[i,] <- range01(TotalScore_scaled[i,])
  }
  colnames(TotalScore_scaled) <- colnames(TotalScore)
  
  # load training data with known subtype assignments
  trainData <- TotalScore_scaled_training_NH
  colnames(trainData)
  table(trainData$subtype)
  TrainingSubtypes <- trainData$subtype
  rownames(trainData) <- paste0("train_", rownames(trainData))
  names(TrainingSubtypes) <- rownames(trainData)
  trainData <- trainData[,match(colnames(TotalScore_scaled), colnames(trainData))]
  #################
  TotalScore_scaled_ALL <- rbind(TotalScore_scaled,trainData)
  DIST_ALL <- stats::dist(TotalScore_scaled_ALL)
  DIST_ALL <- as.data.frame(as.matrix(DIST_ALL))
  
  trainPos <- (nrow(TotalScore_scaled)+1):ncol(DIST_ALL)
  testPos <- 1:nrow(TotalScore_scaled)
  
  # calculate ssGSEA NN results and the final prediction score 
  geneSetPreds_df_NH <- as.data.frame(matrix(NA,ncol = length(sort(unique(TrainingSubtypes)))))
  colnames(geneSetPreds_df_NH) <- sort(unique(TrainingSubtypes))
  i <- 4
  for (i in 1:length(testPos)) {
    NN10 <- colnames(DIST_ALL)[trainPos][order(as.numeric(DIST_ALL[testPos[i],trainPos]))][2:11]
    NN10 <- as.character(TrainingSubtypes[match(NN10, names(TrainingSubtypes))])
    for (j in 1:ncol(geneSetPreds_df_NH)) {
      geneSetPreds_df_NH[i,j] <- length(which(NN10 == colnames(geneSetPreds_df_NH)[j]))
    }
  }
  
  geneSetPreds_df_NH <- geneSetPreds_df_NH/10
  geneSet_preds_NH <- colnames(geneSetPreds_df_NH)[apply(geneSetPreds_df_NH, 1, which.max)]
  colnames(geneSetPreds_df_NH) <- paste0("NN_",colnames(geneSetPreds_df_NH))
  
  utils::head(geneSetPreds_df)
  geneSetPreds_df$NN_Near.haploid <- 0
  
  if (length(which(Prediction20 == "Near.haploid")) > 0 ) {
    ma <- which(Prediction20 == "Near.haploid")
    geneSetPreds_df$NN_Near.haploid[ma] <- geneSetPreds_df_NH$`NN_Near haploid`[ma]
    for (i in 1:length(ma)) {
    ML_KNN[ma] <- mean(c(mat20$ML_ML_Near.haploid[ma[i]], geneSetPreds_df$NN_Near.haploid[ma[i]]))
    }
  }
  
  tier <- c()
  i <- 1
  for (i in 1:nrow(mat20)) {
    ma <- match(Prediction20[i], cutoffs$class)
    if (ML_KNN[i] >= cutoffs$upper[ma]) {
      tier[i] <- "high-confidence"
    }
    if (ML_KNN[i] < cutoffs$upper[ma] & 
        ML_KNN[i] >= cutoffs$lower[ma]) {
      tier[i] <- "candidate"
    }
    if (ML_KNN[i] < cutoffs$lower[ma]) {
      tier[i] <- "unclassified"
    }
  }
table(tier)

################################################################################  
#### blast count prediction ####################################################  
################################################################################
  if(is.na(Counts.file)){
    Counts <- test_data
    cat("test counts loaded...\n")
  }else{
    Counts <- utils::read.csv(Counts.file, sep = sep, stringsAsFactors = F, row.names = 1)
  }
    
  if (length(rownames(Counts)) == length(which(rownames(Counts) == as.character(1:nrow(Counts))))) {
    #    stop("Error: symbol, ensemble or entrez are not provided in the first column")
  }
  ID_conv <- ID_conversion_BC
  # select the genes used for classifier trainig
  ma <- match(ID_conv[,match(ID_class, colnames(ID_conv))], rownames(Counts))
  Counts <- Counts[ma[!is.na(ma)],]
  
  # convert to symbol (classifier was trained on symbols)
  ma <- match(rownames(Counts), ID_conv[,match(ID_class, colnames(ID_conv))])
  Counts <- Counts[!is.na(ma),]
  ma <- match(rownames(Counts), ID_conv[,match(ID_class, colnames(ID_conv))])
  rownames(Counts) <- ID_conv$symbol[ma]
  
  # normalize data and scale between 0 and 1
  Counts.norm <- Counts+1
  Counts.norm <- apply(Counts.norm, 2, log10)
  Counts.norm <- apply(Counts.norm, 2, scale)
  #range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  #for (i in 1:ncol(Counts.norm)) {
  #  Counts.norm[,i] <- range01(Counts.norm[,i])
  #}
  
  # transpose data
  Counts.norm <- as.data.frame(t(Counts.norm))
  colnames(Counts.norm) <- rownames(Counts)
  colnames(Counts.norm) <- make.names(colnames(Counts.norm))
  
  # find genes not provided by user
  ma <- match(ID_conv$symbol, rownames(Counts))
  GenesNoFound <- ID_conv$symbol[is.na(ma)]
  
  # Print number of missing genes
  #cat(paste0(length(GenesNoFound)), " of ", nrow(ID_conv), " genes not found\n")
  if ((length(GenesNoFound)) > 0) {
    #  cat("impute missing genes...\n")
  }
  
  # impute missing genes
  GenesNoFound_df <- matrix(ID_conv$norm_exp[match(GenesNoFound, ID_conv$symbol)], nrow = nrow(Counts.norm), ncol = length(GenesNoFound))
  colnames(GenesNoFound_df) <- GenesNoFound
  Counts.norm <- cbind(Counts.norm, GenesNoFound_df)
  colnames(Counts.norm) <- make.names(colnames(Counts.norm))

  preds <- list()
  for (i in 1:10) {
    preds[[i]] <- stats::predict(BC_model_GMALL[[i]], Counts.norm)
  }
  preds <- do.call("cbind",preds)
  preds1 <- apply(preds, 1, mean)
  preds <- list()
  for (i in 1:10) {
    preds[[i]] <- stats::predict(BC_model_MLL[[i]], Counts.norm)
  }
  preds <- do.call("cbind",preds)
  preds2 <- apply(preds, 1, mean)
  predsBC <- apply(cbind(preds1,preds2), 1, mean)
  
################################################################################  
#### sex prediction ############################################################  
################################################################################  
  cat("Patient's sex prediction...\n")
  # 1. preprocessing ############################################################
  # load count data, where the first column should be gene identifiers
  if(is.na(Counts.file)){
    Counts <- test_data
#    cat("test counts loaded...\n")
  }else{
    Counts <- utils::read.csv(Counts.file, sep = sep, stringsAsFactors = F, row.names = 1)
#    cat("counts loaded...\n")
  }
  
  if (length(rownames(Counts)) == length(which(rownames(Counts) == as.character(1:nrow(Counts))))) {
#    stop("Error: symbol, ensemble or entrez are not provided in the first column")
  }
  ID_conv <- ID_conversion_sex
  # select the genes used for classifier trainig
  ma <- match(ID_conv[,match(ID_class, colnames(ID_conv))], rownames(Counts))
  Counts <- Counts[ma[!is.na(ma)],]
  
  # convert to symbol (classifier was trained on symbols)
  ma <- match(rownames(Counts), ID_conv[,match(ID_class, colnames(ID_conv))])
  Counts <- Counts[!is.na(ma),]
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
  #cat(paste0(length(GenesNoFound)), " of ", nrow(ID_conv), " genes not found\n")
  if ((length(GenesNoFound)) > 0) {
  #  cat("impute missing genes...\n")
  }
  
  # impute missing genes
  GenesNoFound_df <- matrix(ID_conv$norm_exp[match(GenesNoFound, ID_conv$symbol)], nrow = nrow(Counts.norm), ncol = length(GenesNoFound))
  colnames(GenesNoFound_df) <- GenesNoFound
  Counts.norm <- cbind(Counts.norm, GenesNoFound_df)
  colnames(Counts.norm) <- make.names(colnames(Counts.norm))
  
  preds <- list()
  for (i in 1:length(models_L_sex)) {
    preds[[i]] <- as.character(stats::predict(models_L_sex[[i]], Counts.norm))
  }
  
  preds <- do.call("cbind",preds)
  
  utils::head(preds)
  rownames(preds) <- rownames(Counts.norm)
  mat <- matrix(ncol =  2)
  colnames(mat) <- c("F","M")
  mat <- as.data.frame(mat)
  types <- ncol(mat)
  for (i in 1:nrow(preds)) {
    TAB <- table(as.character(preds[i,]))
    mat[i,] <- TAB[match(colnames(mat),names(TAB))]
  }
  mat[is.na(mat)] <- 0
  mat <- mat/10
  MaxSex <- apply(mat, 1, max)
  PredictionSex <- colnames(mat)[apply(mat, 1, function(x) which.max(x))]
  PredictionSex[which(MaxSex<0.5)] <- "Unclassified"

################################################################################  
#### immunophenotype prediction #################################################  
################################################################################  
  preds <- list()
  for (i in 1:length(models_L_Immuno)) {
    preds[[i]] <- as.character(stats::predict(models_L_Immuno[[i]], Counts.norm))
  }
  preds <- do.call("cbind",preds)
  utils::head(preds)
  rownames(preds) <- rownames(Counts.norm)
  mat <- matrix(ncol =  2)
  colnames(mat) <- c("pro_B", "common_B")
  mat <- as.data.frame(mat)
  types <- ncol(mat)
  for (i in 1:nrow(preds)) {
    TAB <- table(as.character(preds[i,]))
    mat[i,] <- TAB[match(colnames(mat),names(TAB))]
  }
  mat[is.na(mat)] <- 0
  mat <- mat/10
  MaxImmuno <- apply(mat, 1, max)
  PredictionImmuno <- colnames(mat)[apply(mat, 1, function(x) which.max(x))]
  PredictionImmuno[which(MaxImmuno<0.5)] <- "Unclassified"  
  
################################################################################
###### progenitor ssGSEA #######################################################
################################################################################
cat("assign putative progenitor...", getwd(),"\n")
    dim(Counts)
  ma <- match(genesMini,rownames(Counts))
  ma <- ma[!is.na(ma)]                                        
  Counts <- Counts[ma,]
  cds_lengthNum <-  cds_length$cds_length[match(rownames(Counts), cds_length$gene)]
  # calculate mean cds length for missing values
  cds_lengthNum[is.na(cds_lengthNum)] <- mean(cds_lengthNum[is.na(cds_lengthNum)], na.rm = T)
  # calculate tpmps itself # https://support.bioconductor.org/p/91218/
  tpms <- tpm3(Counts[,1:ncol(Counts)],cds_lengthNum)
  tpms <- as.data.frame(tpms)
  
  # rank data
  rankData <- singscore::rankGenes(tpms)
  
  # load gene sets used for ssGSEA and transform them to a list
  Genes_up <- Genes_upMini
  Genes_dn <- Genes_dnMini
  scoredf <- suppressWarnings({singscore::simpleScore(rankData, upSet = Genes_up[[1]], downSet = Genes_dn[[1]])})
  utils::head(scoredf)
  TotalScore <- scoredf[,1, drop = FALSE]
  
  for (i in 1:length(Genes_up)) {
    scoredf <- suppressWarnings({singscore::simpleScore(rankData, upSet = Genes_up[[i]], downSet = Genes_dn[[i]])})
    TotalScore[,i] <- scoredf[,1]
  }
  
  # assign pathways names to data
  colnames(TotalScore) <- names(Genes_up)
  TotalScore <- TotalScore[c(1,7,4,5,6,2,3)]
  utils::head(TotalScore)
  
  # scale total enrichment scores
  TotalScore_scaled <- as.data.frame(t(apply(TotalScore, 1, scale)))
  for (i in 1:nrow(TotalScore_scaled)) {
    TotalScore_scaled[i,] <- range01(TotalScore_scaled[i,])
  }
  colnames(TotalScore_scaled) <- colnames(TotalScore)
  utils::head(TotalScore_scaled)
  
################################################################################  
# 4. generate output ###########################################################
################################################################################  
  
output <- cbind(sample = rownames(mat20),
                  Score = ML_KNN,
                  Prediction = Prediction20, 
                  Confidence = tier,
                  BlastCounts = predsBC,
                  Sex = PredictionSex,
                  Score_sex = MaxSex,
                  Immuno = PredictionImmuno,
                  ScoreImmuno = MaxImmuno,
                  TotalScore_scaled,
                  mat20, 
                  geneSetPreds_df)

table(output$Prediction)  
output$Prediction <- cutoffs$subtype[match(output$Prediction,cutoffs$class)]
ML_cols <- grep("ML_", colnames(output))
colnames(output)[ML_cols] <- paste0("ML_",cutoffs$subtype[match(gsub("ML_","", colnames(output)[ML_cols]), cutoffs$class)])
NN_cols <- grep("NN_", colnames(output))
colnames(output)[NN_cols] <- paste0("NN_",cutoffs$subtype[match(gsub("NN_","", colnames(output)[NN_cols]), cutoffs$class)])

cat("predictions saved in:", getwd(),"\n")
  # save predictions
  cat("Writing output file:",paste0(getwd(), "/predictions.tsv"),"...\n")
  utils::write.table(output,"predictions.tsv", sep = "\t", row.names = F)
  return(invisible(output))
}
