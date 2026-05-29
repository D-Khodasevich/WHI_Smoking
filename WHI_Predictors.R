### smoking status predictor
# meth: a DNA methylation matrix with CpGs as row and samples as columns.  
# predictors: list of predictors to calculate
### default setting is to estimate the "packyears", "cessation", and "SMOKE_3cat" predictors) 
# array: default is EPICv2 format (eg. probe names formatted as cg25324105_BC11 rather than cg25324105). 
### If using older arrays ("450K", "EPICv1"), specify this with the optional array argument.  

smoke_status <- function(meth, predictors = c("packyears", "cessation", "SMOKE_3cat"), 
                         array = "EPICv2"){
  
  preds <- as.data.frame(colnames(meth))
  colnames(preds) <- "SampleID"
  
  argue <- setdiff(predictors, c("packyears", "cessation", "SMOKE_3cat"))
  if(length(argue) > 0){
    message(paste0(argue, " not a valid predictor type, 
                   \n Options include packyears, cessation, and SMOKE_3cat"))
    stop()
  }
  
  if("packyears" %in% predictors){
    coefs <- read.csv("https://raw.githubusercontent.com/D-Khodasevich/WHI_Smoking/refs/heads/main/predictors/packyears_coefficients.csv")
    if(array != "EPICv2"){coefs$name <- gsub("\\_.*", "", coefs$name)}
    preds$packyears <- apply_model(coefs, as.matrix(meth))
  }
  if("cessation" %in% predictors){
    coefs <- read.csv("https://raw.githubusercontent.com/D-Khodasevich/WHI_Smoking/refs/heads/main/predictors/cessation_coefficients_nocur.csv")
    if(array != "EPICv2"){coefs$name <- gsub("\\_.*", "", coefs$name)}
    preds$cessation <- apply_model(coefs, as.matrix(meth))
  } 
  if("SMOKE_3cat" %in% predictors){
    coefs <- read.csv("https://raw.githubusercontent.com/D-Khodasevich/WHI_Smoking/refs/heads/main/predictors/3cat_coefficients.csv")
    if(array != "EPICv2"){coefs$name <- gsub("\\_.*", "", coefs$name)}
    never <- subset(coefs, coefs$class == "Never"); never$class <- NULL
    former <- subset(coefs, coefs$class == "Former"); former$class <- NULL
    current <- subset(coefs, coefs$class == "Current"); current$class <- NULL
    
    nev <- apply_model(never, as.matrix(meth))
    form <- apply_model(former, as.matrix(meth))
    cur <- apply_model(current, as.matrix(meth))
    
    pheno <- as.data.frame(cbind(nev, form, cur))
    pheno$SMOKE_3cat <- ifelse(pheno$nev > pheno$form & pheno$nev > pheno$cur, "0", 
                               ifelse(pheno$form > pheno$nev & pheno$form > pheno$cur, "1", 
                                      ifelse(pheno$cur > pheno$form & pheno$cur > pheno$form, "2", 
                                             "Inconclusive")))
    preds$SMOKE_3cat <- pheno$SMOKE_3cat
  } 
  return(preds)
}


### general application function for predictors generated with glmnet
##### coefficients: coefficients file generated using glmnet
##### meth: methylation input data
apply_model <- function(coefficients, meth){
  cpg_list <- coefficients$name
  meth <- meth[c(row.names(meth) %in% cpg_list), ]  # subset DNAm matrix to necessary CpGs
  present_cpgs <- row.names(meth)
  missing <- setdiff(cpg_list, present_cpgs)
  missing <- subset(missing, missing != "(Intercept)")
  
  # mean imputation
  meth <- t(meth)
  meth.na <- sum(is.na(meth))
  if(meth.na > 0){
    message("Missing Data: performing mean imputation to fill in missing data")
    toremove <- c()
    for(i in 1:ncol(meth)){
      if(sum(is.na(meth[, i])) > 0.5*nrow(meth)){
        toremove <- append(toremove, i)
      } else{
        if(sum(is.na(meth[, i])) > 0){
          meth[, i][is.na(meth[, i])] <- mean(meth[, i], na.rm = TRUE)
        }
      }
    }
    if(length(toremove) > 0){meth <- meth[, -toremove]}
  }
  
  int_val <- subset(coefficients, name == "(Intercept)")$coefficient 
  coefficients <- subset(coefficients, name != "(Intercept)")
  meth_sub <- as.data.frame(t(meth))
  meth_sub$name <- row.names(meth_sub)
  
  testdat <- dplyr::left_join(coefficients, meth_sub, by="name")
  
  testb <- testdat[,2]
  testx <- testdat[,3:length(testdat)]
  testval <- testb*testx
  testval <- colSums(testval, na.rm = TRUE)
  testval <- as.data.frame(testval)
  predict_val <- int_val + testval$testval 
  return(predict_val)
}
