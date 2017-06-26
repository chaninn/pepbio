# Rnadom Forest models with 
library(protr)
cancer <- readFASTA("cancer.fasta")
fungus <- readFASTA("fungus.fasta")
virus <- readFASTA("virus.fasta")
bacteria <- readFASTA("bacteria.fasta")
negative <- readFASTA("negative.fasta")


cancer <- cancer[(sapply(cancer, protcheck))]
fungus <- fungus[(sapply(fungus, protcheck))]
bacteria <- bacteria[(sapply(bacteria, protcheck))]
virus <- virus[(sapply(virus, protcheck))]
negative <- negative[(sapply(negative, protcheck))]
## Descriptors generation
cancer_des <- t(sapply(cancer, extractAAC))

fungus_des <- t(sapply(fungus, extractAAC))

bacteria_des <- t(sapply(bacteria, extractAAC))

virus_des <- t(sapply(virus, extractAAC))

negative_des <- t(sapply(negative, extractAAC))


#### label the labels 
cancer_des <- as.data.frame(cancer_des)
cancer_des$Label <- "Cancer"
fungus_des <- as.data.frame(fungus_des)
fungus_des$Label <- "Fungus"
bacteria_des <- as.data.frame(bacteria_des)
bacteria_des$Label <- "Bacteria"
virus_des <- as.data.frame(virus_des)
virus_des$Label <- "Virus"
negative_des <- as.data.frame(negative_des)
negative_des$Label <- "Negative"

cancer <- rbind(cancer_des, negative_des)
cancer$Label <- as.factor(cancer$Label)
fungus <- rbind(fungus_des, negative_des)
fungus$Label <- as.factor(fungus$Label)
bacteria <- rbind(fungus_des, negative_des)
bacteria$Label <- as.factor(bacteria$Label)
virus <- rbind(virus_des, negative_des)
virus$Label <- as.factor(virus$Label)

cancer <- na.omit(cancer)
fungus <- na.omit(fungus)
bacteria <- na.omit(bacteria)
virus <- na.omit(virus)

cancer <- data.frame(cancer)
fungus <- data.frame(fungus)
bacteria <- data.frame(bacteria)
virus <- data.frame(virus)



input <- list(cancer = cancer, fungus = fungus, bacteria = bacteria, 
              virus = virus)
mean_and_sd <- function(x) {
  c(round(mean(x, na.rm = TRUE), digits = 4),
    round(sd(x, na.rm = TRUE), digits = 4))
}

rf_training <- function(x) {
  results <- list(100)
  for (i in 1:100) {
    in_train <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
    train <- x[in_train, ]
    test <- x[-in_train, ]
    model_train <- ranger::ranger(Label~., data = train, write.forest = TRUE, save.memory = TRUE)
    prediction <- predict(model_train, train)
    prediction <- prediction$predictions
    actual <- train$Label
    result <- caret::confusionMatrix(prediction, actual)
    result <- result$table
    result <- as.numeric(result)
    results[[i]] <- result
  }
  return(results)
}

rf_CV <- function(x) {
  library(ranger)
  
  results <- list(100)
  for (i in 1:100) {
    in_train <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
    myData <- x[in_train, ]
    test <- x[-in_train, ]
    
    k = 10
    index <- sample(1:k, nrow(myData), replace = TRUE)
    folds <- 1:k
    myRes <- data.frame()
    for (j in 1:k) {
      training <- subset(myData, index %in% folds[-j])
      testing <- subset(myData, index %in% c(j))
      model <- ranger::ranger(Label~., data = training, write.forest = TRUE, 
                              save.memory = TRUE)
      prediction <- predict(model, testing)
      prediction <- prediction$predictions
      actual <- testing$Label
      data <- data.frame(prediction, actual)
      myRes <- rbind(myRes, data)
    } 
    
    prediction <- myRes$prediction
    actual <- myRes$actual
    result <- caret::confusionMatrix(prediction, actual)
    result <- result$table
    result <- as.numeric(result)
    results[[i]] <- result
  }
  return(results)
}

rf_testing <- function(x) {
  results <- list(100)
  for (i in 1:100) {
    in_train <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
    train <- x[in_train, ]
    test <- x[-in_train, ]
    model_train <- ranger::ranger(Label~., data = train, write.forest = TRUE, save.memory = TRUE)
    prediction <- predict(model_train, test)
    prediction <- prediction$predictions
    actual <- test$Label
    result <- caret::confusionMatrix(prediction, actual)
    result <- result$table
    result <- as.numeric(result)
    results[[i]] <- result
  }
  return(results)
}


assessment <- function(x) {
  results <- data.frame(x)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  return(results_all)
}

rf_modelling <- function(x, type){
  if (type == "Train") {
    ok <- rf_training(x)
    results <- assessment(ok)
  } else if (type == "CV") {
    ok <- rf_CV(x)
    results <- assessment(ok)
    
  } else if (type == "Test") {
    ok <- rf_testing(x)
    results <- assessment(ok)
  }
  return(results)
}

results_RF <- function(x) {
  training <- rf_modelling(x, type = "Train")
  cross_validation <- rf_modelling(x, type = "CV")
  testing <- rf_modelling(x, type = "Test")
  results <- data.frame(Training = training, Cross_VAlidation = cross_validation, Testing = testing)
  colnames(results) <- c("Training set", "Internal Validatin", "External Set")
  return(results)
}

result_RF_performance <- lapply(input, function(x) {
  models <- results_RF(x)
  return(models)
})



####################HDP Binary models

library(protr)
library(protr)
cancer <- readFASTA("cancer.fasta")
fungus <- readFASTA("fungus.fasta")
virus <- readFASTA("virus.fasta")
bacteria <- readFASTA("bacteria.fasta")
negative <- readFASTA("negative.fasta")


cancer <- cancer[(sapply(cancer, protcheck))]
fungus <- fungus[(sapply(fungus, protcheck))]
bacteria <- bacteria[(sapply(bacteria, protcheck))]
virus <- virus[(sapply(virus, protcheck))]
negative <- negative[(sapply(negative, protcheck))]
## Descriptors generation
cancer_des <- t(sapply(cancer, extractAAC))

fungus_des <- t(sapply(fungus, extractAAC))

bacteria_des <- t(sapply(bacteria, extractAAC))

virus_des <- t(sapply(virus, extractAAC))

negative_des <- t(sapply(negative, extractAAC))

#### label the labels 
cancer_des <- as.data.frame(cancer_des)
fungus_des <- as.data.frame(fungus_des)
bacteria_des <- as.data.frame(bacteria_des)
virus_des <- as.data.frame(virus_des)
hdp <- rbind(cancer_des, fungus_des, bacteria_des, virus_des)
hdp$Label <- "Positive"
hdp <- as.data.frame(hdp)
negative_des <- as.data.frame(negative_des)
negative_des$Label <- "Negative"

data <- rbind(hdp, negative_des)

data <- data.frame(data)
data <- na.omit(data)
data$Label <- as.factor(data$Label)


results_RF_Positive_Negative <- function(x) {
  training <- rf_modelling(x, type = "Train")
  cross_validation <- rf_modelling(x, type = "CV")
  testing <- rf_modelling(x, type = "Test")
  results <- data.frame(Training = training, Cross_VAlidation = cross_validation, Testing = testing)
  colnames(results) <- c("Training set", "Internal Validatin", "External Set")
  return(results)
}

results_positive_negative <- results_RF_Positive_Negative(data)


#################Multiple class

library(protr)
cancer <- readFASTA("cancer.fasta")
fungus <- readFASTA("fungus.fasta")
bacteria <- readFASTA("bacteria.fasta")
virus <- readFASTA("virus.fasta")
negative <- protr::readFASTA("negative.fasta")

### removed wired protein
cancer <- cancer[(sapply(cancer, protcheck))]
fungus <- fungus[(sapply(fungus, protcheck))]
bacteria <- bacteria[(sapply(bacteria, protcheck))]
virus <- virus[(sapply(virus, protcheck))]
negative <- negative[(sapply(negative, protcheck))]


### generating the Descriptors
cancer_des <- t(sapply(cancer, extractAAC))

fungus_des <- t(sapply(fungus, extractAAC))

bacteria_des <- t(sapply(bacteria, extractAAC))

virus_des <- t(sapply(virus, extractAAC))

negative_des <- t(sapply(negative, extractAAC))

#### label the labels 
cancer_des <- as.data.frame(cancer_des)
cancer_des$Label <- "Cancer"
fungus_des <- as.data.frame(fungus_des)
fungus_des$Label <- "Fungus"
bacteria_des <- as.data.frame(bacteria_des)
bacteria_des$Label <- "Bacteria"
virus_des <- as.data.frame(virus_des)
virus_des$Label <- "Virus"
combine_data <- rbind(cancer_des, fungus_des,
                      bacteria_des, virus_des)
combine_data$Label <- as.factor(combine_data$Label)

combine_data <- data.frame(combine_data)
combine_data <- na.omit(combine_data)
RF_training <- function(x, Label){
  
  
  
  ok <- list(100)
  for (i in 1:100) { 
    in_train <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
    train <- x[in_train, ]
    test <- x[-in_train, ]
    rm(test)
    model_train <- ranger::ranger(Label~., data = train, write.forest = TRUE, save.memory = TRUE)
    prediction <- predict(model_train, train)
    prediction <- prediction$predictions
    actual <- train$Label
    result <- caret::confusionMatrix(prediction, actual)
    result <- result$table
    results <- as.numeric(result)
    rm(model_train)
    rm(prediction)
    rm(actual)
    results <- as.numeric(results)
    if (Label == "Bacteria") {
      ok[[i]] <- cbind(results[[1]], (results[[5]] + results[[9]] + results[[13]]),
                       (results[[2]] + results[[3]] + results[[4]]), (results[[6]] + results[[11]] + results[[16]]))
      
      
    }  else if (Label == "Cancer") {
      
      
      ok[[i]] <- cbind(results[[6]], (results[[2]] + results[[10]] + results[[14]]), 
                       (results[[5]] + results[[7]] + results[[8]]), (results[[1]] + results[[11]] + results[[16]]))
      
      
    }  else if (Label == "Fungus") {
      
      
      ok[[i]] <- cbind(results[[11]], (results[[3]] + results[[7]] + results[[15]]),
                       (results[[9]] + results[[10]] + results[[12]]), (results[[1]] + results[[6]] + results[[16]]))
      
      
    }  else if (Label == "Virus") {
      
      
      ok[[i]] <- cbind(results[[16]], (results[[4]] + results[[8]] + results[[12]]), 
                       (results[[13]] + results[[14]] + results[[15]]), (results[[1]] + results[[6]] + results[[11]]))    
    } }
  return(ok)
  stopCluster(cl) } 


mean_and_sd <- function(x) {
  c(round(mean(x, na.rm = TRUE), digits = 4),
    round(sd(x, na.rm = TRUE), digits = 4))
}

assessment <- function(x) {
  results <- data.frame(x)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  return(results_all)
}

rf_training_all <- function(x, type){
  if (type == "Bacteria") {
    ok <- RF_training(x, Label = "Bacteria")
    results <- assessment(ok)
  } else if (type == "Cancer") {
    ok <- RF_training(x, Label = "Cancer")
    results <- assessment(ok)
    
  } else if (type == "Fungus") {
    ok <- RF_training(x, Label = "Fungus")
    results <- assessment(ok)
  } else if (type == "Virus") {
    ok <- RF_training(x, Label = "Virus")
    results <- assessment(ok)
  }
  return(results)
}


RF_training_all <- function(x) {
  bacteria <- rf_training_all(x, type = "Bacteria")
  cancer <- rf_training_all(x, type = "Cancer")
  fungus <- rf_training_all(x, type = "Fungus")
  virus <- rf_training_all(x, type = "Virus")
  results_all <- cbind(bacteria, cancer, fungus, virus)
  rm(bacteria)
  rm(cancer)
  rm(fungus)
  rm(virus)
  total <- apply(results_all, 1, mean)
  results_all_mean <- cbind(results_all, total)
  rm(results_all)
  rm(total)
  colnames(results_all_mean) <- c("Bacteria", "Cancer", "Fungus", "Virus", "Overall")
  return(results_all_mean)
}





###############Multiple models
RF_10_CV <- function(x, Label){
  ok <- list(100)
  for (i in 1:100) { 
    in_train <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
    myData <- x[in_train, ]
    test <- x[-in_train, ]
    
    k = 10
    index <- sample(1:k, nrow(myData), replace = TRUE)
    folds <- 1:k
    myRes <- data.frame()
    for (j in 1:k) {
      training <- subset(myData, index %in% folds[-j])
      testing <- subset(myData, index %in% c(j))
      model <- ranger::ranger(Label~., data = training, write.forest = TRUE, 
                              save.memory = TRUE)
      prediction <- predict(model, testing)
      prediction <- prediction$predictions
      actual <- testing$Label
      data <- data.frame(prediction, actual)
      myRes <- rbind(myRes, data)
    } 
    
    prediction <- myRes$prediction
    actual <- myRes$actual
    result <- caret::confusionMatrix(prediction, actual)
    result <- result$table
    results <- as.numeric(result)
    
    if (Label == "Bacteria") {
      
      ok[[i]] <- cbind(results[[1]], (results[[5]] + results[[9]] + results[[13]]),
                       (results[[2]] + results[[3]] + results[[4]]), (results[[6]] + results[[11]] + results[[16]]))
    } else if (Label == "Cancer") {
      ok[[i]] <- cbind(results[[6]], (results[[2]] + results[[10]] + results[[14]]), 
                       (results[[5]] + results[[7]] + results[[8]]), (results[[1]] + results[[11]] + results[[16]]))
    } else if (Label == "Fungus") {
      ok[[i]] <- cbind(results[[11]], (results[[3]] + results[[7]] + results[[15]]),
                       (results[[9]] + results[[10]] + results[[12]]), (results[[1]] + results[[6]] + results[[16]]))
    } else if (Label == "Virus") {
      ok[[i]] <- cbind(results[[16]], (results[[4]] + results[[8]] + results[[12]]), 
                       (results[[13]] + results[[14]] + results[[15]]), (results[[1]] + results[[6]] + results[[11]]))
    } }
  
  return(ok)
  stopCluster(cl)
} 

mean_and_sd <- function(x) {
  c(round(mean(x, na.rm = TRUE), digits = 4),
    round(sd(x, na.rm = TRUE), digits = 4))
}

assessment <- function(x) {
  results <- data.frame(x)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  return(results_all)
}

rf_10_CV_all <- function(x, type){
  if (type == "Bacteria") {
    ok <- RF_10_CV(x, Label = "Bacteria")
    results <- assessment(ok)
  } else if (type == "Cancer") {
    ok <- RF_10_CV(x, Label = "Cancer")
    results <- assessment(ok)
    
  } else if (type == "Fungus") {
    ok <- RF_10_CV(x, Label = "Fungus")
    results <- assessment(ok)
  } else if (type == "Virus") {
    ok <- RF_10_CV(x, Label = "Virus")
    results <- assessment(ok)
  }
  return(results)
}

RF_10_CV_all <- function(x) {
  bacteria <- rf_10_CV_all(x, type = "Bacteria")
  cancer <- rf_10_CV_all(x, type = "Cancer")
  fungus <- rf_10_CV_all(x, type = "Fungus")
  virus <- rf_10_CV_all(x, type = "Virus")
  results_all <- cbind(bacteria, cancer, fungus, virus)
  rm(bacteria)
  rm(cancer)
  rm(fungus)
  rm(virus)
  total <- apply(results_all, 1, mean)
  results_all_mean <- cbind(results_all, total)
  rm(results_all)
  rm(total)
  colnames(results_all_mean) <- c("Bacteria", "Cancer", "Fungus", "Virus", "Overall")
  return(results_all_mean)
}
######################Multiple Models


library(protr)
cancer <- readFASTA("cancer.fasta")
fungus <- readFASTA("fungus.fasta")
bacteria <- readFASTA("bacteria.fasta")
virus <- readFASTA("virus.fasta")
negative <- protr::readFASTA("negative.fasta")

### removed wired protein
cancer <- cancer[(sapply(cancer, protcheck))]
fungus <- fungus[(sapply(fungus, protcheck))]
bacteria <- bacteria[(sapply(bacteria, protcheck))]
virus <- virus[(sapply(virus, protcheck))]
negative <- negative[(sapply(negative, protcheck))]




### generating the Descriptors
cancer_des <- t(sapply(cancer, extractAAC))

fungus_des <- t(sapply(fungus, extractAAC))

bacteria_des <- t(sapply(bacteria, extractAAC))

virus_des <- t(sapply(virus, extractAAC))

negative_des <- t(sapply(negative, extractAAC))

#### label the labels 
cancer_des <- as.data.frame(cancer_des)
cancer_des$Label <- "Cancer"
fungus_des <- as.data.frame(fungus_des)
fungus_des$Label <- "Fungus"
bacteria_des <- as.data.frame(bacteria_des)
bacteria_des$Label <- "Bacteria"
virus_des <- as.data.frame(virus_des)
virus_des$Label <- "Virus"
combine_data <- rbind(cancer_des, fungus_des,
                      bacteria_des, virus_des)
combine_data$Label <- as.factor(combine_data$Label)
combine_data <- na.omit(combine_data)
combine_data <- data.frame(combine_data)

RF_testing <- function(x, Label){
  
  
  
  ok <- list(100)
  for (i in 1:100) { 
    in_train <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
    train <- x[in_train, ]
    test <- x[-in_train, ]
    model_train <- ranger::ranger(Label~., data = train, write.forest = TRUE)
    prediction <- predict(model_train, test)
    prediction <- prediction$predictions
    actual <- test$Label
    result <- caret::confusionMatrix(prediction, actual)
    result <- result$table
    results <- as.numeric(result)
    rm(train)
    rm(test)
    rm(model_train)
    rm(prediction)
    rm(actual)
    results <- as.numeric(results)
    if (Label == "Bacteria") {
      ok[[i]] <- cbind(results[[1]], (results[[5]] + results[[9]] + results[[13]]),
                       (results[[2]] + results[[3]] + results[[4]]), (results[[6]] + results[[11]] + results[[16]]))
      
      
    }  else if (Label == "Cancer") {
      
      
      ok[[i]] <- cbind(results[[6]], (results[[2]] + results[[10]] + results[[14]]), 
                       (results[[5]] + results[[7]] + results[[8]]), (results[[1]] + results[[11]] + results[[16]]))
      
      
    }  else if (Label == "Fungus") {
      
      
      ok[[i]] <- cbind(results[[11]], (results[[3]] + results[[7]] + results[[15]]),
                       (results[[9]] + results[[10]] + results[[12]]), (results[[1]] + results[[6]] + results[[16]]))
      
      
    }  else if (Label == "Virus") {
      
      
      ok[[i]] <- cbind(results[[16]], (results[[4]] + results[[8]] + results[[12]]), 
                       (results[[13]] + results[[14]] + results[[15]]), (results[[1]] + results[[6]] + results[[11]]))    
    } }
  return(ok)
  stopCluster(cl) } 


mean_and_sd <- function(x) {
  c(round(mean(x, na.rm = TRUE), digits = 4),
    round(sd(x, na.rm = TRUE), digits = 4))
}




assessment <- function(x) {
  result <- data.frame(x)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by= 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, result[TP], result[FN], result[FP], result[TN])
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  return(results_all)
}

rf_testing_all <- function(x, type){
  if (type == "Bacteria") {
    ok <- RF_testing(x, Label = "Bacteria")
    results <- assessment(ok)
  } else if (type == "Cancer") {
    ok <- RF_testing(x, Label = "Cancer")
    results <- assessment(ok)
    
  } else if (type == "Fungus") {
    ok <- RF_testing(x, Label = "Fungus")
    results <- assessment(ok)
  } else if (type == "Virus") {
    ok <- RF_testing(x, Label = "Virus")
    results <- assessment(ok)
  }
  return(results)
}

RF_testing_all <- function(x) {
  bacteria <- rf_testing_all(x, type = "Bacteria")
  cancer <- rf_testing_all(x, type = "Cancer")
  fungus <- rf_testing_all(x, type = "Fungus")
  virus <- rf_testing_all(x, type = "Virus")
  results_all <- cbind(bacteria, cancer, fungus, virus)
  rm(bacteria)
  rm(cancer)
  rm(fungus)
  rm(virus)
  total <- apply(results_all, 1, mean, na.rm = TRUE)
  results_all_mean <- cbind(results_all, total)
  rm(results_all)
  rm(total)
  colnames(results_all_mean) <- c("Bacteria", "Cancer", "Fungus", "Virus", "Overall")
  return(results_all_mean)
}

results_testing_all <- RF_testing_all(combine_data)
results_training_all <- RF_10_CV_all(combine_data)
results_CV_all <- RF_10_CV_all(combine_data)