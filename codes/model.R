#get the external testing data from dover anlayzer
library(Biostrings)
cancer <- readAAStringSet("cancer.fasta")
dover <- readAAStringSet("AMPs (nrdb-1.00).fasta")
fungus <- readAAStringSet("fungus.fasta")
virus <- readAAStringSet("virus.fasta")  
bacteria <- readAAStringSet("bacteria.fasta")  
seq_cancer_name <- names(cancer)
seq_cancer <- paste(cancer)
seq_dover_name <- names(dover)
seq_dover <- paste(dover)
seq_fungus_name <- names(fungus)
seq_fungus <- paste(fungus)
seq_virus_name <- names(virus)
seq_virus <- paste(virus)
seq_bacteria_name <- names(bacteria)
seq_bacteria <- paste(bacteria)
  
pepBio <- c(seq_cancer, seq_fungus, seq_virus, seq_bacteria)
df_dover <- data.frame(names = seq_dover_name, sequence = seq_dover)
validation_set <- df_dover[!df_dover$sequence %in% pepBio, ]
names(validation_set) <- c("names", "sequences")
library(seqRFLP)
validation_set_fasta <- dataframe2fas(validation_set, file = "validation_set_fasta.fasta")

  
library(protr)  
cancer <- readFASTA("cancer.fasta")
fungus <- readFASTA("fungus.fasta")
bacteria <- readFASTA("bacteria.fasta")
virus <- readFASTA("virus.fasta")
negative <- protr::readFASTA("negative.fasta")
validation <- protr::readFASTA("validation_set_fasta.fasta")

### Remove problematic peptide
cancer <- cancer[(sapply(cancer, protcheck))]
fungus <- fungus[(sapply(fungus, protcheck))]
bacteria <- bacteria[(sapply(bacteria, protcheck))]
virus <- virus[(sapply(virus, protcheck))]
negative <- negative[(sapply(negative, protcheck))]
validation <- validation[(sapply(validation, protcheck))]

### Calculate descriptors
cancer_des <- t(sapply(cancer, extractAAC))
fungus_des <- t(sapply(fungus, extractAAC))
bacteria_des <- t(sapply(bacteria, extractAAC))
virus_des <- t(sapply(virus, extractAAC))
negative_des <- t(sapply(negative, extractAAC))
validation_des <- t(sapply(validation, extractAAC))

cancer_des <- as.data.frame(cancer_des)
fungus_des <- as.data.frame(fungus_des)
bacteria_des <- as.data.frame(bacteria_des)
virus_des <- as.data.frame(virus_des)
negative_des <- as.data.frame(negative_des)
validation_des <- as.data.frame(validation_des)
#### label the labels 
hdp <- rbind(cancer_des, fungus_des, bacteria_des, virus_des)
hdp$Label <- c("Positive")
hdp <- as.data.frame(hdp)
negative_des <- as.data.frame(negative_des)
negative_des$Label <- c("Negative")

validation_des$Label <- c("Positive")




data <- rbind(hdp, negative_des)
data$Label <- as.factor(data$Label)
validation <- validation_des
validation$Label <- as.factor(validation$Label)

cancer <- rbind(cancer_des, negative_des)
cancer$Label <- as.factor(cancer$Label)
fungus <- rbind(fungus_des, negative_des)
fungus$Label <- as.factor(fungus$Label)
bacteria <- rbind(fungus_des, negative_des)
bacteria$Label <- as.factor(bacteria$Label)
virus <- rbind(virus_des, negative_des)
virus$Label <- as.factor(virus$Label)

input <- list(cancer = cancer, fungus = fungus, bacteria = bacteria, 
              virus = virus)


J48_training <- function(x) {
  
  results <- list(100)
  for (i in 1:100) {
    in_train <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
    train <- x[in_train, ]
    test <- x[-in_train, ]
    model_train <- RWeka::J48(Label~., data = train)
    summary <- summary(model_train)
    confusionmatrix <- summary$confusionMatrix
    results[[i]] <- as.numeric(confusionmatrix)
  }
  return(results)
}
set.seed(1)
in_train <- caret::createDataPartition(data$Label, p = 0.8, list = FALSE)
train <- data[in_train, ]
test <- data[-in_train, ]
write.csv(train[, ], file = "Train.csv", row.names = FALSE)
write.csv(test[, ], file = "Test.csv", row.names = FALSE)


model_train <- RWeka::J48(Label~., data = data)
prediction <- predict(model_train, validation)
actual <- validation$Label
actual <- factor(actual, levels = c("Negative", "Positive"))
value <- data.frame(obs = actual, pred = prediction)
result <- kebabs::evaluatePrediction(prediction, actual, print = FALSE)
result_2 <- caret::multiClassSummary(value, lev = c("Negative", "Positive"))
caret::twoClassSummary(value)

reviewer_descriptors <- function(x) {
  c(protr::extractAAC(x), protr::extractDC(x))
}



library(protr)
cancer <- readFASTA("cancer.fasta")
fungus <- readFASTA("fungus.fasta")
bacteria <- readFASTA("bacteria.fasta")
virus <- readFASTA("virus.fasta")
negative <- protr::readFASTA("negative.fasta")

### Remove problematic peptide
cancer <- cancer[(sapply(cancer, protcheck))]
fungus <- fungus[(sapply(fungus, protcheck))]
bacteria <- bacteria[(sapply(bacteria, protcheck))]
virus <- virus[(sapply(virus, protcheck))]
negative <- negative[(sapply(negative, protcheck))]


### Calculate descriptors
cancer_des <- t(sapply(cancer, extractCTDC))
fungus_des <- t(sapply(fungus, extractCTDC))
bacteria_des <- t(sapply(bacteria, extractCTDC))
virus_des <- t(sapply(virus, extractCTDC))
negative_des <- t(sapply(negative, extractCTDC))
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

input <- list(cancer = cancer, fungus = fungus, bacteria = bacteria, 
              virus = virus)

#### Training results using J48
J48_training <- function(x) {
  
  results <- list(100)
  for (i in 1:100) {
    in_train <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
    train <- x[in_train, ]
    test <- x[-in_train, ]
    model_train <- RWeka::J48(Label~., data = train)
    summary <- summary(model_train)
    confusionmatrix <- summary$confusionMatrix
    results[[i]] <- as.numeric(confusionmatrix)
  }
  return(results)
}

mean_and_sd <- function(x) {
  c(round(mean(x, na.rm = TRUE), digits = 4),
    round(sd(x, na.rm = TRUE), digits = 4))
}

J48_train <- function(x) {
  ok <- J48_training(x)
  results <- data.frame(ok)
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

#### 10-fold results using J48
J48_10fold <- function(x) {
  results <- list(100)
  for (i in 1:100) {
    in_train <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
    train <- x[in_train, ]
    test <- x[-in_train, ]
    model_train <- RWeka::J48(Label~., data = train)
    eval_j48 <- RWeka::evaluate_Weka_classifier(model_train, numFolds = 10, complexity = FALSE, seed = 1, class = TRUE)
    confusionmatrix <- eval_j48$confusionMatrix
    results[[i]] <- as.numeric(confusionmatrix)
  }
  return(results)
}

J48_cross_validation <- function(x) {
  ok <- J48_10fold(x)
  results <- data.frame(ok)
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

### Testing results using J48
J48_testing <- function(x) {
  results <- list(100)
  for (i in 1:100) {
    in_train <- caret::createDataPartition(x$Label, p = 0.80, list = FALSE)
    train <- x[in_train, ]
    test <- x[-in_train, ]
    model_train <- RWeka::J48(Label~., data = train)
    eval_external <- RWeka::evaluate_Weka_classifier(model_train, newdata = test, numFolds = 0, complexity = FALSE, seed = 1, class = TRUE)
    confusionmatrix <- eval_external$confusionMatrix
    results[[i]] <- as.numeric(confusionmatrix)
  }
  return(results)
}


J48_external <- function(x) {
  ok <- J48_testing(x)
  results <- data.frame(ok)
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

results_J48 <- function(x) {
  training <- J48_train(x)
  cross_validation <- J48_cross_validation(x)
  testing <- J48_external(x)
  results <- data.frame(Training = training, Cross_Validation = cross_validation, Testing = testing)
  colnames(results) <- c("Training set", "Internal validation", "External set")
  return(results)
}
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(doSNOW))
cl <- makeCluster(7)
registerDoSNOW(cl)
clusterExport(cl = cl, ls())

result_J48_performance <- parLapply(cl = cl, input, function(x) {
  models <- suppressWarnings(results_J48(x))
  return(models)
})
print(result_J48_performance)
stopCluster(cl)


##### Euclidean Distance

library(fields)
out <- rdist(train[-21], test[-21])
results <- data.frame()
for (i in 1:ncol(out)) {
  EuclideanDistance <- sum(out[, i])
  mean_distance <- mean(out[, i][-i])
  near_5_neighbours_position <- order(out[, i])[1:6]
  near_5_neighbours_value <- out[, i][order(out[, i])[1:6]]
  mean_near_5_neightoubrs_value <- mean(near_5_neighbours_value[-1], na.rm = TRUE)
  result <- c(EuclideanDistance, mean_distance, near_5_neighbours_value, mean_near_5_neightoubrs_value)
  result <- t(result)
  result <- data.frame(result)
  names(result) <- c("EuclideanDistance", "Mean_Distance", "Original", "NN_1", "NN_2",
                     "NN_3", "NN_4", "NN_5", "Mean_NN")
  
  results <- rbind(results, result)
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
Norm_MeanDistance <- range01(results$Mean_NN)
great <- cbind(results, Norm_MeanDistance, df$CompNo)
write.csv(great, file = "EuclideanDistance/EuclideanDistance/Euclidean_distance_NN.csv", row.names = FALSE)



library(RWeka)
fit <- J48(Label ~., data = train)
prediction <- predict(fit, test, type = "class")

gg_df <- data.frame(predictions = prediction, NN_5 = Norm_MeanDistance, Activity = test$Label)
ggplot(gg_df, aes(x = NN_5, y = probability_score)) + geom_point(aes(color = Activity)) +
  ggtitle("Positive") + xlab("Normalized 5 Nearest Neighbour Eucliean Distance") + 
  ylab("Prediction Probability Score")

#bin_NN <- cut(gg_df$NN_5, breaks = quantile(Norm_MeanDistance), include.lowest = TRUE, labels = FALSE)
#gg_df <- cbind(gg_df, bin_NN)
positive <- subset(gg_df, Activity == "Positive")
bin_NN <- cut(positive$NN_5, breaks = quantile(positive$NN_5), include.lowest = TRUE, labels = FALSE)
positive <- cbind(positive, bin_NN)
negative <- subset(gg_df, Activity == "Negative")
bin_NN <- cut(negative$NN_5, breaks = quantile(negative$NN_5), include.lowest = TRUE, labels = FALSE)
negative <- cbind(negative, bin_NN)

library(easyGgplot2)

ggplot2.dotplot(data = positive, xName ="bin_NN", yName = "NN_5", groupName = "predictions")


###### individual pie chart
a_positive <- subset(positive, bin_NN == "1")
a_positive <- table(a_positive$predictions)
a_positive <- data.frame(a_positive)
names(a_positive) <- c("Activity", "Count")
a_positive$lab <- as.character(round(100 * a_positive$Count / sum(a_positive$Count), 1))
a_positive_plot <- ggplot(a_positive, aes(x = 1, y = Count, fill = Activity)) +
  geom_bar(stat = "identity") +  #geom_text(aes(label = lab), vjust = 1, size = 5) +
  ggtitle("Accuracy = 100.0% (N = 421)") + coord_polar(theta = 'y') + theme(text = element_text(size = 22),
                                                legend.position = ("none"),
                                                plot.title = element_text(vjust = 10),
                                                 axis.ticks = element_blank(),
                                                 axis.text = element_blank(),
                                                 axis.title = element_blank())


b_positive <- subset(positive, bin_NN == "2")
b_positive <- table(b_positive$predictions)
b_positive <- data.frame(b_positive)
names(b_positive) <- c("Activity", "Count")
b_positive$lab <- as.character(round(100 * b_positive$Count / sum(b_positive$Count), 1))
b_positive_plot <- ggplot(b_positive, aes(x = 1, y = Count, fill = Activity)) +
  geom_bar(stat = "identity", position = "stack") + # geom_text(aes(label = lab), vjust = 1, size = 5) +
  ggtitle("Accuracy = 100.0% (N = 420)") + coord_polar(theta = 'y') + theme(text = element_text(size = 22),
                                                                            legend.position = ("none"),
                                               axis.ticks = element_blank(),
                                               axis.text = element_blank(),
                                               axis.title = element_blank())


c_positive <- subset(positive, bin_NN == "3")
c_positive <- table(c_positive$predictions)
c_positive <- data.frame(c_positive)
names(c_positive) <- c("Activity", "Count")
c_positive$lab <- as.character(round(100 * c_positive$Count / sum(c_positive$Count), 1))
c_positive_plot <- ggplot(c_positive, aes(x = 1, y = Count, fill = Activity)) +
  geom_bar(stat = "identity") +  #geom_text(aes(label = lab), vjust = 1) +
  ggtitle("Accuracy = 99.8% (N = 420)") + coord_polar(theta = 'y') + theme(text = element_text(size = 22),
                                                                           legend.position = ("none"),
                                                 axis.ticks = element_blank(),
                                                 axis.text = element_blank(),
                                                 axis.title = element_blank())

d_positive <- subset(positive, bin_NN == "4")
d_positive <- table(d_positive$predictions)
d_positive <- data.frame(d_positive)
names(d_positive) <- c("Activity", "Count")
d_positive$lab <- as.character(round(100 * d_positive$Count / sum(d_positive$Count), 1))
d_positive_plot <- ggplot(d_positive, aes(x = 1, y = Count, fill = Activity)) +
  geom_bar(stat = "identity") + # geom_text(aes(label = lab), vjust = 1) +
  ggtitle("Accuracy = 98.3% (N = 421)") + coord_polar(theta = 'y') + theme(text = element_text(size = 22),
                                                                           legend.position = ("none"),
                                                 axis.ticks = element_blank(),
                                                 axis.text = element_blank(),
                                                 axis.title = element_blank()) 

library(cowplot)
positive_plot <- plot_grid(a_positive_plot, b_positive_plot, c_positive_plot, d_positive_plot, labels = c("Q1 (0.0 - 0.18)", "Q2 (0.18-0.27)", "Q3 (0.27-0.36)", "Q4 (0.36-1.0)"),
          nrow = 1)

#### Negative 
a_negative <- subset(negative, bin_NN == "1")
a_negative <- table(a_negative$predictions)
a_negative <- data.frame(a_negative)
names(a_negative) <- c("Activity", "Count")
a_negative$lab <- as.character(round(100 * a_negative$Count / sum(a_negative$Count), 1))
a_negative_plot <- ggplot(a_negative, aes(x = 1, y = Count, fill = Activity)) +
  geom_bar(stat = "identity") +  #geom_text(aes(label = lab), vjust = 1, size = 5) +
  ggtitle("Accuracy = 100.0% (N = 86)") + coord_polar(theta = 'y') + theme(text = element_text(size = 22),
                                                                            legend.position = ("none"),
                                                                            axis.ticks = element_blank(),
                                                                            axis.text = element_blank(),
                                                                            axis.title = element_blank())



b_negative <- subset(negative, bin_NN == "2")
b_negative <- table(b_negative$predictions)
b_negative <- data.frame(b_negative)
names(b_negative) <- c("Activity", "Count")
b_negative$lab <- as.character(round(100 * b_negative$Count / sum(b_negative$Count), 1))
b_negative_plot <- ggplot(b_negative, aes(x = 1, y = Count, fill = Activity)) +
  geom_bar(stat = "identity") +  #geom_text(aes(label = lab), vjust = 1, size = 5) +
  ggtitle("Accuracy = 97.6% (N = 85)") + coord_polar(theta = 'y') + theme(text = element_text(size = 22),
                                                                            legend.position = ("none"),
                                                                            axis.ticks = element_blank(),
                                                                            axis.text = element_blank(),
                                                                            axis.title = element_blank())

c_negative <- subset(negative, bin_NN == "3")
c_negative <- table(c_negative$predictions)
c_negative <- data.frame(c_negative)
names(c_negative) <- c("Activity", "Count")
c_negative$lab <- as.character(round(100 * c_negative$Count / sum(c_negative$Count), 1))
c_negative_plot <- ggplot(c_negative, aes(x = 1, y = Count, fill = Activity)) +
  geom_bar(stat = "identity") +  #geom_text(aes(label = lab), vjust = 1, size = 5) +
  ggtitle("Accuracy = 98.8% (N = 85)") + coord_polar(theta = 'y') + theme(text = element_text(size = 22),
                                                                            legend.position = ("none"),
                                                                            axis.ticks = element_blank(),
                                                                            axis.text = element_blank(),
                                                                            axis.title = element_blank())
d_negative <- subset(negative, bin_NN == "4")
d_negative <- table(d_negative$predictions)
d_negative <- data.frame(d_negative)
names(d_negative) <- c("Activity", "Count")
d_negative$lab <- as.character(round(100 * d_negative$Count / sum(d_negative$Count), 1))
d_negative_plot <- ggplot(d_negative, aes(x = 1, y = Count, fill = Activity)) +
  geom_bar(stat = "identity") +  #geom_text(aes(label = lab), vjust = 1, size = 5) +
  ggtitle("Accuracy = 80.2% (N = 86)") + coord_polar(theta = 'y') + theme(text = element_text(size = 22),
                                                                            legend.position = ("none"),
                                                                            axis.ticks = element_blank(),
                                                                            axis.text = element_blank(),
                                                                            axis.title = element_blank())


library(cowplot)
negative_plot <- plot_grid(a_negative_plot, b_negative_plot, c_negative_plot, d_negative_plot, labels = c("Q1 (0.0 - 0.18)", "Q2 (0.18-0.27)", "Q3 (0.27-0.36)", "Q4 (0.36-1.0)"),
          nrow = 1)
plot_grid(positive_plot, negative_plot, nrow = 2)


library(ggplot2)
a <- ggplot(positive, aes())
#ggplot(gg_df, aes(x = bin_NN, y = probability_score)) + geom_point(stat = "identity", aes(color = Activity)) +
#  scale_x_discrete(breaks = c("1", "2", "3", "4"), labels = c("0.25", "0.50", "0.75", "1.0"))
#  scale_x_discrete(labels = c("1" = "0-0.25", "2" = "0.25-0.50", "3" = "0.50-0.75", "4" = "0.75-1.0"))

################################Feature importances 
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


library(ggplot2)
### BActeria T, Q, N, E, P
### CAncer F, G, W, P, V
### Fungus T, L, P, C, N
### Virus W, T, F, V, R
a <- summarySE(virus, measurevar = c("W"), groupvars = "Label")
MinMeanSEMMax <- function(x) {
  v <- c(min(x), mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)), max(x))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v
}

a <- ggplot(virus, aes(factor(Label), W)) +
  geom_errorbar(data = a, aes(ymin=W-se, ymax = W+se), width = 0.5) +
  stat_summary(fun.data = MinMeanSEMMax, geom = "boxplot", colour = "red", aes(fill = factor(Label)), alpha = 0.7) +
  
  theme(legend.position = ("none"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1)) +
  coord_cartesian(ylim = c(-1, 1)) + labs(y = "Percentage of Amino Acid", x = " ") + ggtitle("W")

b <- ggplot(virus, aes(factor(Label), T)) +
  geom_boxplot(aes(fill = factor(Label)), alpha = 0.7) +
  theme(legend.position = ("none"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1)) +
  coord_cartesian(ylim = c(0, 1)) + labs(y = "Percentage of Amino Acid", x = " ") + ggtitle("T")

c <- ggplot(virus, aes(factor(Label), F)) +
  geom_boxplot(aes(fill = factor(Label)), alpha  = 0.7) +
  theme(legend.position = ("none"), 
        panel.border  = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1)) +
  coord_cartesian(ylim = c(0, 1)) + labs(y = "Percentage of Amino Acid", x = " ") + ggtitle("F")

d <- ggplot(virus, aes(factor(Label), V)) +
  geom_boxplot(aes(fill = factor(Label)), alpha = 0.7) +
  theme(legend.position = ("none"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1)) +
  coord_cartesian(ylim= c(0, 1)) + labs(y = "Percentage of Amino ACid", x = " ") + ggtitle("V")

e <- ggplot(virus, aes(factor(Label), R)) +
  geom_boxplot(aes(fill = factor(Label)), alpha = 0.7) + 
  theme(legend.position = ("none"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1)) +
  coord_cartesian(ylim = c(0, 1)) + labs(y = "Percentage of Amino Acid", x = " ") + ggtitle("R")

plot_grid(a, b, c, d, e)
###########################Random forest
