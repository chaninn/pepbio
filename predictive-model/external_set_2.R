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



################################# AAC 


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

library(RWeka)

fit_J48 <- J48(Label ~., data)
prediction <- predict(fit_J48, validation)
dat <- data.frame(prediction, validation$Label)
dat <- data.frame(pred = prediction, truth = validation$Label)
dat$truth <- factor(dat$truth, levels = c("Negative", "Positive"))
xtab <- table(dat$pred, dat$truth)
caret::confusionMatrix(xtab, positive = "Positive")

library(ranger)
fit_rf <- ranger::ranger(Label ~., data)
prediction <- predict(fit_rf, validation)
prediction <- prediction$predictions
dat <- data.frame(prediction, validation$Label)
dat <- data.frame(pred = prediction, truth = validation$Label)
dat$truth <- factor(dat$truth, levels = c("Negative", "Positive"))
xtab <- table(dat$pred, dat$truth)
caret::confusionMatrix(xtab, positive = "Positive")

####################################DPC

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
cancer_des <- t(sapply(cancer, extractDC))
fungus_des <- t(sapply(fungus, extractDC))
bacteria_des <- t(sapply(bacteria, extractDC))
virus_des <- t(sapply(virus, extractDC))
negative_des <- t(sapply(negative, extractDC))
validation_des <- t(sapply(validation, extractDC))

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

data <- data.frame(data)
data <- na.omit(data)
validation <- data.frame(validation)
validation <- na.omit(validation)

library(RWeka)

fit_J48 <- J48(Label ~., data)
prediction <- predict(fit_J48, validation)
dat <- data.frame(prediction, validation$Label)
dat <- data.frame(pred = prediction, truth = validation$Label)
dat$truth <- factor(dat$truth, levels = c("Negative", "Positive"))
xtab <- table(dat$pred, dat$truth)
caret::confusionMatrix(xtab, positive = "Positive")

library(ranger)
fit_rf <- ranger::ranger(Label ~., data)
prediction <- predict(fit_rf, validation)
prediction <- prediction$predictions
dat <- data.frame(prediction, validation$Label)
dat <- data.frame(pred = prediction, truth = validation$Label)
dat$truth <- factor(dat$truth, levels = c("Negative", "Positive"))
xtab <- table(dat$pred, dat$truth)
caret::confusionMatrix(xtab, positive = "Positive")



#################################### Composition Class

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
cancer_des <- t(sapply(cancer, extractCTDC))
fungus_des <- t(sapply(fungus, extractCTDC))
bacteria_des <- t(sapply(bacteria, extractCTDC))
virus_des <- t(sapply(virus, extractCTDC))
negative_des <- t(sapply(negative, extractCTDC))
validation_des <- t(sapply(validation, extractCTDC))

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

data <- data.frame(data)
data <- na.omit(data)
validation <- data.frame(validation)
validation <- na.omit(validation)

library(RWeka)

fit_J48 <- J48(Label ~., data)
prediction <- predict(fit_J48, validation)
dat <- data.frame(prediction, validation$Label)
dat <- data.frame(pred = prediction, truth = validation$Label)
dat$truth <- factor(dat$truth, levels = c("Negative", "Positive"))
xtab <- table(dat$pred, dat$truth)
caret::confusionMatrix(xtab, positive = "Positive")

library(ranger)
fit_rf <- ranger::ranger(Label ~., data)
prediction <- predict(fit_rf, validation)
prediction <- prediction$predictions
dat <- data.frame(prediction, validation$Label)
dat <- data.frame(pred = prediction, truth = validation$Label)
dat$truth <- factor(dat$truth, levels = c("Negative", "Positive"))
xtab <- table(dat$pred, dat$truth)
caret::confusionMatrix(xtab, positive = "Positive")


#################################### ACC + DPC

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

AAC_DPC <- function(x) {
  c(extractAAC(x), extractDC(x))
}
### Calculate descriptors
cancer_des <- t(sapply(cancer, AAC_DPC))
fungus_des <- t(sapply(fungus, AAC_DPC))
bacteria_des <- t(sapply(bacteria, AAC_DPC))
virus_des <- t(sapply(virus, AAC_DPC))
negative_des <- t(sapply(negative, AAC_DPC))
validation_des <- t(sapply(validation, AAC_DPC))

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

data <- data.frame(data)
data <- na.omit(data)
validation <- data.frame(validation)
validation <- na.omit(validation)

library(RWeka)

fit_J48 <- J48(Label ~., data)
prediction <- predict(fit_J48, validation)
dat <- data.frame(prediction, validation$Label)
dat <- data.frame(pred = prediction, truth = validation$Label)
dat$truth <- factor(dat$truth, levels = c("Negative", "Positive"))
xtab <- table(dat$pred, dat$truth)
caret::confusionMatrix(xtab, positive = "Positive")

library(ranger)
fit_rf <- ranger::ranger(Label ~., data)
prediction <- predict(fit_rf, validation)
prediction <- prediction$predictions
dat <- data.frame(prediction, validation$Label)
dat <- data.frame(pred = prediction, truth = validation$Label)
dat$truth <- factor(dat$truth, levels = c("Negative", "Positive"))
xtab <- table(dat$pred, dat$truth)
caret::confusionMatrix(xtab, positive = "Positive")




