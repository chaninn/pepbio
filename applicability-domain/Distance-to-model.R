
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

set.seed(1)
in_train <- caret::createDataPartition(data$Label, p = 0.8, list = FALSE)
train <- data[in_train, ]
test <- data[-in_train, ]

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
  names(result) <- c("EuclideanDistance", "Mean_Distance", "Position", "NN_1", "NN_2",
                     "NN_3", "NN_4", "NN_5", "Mean_NN")
  
  results <- rbind(results, result)
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
Norm_MeanDistance <- range01(results$Mean_NN)
library(RWeka)
fit <- J48(Label ~., data = train)
prediction <- predict(fit, test, type = "class")

gg_df <- data.frame(predictions = prediction, NN_5 = Norm_MeanDistance, Activity = test$Label)

positive <- subset(gg_df, Activity == "Positive")
bin_NN <- cut(positive$NN_5, breaks = quantile(positive$NN_5), include.lowest = TRUE, labels = FALSE)
positive <- cbind(positive, bin_NN)
negative <- subset(gg_df, Activity == "Negative")
bin_NN <- cut(negative$NN_5, breaks = quantile(negative$NN_5), include.lowest = TRUE, labels = FALSE)
negative <- cbind(negative, bin_NN)

View(positive)


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
