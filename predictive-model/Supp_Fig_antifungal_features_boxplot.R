library(protr) 

### Load the data
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
cancer_des <- t(sapply(cancer, extractAAC))
fungus_des <- t(sapply(fungus, extractAAC))
bacteria_des <- t(sapply(bacteria, extractAAC))
virus_des <- t(sapply(virus, extractAAC))
negative_des <- t(sapply(negative, extractAAC))


cancer_des <- as.data.frame(cancer_des)
cancer_des$Label <- "Anticancer"
fungus_des <- as.data.frame(fungus_des)
fungus_des$Label <- "Antifungal"
bacteria_des <- as.data.frame(bacteria_des)
bacteria_des$Label <- "Antibacterial"
virus_des <- as.data.frame(virus_des)
virus_des$Label <- "Antiviral"
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


### Making the box plots
library(ggplot2)

HDP_class <- "Antifungal"

A <- ggplot(virus, aes(factor(Label), A)) +
  geom_boxplot(aes(fill = factor(Label)), alpha = 0.7) + scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  theme(legend.position = ("none"), plot.margin = unit(c(0, 0, -0.5, 0), "cm"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1), axis.text.x=element_blank(), axis.text.y=element_blank() ) +
  coord_cartesian(ylim = c(0, 1)) + labs(x = " ", y = " ") + annotate(geom = "text", x = 1.5, y = 0.97, fontface = "bold", size = 6, label = "A")
  # + labs(y = "Percentage of Amino Acid", x = " ")
  

C <- ggplot(virus, aes(factor(Label), C)) +
  geom_boxplot(aes(fill = factor(Label)), alpha = 0.7) + scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  theme(legend.position = ("none"), plot.margin = unit(c(0, 0, -0.5, 0), "cm"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1), axis.text.x=element_blank(), axis.text.y=element_blank() ) +
  coord_cartesian(ylim = c(0, 1)) + labs(x = " ", y = " ") + annotate(geom = "text", x = 1.5, y = 0.97, fontface = "bold", size = 6, label = "C")

D <- ggplot(virus, aes(factor(Label), D)) +
  geom_boxplot(aes(fill = factor(Label)), alpha = 0.7) + scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  theme(legend.position = ("none"), plot.margin = unit(c(0, 0, -0.5, 0), "cm"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1) ) +
  coord_cartesian(ylim = c(0, 1)) + labs(x = " ", y = " ") + annotate(geom = "text", x = 1.5, y = 0.97, fontface = "bold", size = 6, label = "D")

E <- ggplot(virus, aes(factor(Label), E)) +
  geom_boxplot(aes(fill = factor(Label)), alpha = 0.7) + scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  theme(legend.position = ("none"), plot.margin = unit(c(0, 0, -0.5, 0), "cm"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1), axis.text.x=element_blank(), axis.text.y=element_blank() ) +
  coord_cartesian(ylim = c(0, 1)) + labs(x = " ", y = " ") + annotate(geom = "text", x = 1.5, y = 0.97, fontface = "bold", size = 6, label = "E")

F <- ggplot(virus, aes(factor(Label), F)) +
  geom_boxplot(aes(fill = factor(Label)), alpha  = 0.7) + scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  theme(legend.position = ("none"), plot.margin = unit(c(0, 0, -0.5, 0), "cm"), 
        panel.border  = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1), axis.text.x=element_blank(), axis.text.y=element_blank() ) +
  coord_cartesian(ylim = c(0, 1)) + labs(x = " ", y = " ") + annotate(geom = "text", x = 1.5, y = 0.97, fontface = "bold", size = 6, label = "F")

G <- ggplot(virus, aes(factor(Label), G)) +
  geom_boxplot(aes(fill = factor(Label)), alpha = 0.7) + scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  theme(legend.position = ("none"), plot.margin = unit(c(0, 0, -0.5, 0), "cm"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1), axis.text.x=element_blank(), axis.text.y=element_blank() ) +
  coord_cartesian(ylim = c(0, 1)) + labs(x = " ", y = " ") + annotate(geom = "text", x = 1.5, y = 0.97, fontface = "bold", size = 6, label = "G")

H <- ggplot(virus, aes(factor(Label), H)) +
  geom_boxplot(aes(fill = factor(Label)), alpha = 0.7) + scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  theme(legend.position = ("none"), plot.margin = unit(c(0, 0, -0.5, 0), "cm"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1), axis.text.x=element_blank() ) +
  coord_cartesian(ylim = c(0, 1)) + labs(x = " ", y = " ") + annotate(geom = "text", x = 1.5, y = 0.97, fontface = "bold", size = 6, label = "H")

I <- ggplot(virus, aes(factor(Label), I)) +
  geom_boxplot(aes(fill = factor(Label)), alpha = 0.7) + scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  theme(legend.position = ("none"), plot.margin = unit(c(0, 0, -0.5, 0), "cm"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1), axis.text.y=element_blank() ) +
  coord_cartesian(ylim= c(0, 1)) + labs(x = " ", y = " ") + annotate(geom = "text", x = 1.5, y = 0.97, fontface = "bold", size = 6, label = "I")

K <- ggplot(virus, aes(factor(Label), K)) +
  geom_boxplot(aes(fill = factor(Label)), alpha = 0.7) + scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  theme(legend.position = ("none"), plot.margin = unit(c(0, 0, -0.5, 0), "cm"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1), axis.text.x=element_blank(), axis.text.y=element_blank() ) +
  coord_cartesian(ylim = c(0, 1)) + labs(x = " ", y = " ") + annotate(geom = "text", x = 1.5, y = 0.97, fontface = "bold", size = 6, label = "K")

L <- ggplot(virus, aes(factor(Label), L)) +
  geom_boxplot(aes(fill = factor(Label)), alpha = 0.7) + scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  theme(legend.position = ("none"), plot.margin = unit(c(0, 0, -0.5, 0), "cm"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1), axis.text.x=element_blank(), axis.text.y=element_blank() ) +
  coord_cartesian(ylim = c(0, 1)) + labs(x = " ", y = " ") + annotate(geom = "text", x = 1.5, y = 0.97, fontface = "bold", size = 6, label = "L")

M <- ggplot(virus, aes(factor(Label), M)) +
  geom_boxplot(aes(fill = factor(Label)), alpha = 0.7) + scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  theme(legend.position = ("none"), plot.margin = unit(c(0, 0, -0.5, 0), "cm"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1), axis.text.y=element_blank() ) +
  coord_cartesian(ylim = c(0, 1)) + labs(x = " ", y = " ") + annotate(geom = "text", x = 1.5, y = 0.97, fontface = "bold", size = 6, label = "M")

N <- ggplot(virus, aes(factor(Label), N)) +
  geom_boxplot(aes(fill = factor(Label)), alpha = 0.7) + scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  theme(legend.position = ("none"), plot.margin = unit(c(0, 0, -0.5, 0), "cm"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1), axis.text.x=element_blank(), axis.text.y=element_blank() ) +
  coord_cartesian(ylim = c(0, 1)) + labs(x = " ", y = " ") + annotate(geom = "text", x = 1.5, y = 0.97, fontface = "bold", size = 6, label = "N")

P <- ggplot(virus, aes(factor(Label), P)) +
  geom_boxplot(aes(fill = factor(Label)), alpha = 0.7) + scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  theme(legend.position = ("none"), plot.margin = unit(c(0, 0, -0.5, 0), "cm"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1), axis.text.x=element_blank(), axis.text.y=element_blank() ) +
  coord_cartesian(ylim = c(0, 1)) + labs(x = " ", y = " ") + annotate(geom = "text", x = 1.5, y = 0.97, fontface = "bold", size = 6, label = "P")

Q <- ggplot(virus, aes(factor(Label), Q)) +
  geom_boxplot(aes(fill = factor(Label)), alpha = 0.7) + scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  theme(legend.position = ("none"), plot.margin = unit(c(0, 0, -0.5, 0), "cm"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1), axis.text.x=element_blank() ) +
  coord_cartesian(ylim = c(0, 1)) + labs(x = " ", y = " ") + annotate(geom = "text", x = 1.5, y = 0.97, fontface = "bold", size = 6, label = "Q")

R <- ggplot(virus, aes(factor(Label), R)) +
  geom_boxplot(aes(fill = factor(Label)), alpha = 0.7) + scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  theme(legend.position = ("none"), plot.margin = unit(c(0, 0, -0.5, 0), "cm"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1), axis.text.x=element_blank(), axis.text.y=element_blank() ) +
  coord_cartesian(ylim = c(0, 1)) + labs(x = " ", y = " ") + annotate(geom = "text", x = 1.5, y = 0.97, fontface = "bold", size = 6, label = "R")

S <- ggplot(virus, aes(factor(Label), S)) +
  geom_boxplot(aes(fill = factor(Label)), alpha = 0.7) + scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  theme(legend.position = ("none"), plot.margin = unit(c(0, 0, -0.5, 0), "cm"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1), axis.text.y=element_blank() ) +
  coord_cartesian(ylim = c(0, 1)) + labs(x = " ", y = " ") + annotate(geom = "text", x = 1.5, y = 0.97, fontface = "bold", size = 6, label = "S")

T <- ggplot(virus, aes(factor(Label), T)) +
  geom_boxplot(aes(fill = factor(Label)), alpha = 0.7) + scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  theme(legend.position = ("none"), plot.margin = unit(c(0, 0, -0.5, 0), "cm"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1), axis.text.x=element_blank() ) +
  coord_cartesian(ylim = c(0, 1)) + labs(x = " ", y = " ") + annotate(geom = "text", x = 1.5, y = 0.97, fontface = "bold", size = 6, label = "T")

V <- ggplot(virus, aes(factor(Label), V)) +
  geom_boxplot(aes(fill = factor(Label)), alpha = 0.7) + scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  theme(legend.position = ("none"), plot.margin = unit(c(0, 0, -0.5, 0), "cm"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1), axis.text.x=element_blank(), axis.text.y=element_blank() ) +
  coord_cartesian(ylim = c(0, 1)) + labs(x = " ", y = " ") + annotate(geom = "text", x = 1.5, y = 0.97, fontface = "bold", size = 6, label = "V")

W <- ggplot(virus, aes(factor(Label), W)) +
  geom_boxplot(aes(fill = factor(Label)), alpha = 0.7) + scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  theme(legend.position = ("none"), plot.margin = unit(c(0, 0, -0.5, 0), "cm"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1), axis.text.x=element_blank(), axis.text.y=element_blank() ) +
  coord_cartesian(ylim = c(0, 1)) + labs(x = " ", y = " ") + annotate(geom = "text", x = 1.5, y = 0.97, fontface = "bold", size = 6, label = "W")

Y <- ggplot(virus, aes(factor(Label), Y)) +
  geom_boxplot(aes(fill = factor(Label)), alpha = 0.7) + scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  theme(legend.position = ("none"), plot.margin = unit(c(0, 0, -0.5, 0), "cm"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1), axis.text.y=element_blank() ) +
  coord_cartesian(ylim = c(0, 1)) + labs(x = " ", y = " ") + annotate(geom = "text", x = 1.5, y = 0.97, fontface = "bold", size = 6, label = "Y")

# + theme(axis.text.x = element_blank(), axis.text.y = element_blank())

### Fusing the plots

### Antibacterial T, Q, N, E, P, H, V, I, G, F, D, W, L, C, Y, S, A, R, M, K
### Anticancer    F, G, W, P, V, T, C, L, Y, K, Q, S, R, N, M, I, H, E, D, A
### Antifungal    T, L, P, C, N, Q, V, F, G, W, H, K, R, E, A, D, Y, S, M, I
### Antivirus     W, T, F, V, R, M, C, Y, K, L, I, G, P, S, A, Q, N, H, E, D

library(cowplot)

antifungal_plot <- plot_grid(T, L, P, C, N, Q, V, F, G, W, H, K, R, E, A, D, Y, S, M, I,
          align = 'h',
          rel_widths = c(1.14, 1, 1, 1, 1, 1.14, 1, 1, 1, 1, 1.14, 1, 1, 1, 1, 1.14, 1, 1, 1, 1),
          ncol = 5)

save_plot("Fig-antifungal-features.pdf", antifungal_plot, 
          base_height=10, base_width=12)
