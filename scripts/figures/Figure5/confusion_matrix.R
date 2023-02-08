library(ggplot2)
library(ggpubr)

outdir <- "/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/integration"

### Make confusion matrices
# Positive samples
data_pos_sample <- data.frame(A = c("Positive", "Positive"),
                              B = c("Negative", "Positive"),
                              Y = c(7, 31),
                              color = c("bad", "good"))
sum_pos_sample <- sum(data_pos_sample$Y)

# Positive patients
data_pos_patient <- data.frame(A = c("Positive", "Positive"),
                               B = c("Negative", "Positive"),
                               Y = c(3, 22),
                               color = c("bad", "good"))
sum_pos_patient <- sum(data_pos_patient$Y)

# Negative samples
data_neg_sample <- data.frame(A = c("Negative", "Negative", "Positive", "Positive"),
                              B = c("Negative", "Positive", "Negative", "Positive"),
                              Y = c(72, 19, 0, 40),
                              color = c("good", "bad", "bad", "good"))
sum_neg_sample <- sum(data_neg_sample$Y)

# Negative patients
data_neg_patient <- data.frame(A = c("Negative", "Negative", "Positive", "Positive"),
                               B = c("Negative", "Positive", "Negative", "Positive"),
                               Y = c(33, 15, 0, 25),
                               color = c("good", "bad", "bad", "good"))
sum_neg_patient <- sum(data_neg_patient$Y)

### Set theme
theme <- theme(plot.title = element_text(hjust = 0.5, size = 13), 
               axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.position = "none",
               legend.key = element_rect(fill = "white"),
               legend.text = element_text(size = 12),
               axis.text = element_text(size = 13),
               axis.title = element_text(size = 13))

### Plot cancer positives
mat_pos <- ggplot(data_pos_sample, aes(A, B)) +
  geom_tile(aes(fill = color), color = "black", size = 1, alpha = 0.25) +
  geom_text(aes(label = paste0(Y, "\n", round(Y/sum_pos_sample, 4)*100, "%")), vjust = 0.5, size = 4) +
  scale_fill_manual(values = c(good = "#E41A1C", bad = "white")) +
  xlab("Clinical Finding") +
  ylab("ctDNA Finding") +
  ggtitle("LFS Cancer\nPositive Samples") +
  theme +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))
mat_pos

mat_pos2 <- ggplot(data_pos_patient, aes(A, B)) +
  geom_tile(aes(fill = color), color = "black", size = 1, alpha = 0.25) +
  geom_text(aes(label = paste0(Y, "\n", round(Y/sum_pos_patient, 4)*100, "%")), vjust = 0.5, size = 4) +
  scale_fill_manual(values = c(good = "#E41A1C", bad = "white")) +
  xlab("Clinical Finding") +
  ylab("ctDNA Finding") +
  ggtitle("LFS Cancer\nPositive Patients") +
  theme +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))
mat_pos2

fig <- ggarrange(mat_pos, mat_pos2, align = "hv")
fig
ggsave(file.path(outdir, "confusion_matrix_pos.pdf"), fig, width = 5.35, height = 2)

### Plot cancer negatives
mat_neg <- ggplot(data_neg_sample, aes(A, B)) +
  geom_tile(aes(fill = color), color = "black", size = 1, alpha = 0.25) +
  geom_text(aes(label = paste0(Y, "\n", round(Y/sum_neg_sample, 4)*100, "%")), vjust = 0.5, size = 4) +
  scale_fill_manual(values = c(good = "#E41A1C", bad = "white")) +
  xlab(paste0("Clinical Finding\n\nPPV = ", round(data_neg_sample$Y[4]/sum(data_neg_sample$Y[c(2,4)])*100, 2), "%\n",
              "NPV = ", round(data_neg_sample$Y[1]/sum(data_neg_sample$Y[c(1,3)])*100, 2), "%")) +
  ylab("ctDNA Finding") +
  ggtitle("LFS Cancer\nNegative Samples") +
  theme +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))
mat_neg

mat_neg2 <- ggplot(data_neg_patient, aes(A, B)) +
  geom_tile(aes(fill = color), color = "black", size = 1, alpha = 0.25) +
  geom_text(aes(label = paste0(Y, "\n", round(Y/sum_neg_patient, 4)*100, "%")), vjust = 0.5, size = 4) +
  scale_fill_manual(values = c(good = "#E41A1C", bad = "white")) +
  xlab(paste0("Clinical Finding\n\nPPV = ", round(data_neg_patient$Y[4]/sum(data_neg_patient$Y[c(2,4)])*100, 2), "%\n",
              "NPV = ", round(data_neg_patient$Y[1]/sum(data_neg_patient$Y[c(1,3)])*100, 2), "%")) +
  ylab("ctDNA Finding") +
  ggtitle("LFS Cancer\nNegative Patients") +
  theme +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))
mat_neg2

fig <- grid.arrange(arrangeGrob(mat_neg + theme(plot.margin=unit(c(0.2,0.2,0.2,0.2), "cm")), 
                                mat_neg2 + theme(axis.title.y = element_blank(),
                                                 axis.text.y = element_blank(),
                                                 plot.margin=unit(c(0.2,0.2,0.2,0.2), "cm")),
                                widths = c(1.5, 1),
                                nrow = 1))
fig
ggsave(file.path(outdir, "confusion_matrix_neg.pdf"), fig, width = 5, height = 4)


