#####################
###  CHROMOSOMES  ###
#####################

# input total number of fragments
frag_num_chromo <- 86570

# fragments which are not 100% chromosome (plasmid > 0%)
data_chromo <- read.table("stats-chromosome.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
frag_gt_0_chromo <- nrow(data_chromo)

# misclassified fragments (plasmid > 50%)
data_cut <- data_chromo[data_chromo$Plasmid / data_chromo$Total > 0.5,]
frag_gt_50 <- nrow(data_cut)


# pie chart 1
perc1 <- round(frag_gt_0_chromo / frag_num_chromo * 100, 2)
perc2 <- 100 - perc1
lbls <- c(paste(perc1, "%"), paste(perc2, "%"))
colors <- c("red", "green")
pie(c(perc1, perc2), labels = lbls, main="Correctly Classified Chromosome Fragments", col=colors)
legend("topleft", c("yes", "no"), cex = 0.8, fill = colors, title=expression(bold(">0% plasmid")))

# pie chart 2
perc1 <- round(frag_gt_50 / frag_gt_0_chromo * 100, 2)
perc2 <- 100 - perc1
lbls <- c(paste(perc1, "%"), paste(perc2, "%"))
colors <- c("red", rgb(51, 0, 102, max=255))
pie(c(perc1, perc2), labels = lbls, main="Potentially Ambiguous Chromosome Fragments", col=colors)
legend(-1.10, 1.08, c("yes", "no"), cex = 0.8, fill = colors, title=expression(bold(">50% plasmid")))

# histogram 1
hist(data_chromo$Plasmid / data_chromo$Total * 100, 
     main="Potentially Ambiguous Chromosome Fragments", 
     xlab="Plasmid (%)", 
     col="red",
     breaks=10)

# histogram 2
hist(data_chromo$None / data_chromo$Total * 100, 
     main="Potentially Ambiguous Chromosome Fragments", 
     xlab="None (%)", 
     col="red",
     breaks=10)



##################
###  PLASMIDS  ###
##################

# input total number of fragments
frag_num_pla <- 84019

# fragments which are not 100% plasmid (chromosome > 0%)
data_pla <- read.table("stats-plasmid.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
frag_gt_0_pla <- nrow(data_pla)

# misclassified fragments (chromosome > 50%)
data_cut <- data_pla[data_pla$Chromosome / data_pla$Total > 0.5,]
frag_gt_50 <- nrow(data_cut)


# pie chart 1
perc1 <- round(frag_gt_0_pla / frag_num_pla * 100, 2)
perc2 <- 100 - perc1
lbls <- c(paste(perc1, "%"), paste(perc2, "%"))
colors <- c("red", "green")
pie(c(perc1, perc2), labels = lbls, main="Correctly Classified Plasmid Fragments", col=colors)
legend("topleft", c("yes", "no"), cex = 0.8, fill = colors, title=expression(bold(">0% chromosome")))

# pie chart 2
perc1 <- round(frag_gt_50 / frag_gt_0_pla * 100, 2)
perc2 <- 100 - perc1
lbls <- c(paste(perc1, "%"), paste(perc2, "%"))
colors <- c("red", rgb(51, 0, 102, max=255))
pie(c(perc1, perc2), labels = lbls, main="Potentially Ambiguous Plasmid Fragments", col=colors)
legend("topleft", c("yes", "no"), cex = 0.8, fill = colors, title=expression(bold(">50% chromosome")))

# histogram 1
hist(data_pla$Chromosome / data_pla$Total * 100, 
     main="Potentially Ambiguous Plasmid Fragments", 
     xlab="Chromosome (%)", 
     col="red",
     breaks=10)

# histogram 2
hist(data_pla$None / data_pla$Total * 100, 
     main="Potentially Ambiguous Plasmid Fragments", 
     xlab="None (%)", 
     col="red",
     breaks=10)



###################
###  ROC CURVE  ###
###################

# true and false negative - chromosomes
x <- c()
y <- c()
i = 0

while (i < 1.05) { 
  data_cut <- data_chromo[data_chromo$Chromosome / data_chromo$Total >= i,]
  num_cut <- nrow(data_cut)
  perc1 <- round((num_cut + frag_num_chromo-frag_gt_0_chromo) / frag_num_chromo * 100, 2)
  x <- c(x, i*100)
  y <- c(y, perc1)
  i <- i + 0.05
}

plot(x, y, ylim=c(0,100), xlab='Xc', ylab='% classified as chromosome', type = "o")
lines(x, y, type = "o", col = "blue")

TN <- y


x <- c()
y <- c()
i = 0

while (i < 1.05) { 
  data_cut <- data_pla[data_pla$Chromosome / data_pla$Total >= i,]
  num_cut <- nrow(data_cut)
  perc1 <- round(num_cut / frag_num_pla * 100, 2)
  x <- c(x, i*100)
  y <- c(y, perc1)
  i <- i + 0.05
}

lines(x, y, type = "o", col = "red")

FN <- y


# true and false positive - plasmids
x <- c()
y <- c()
i = 0

while (i < 1.05) { 
  data_cut <- data_pla[data_pla$Plasmid / data_pla$Total >= i,]
  num_cut <- nrow(data_cut)
  perc1 <- round((num_cut + frag_num_pla-frag_gt_0_pla) / frag_num_pla * 100, 2)
  x <- c(x, i*100)
  y <- c(y, perc1)
  i <- i + 0.05
}

plot(x, y, ylim=c(0,100), xlab='Xp', ylab='% classified as plasmid', type = "o")
lines(x, y, type = "o", col = "red")

TP <- y


x <- c()
y <- c()
i = 0

while (i < 1.05) { 
  data_cut <- data_chromo[data_chromo$Plasmid / data_chromo$Total >= i,]
  num_cut <- nrow(data_cut)
  perc1 <- round(num_cut / frag_num_chromo * 100, 2)
  x <- c(x, i*100)
  y <- c(y, perc1)
  i <- i + 0.05
}

lines(x, y, type = "o", col = "blue")

FP <- y


# ROC curve (plasmids being true positive)
sensitivity <- 1 - TP / (TP + FN)
specificity <- 1 - TN / (TN + FP)

plot(specificity, sensitivity, type = "n", main="ROC Curve for Plasmids", xlab="True Positive Rate", ylab="True Negative Rate")
lines(specificity, sensitivity, col = "purple", lwd=5)

