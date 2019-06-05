# input total number of fragments
frag_num <- 84019

# fragments which are not 100% plasmid (chromosome > 0%)
data <- read.table("stats-plasmid.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
frag_gt_0 <- nrow(data)

# misclassified fragments (chromosome > 50%)
data_cut <- data[data$Chromosome / data$Total > 0.5,]
frag_gt_50 <- nrow(data_cut)


# graph 1
slices <- c(frag_num - frag_gt_0, frag_gt_0)
perc1 <- round((frag_num - frag_gt_0) / frag_num * 100, 2)
perc2 <- 100 - perc1
lbls <- c(paste(perc1, "%"), paste(perc2, "%"))
colors <- c("green", "red")
pie(slices, labels = lbls, main="Classification of Plasmids", col=colors)
legend("topleft", c("yes", "no"), cex = 0.8, fill = colors, title=expression(bold("100% plasmid")))

# graph 2
slices <- c(frag_gt_50, frag_gt_0 - frag_gt_50)
perc1 <- round(frag_gt_50 / frag_gt_0 * 100, 2)
perc2 <- 100 - perc1
lbls <- c(paste(perc1, "%"), paste(perc2, "%"))
colors <- c("red", rgb(51, 0, 102, max=255))
pie(slices, labels = lbls, main="Misclassified Plasmids", col=colors)
legend("topleft", c("yes", "no"), cex = 0.8, fill = colors, title=expression(bold(">50% chromosome")))

# graph 3
hist(data$Chromosome / data$Total * 100, 
     main="Misclassified Plasmids", 
     xlab="Chromosome (%)", 
     col="red",
     breaks=10)

# graph 4
hist(data$None / data$Total * 100, 
     main="Misclassified Plasmids", 
     xlab="None (%)", 
     col="red",
     breaks=10)
