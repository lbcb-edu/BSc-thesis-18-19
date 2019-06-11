library(plyr)
library(slcp)
source("doBlast.R")

ref_fna_file <- "Gamma_plasmids.fna"
ref_meta_file <- "Gamma_plasmids_meta.txt"
ref_meta <- read.table(ref_meta_file, sep="\t", quote="", header=T, comment.char="", stringsAsFactors=F)
rownames(ref_meta) <- ref_meta$Accession
download_tempdir <- tempdir()

# location where fragmented chromosomes/plasmids are stored
files <- list.files(path="./database/class-chromosome", pattern="*.fna", full.names=TRUE, recursive=FALSE)

out <- data.frame(Contig=character(),
                  Total=integer(), 
                  Plasmid=integer(), 
                  Chromosome=integer(),
                  None=integer()) 

for (f in 1:length(files)) { 
  name <- sapply(strsplit(files[f], "/"), tail, 1)
  name <- substr(name, 1, nchar(name)-4)
  
  blastn_data <- doBlast2(files[f], name=name, database=ref_fna_file)
  if (nrow(blastn_data) == 0) next
  summary_class <- ddply(blastn_data, .(Query, Contig), contigSummary, .progress="text")
  p_contigs <- subset(summary_class[order(summary_class$Length, decreasing=T),], Plasmid > 0)$Contig
  if (length(p_contigs) == 0) next
  
  for (i in 1:length(p_contigs)) { 
    x <- contigDetail(subset(blastn_data, Contig==p_contigs[i]))
    full_IR <- IRanges(start = c(start(x[[1]]), start(x[[2]]), start(x[[3]])), end = c(end(x[[1]]), end(x[[2]]), end(x[[3]])))
    segments <- disjoin(full_IR)
    totalnt <- max(end(x[[1]]))
    plasmidnt <- length(IRanges::which(x[[1]]=="Plasmid"))
    chromosoment <- length(IRanges::which(x[[1]]=="Chromosome"))
    nonent <- length(IRanges::which(x[[2]]==0))
    
    tmp <- data.frame(p_contigs[i], totalnt, plasmidnt, chromosoment, nonent)
    names(tmp) <- c('Contig', 'Total', 'Plasmid', 'Chromosome', 'None')
    out <- rbind(out, tmp)
  }
}

# location where the classification statistics are stored
write.table(out, "./database/stats/stats-chromosome.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
