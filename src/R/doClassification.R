library(IRanges)
library(plyr)
library(ggplot2)
library(ape)
# need to figure out which actually are still needed
library(Biostrings)
library(biomartr)
library(slcp)
source("myDoBlast.R")

ref_fna_file <- "./Gamma_plasmids.fna"
ref_meta_file <- "./Gamma_plasmids_meta.txt"
ref_meta <- read.table(ref_meta_file, sep="\t", quote="", header=T, comment.char="", stringsAsFactors=F)
rownames(ref_meta) <- ref_meta$Accession


download_tempdir <- tempdir()
##### class-... #####
files <- list.files(path="./database/class-plasmid", pattern="*.fna", full.names=TRUE, recursive=FALSE)

out <- data.frame(Contig=character(),
                  Total=integer(), 
                  Plasmid=integer(), 
                  Chromosome=integer(),
                  None=integer()) 

for (f in 1:length(files)) { 
  name <- sapply(strsplit(files[f], "/"), `[`, 4)
  name <- substr(name, 1, nchar(name)-4)

  ##### cutoff as last argument (default 100) #####
  blastn_data <- myDoBlast(files[f], name=name, database=ref_fna_file, 99)
  if (nrow(blastn_data) == 0) next
  summary_class <- ddply(blastn_data, .(Query, Contig), contigSummary, .progress="text")
  ##### Plasmid > 0.5 | Plasmid > 0 | Chromosome > 0.5 | Chromosome > 0 #####
  p_contigs <- subset(summary_class[order(summary_class$Length, decreasing=T),], Plasmid > 0)$Contig
  if (length(p_contigs) == 0) next
  
  for (i in 1:length(p_contigs)) { 
    x <- contigDetail(subset(blastn_data, Contig==p_contigs[i]))
    full_IR <- IRanges(start = c(start(x[[1]]), start(x[[2]]), start(x[[3]])), end = c(end(x[[1]]), end(x[[2]]), end(x[[3]])))
    segments <- disjoin(full_IR)
    totalnt = max(end(x[[1]]))
    plasmidnt <- length(IRanges::which(x[[1]]=="Plasmid"))
    chromosoment <- length(IRanges::which(x[[1]]=="Chromosome"))
    nonent <- length(IRanges::which(x[[2]]==0))
    
    tmp <- data.frame(p_contigs[i], totalnt, plasmidnt, chromosoment, nonent)
    names(tmp) <- c('Contig', 'Total', 'Plasmid', 'Chromosome', 'None')
    out <- rbind(out, tmp)
  }
  print(f)
}
##### stats-... #####
write.table(out, "./database/stats/stats-plasmid-p-99.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)



# predict plasmid contigs and their best hit
query_fasta <- readDNAStringSet(file)
blastn_data <- myDoBlast(file, name=name, database=ref_fna_file)
# TODO: problem here because need to capture sequence names during db generation
summary_class <- ddply(blastn_data, .(Query, Contig), contigSummary, .progress="text")
#summary_classification(blastn_data, ref_meta)
stackPlot(summary_class, col=c("gray", "red", "blue"))  ##graph1
p_contigs <- subset(summary_class[order(summary_class$Length, decreasing=T),], Plasmid > 0.5)$Contig

#print statistics
for (i in 1:length(p_contigs)) { 
  print(contigDetail(subset(blastn_data, Contig==p_contigs[i])))
}

plotDetail(contigDetail(subset(blastn_data, Contig==p_contigs[1])))  ##graph2

# analyze against a single plasmid or assembly
blastn <- myDoBlast(file, name=name, database=ref_fna_file)
subject_meta <- ddply(blastn, .(Subject), summarize, Length=max(Slen), Type="None")
rownames(subject_meta) <- subject_meta$Subject
# go fix up the real types
subject_meta$Type <- c("Chromosome", "Plasmid", "None", "None", "None")
ref_meta <- subject_meta
summary_class <- ddply(blastn, .(Query, Contig), contigSummary, .progress="text")
plotDetail(subjectDetail(subset(blastn, Subject=="LN_...")), classification=summary_class)  ##graph3

# analyze against a plasmid in Genbank
# for ex CP024856.1
accession <- "LN_..."
plasmid_blastn <- doBlast(file, name=name, database=accession)
plotDetail(subjectDetail(plasmid_blastn))  ##graph4

# need to add gene annotation, from gff file, ? prokka to automate?
# need to add mouseover information
# dynamic zoom? shiny, D3