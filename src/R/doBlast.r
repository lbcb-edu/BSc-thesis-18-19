# database could have been passed a filename or an accession number
doBlast2 <- function(file, name="", database=ref_fna_file, cutoff=100) {
  
  # contains information about strains which consist of chromosome and plasmids
  ref_meta_file <- "./Gamma_plasmids_meta.txt"
  ref_meta <- read.table(ref_meta_file, sep="\t", quote="", header=T, comment.char="", stringsAsFactors=F)
  rownames(ref_meta) <- ref_meta$Accession

  if(!file.exists(database)) {
    # the only other option right now is we assume it's a Genbank accession
    tempfilename <- paste(download_tempdir, "/", database, ".fasta", sep="")
    write.dna(read.GenBank(database), file=tempfilename, format="fasta", nbcol=1, colw=60)
    if(!file.exists(tempfilename)) {
      return(NA)
    } else {
      db_file <- tempfilename
    }
  } else {
    # we have a file that was passed in
    db_file <- database
  }
  # check for blast indices
  if(!file.exists(paste(db_file, ".nin", sep="")) &
     !file.exists(paste(db_file, ".00.nin", sep=""))) {
    system(paste("makeblastdb -dbtype nucl -in", db_file, sep=" "))
  }
  evalue_cutoff <- 1e-10
  command <- paste(sep=" ",
                   "blastn -db", db_file,
                   "-query", file,
                   "-out - -evalue", evalue_cutoff,
                   "-dust no",
                   '-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"',
                   "-max_target_seqs 100")
  blastcon <- pipe(command, open="r")
  blastdata <- read.table(blastcon, header=F, stringsAsFactors=F, sep="\t", quote="")
  close(blastcon)
  
  colnames(blastdata) <- c("Contig", "Subject", "ID", "AlignLen", "Mismatch", "Gaps", "Qstart", "Qend", "Sstart", "Send", "Evalue", "Bitscore", "Qlen", "Slen")
  if (name == "") {
    blastdata$Query <- basename(file)
  } else {
    blastdata$Query <- name
  }
  
  blastdata <- blastdata[,c(ncol(blastdata),1:(ncol(blastdata)-1))]
  
  # filter out the self-hit
  blastdata <- blastdata[blastdata$Query != blastdata$Subject,]
  
  # default: filter out the chromosomes with identity higher than 'cutoff'
  blastdata <- blastdata[ref_meta[blastdata$Subject, 4] == "Chromosome" & blastdata$ID <= cutoff | ref_meta[blastdata$Subject, 4] == "Plasmid",]
  # filter out the plasmids with identity higher than 'cutoff'
  # blastdata <- blastdata[ref_meta[blastdata$Subject, 4] == "Plasmid" & blastdata$ID <= cutoff | ref_meta[blastdata$Subject, 4] == "Chromosome",]
  
  rownames(blastdata) <- seq(length=nrow(blastdata))

  return(blastdata)
}


# option using command line
args <- commandArgs(TRUE)
name <- sapply(strsplit(args[1], "/"), tail, 1)
name <- substr(name, 1, nchar(name)-4)

if (length(args) == 2) {
  print(doBlast2(file=args[1], name, database=args[2]))
  
} else if (length(args) == 3) {
  print(doBlast2(file=args[1], name, database=args[2], cutoff=as.numeric(args[3])))
  
} else {
  print("Incorrect number of arguments - must be 2 or 3!")
}
