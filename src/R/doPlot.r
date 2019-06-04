#!/usr/bin/Rscript

args <- commandArgs(TRUE)

title = ifelse(args[2]=="p", "Plasmid", "Chromosome")
if(args[1]=="s") {
	title = paste(title, args[3], sep=" ")
}

data <- read.table("temp.txt", header=TRUE, sep=" ", stringsAsFactors=FALSE)
attach(data)

# optional - if working with X11 forwarding
X11()

switch(args[4], 
"s"={
	plot(Identity, log10(Hits), ylab="Hit Length (log10)", main=title, pch=18, col=ifelse(Type=="Plasmid", "Red", "Blue"))
	legend("topleft", inset=.02, legend=c("Plasmid", "Chromosome"), col=c("red", "blue"), lty=c(0,0), pch=c(18,18))
},
"m"={
  library(ggplot2)
  library(ggExtra)
  
  Type = ifelse(Type=="Plasmid", "Plasmid", "Chromosome")
  
  p <- ggplot(data, aes(Identity, log10(Hits), col=Type)) +
    ggtitle(title) +
    ylab("Hit Length (log10)") +
    geom_point() +
    theme(legend.position = c(0, 1), legend.justification = c(0, 1), plot.title = element_text(hjust=0.5, size=12, face="bold")) +
    scale_color_manual(values=c("Blue", "Red"))
  
  ggMarginal(p, type="histogram", groupFill=T)
}
)

message("Press Return To Continue")
invisible(readLines("stdin", n=1))

