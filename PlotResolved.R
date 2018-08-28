#setwd("/net/eichler/vol19/projects/CHM1_project/nobackups/denovo_assemblies/mhap/GenBank")

args <- commandArgs(trailingOnly = TRUE)

#args <- c("MHAP.published.GCA_000772585.3.fasta.resolved", "MHAP.published.GCA_000772585.3.fasta.unresolved", "MHAP", "test.pdf")

#args <- c("CHM1P6_v030_070102015_pa.consensus.fasta.resolved", "CHM1P6_v030_070102015_pa.consensus.fasta.unresolved", "CHM1.v27", "test.pdf")
res <- read.table(args[1], header=F)
unr <- read.table(args[2], header=F)



name <- args[3]
plotfile1 <- sprintf("%s-SegdupResolutionByDensity.pdf", args[4])
plotfile2 <- sprintf("%s-SegdupResolutionByPoint.pdf", args[4])

res$l <- res$V3 - res$V2
res$logl <- log(res$l)
unr$l <- unr$V3 - unr$V2
unr$logl <- log(unr$l)
minq <- .00001
res$q <- round(-10*log10(1-res$V4/100 + minq))
unr$q <- round(-10*log10(1-unr$V4/100 + minq))

#res$x <- res$l
#res$y <- res$q

ml <- max(res$l,unr$l)

maxX <- 2*10^floor(log10(ml))

#minI <- min(res$V4, unr$V4)
#maxI <- max(res$V4, unr$V4)

minI <- min(res$q, unr$q)
maxI <- max(res$q, unr$q)

#r <- lapply(seq(minI,maxI+0.01, by=0.01), function(x) which(res$V4 == x))
#u <- lapply(seq(minI,maxI+0.01, by=0.01), function(x) which(unr$V4 == x))
require(RColorBrewer)
library(ggplot2)
p <- brewer.pal(3, "Set1")

pdf(plotfile1, family="Helvetica")

green <- paste(p[3], "55", sep="")
black <- "#00000099"
red <- "#ff000099"
#red <- paste(p[1], "55", sep="")
#pdf("MHAP-GenBank.SegDupResolutionByDensity.pdf")

#ggplot(rbind(data.frame(res, group="Resolved"), data.frame(unr, group="Unresolved")), aes(x=l,y=q)) + geom_density2d(aes(colour=group)) + scale_colour_brewer(type="sequential", palette="Set1") + theme_minimal() + xlab("Length of duplication") + ylab("Percent identity (PHRED)") + ggtitle(args[3]) + ylim(5,30) + scale_x_log10()
ggplot(rbind(data.frame(res, group="Resolved"), data.frame(unr, group="Unresolved")), aes(x=l,y=V4)) + geom_density2d(aes(colour=group)) + scale_colour_brewer(type="sequential", palette="Set1") + theme_minimal() + xlab("Length of duplication") + ylab("SegDup percent identity") + ggtitle(args[3]) + scale_x_log10()
dev.off()
#x11(type="dbcairo")

xScale <- c("10","100","1k","10k","100k","1000K","10M","100M","1G","10G","100G")

#setwd("/net/eichler/vol19/projects/CHM1_project/nobackups/denovo_assemblies/falcon.v27/SegdupAnalysis")

#pdf("MHAP-GenBank.SegDupResolutionByPoint.pdf")
pdf(plotfile2, family="Helvetica")
#x11(type="dbcairo")
minIdent = 90
maxIdent = max(res$V4, unr$V4)
options(scipen=100)
minX <- 1000
xlabelIndices <- log10(minX):log10(maxX)


plot(c(), xlim=c(minX,maxX), ylim=c(minIdent, maxIdent), xaxt='n', ylab="Percent sequence identity ", xlab="Segmental duplication length (bp)", main=name, log="x", las=1)

axis(1,at=10^xlabelIndices, labels=xScale[xlabelIndices])


identRange <- ceiling(maxIdent) - floor(minIdent)

#
#

plotp <- function(r, bin, offset, color) {
  i <- which(r$V4 > (bin) & r$V4 < (bin+step));
  points(r$l[i], offset+ bin+((r$V4[i]-bin)*(step*0.5)), col=color, pch=16,cex=0.5);
}

bin <- 90
step=1
lapply(90:99, function(bin) plotp(res, bin, 0, black))
lapply(90:99, function(bin) plotp(unr, bin, 0.5, red))


legend("bottomright", legend=c("resolved", "unresolved"), col=c("black", "red"), pch=19)
#
dev.off()

