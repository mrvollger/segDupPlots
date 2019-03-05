#!/usr/bin/env Rscript
#.libPaths(c(.libPaths(), "/net/eichler/vol2/eee_shared/modules/anaconda3/envs/py3.20161130/lib/R/library"))
#install.packages("devtools")
#devtools::install_github("daattali/ggExtra")
library(ggExtra)
library(ggplot2)
library(plyr)
require(gridExtra)
require(scales)
library(RColorBrewer)
library(reshape2)
#install.packages("tidyr")
library(tidyr)
library(dplyr)
library(stringr)
library(data.table)
#install.packages("argparse")
library(argparse)


#suppressPackageStartupMessages(library("argparse"))
# create defualt files to run
genome = "Mitchell_CHM1"
genome = "Mitchell_CHM1_V3"
genome = "CHM13"
genome = "Yoruban_feb_2018"
genome = "AK1"
genome = "NA12878"
genome = "HX1"
genome = "Mitchell_CHM1_V2"


euchro <- Sys.glob(sprintf("~/Desktop/work/assemblies/%s/*/Segdups/*.euchromatic.all", genome) )[1]
all <- Sys.glob(sprintf("~/Desktop/work/assemblies/%s/*/Segdups/*.all", genome) )[1]
dest <- sprintf("~/Desktop/data/genomeWide/%s/plots/SegDupPlots/", genome)
gene <- Sys.glob(sprintf("~/Desktop/work/assemblies/%s/*/Segdups/*.unresolved.gene.intersection", genome) )[1]

gene

# create parser object
parser <- ArgumentParser()
parser$add_argument("-a", "--all", default=all, help=" list of all seg dups")
parser$add_argument("-g", "--gene", default=gene, help=" list of all seg dups")
parser$add_argument("-e", "--euchro", default=euchro, help=" list of euchro segdups")
parser$add_argument("-d", "--dest", default=dest, help="destination directory")

args <- parser$parse_args()
args
# make the output dir if not already there
system(paste0("mkdir -p ", args$dest))

#
# ploting setup
#
black <- "#000000"
red <- "#FF0000"
green <- "#228B22"
myColors = c(green, black, red)
names(myColors) <- levels(as.factor(c("unresolved", "resolved", "all")))
colScale <- scale_colour_manual(name = "Status",values = myColors)
myColors
xkb=c("1", "10", "100", "1,000")
xbreaks=c(1000, 10000, 100000, 1000000)
ybreaks=seq(90,100,1)
theme_set(theme_gray(base_size = 18))
h=20
w=30
myTheme =theme(plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
               text = element_text(),
               panel.background = element_blank(),
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               plot.background = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(1)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(), 
               axis.line = element_line(colour="black"),
               axis.ticks = element_line(),
               #panel.grid.major = element_line(colour="#f0f0f0"),
               #panel.grid.minor = element_blank(),
               plot.margin=unit(c(10,5,5,5),"mm"),
               #strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
)
mysave <- function(name, p){
  if(args$dest != ""){
    args$dest = paste0(args$dest, "/")
    ggsave(paste0(args$dest, name), plot = p,  width = w, height = h, units = "cm")
  }
}

#
# data setup
#
print("reading data")
segdups = read.table(args$all)
euchro = read.table(args$euchro)
ldfs = list(segdups, euchro)
dfs = lapply(ldfs, function(df) {
  colnames(df) = c("chr", "start", "end", "perID", "dist")
  df$length = df$end - df$start
  df$Status = "unresolved"
  df$perID = df$perID * 100
  df$perID[df$perID == 100] = 99.99999999
  df = df[df$perID >= 90, ]
  return( df )
})
segdups = dfs[[1]]
euchro = dfs[[2]]
print("done reading data")

#
# creating plots 
#
segs = c()
vals = seq(0, 100000,10000)
for( minDist in vals){
  #print(minDist)
  segdups$Status = "unresolved"
  segdups$Status[segdups$dist >= minDist] = "resolved"
  segdups$dec=segdups$perID - floor(segdups$perID)
  segdups$dec[segdups$Status == "resolved"] = rescale( segdups$dec[segdups$Status == "resolved"] , to=c(0,1/2), from=c(0,1) )
  segdups$dec[segdups$Status == "unresolved"] = rescale( segdups$dec[segdups$Status == "unresolved"] , to=c(1/2,1), from=c(0,1) )
  segdups$scaled = floor(segdups$perID) + segdups$dec

  p1 <- ggplot(segdups, aes(x=length, y=scaled, color=Status) ) + geom_point(size=2) + colScale +
    scale_x_continuous(trans='log10',labels = xkb, breaks = xbreaks) +
    scale_y_continuous(breaks = ybreaks) + 
    xlab("Segmental duplication length (kbp)") + ylab("Percent sequence identity") + myTheme + theme(legend.position="none")
  p1=ggMarginal(p1, type="histogram", alpha=.8, groupFill=T, yparams = list(bins=10), xparams = list(bins=15), size=6);
  mysave(paste0(minDist, ".pdf"), p1)
  
}
#system("convert ~/Desktop/gif/*.png -delay 3 -loop 0 ~/Desktop/gif/resolved.gif")

myxlab = "Minimal extension from duplication boundary (kbp)"

Extra = seq(0,1000000,1000)
PercentResolved = c()
for(minDist in Extra ){
  euchro$Status = "unresolved"
  euchro$Status[euchro$dist >= minDist] = "resolved"
  un = sum(euchro$length[euchro$Status == "unresolved" ])
  r = sum(euchro$length[euchro$Status == "resolved" ])
  PercentResolved = c(PercentResolved, r*100.0/(un+r))
}
perE = data.frame(Extra,PercentResolved)
p = ggplot(perE) + geom_point(aes(Extra/1000, PercentResolved)) + myTheme +
  ylab("Percent of euchromatic segmental duplications bases resolved") +
  xlab(myxlab) +  scale_x_continuous(labels = comma)
p
mysave("EchroZeroToOneMillion.pdf", p)

p2 = p + xlim(0,250)
p2
mysave("EuchroZeroTo250k.pdf", p2)



PercentResolved = c()
for(minDist in Extra ){
  segdups$Status = "unresolved"
  segdups$Status[segdups$dist >= minDist] = "resolved"
  un = sum(segdups$length[segdups$Status == "unresolved" ])
  r = sum(segdups$length[segdups$Status == "resolved" ])
  PercentResolved = c(PercentResolved, r*100.0/(un+r))
}
perE = data.frame(Extra,PercentResolved)
maxy = perE[perE$Extra==50000, "PercentResolved"]
pos = 3
p = ggplot(perE) + geom_point(aes(Extra/1000, PercentResolved), color = "darkblue") + myTheme +
  ylab("Percent of segmental duplications resolved") +
  xlab(myxlab) +  scale_x_continuous(labels = comma) +
  geom_segment( aes(x=50, xend=50, y = 0, yend = maxy), linetype="dashed", color = "red", size=.5)
mysave("zeroToOneMillion.pdf", p)

mymax = 250
p2 = p + xlim(0,mymax)  #+
  #geom_segment(aes(x=50, xend=0    , y=maxy/pos, yend=maxy/pos), size = 1.5, color = "darkred", arrow = arrow(length = unit(0.5, "cm"))) +
  #geom_segment(aes(x=50, xend=mymax, y=maxy/pos, yend=maxy/pos), size = 1.5, arrow = arrow(length = unit(0.5, "cm"))) +
  #annotate("text", label = "Unresolved", x =(0+50)/2, y = maxy/pos*1.1, size = 5, color = "darkred") +
  #annotate("text", label = "Resolved", x =(mymax+50)/2, y = maxy/pos*1.1, size = 5)

p2; mysave("zeroTo250k.pdf", p2)





#
# read in genes
#
print(args$gene)
if(!(is.na(args$gene)) & (args$gene != "NA")){

genes = fread(args$gene)
colnames(genes) = c("chr", "start", "end", "gene", "chr2", "start2", "end2", "perID")
genes$perID = genes$perID*100
genes = genes[order(-genes$perID),]
dups = duplicated(genes, by=c("gene"))
genes = genes[!dups]
print(dim(genes))
perid = genes$perID
status = rep("gene", length(perid))
gene_ecdf = data.table(perid, status)
#
#
#

segdups$Status == "unresolved"
segdups[segdups$dist >= 50000,]$Status = "resolved"
divider  = 100
perid = rep(0.0, sum(segdups$length)/divider )
status = rep("Status", sum(segdups$length)/divider )
curidx=1
for(i in 1:dim(segdups)[1]){
  row = segdups[i,]
  mynum = row$length/divider
  perid[curidx:(curidx+mynum-1)] =  row$perID
  status[curidx:(curidx+mynum-1)] = row$Status 
  curidx = curidx + mynum
}
df_ecdf = data.table(perid,status)
df_ecdf = df_ecdf[df_ecdf$perid != 0.0]
df_ecdf
curidx
length(df_ecdf)
summary(df_ecdf)
df_ecdf = df_ecdf[df_ecdf$status == "unresolved"]
df_ecdf = rbind(df_ecdf, gene_ecdf)
col2 = c("blue", "red")
col2 = c("#CE0000", "#ff9a00")
names(col2) = c("unresolved", "gene")

y1 = "Mb of unresolved SDs"
y2 = "Number of genes in unresolved SDs"
breaks = c(0, 0.25, .5, .75, 1.0)
y1lables = round( breaks * length(df_ecdf[df_ecdf$status=="unresolved"]$perid)/10000,2 )
y2lables = round( breaks * length(df_ecdf[df_ecdf$status=="gene"]$perid),0 )
p3.0 = ggplot(df_ecdf, aes(perid, color = status)) + stat_ecdf(geom = "line", size =2, alpha = 0.65) + 
  ylab(y1) +
  xlab("Percent identity of duplication block") +
  scale_x_continuous(limits=c(min(df_ecdf$perid),100), expand = c(0.0, 0.1)) +
  theme(axis.text.y.right = element_text(color = col2[2]))+
  theme(axis.text.y = element_text(color = col2[1]))+
  scale_color_manual(values=col2) +
  theme(legend.position="none") + 
  scale_y_continuous(y1, breaks = breaks, labels = y1lables,
                     sec.axis = sec_axis(~.*1, name = y2, breaks = breaks, labels = y2lables)) +
  myTheme 

p3.0; mysave("MBAndGenesUnresolvedByPerID.pdf", p3.0)




#p3.1 = p3.0 + coord_cartesian(xlim=c(99.8, 100) )
#p3.1

df_gene = df_ecdf[df_ecdf$status=="gene",] 
y2lables = c(0,500, 1000, 1500, 2000, 2500)
perc.rank <- ecdf(1:nrow(df_gene))
breaks = perc.rank(y2lables)
p4.0 = ggplot(df_gene, aes(perid)) + stat_ecdf(geom = "line", size =2, alpha = 0.65) + 
  ylab(y2) +
  xlab("Percent identity of duplication block") +
  scale_x_continuous(limits=c(min(df_ecdf$perid),100), expand = c(0.0, 0.1)) +
  #theme(axis.text.y.right = element_text(color = col2[2]))+
  #theme(axis.text.y = element_text(color = col2[1]))+
  scale_color_manual(values=col2) +
  theme(legend.position="none") + 
  scale_y_continuous(y2, breaks = breaks, labels = y2lables )+
                     #, sec.axis = sec_axis(~.*1, name = y2, breaks = breaks, labels = y2lables)) +
  myTheme 

p4.0; mysave("MBunresolvedByPerID.pdf", p4.0)

}





p1


p7=ggMarginal(p1, type="histogram", alpha=.8, groupFill=T, yparams = list(bins=10), xparams = list(bins=15), size=6); p7
