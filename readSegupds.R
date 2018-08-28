#!/usr/bin/env Rscript
#.libPaths(c(.libPaths(), "/net/eichler/vol2/eee_shared/modules/anaconda3/envs/py3.20161130/lib/R/library"))
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


ls




sd = fread("~/Desktop/work/assemblies/hg38/ucsc.unmerged.sorted.segdups.bed")
nr = fread("~/Desktop/work/assemblies/hg38/ucsc.collapsed.sorted.segdups")
meansd = fread("~/Desktop/work/assemblies/hg38/ucsc.merged.mean.segdups.bed")
maxsd = fread("~/Desktop/work/assemblies/hg38/ucsc.merged.max.segdups.bed")

maxsd = fread("~/Desktop/work/assemblies/Mitchell_CHM1_V3/internalassmbly/Segdups/asm.unresolved")
nr = fread("~/Desktop/work/assemblies/Mitchell_CHM1_V3/internalassmbly/Segdups/asm.col.unresolved")


df = nr
colnames(df) = c("chr", "start","end","id")
df$length = df$end -df$start
max(df$length)/10^6

ggplot(df) + geom_histogram(aes(log10(length)), bins =100)
