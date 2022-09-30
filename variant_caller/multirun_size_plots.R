#!/usr/bin/env Rscript --vanilla
# R script to process output of parse_file_size.sh

setwd("~/bioss/serverless_genomics/variant_caller/")

library(tidyverse)
library(ggplot2)
library(plyr)
library(plotly)
library(rgl)
library(viridis)

dir="varcall_out/"

# GRAPH WITH FIXED SUBSET (FASTA OR FASTQ) ----
filename<-"csv_sizes_1000fq_subset.txt"

sizes <- read.table(paste0(dir,filename), header=F, fill= TRUE, sep= "\t", stringsAsFactors = FALSE)
colnames(sizes) <- c("FILE", "FASTQ", "FASTA", "TARGET_FILE_TYPE", "TARGET_FILE_SIZE")
file_type=sizes$TARGET_FILE_TYPE[1]
# convert file size to Mb
sizes$TARGET_FILE_SIZE <- sizes$TARGET_FILE_SIZE / 1000000

sizes_part <- sizes %>%
  filter(FILE=="varcall_8_1_SRR15068323_1000Kfq_100Mfa__rep1.log")

# check if fixed fastq file size
if (length(unique(sizes$FASTQ))==1) {
  fastq_size=sizes$FASTQ[1]
  p <- sizes %>% ggplot(aes(x = FASTA, y = TARGET_FILE_SIZE, colour=FILE)) + 
    geom_point(alpha=0.3) +
    labs(title=paste0(file_type, " file sizes for ",fastq_size,"Kfa fastq file"),
                      x ="fasta file size Mb", y = "csv file size - Mb") 
  p
  ggplotly(p)
  p <- sizes_part %>% ggplot(aes(x = FASTA, y = TARGET_FILE_SIZE, colour=FILE)) + 
    geom_point(alpha=0.3) +
    labs(title=paste0(file_type, " file sizes for ",fastq_size,"Kfa fastq file"),
         x ="fasta file size Mb", y = "csv file size - Mb") 
  p
}

# GRAPH WITH ALL DATA
filename<-"csv_sizes_fq_subset.txt"
sizes_all <- read.table(paste0(dir,filename), header=F, fill= TRUE, sep= "\t", stringsAsFactors = FALSE)
colnames(sizes_all) <- c("FILE", "FASTQ", "FASTA", "TARGET_FILE_TYPE", "TARGET_FILE_SIZE")
file_type=sizes_all$TARGET_FILE_TYPE[1]
# convert file size to Mb
sizes_all$TARGET_FILE_SIZE <- sizes_all$TARGET_FILE_SIZE / 1000000

#colours for plot
# Add a new column with color
mycolors <- viridis(length(unique(sizes_all$FILE)))
sizes_all$color <- mycolors[ as.numeric(as.factor(sizes_all$FILE)) ]

plot3d( 
  x=sizes_all$FASTQ, y=sizes_all$FASTA, z=sizes_all$TARGET_FILE_SIZE, 
  col = sizes_all$color, 
  type = 's', 
  radius = .1,
  xlab="fastq number of reads", ylab="fasta file size Mb", zlab=paste0(file_type, " size Mb"))


data <- iris

# Add a new column with color
mycolors <- c('royalblue1', 'darkcyan', 'oldlace')
data$color <- mycolors[ as.numeric(data$Species) ]

# Plot
plot3d( 
  x=data$`Sepal.Length`, y=data$`Sepal.Width`, z=data$`Petal.Length`, 
  col = data$color, 
  type = 's', 
  radius = .1,
  xlab="Sepal Length", ylab="Sepal Width", zlab="Petal Length") 
