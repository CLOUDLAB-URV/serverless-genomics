#!/usr/bin/env Rscript --vanilla
# R script called by redis_plots.sh

setwd("~/bioss/serverless_genomics/variant_caller/")
args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(ggplot2)
library(plyr)
library(plotly)

print_png <- function(file_name, suffix, p) {
  ggsave(filename = paste0(file_name, suffix, ".png"), 
         device = "png", 
         plot = p,
         width = 10,
         height = 4.5,
         dpi = 1200
  )
}


debug <- FALSE 

# PARSE COMMAND LINE----
if (debug==FALSE) {
  redis_start_file <- args[1]
  failed_functions_file <- args[2]
  redis_exec_file <- args[3]
} else {
  dir="varcall_out/"
  name="varcall_56_1_ERR9729866_1000Kfq_100Mfa_500itn_rep1.log_functions_sorted.log"
  redis_start_file <- paste0(dir,name,"_redis_start_times.txt")
  failed_functions_file <- paste0(dir,name,"_redis_failed_functions.txt")
  redis_exec_file <- paste0(dir,name,"_redis_execution_time.txt")
}
  
# REDIS START TIME----

redis_start <- read.table(redis_start_file, header=F, fill= TRUE, sep= "\t", stringsAsFactors = FALSE)
colnames(redis_start) <- c("FUNCTION", "TIME")

failed_functions_file<- paste0(dir,"varcall_57_1_ERR9729866_1300Kfq_100Mfa_3001itn_rep1.log_functions_sorted.log_redis_failed_functions.txt")
failed_functions <- read.table(failed_functions_file, header=F, fill= TRUE, sep= "\t", stringsAsFactors = FALSE)
tryCatch(read.table(x, header = TRUE, sep = '|'), error=function(e) NULL)
colnames(failed_functions) <- c("FUNCTION")


if fa
# filter dataframe to get data to be highlighted (failed functions)
highlight_redis <- subset(redis_start, FUNCTION %in% failed_functions$FUNCTION)

p <- redis_start %>% ggplot(aes(x = FUNCTION, y = TIME)) + 
  geom_point(alpha=0.3) +
  geom_point(data=highlight_redis, 
             aes(x = FUNCTION, y = TIME), 
             color='red',
             size=3) +
  labs(title="Redis start time - map function",
       x ="", y = "") 
#p
#ggplotly(p)
print_png(redis_start_file, "", p)

# REDIS EXECUTION TIME----
redis_exec <- read.table(redis_exec_file, header=F, fill= TRUE, sep= "\t", stringsAsFactors = FALSE)
colnames(redis_exec) <- c("FUNCTION", "TIME")
p <- redis_exec %>% ggplot(aes(x = FUNCTION, y = TIME)) + 
  geom_point(alpha=0.3) +
  labs(title="Redis execution time - map function",
       x ="", y = "seconds") 
p
print_png(redis_exec_file, "", p)


