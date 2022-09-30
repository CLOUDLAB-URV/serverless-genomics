#!/usr/bin/env Rscript --vanilla

#setwd("~/bioss/serverless_genomics/variant_caller/")
args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(ggplot2)
library(plyr)
count <- dplyr::count

#FUNCTIONS----
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

print_png <- function(file_name, suffix, p) {
  ggsave(filename = paste0(file_name, suffix, ".png"), 
         device = "png", 
         plot = p,
         width = 10,
         height = 4.5,
         dpi = 1200
  )
}

## http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#Helper%20functions
## Summarizes data. 
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {

  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}



create_time_plot <- function(times, type, suffix, debug) {
  if (type=="Average") {
    p <- ggplot(times) + 
    geom_point(mapping = aes(x = COMMAND, y = DATA, colour = STAGE), size=4) +
    geom_errorbar(aes(x = COMMAND, ymin=DATA-se, ymax=DATA+se), width=.1, color="#808080") +
    labs(title=paste0("Variant Caller Execution Time -",type),
         subtitle = paste0("using ",map_functions, " map functions, and ",reduce_functions," reduce functions"),
         caption = command,
         x ="", y = "seconds") +
    scale_x_discrete(labels=function(x){sub("\\s", "\n", x)})+
    theme_bw() + 
    theme(panel.border = element_blank(), 
          #panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))+
    theme(legend.position="none") +
    scale_color_manual(values = c("PP" = "#D64324", #red
                                  "A" = "#F59513", # orange..
                                  "C" = "#228B22")) # green 2
  }
  else {
    p <- ggplot(times) + 
      geom_point(mapping = aes(x = COMMAND, y = DATA, colour = STAGE), size=4) +
      #geom_errorbar(aes(x = COMMAND, ymin=DATA-se, ymax=DATA+se), width=.1, color="#808080") +
      labs(title=paste0("Variant Caller Execution Time -",type),
           subtitle = paste0("using ",map_functions, " map functions, and ",reduce_functions," reduce functions"),
           caption = command,
           x ="", y = "seconds") +
      scale_x_discrete(labels=function(x){sub("\\s", "\n", x)})+
      theme_bw() + 
      theme(panel.border = element_blank(), 
            #panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"))+
      theme(legend.position="none") +
      scale_color_manual(values = c("PP" = "#D64324", #red
                                    "A" = "#F59513", # orange..
                                    "C" = "#228B22")) # green 2
    
  }
  if (debug == FALSE) {
    print_png(times_file, paste0("_",suffix), p)
    #return(p)
  }
  else {
  return(p)
  }
}


debug <- FALSE 

# PARSE COMMAND LINE----
if (debug==FALSE) {
  times_file <- args[1]
  sizes_file <- args[2]
  command <- args[3]
} else {
  directory="varcall_out/"
  name="varcall_93_1_ERR9856489_1300Kfq_100Mfa__"
  name="varcall_92_2_ERR9856489_1300Kfq_100Mfa_18itn_"
  times_file <- paste0(directory,name, "exec_time.log")
  sizes_file <- paste0(directory,name,"file_size.log")
  command <- "python varcall_lithops_demo_v7.py -fq ERR9729866 -nfq 1300000 -nfa 100000000 -itn 1000 -cf 10000 -fa hg19.fa -cl aws -b cloudbutton-variant-caller-input -fb ayman-lithops-meta-cloudbutton-hutton -ds SRA -ofa 300 -rl 152 -t 0 -ff csv -s3w False -rt lumimar/hutton-genomics-v03:18 -rtm 4096 -rtr 4096 -bs 75% -ftm 900 -ftr 900 -sk False -lb manual -ip 54.146.89.181 -id i-0cee52f66655d990b -rg us-east-1"
  debug <- TRUE
}

# add newlines to command every x characters
command = gsub("(.{100,}?)\\s", "\\1\n", command)



print(paste0("process_execution_times.R","\n",".....processing execution times for ", times_file))


# PROCESS TIME DATA----

# 1. format table ----
times <- read.table(times_file, header=F, fill= TRUE, sep= "\t", stringsAsFactors = FALSE)
colnames(times) <- c("ITERATION", "STAGE", "FUNCTION_N", "COMMAND", "DATA_INFO", "DATA", "DATA_UNITS")
# remove leading and trailing white space
times[1:7] <- lapply(times[1:7], function(x) trimws(x, which = c("both")))
times$DATA <- as.numeric(times$DATA)
# replace map stage number with A alone
times$STAGE <- gsub('[0-9]*','',times$STAGE)
# create levels for ggplot
times$STAGE <- factor(times$STAGE,levels = c("PP", "A", "B", "C", "D", "MR"))
times <- times[order(times$STAGE), ]
times$COMMAND <- factor(times$COMMAND,levels = unique(times$COMMAND))

# 2. check table contents----
count(times, STAGE)
count(times, DATA_INFO)
count(times, COMMAND)
# check all
times %>% 
  count(STAGE, COMMAND, sort=TRUE) 
times %>% 
  count(COMMAND, DATA_INFO, sort=TRUE) 
# total length
start_end_times <- times %>% 
  filter(!(grepl("execution", DATA_INFO)))
min_time <- min(start_end_times$DATA)
max_time <- max(start_end_times$DATA)
run_length=(max_time-min_time)/60

# check pre-processing
times %>% 
  filter(STAGE=="PP") %>%
  count(COMMAND, DATA_INFO, sort=TRUE) 
# check map function times
map_counts <- times %>% 
  filter(STAGE=="A") %>%
  count(COMMAND, DATA_INFO, sort=TRUE) 
# check reduce
reduce1 <- times %>% 
  filter(STAGE=="B")
reduce3 <- times %>% 
  filter(STAGE=="D")
map_reduce <- times %>% 
  filter(STAGE=="MR")
# check reduce function times
reduce2_counts <- times %>% 
  filter(STAGE=="C") %>%
  count(COMMAND, DATA_INFO, sort=TRUE) 


# number of functions 
map_functions <- map_counts[1,3]
print(paste0("number of map functions: ", map_functions))
reduce_functions <- reduce_counts[1, 3]
print(paste0("number of reduce functions: ", reduce_functions))


# 3. prepare times tables----

# overall execution times for preprocessing, map and reduce
times_main <- times %>%
  filter((grepl("PP", STAGE) & grepl("execution", DATA_INFO))| grepl("total", DATA_INFO) )
# check
times_main %>% 
  count(COMMAND, DATA_INFO, sort=TRUE) 

# times for sub-stages within functions
times_substages <- times %>%
  filter(!(grepl("total", DATA_INFO)) & (grepl("A", STAGE) | grepl("C", STAGE)) & grepl("execution", DATA_INFO))

times_substages %>% 
  count(COMMAND, DATA_INFO, sort=TRUE) 

# summary tables (average)
times_main_summary <- summarySE(times_main, measurevar="DATA", groupvars=c("STAGE", "COMMAND"))
times_substages_summary <- summarySE(times_substages, measurevar="DATA", groupvars=c("STAGE", "COMMAND"))

# table with overall execution time for preprocessing - map and reduce
times %>% 
  filter(STAGE=="PP" & grepl("start", DATA_INFO)) %>%
  min(data)

# 4. time plots----

# 4.1. summary plot: time sub stages 
create_time_plot(times_substages_summary, "Average", "substages_average", debug)

# 4.2. summary plot: times main
create_time_plot(times_main_summary, "Average", "main_average", debug)

# 4.3. alldata plot: time sub stages
create_time_plot(times_substages, "All data", "substages_alldata", debug)

# 4.4 alldata plot: times main
create_time_plot(times_main, "All data", "main_alldata", debug)






~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PROCESS FILE SIZE DATA----
print(paste0("process_execution_times.R","\n",".....processing execution file sizes for ", sizes_file))
sizes <- read.table(sizes_file, header=F, fill= TRUE, sep= "\t", stringsAsFactors = FALSE)
colnames(sizes) <- c("ITERATION", "STAGE","FUNCTION_N", "FILE_TYPE", "FILE_NAME", "SIZE")
# replace function_n column by splitting stage into stage and function_n

# convert SIZE to Mb
sizes$SIZE <- sizes$SIZE/1000000
sizes$FILE_TYPE <- factor(sizes$FILE_TYPE,levels = unique(sizes$FILE_TYPE))
sizes_summary <- summarySE(sizes, measurevar="SIZE", groupvars=c("STAGE", "FILE_TYPE"))


p_sizes <- ggplot(sizes_summary) + 
  geom_point(mapping = aes(x = FILE_TYPE, y = SIZE, colour = STAGE), size=4) +
  geom_errorbar(aes(x = FILE_TYPE, ymin=SIZE-se, ymax=SIZE+se), width=.1, color="#808080") +
  labs(title="Variant Caller - File sizes",
       subtitle = paste0(iterations," program iterations, using ",no_functions," functions each"),
       caption = command,
       x ="", y = "Mb") +
  scale_x_discrete(labels=function(x){sub("\\s", "\n", x)})+
  theme_bw() + 
  theme(panel.border = element_blank(), 
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  theme(legend.position="none") +
  scale_color_manual(values = c("PP" = "#D64324", #red
                                A_id = "#F59513", # orange..
                                "B" = "#A0D624", # green 1
                                "C" = "#228B22")) # green 2
p_sizes
print_png(sizes_file, "_summary", p_sizes)

p_sizes_all <- ggplot(sizes) + 
  geom_point(mapping = aes(x = FILE_TYPE, y = SIZE, colour = STAGE), size=4) +
  labs(title="Variant Caller - File sizes - all",
       subtitle = paste0(iterations," program iterations, using ",no_functions," functions each"),
       caption = command,
       x ="", y = "Mb") +
  scale_x_discrete(labels=function(x){sub("\\s", "\n", x)})+
  theme_bw() + 
  theme(panel.border = element_blank(), 
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  theme(legend.position="none") +
  scale_color_manual(values = c("PP" = "#D64324", #red
                                A_id = "#F59513", # orange..
                                "B" = "#A0D624", # green 1
                                "C" = "#228B22")) # green 2
p_sizes_all
print_png(sizes_file, "_all", p_sizes_all)
# p_sizes_all_jitter<- p_sizes_all + geom_jitter()
# print_png(sizes_file, "_all_jitter", p_sizes_all_jitter)

# THROUGHPUT ANALYSIS----

map_substages <- times_summary %>%
  filter(STAGE=="M") 
map_time <- sum(map_substages$DATA)

fasta_input_size <- sizes_summary$SIZE[sizes_summary$FILE_TYPE=="fasta"]
map_input_files <- sizes_summary 
#map_input_size <- 
  

# "leftover" code ----
#"B" = "#A0D624", # green 1
#print_png(paste0(out_dir,times_file), p)
# # jitter version
# p_times_main_jitter<- p_times_main + geom_jitter()
# print_png(times_file, "_all_jitter", p_times_main_jitter)

# jitter test
# mpg
# p <- ggplot(mpg, aes(cyl, hwy))
# print_png("varcall_out/test", "_1", p)
# p2 <- p + geom_point()
# print_png("varcall_out/test", "_2", p2)
# jitter_test <- p2 + geom_jitter()
# print_png("varcall_out/jitter_test", "_1", jitter_test)
