# Load libraries 
library(data.table)
library(tidyverse)
library(plyr)
library(dplyr)
library(ggpubr)
library(PMCMR)
library(stringr)
library(readxl)
library(readr)
library(pracma)
library(bayestestR)
library(zoo)

# Should define here how large the rollmean and amount of iterations. Also the 
# threshold for peak detection. 


# Clear working space
rm(list = ls())
setwd('~/Python/Raw data for article/Ulgu cells')

savepath = paste0(getwd(), '/graphs')
dir.create(savepath)

# import data and organize data
d = read.delim("Summary final.csv", header = TRUE, sep = (','))
d$name = sub('tiff', 'tif', d$Name)
colnames(d) = c('x', 'Name', 'length', 'Angle', 'Organization', 'alignment')

# Retrieve framenr from imagename 
d$frame = as.numeric(substr(d$Name, nchar(d$Name)-7, nchar(d$Name)-4))
d$Name = substr(d$Name, 1, nchar(d$Name)-8)
cells <-unique(d$Name)

short = data.frame()
for(i in cells){
  temp = subset(d, d$Name == i)
  
  setorder(temp,frame)
  rownames(temp) <- 1:nrow(temp)
  temp$length2 = temp$length*-1
  temp$length2 = rollmean(temp$length2, 3, na.pad=TRUE)
  temp$length2 = rollmean(temp$length2, 3, na.pad=TRUE)

  temp$diff = c(diff(temp$length2), NA)
  
  peaks = data.frame(findpeaks(temp$length2, threshold = mean(temp$length2[abs(temp$diff) <= 0.001], na.rm = T)+0.02,  minpeakheight = mean(temp$length2[abs(temp$diff) <= 0.001], na.rm = T)+0.02, , minpeakdistance =  25))
  
  pkx = c()
  xx = c()
  
  mean=mean(temp$length2[abs(temp$diff) <= 0.001 & temp$length2 <= mean(temp$length2[abs(temp$diff) <= 0.001], na.rm = T)], na.rm = T)
  
  for (j in 1:nrow(peaks)){
    peaks$delta[j] = abs((temp$length2[peaks[ j , 2]]-mean)/mean)*100
  }
  
  peaks = subset(peaks, peaks$delta >= 1.5)
  
  for (j in 1:nrow(peaks)){
    pkx[j] = abs((temp$length2[peaks[ j , 2]]-mean)/mean)*100
  }
  
  
  meanorg=mean(temp$Organization[abs(temp$diff) <= 0.001 & temp$length2 <= mean(temp$length2[abs(temp$diff) <= 0.001], na.rm = T)], na.rm = T)
  
  for (j in 1:nrow(peaks)){
    xx[j] = temp$frame[peaks[j, 4]]-temp$frame[peaks[j, 2]]
  }
  
  sum = data.frame( 'Name' = temp[1,2], 
                    'interval' = mean(diff(sort(temp$frame[peaks[,2]]))),
                    'short' = max(pkx),
                    'var' = sd(pkx),
                    'baseline_length' = mean,
                    'contractions' = nrow(peaks),
                    'baseline_org' = meanorg,
                    'duration' = mean(xx)
  )
  short = rbind(short, sum)  
  
  wat = data.frame(wat2 = temp$frame[peaks[ 1:nrow(peaks) , 2]], wat3 = temp$frame[peaks[ 1:nrow(peaks) , 3]], wat4 = temp$frame[peaks[ 1:nrow(peaks) , 4]])

  #plot the sarcomere organization score over different conditions

 ggplot(temp, aes(x = frame, y = length2*-1))+
    geom_point()+
    geom_vline(data = wat, aes(xintercept = wat2))+
    geom_vline(data = wat, aes(xintercept = wat3), color = 'red')+
    geom_vline(data = wat, aes(xintercept = wat4), color = 'green')+
    geom_line(aes(y = abs(mean(length2[abs(diff) <= 0.001], na.rm = T)+0.03)))+
    geom_line(aes(y = abs(mean)))+
    geom_line(aes(y = abs(mean(length2[peaks[ 1:nrow(peaks) , 2]]))))+
    geom_line(aes(y = abs(meanorg)))+
    ggtitle(temp$Name[1])+
    ylab("Sarcomere length [um]")+
    xlab('Time [s]')+
    theme_bw(base_size = 20)+
    scale_y_continuous(limits = c(1.2, 2.1))
    ggsave(paste0(savepath, sub('\\\\', '/', temp$Name[1]), ".jpeg"),  device = "jpeg", plot = last_plot())
}

write.table(short, (paste(savepath, 'data_summary.txt', sep = "/")), sep = "\t", quote = FALSE)


