#libraries

library(RColorBrewer)
library(reshape2)
library(abind)
library(gplots)
library(EnvStats)
library(xtable)
library(corrgram)

#import data

panel1 <- read.csv("panel1.csv")
panel2 <- read.csv("panel2.csv")
panel3 <- read.csv("panel3.csv")
panel4 <- read.csv("panel4.csv")
panel5 <- read.csv("panel5.csv")

#merge data

df <- rbind(panel1, panel2, panel3, panel4, panel5)

write.csv(df, file = "workingfile.csv") #the file

#filename generation

count <<- 0
gen <- function(x){
	count <<- count + 1
	countf <- formatC(count, width = 2, format = "d", flag = "0")
	paste0(countf,x,".png")
}

#extract sample number info

pos <- regexpr("CTCL[[:digit:]]{3}[A-Z]{2}[ ]?[[:digit:]]{2}[A-Z]?",df$Data.Set)
samplestring <- substr(df$Data.Set,pos,pos+attributes(pos)[[1]]-1)
pos2 <- regexpr("CTCL[[:digit:]]{3}[A-Z]{2}",samplestring)
patientnumber <- substr(samplestring,pos2,pos2+attributes(pos2)[[1]]-1)
samplenumber <- gsub(" ","",samplestring, fixed=TRUE)

levels(as.factor(patientnumber))
levels(as.factor(samplenumber))


