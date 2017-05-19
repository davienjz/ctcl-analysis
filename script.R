#libraries

library(RColorBrewer)
library(reshape2)
library(abind)
library(gplots)
library(EnvStats)
library(xtable)
library(corrgram)

### import data

panel1 <- read.csv("panel1.csv")
panel2 <- read.csv("panel2.csv")
panel3 <- read.csv("panel3.csv")
panel4 <- read.csv("panel4.csv")
panel5 <- read.csv("panel5.csv")

### merge data

df <- rbind(panel1, panel2, panel3, panel4, panel5)

write.csv(df, file = "workingfile.csv")

### filename generation

count <<- 0
gen <- function(x){
	count <<- count + 1
	countf <- formatC(count, width = 2, format = "d", flag = "0")
	paste0(countf,x,".png")
}

### extract sample number info

# find CTCL or CTL (typo) plus rest of string 
pos <- regexpr(
							 "CT[C]?L[[:digit:]]{3}[A-Z]{2}[ ]?[[:digit:]]{2}[A-Z]{0,2}",
							 df$Data.Set,
							 perl = TRUE)

samplestring_unfixed <- substr(df$Data.Set,pos,pos+attributes(pos)[[1]]-1)

# remove CTCL or CTL
pos_fix <- regexpr("(?<=L)[[:digit:]]{3}[A-Z]{2}[ ]?[[:digit:]]{2}[A-Z]{0,2}",
									 samplestring_unfixed,
									 perl = TRUE)

samplestring <- substr(samplestring_unfixed,pos_fix,pos_fix+attributes(pos_fix)[[1]]-1)

# add CTCL back
samplestring <- paste0("CTCL",samplestring)

pos2 <- regexpr(
								"CT[C]?L[[:digit:]]{3}[A-Z]{0,2}",
								samplestring,
								perl = TRUE)

# find patient number
patientnumber <- substr(samplestring,pos2,pos2+attributes(pos2)[[1]]-1)

# find sample number
samplenumber <- gsub(" ","",samplestring, fixed=TRUE)

# check results
levels(as.factor(patientnumber))
levels(as.factor(samplenumber))

df <- cbind(patientnumber,samplenumber,df)

df[df$patientnumber=="CTCL013DB",]

df_facs <- read.csv("sample_facs_data.csv")
df_samples <- read.csv("samples.csv")

### merge into facs data and samples sheet
df3 <- merge(df2,df_facs,by="galliosfilenumber")
df4 <- merge(df3,df_samples,by="samplenumber")

#drop unnecessary columns

drop <- c("notes.x","notes.y","location.1","run","run.","clonotypic","vb","pe_fitc_ab",
					"galliosfilenumber","galliospatientnumber","Data.Set")
df5 <- df4[,!(names(df4) %in% drop)]

write.csv(df5, file = "workingfile5.csv")

### 

