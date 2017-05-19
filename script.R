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
galliospatientnumber <- substr(samplestring,pos2,pos2+attributes(pos2)[[1]]-1)

# find sample number
galliosfilenumber <- gsub(" ","",samplestring, fixed=TRUE)

# check results
levels(as.factor(galliospatientnumber))
levels(as.factor(galliosfilenumber))

### add patient number and file number in
df2 <- cbind(galliospatientnumber,galliosfilenumber,df)

df_facs <- read.csv("sample_facs_data.csv")
df_samples <- read.csv("samples.csv")

### merge into facs data and samples sheet
df3 <- merge(df2,df_facs,by="galliosfilenumber")
df4 <- merge(df3,df_samples,by="samplenumber")

### tidy dataframe

#make numeric
numerics <- c("X.Total","X.Gated","X.Med","X.AMean","X.Mode","X.Stdev","X.CV","HP.X.CV","X.Min","X.Max","X.GMean","Y.Med","Y.AMean","Y.Mode","Y.Stdev","Y.CV","HP.Y.CV","Y.Min","Y.Max","Y.GMean")

numerics_select <- names(df4) %in% numerics

df4[,numerics_select] <- lapply(df4[,numerics_select], function(column) as.numeric(column))

#drop unnecessary columns
drop <- c("notes.x","notes.y","location.1","run","run.","clonotypic","vb","pe_fitc_ab",
					"galliosfilenumber","galliospatientnumber","Data.Set")
df5 <- df4[,!(names(df4) %in% drop)]

write.csv(df5, file = "workingfile5.csv")

#drop other factors, by factoring over length of dataframe
df5[] <- lapply(df5, function(column) if(is.factor(column)) factor(column) else column)

#remove square brackets
df5$population <- as.factor(gsub("\\[|\\]","",df5$Input.Gate))
df5$expression <- df5$X.Parameter
df5$gmean <- df5$X.GMean

### find gated data
##find t reg gates
#make space
df5$expression <- factor(df5$expression,c(levels(df5$expression), "treg", "foxp3_pos"))

#find ++ gates (t regs or foxp3+ cd25+)
select <- regexpr("\\+\\+",df5$Gate) != -1
df5[select,"expression"] <- "treg"
df5[select,"pcgate"] <- df5[select,"X.Gated"]
df5[select,c("expression","pcgate")]

#find +- gates (foxp3+ cd25-)
select2 <- regexpr("\\+\\-",df5$Gate) != -1
df5[select2,"expression"] <- "foxp3_pos"
df5[select2,"pcgate"] <- df5[select2,"X.Gated"]
df5[select2,c("expression","pcgate")]

##find panel 5 gates
panel5expression <- c("ifngamma","il4","il10","il17a")
select3a <- df5$Protocol == "panel5"
select3b <- df5$X.Parameter %in% panel5expression
select3c <- df5$Y.Parameter == "Count"
select3d <- df5$Gate != "All"
select3 <- select3a & select3b & select3c & select3d

df5[select3,"pcgate"] <- df5[select3,"X.Gated"]

#check data
table(df5$samplenumber,df5$clone)
table(df5$samplenumber,df5$storage)
table(df5$samplenumber,df5$location)

### subset down dataframe to make a more workable df

df6 <- df5

write.csv(df6, file = "workingfile6.csv")

columns <- c("samplenumber","studynumber","run_date","clone","vbeta","vb_cd7","vb_all","cd2","sampletype","location","storage","population","expression","gmean","pcgate")

columns_select<- names(df6) %in% columns

rows1 <- c("pd1",
					 "pdl1",
					 "pdl2",
					 "tigit",
					 "tim3",
					 "gal9",
					 "lag3",
					 "hladr",
					 "cd25",
					 "foxp3",
					 "fas",
					 "fasl",
					 "mhc1",
					 "ifngamma",
					 "il4",
					 "il10",
					 "il17a",
					 "treg",
					 "foxp3_pos"
					 )

#drop panel3 il10 and tfgb1 to prevent confusion with panel5
rows2_select <- !(df6$Protocol == "panel3" & df6$X.Parameter %in% c("il10","tgfb1"))
mean(rows2_select)

rows1_select <- df6$expression %in% rows1 
mean(rows1_select)

df7 <- df6[rows1_select & rows2_select, columns_select]

#drop other factors, by factoring over length of dataframe
df7[] <- lapply(df7, function(column) if(is.factor(column)) factor(column) else column)

### analyse clonal data that has geometric means
#drop pcgate and take clonal
df8 <- df7[df7$clone == TRUE & !is.na(df7$gmean),!(names(df7) %in% c("pcgate"))]

#melt dataframe
df9 <- melt(df8,c("samplenumber","population","expression"),c("gmean"))

#cast to 3D array
df10 <- acast(df9[,-"Gmean"], samplenumber ~ population ~ expression)

