#package installation
install.packages("RColorBrewer")
install.packages("reshape2")
install.packages("abind")
install.packages("gplots")
install.packages("EnvStats")
install.packages("xtable")
install.packages("corrgram")


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

#convert all empty spaces and N/A to NA AND make numeric

df4[df4 == "" | df4 == "N/A"] <- NA

numerics <- c("X.Total","X.Gated","X.Med","X.AMean","X.Mode", "X.Stdev",
              "X.CV","HP.X.CV","X.Min","X.Max","X.GMean","Y.Med","Y.AMean","Y.Mode",
              "Y.Stdev","Y.CV","HP.Y.CV","Y.Min","Y.Max","Y.GMean") 
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
df5$gated <- df5$X.Gated

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
select3a <- df5$ï..Protocol == "panel5"
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

columns <- c("samplenumber",
             "studynumber",
             "run_date",
             "clone",
             "vbeta",
             "vb_cd7",
             "vb_all",
             "cd2",
             "sampletype",
             "location",
             "storage",
             "population",
             "expression",
             "gmean",
             "pcgate",
             "Gate",
             "Y.Parameter",
             "Protocol",
             "gated")

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
rows2_select <- !(df6$ï..Protocol == "panel3" & df6$X.Parameter %in% c("il10","tgfb1"))
mean(rows2_select)

rows1_select <- df6$expression %in% rows1 
mean(rows1_select)

df7 <- df6[rows1_select & rows2_select, columns_select]

#drop other factors, by factoring over length of dataframe
df7[] <- lapply(df7, function(column) if(is.factor(column)) factor(column) else column)

write.csv(df7, file = "workingfile7.csv")

###analysis of panel 5
#subset panel 5
dfpan5 <- subset(df7, expression == "ifngamma"|expression == "il4"|expression == "il10"|expression == "il17a")
dfpan5

write.csv(dfpan5, file = "workingfilepan5.csv")

##drop gates that are not 'all'
dfpan5a <- dfpan5[!dfpan5$Gate == "All",]

##Calculate iMFI
#convert gated percentage to decimal
dfpan5a$gated <- dfpan5a$gated/100

#calculate iMFI

dfpan5a$iMFI <- dfpan5a$gmean * dfpan5a$gated
write.csv(dfpan5a, file = "workingfileiMFI.csv")

#subset clonal
dfpan5b <- subset(dfpan5a, clone == TRUE & !is.na(dfpan5a$iMFI))
write.csv(dfpan5b, file = "clonalpanel5.csv")

#melt dataframe
dfpan5c <- melt(dfpan5b,c("samplenumber","population","expression"),c("iMFI"))

#cast to 3D array
dfpan5d <- acast(dfpan5c, samplenumber ~ population + expression)

#names
justfornames <- acast(dfpan5c, samplenumber ~ population ~ expression)

naming <- justfornames[, , c("ifngamma",
                               "il4",
                               "il17a",
                               "il10")]

heatmap.2(dfpan5d, 
          trace = "none",
          Colv = NA,
          dendrogram = "row",
          key = TRUE,
          symkey = FALSE,
          key.title = NA,
          margins = c(8.7,7),
          colsep = c(4,8, 12, 16, 20, 24, 28))


### analyse clonal data that has geometric means for panels 1-4
#drop pcgate and take clonal
df8 <- df7[df7$clone == TRUE & !is.na(df7$gmean),!(names(df7) %in% c("pcgate"))]

#drop gates that are not 'all'
df8b <- df8[!df8$Gate == "All",]

#just take panel3 histogram data
exclude <- df8b$Y.Parameter != "Count" & df8b$Protocol == "panel3"
df8b[exclude,]
df8c <- df8b[!exclude,]

#drop panel5
df8d <- df8c[!df8c$Protocol == "panel5",]
df8c[df8c$Protocol == "Panel5",]

#melt dataframe
df9 <- melt(df8d,c("samplenumber","population","expression"),c("gmean"))

#cast to 3D array
df10 <- acast(df9, samplenumber ~ population ~ expression)

#calculate fold change
cd4tilfc <- asinh(df10[,"til_cd4",]) - asinh(df10[,"pb_cd4",])

tumourfc <- asinh(df10[,"tumour",]) - asinh(df10[,"pb_cd4",])

cd8tilfc <- asinh(df10[,"til_cd8",]) - asinh(df10[,"pb_cd8",])

df11 <- abind(df10, "cd4tilfc" = cd4tilfc, "tumourfc" = tumourfc, "cd8tilfc" = cd8tilfc, along = 2)

names(dimnames(df11)) <- c("samplenumber","population","expression")

df12 <- melt(df11)

df13 <- df12[df12$population %in% c("tumourfc","cd4tilfc","cd8tilfc"),]

#re-order by populations
df14 <- acast(df13,samplenumber ~ population + expression)
df14names <- acast(df13, samplenumber ~ population ~ expression)

#clinical groups
clinical <- match(rownames(df14), df7$samplenumber)

#colours
col_breaks <- c(seq(-1,-0.01,length=100),0,seq(0.01,1,length = 100),seq(1.01,2,length=100))
my_palette <- colorRampPalette(c("#3540FF","black","#D42C2C","#FF3535"))(n = length(col_breaks)-1)

col1 <- palette(brewer.pal(8, "Pastel2"))[as.numeric(factor(df7$sampletype[clinical]))+5]

col2 <- palette(brewer.pal(8, "Pastel2"))[as.numeric(factor(rep(dimnames(df14names)[[2]],each = length(dimnames(df14names)[[3]]))))]

col3 <- palette(brewer.pal(8, "Pastel2"))[as.numeric(factor(dimnames(df14names)[[2]]))]

col4 <- palette(brewer.pal(8, "Pastel2"))[seq_along(levels(factor(df7$sampletype[clinical])))+5]

#heatmap
heatmap.2(df14,
					breaks = col_breaks,
					col = my_palette,
					trace = "none",
					Colv = NA,
					RowSideColors = col1,
					ColSideColors = col2
)




