#package installation
#install.packages("RColorBrewer")
#install.packages("reshape2")
#install.packages("abind")
#install.packages("gplots")
#install.packages("EnvStats")
#install.packages("xtable")
#install.packages("corrgram")

#libraries
library(RColorBrewer)
library(reshape2)
library(abind)
library(gplots)
library(EnvStats)
library(xtable)
library(corrgram)
library(tsne)
library(quantmod)
library(xts)
library(igraph)

### import functions
source("functions.R")

### import data
panel1 <- read.csv("data/panel1.csv")
panel2 <- read.csv("data/panel2.csv")
panel3 <- read.csv("data/panel3.csv")
panel4 <- read.csv("data/panel4.csv")
panel5 <- read.csv("data/panel5.csv")

df_facs <- read.csv("data/sample_facs_data.csv")
df_samples <- read.csv("data/samples.csv")

### merge data
df <- rbind(panel1, panel2, panel3, panel4, panel5)

writeCsv(df)

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

writeCsv(df5)

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

writeCsv(df6)

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
rows2_select <- !(df6$Protocol == "panel3" & df6$X.Parameter %in% c("il10","tgfb1"))
mean(rows2_select)

rows1_select <- df6$expression %in% rows1 
mean(rows1_select)

df7 <- df6[rows1_select & rows2_select, columns_select]

#drop other factors, by factoring over length of dataframe
df7[] <- lapply(df7, function(column) if(is.factor(column)) factor(column) else column)

writeCsv(df7)

###analysis of panel 5
#subset panel 5
dfpan5 <- subset(df7, expression == "ifngamma"|expression == "il4"|expression == "il10"|expression == "il17a")
dfpan5

writeCsv(dfpan5)

##drop gates that are not 'all'
dfpan5a <- dfpan5[!dfpan5$Gate == "All",]

##Calculate iMFI
#convert gated percentage to decimal
dfpan5a$gated <- dfpan5a$gated/100

#calculate iMFI

dfpan5a$iMFI <- dfpan5a$gmean * dfpan5a$gated

#subset clonal
dfpan5b <- subset(dfpan5a, clone == TRUE & !is.na(dfpan5a$iMFI))

#melt dataframe
dfpan5c <- melt(dfpan5b,c("samplenumber","population","expression"),c("iMFI"))
dfpan5c

###take tils and tumour
tiltum <- subset(dfpan5c, population == "tumour" | population =="til_cd8" | population == "til_cd4")
tiltum

#cast to 3D array
dfpan5d <- acast(tiltum, samplenumber ~ population + expression)
dfpan5d

##organise
dfpan5d[, c("til_cd4_ifngamma",
        "til_cd4_il4",
        "til_cd4_il17a",
        "til_cd4_il10",
        "tumour_ifngamma",
        "tumour_il4",
        "tumour_il17a",
        "tumour_il10",
        "til_cd8_ifngamma",
        "til_cd8_il4",
        "til_cd8_il17a",
        "til_cd8_il10")]

#names
justfornames <- acast(dfpan5c, samplenumber ~ population ~ expression)

naming <- justfornames[, , c("ifngamma",
                               "il4",
                               "il17a",
                               "il10")]

##colours
#set colours
col_breaks <- c(seq(0.01,2,length = 100),
                seq(2.01,4,length=100),
                seq(4.01,6,length=100),
                seq(6.01, 8, length =100)) #sets where different colours start

palette <- colorRampPalette(c("black","darkred", 
                              "#D42C2C",
                              "#FF3535"
                              ))(n = length(col_breaks)-1)

col2 <- palette(brewer.pal(8, "Pastel2"))[as.numeric(factor(rep(dimnames(naming)[[2]])))]
palette(brewer.pal(8, "Pastel2"))[as.numeric(factor(rep(dimnames(naming)[[2]])))]



heatmap.2(dfpan5d, 
          trace = "none",
          Colv = NA,
          dendrogram = "row",
          key = TRUE,
          symkey = FALSE,
          key.title = NA,
          margins = c(8.7,7),
          col = palette,
          breaks = col_breaks,
          colsep = c(4,8),
          labCol = rep(dimnames(naming)[[3]], 
                       length(dimnames(naming)[[3]])))

text(0.3266486, 0.9, labels = "CD4 TIL")
text(0.5275845, 0.9, labels = "TUMOUR" )
text(0.7384352, 0.9, labels = "CD8 TIL")


### analyse clonal data that has geometric means for panels 1-4
#drop pcgate and take clonal
df8 <- df7[df7$clone == TRUE & !is.na(df7$gmean),!(names(df7) %in% c("pcgate"))]

#take gates that are 'All' (crude gmean analysis)
df8b <- df8[df8$Gate == "All",]

#just take panel3 histogram data
exclude <- df8b$Y.Parameter != "Count" & df8b$Protocol == "panel3"
df8b[exclude,]
df8c <- df8b[!exclude,]

# -> fork df8c for striplots

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

#fold change of fold change

df11 <- abind(df10, "cd4tilfc" = cd4tilfc, "tumourfc" = tumourfc, "cd8tilfc" = cd8tilfc, along = 2)

names(dimnames(df11)) <- c("samplenumber","population","expression")

df12 <- melt(df11)

df13 <- df12[df12$population %in% c("tumourfc","cd4tilfc","cd8tilfc"),]

#re-order by populations
df14 <- acast(df13,samplenumber ~ population + expression)
df14name1 <- dcast(df13,samplenumber ~ population + expression)
df14names <- acast(df13, samplenumber ~ population ~ expression)

#clinical groups
clinical <- match(df14name1$samplenumber, df7$samplenumber)

#colours
col_breaks <- c(seq(-1.5,-0.01,length=200),0,seq(0.01,2.5,length = 200),seq(2.51,5,length=200))
my_palette <- colorRampPalette(c("#3540FF","black","#D42C2C","#FF3535"))(n = length(col_breaks)-1)

col1 <- palette(brewer.pal(8, "Pastel2"))[as.numeric(factor(df7$sampletype[clinical]))+4]

col2 <- palette(brewer.pal(8, "Pastel2"))[as.numeric(factor(rep(dimnames(df14names)[[2]],each = length(dimnames(df14names)[[3]]))))]

col3 <- palette(brewer.pal(8, "Pastel2"))[as.numeric(factor(dimnames(df14names)[[2]]))]

col4 <- palette(brewer.pal(8, "Pastel2"))[seq_along(levels(factor(df7$sampletype[clinical])))+4]

#heatmap
heatmap.2(df14,
					breaks = col_breaks,
					col = my_palette,
					trace = "none",
					Colv = NA,
					RowSideColors = col1,
					ColSideColors = col2,
					symm = F,
					symkey = F,
					symbreaks = FALSE,
					scale = "none",
					margins = c(12,14),
					dendrogram = "row",
					labCol = rep(dimnames(df14names)[[3]], length(dimnames(df14names)[[2]])),
					cexCol = 1.5,
					cexRow = 1.5,
					keysize = 1,
					key.title = "Key",
					key.xlab = "Fold Change"
)

#locator()

leg <- legend(x = "topright",
			 title = "Cell population",
			 legend = c("CD4 TILs","Tumour Cells","CD8 TILs"),
			 fill = c(col3),
			 cex = 1.25
			 )
leg

legend(x = leg$rect$left - 1.25 * leg$rect$w,
			 y = leg$rect$top,
			 title = "Stage of disease",
			 legend = levels(factor(df7$sampletype[clinical])),
			 fill = col4,
			 cex = 1.25
			 )

### heatmap order by patient/pop
#re-order by populations
dh14 <- acast(df13,samplenumber + population ~ expression)
dh14name1 <- dcast(df13,samplenumber + population ~ expression)
dh14names <- acast(df13, samplenumber ~ population ~ expression)

#clinical groups
clinical <- match(dh14name1$samplenumber, df7$samplenumber)

#colours
col_breaks <- c(seq(-1.5,-0.01,length=200),0,seq(0.01,2.5,length = 200),seq(2.51,5,length=200))
my_palette <- colorRampPalette(c("#3540FF","black","#D42C2C","#FF3535"))(n = length(col_breaks)-1)

col1 <- palette(brewer.pal(3, "Set1"))[as.numeric(factor(dh14name1$population))]

col2 <- palette(brewer.pal(8, "Pastel2"))[as.numeric(factor(rep(dimnames(dh14names)[[2]],each = length(dimnames(dh14names)[[3]]))))]

col3 <- palette(brewer.pal(8, "Pastel2"))[as.numeric(factor(dimnames(dh14names)[[2]]))]

col4 <- palette(brewer.pal(8, "Pastel2"))[seq_along(levels(factor(df7$sampletype[clinical])))+4]

#heatmap
heatmap.2(dh14,
					breaks = col_breaks,
					col = my_palette,
					trace = "none",
					Colv = NA,
					RowSideColors = col1,
					#ColSideColors = col2,
					symm = F,
					symkey = F,
					symbreaks = FALSE,
					scale = "none",
					margins = c(12,14),
					dendrogram = "row",
					labCol = rep(dimnames(dh14names)[[3]], length(dimnames(dh14names)[[2]])),
					cexCol = 1.5,
					cexRow = 0.8,
					keysize = 1,
					key.title = "Key",
					key.xlab = "Fold Change"
)

#locator()

leg <- legend(x = "topright",
			 title = "Cell population",
			 legend = c("CD4 TILs","Tumour Cells","CD8 TILs"),
			 fill = c(col3),
			 cex = 1.25
			 )

legend(x = leg$rect$left - 1.25 * leg$rect$w,
			 y = leg$rect$top,
			 title = "Stage of disease",
			 legend = levels(factor(df7$sampletype[clinical])),
			 fill = col4,
			 cex = 1.25
			 )

### strip plots

df8c$population

stripPlot(df8c,c("tumour","pb_cd4"))
stripPlot(df8c,c("til_cd4","pb_cd4"))
stripPlot(df8c,c("til_cd8","pb_cd8"))

### PCA of populations+patients

dpca1 <- t(dh14[complete.cases(dh14),])
rowMeans(dpca1)
s <- svd(dpca1-rowMeans(dpca1))

pc1 <- s$d[1]*s$v[,1]
pc2 <- s$d[2]*s$v[,2]

dh14name1$samplenumber

group <- factor(dh14name1$population)
color1 <- brewer.pal(3,"Set1")

png(gen("pca"))

plot(pc1,
		 pc2,
		 xlab = "PC1",
		 ylab = "PC2",
		 col = c("green","red","blue"),
		 type = "n"
		 )

text(pc1,
		 pc2,
		 labels = dh14name1$samplenumber,
	 col = c("green","red","blue"),
	 )

dev.off()


### plots of PCA variable colour

pca_df <- df13
pca_formula <- samplenumber + population ~ expression
pca_colour <- "pdl2"

plotPCA(pca_df, pca_formula, pca_colour)

for (pca_colour in levels(pca_df$expression)){
	plotPCA(pca_df, pca_formula, pca_colour)
}

### tSNE of populations+patients

#uncomment to run tSNE
#output <- tsne(t(dpca1), perplexity = 15, max_iter = 3000)

#str(output)

#png(gen("tsne"))

#plot(output[,1],
#		 output[,2],
#		 type = "n",
#		 col = c("green","red","blue"),
#		 pch = 19
#		 )

#text(output[,1],
#		 output[,2],
#		 labels = dh14name1$samplenumber,
#		 col = c("green","red","blue"),
#		 )

#dev.off()

### correlation analysis


dc14 <- acast(df13, samplenumber ~ population + expression)
dc15 <- dc14[complete.cases(dc14),]
dc16 <- cov(dc15)

max(dc16)
min(dc16)
cov_breaks <- c(seq(-0.3,-0.01,length = 100),0,seq(0.01,0.8,length = 100))
cov_palette <- colorRampPalette(c("#3540FF","black","#D42C2C"))(n = length(cov_breaks)-1)

heatmap.2(dc16,
					margins = c(10,10),
					trace = "none",
					col = cov_palette,
					breaks = cov_breaks
					)

#get smaller labels
df13
dn14 <- acast(df13,samplenumber ~ expression ~ population)
dn15 <- dn14[complete.cases(dc14),,]
dimnames(dn15)[[3]] <- c("cd4","tmr","cd8")
dn16 <- melt(dn15)
colnames(dn16) <- c("samplenumber","expression","population","value")
dn17 <- acast(dn16, samplenumber ~ expression + population)
dim(dn15)

#define colours and breaks
net_breaks <- c(seq(-1,-0.01,length=100),0,seq(0.01,1,length = 100)) 
net_palette <- colorRampPalette(c("red","white","#006104"))(n = length(net_breaks)-1) 

#define threshold for network
thres_pos <- 0.7
thres_neg <- -0.2

#get correlation matrix, taken upper triangle and those above threshold
dn18<- cor(dn17, method="spearman")
#dn18["pd1_tmr","pdl1_tmr"] <- -0.9 # test correct weights
dn18[ lower.tri(dn18, diag=TRUE) ] <- 0
net_deselect <- thres_neg < dn18 & dn18 < thres_pos
dn18[net_deselect] <- 0

#invert matrix to arrange values in same order as E(graph)
dn19 <- t(dn18)
net_weights <- dn19[dn19!=0]
writeCsv(dn19)

#create graph
graph <- graph.adjacency(dn18, weighted=TRUE, mode="upper", diag = FALSE)

#assign weights graph
E(graph)$weight <- net_weights
E(graph)

#find colours for weights
colours <- cut(E(graph)$weight, breaks = net_breaks, labels = net_palette)
colours <- as.character(colours)

#find absolute weights
E(graph)$weight

#assign colors to graph
E(graph)$color <- colours
V(graph)$color <- rep(c("#B0DCFF","#FFB0B1","#D2FFB0"),dim(dn15)[2])

#assign width of lines to graph
E(graph)$width <- 2

#assign graph layout
graph$layout <- layout.fruchterman.reingold

#delete vertices which aren't connected
graph <- delete.vertices(graph,which(degree(graph) < 1))

plot(graph)
