### filename generation
count <<- 0
gen <- function(x){
	count <<- count + 1
	countf <- formatC(count, width = 2, format = "d", flag = "0")
	paste0("./figures/",countf,x,".png")
}

plotPCA <- function(pca_df, pca_formula, pca_colour){
	#takes melted dataframe, formula and dimension to colour by

	pca_df1 <- acast(pca_df,pca_formula)
	pca_df_name1 <- dcast(pca_df,pca_formula)
	pca_df2 <- pca_df1[complete.cases(pca_df1),]
	pca_df_name2 <- pca_df_name1[complete.cases(pca_df_name1),]

	pca_col_breaks <- c(seq(-1.5,-0.01,length=200),0,seq(0.01,2.5,length = 200),seq(2.51,5,length=200))

	pca_my_palette <- colorRampPalette(c("#3540FF","black","#D42C2C","#FF3535"))(n = length(col_breaks)-1)

	groups <- cut(pca_df_name2[,pca_colour], breaks = pca_col_breaks)
	colours <- pca_my_palette[groups]

	pca_s <- svd(t(pca_df2))

	pca1 <- pca_s$d[1]*pca_s$v[,1]
	pca2 <- pca_s$d[2]*pca_s$v[,2]

	png(gen(paste0("pca",pca_colour)))

	plot(pca1,
			 pca2,
			 main = paste0("PCA analysis - ",pca_colour),
			 type = "n"
			 )

	text(pca1,
			 pca2,
			 labels = pca_df_name2$population,
			 col = colours
			 )

	dev.off()
}


stripPlot <- function(df_a, compare){
	#takes a melted datafram and a vector of two parameters to compare
	df_b <- df_a[df_a$population %in% compare,]
	df_c <- df_b
	df_c$gmean <- asinh(df_b$gmean)
	names <- levels(factor(df_c$expression))
	levels(factor(df_c$expression))
	names2 <- rep(names ,each = 2)
	
	stripChart(gmean ~ population + expression,
					 df_c,
					 col = c("red","blue"),
					 vertical = T,
					 method = "jitter",
					 p.value = TRUE,
					 cex = 0.63,
					 group.names = names2,
					 ylab = "MFI"
					 )

	legend(
				 "topright",
				 inset = 0.05,
				 c(compare[1],compare[2]),
				 fill = c("red","blue")
				 )
	df_d <- acast(df_c, samplenumber ~ population ~ expression, value.var = "gmean")
	
	stats <- apply(df_d, 3, function(x){
							 statistic <- t.test(x[,2],x[,1],paired = TRUE)[[5]]
							 pvalue <- t.test(x[,2],x[,1],paired = TRUE)[[3]]
							 return(c(statistic,pvalue))
					 })

	return(stats)
}

writeCsv <- function(df_csv){
	fn <- deparse(substitute(df_csv))
	write.csv(df_csv, file = paste0("output/",fn,".csv"))
}
