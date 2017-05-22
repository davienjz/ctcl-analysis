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
