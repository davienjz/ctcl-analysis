#load in the data
df_div <- read.csv("data/simpsons_data.csv")

### extract sample number info

# find CTCL plus rest of string 
pos <- regexpr(
  "CTCL[[:digit:]]{3}[A-Z]{2}[ ]?[[:digit:]]{2}[A-Z]{0,2}",
  df_div$Data.Set,
  perl = TRUE)

samplestring <- substr(df_div$Data.Set,pos,pos+attributes(pos)[[1]]-1)

#find the sample number
diversitynumber <- gsub(" ","",samplestring, fixed=TRUE)
levels(as.factor(diversitynumber))
