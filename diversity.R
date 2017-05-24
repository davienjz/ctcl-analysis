#load in the data

df_div <- read.csv("data/simpsons_data.csv", fileEncoding = "UTF-8-BOM")

### extract sample number info

# find CTCL plus rest of string 

pos <- regexpr(
  "CTCL[[:digit:]]{3}[A-Z]{2}[ ]?[[:digit:]]{2}[A-S | U-Z]{0,1}",
  df_div$Data.Set,
  perl = TRUE) #don't include 'T', some strings are continuous with TCRVB e.g. CTCL008PB01TCRVB --> CTCL008PB01T

samplestring <- substr(df_div$Data.Set,pos,pos+attributes(pos)[[1]]-1)

#find the sample number
diversitynumber <- gsub(" ","",samplestring, fixed=TRUE)
diversitynum <- levels(as.factor(diversitynumber))

#load linking file
link <- read.csv("data/diversity_link.csv")

#bind columns
diversity <- cbind(df_div, diversitynumber)

#merge files
diversity_df <- merge(diversity,link,by="diversitynumber")
