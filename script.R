#Libraries

library(RColor
#Import data

panel1 <- read.csv("panel1.csv")
panel2 <- read.csv("panel2.csv")
panel3 <- read.csv("panel3.csv")
panel4 <- read.csv("panel4.csv")
panel5 <- read.csv("panel5.csv")

#Merge data

df <- rbind(panel1, panel2, panel3, panel4, panel5)

write.csv(df, file = "workingfile.csv")


