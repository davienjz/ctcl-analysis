library(openCyto)
library(data.table)

#Set-up
base_dir = ".."
lmd_dir_p1 = "panel_1_lmd/"
lmd_dir_p2 = "panel_2_lmd/"
lmd_dir_p3 = "panel_3_lmd/"
lmd_dir_p4 = "panel_4_lmd/"
gates_dir = "gating/gates/"

#1) Name channels
p1_channel_names <-  c("FL1", "FL2", "dump", "tigit", "pd1", "pdl2", "cd4", "cd3", "pdl1", "cd8")
p2_channel_names <-  c("FL1", "FL2", "dump", "gal9", "hladr", "lag3", "cd4", "cd3", "tim3", "cd8")
p3_channel_names <-  c("FL1", "FL2", "dump", "il10", "cd25", "tgfb1", "cd4", "cd3", "foxp3", "cd8")
p4_channel_names <-  c("FL1", "FL2", "dump", "blank", "mhci", "fasl", "cd4", "cd3", "fas", "cd8")

#2) Load in flow datasets using the lmd directories defined above
f_frames_p1 <- load_compensate_transform(paste(base_dir, lmd_dir_p1, sep="/"), p1_channel_names)
f_frames_p2 <- load_compensate_transform(paste(base_dir, lmd_dir_p2, sep="/"), p2_channel_names)
f_frames_p3 <- load_compensate_transform(paste(base_dir, lmd_dir_p3, sep="/"), p3_channel_names)
f_frames_p4 <- load_compensate_transform(paste(base_dir, lmd_dir_p4, sep="/"), p4_channel_names)

#3) Merge flow files. Keep files with >= 400 cells
merged_fset_1 <- merge_flow_frames(f_frames_p1, 400)
merged_fset_2 <- merge_flow_frames(f_frames_p2, 400)
merged_fset_3 <- merge_flow_frames(f_frames_p3, 400)
merged_fset_4 <- merge_flow_frames(f_frames_p4, 400)

#4) Load in gating strategies from csv
all_gating <- gatingTemplate(paste(base_dir, gates_dir ,"duncan_gating.csv", sep="/"), autostart = 1L)
panel_3_extra_gating <- gatingTemplate(paste(base_dir, gates_dir, "panel_3_extra_gating.csv",sep="/"), autostart = 1L)

#5) Convert flow frames to gating sets
p1_frames <- GatingSet(merged_fset_1)
p2_frames <- GatingSet(merged_fset_2)
p3_frames <- GatingSet(merged_fset_3)
p4_frames <- GatingSet(merged_fset_4)

#6) Load in metadata
facs_data <- data.frame(read_csv(paste(base_dir, "data", "sample_facs_data.csv", sep="/")))

#7) Set panel 1 phenotype data
pData(p1_frames)$name <- row.names(pData(p1_frames))
pData(p1_frames)$tissue <- "none"
pData(p1_frames)$patient_id <- "none"
pData(p1_frames)$vb_cd7 <- "vbeta"
pData(p1_frames)$clone_quadrant <- "tr"
pData(p1_frames)$type <- "none"
for(i in 1:nrow(pData(p1_frames)))
{
  patient_id = strsplit(pData(p1_frames)[i,"name"],"_")[[1]][1]
  pData(p1_frames)[i,"patient_id"] = patient_id
  pData(p1_frames)[i,"tissue"] = strsplit(pData(p1_frames)[i,"name"],"_")[[1]][2]
  
  facs_point = which(facs_data$patientnumber == patient_id)
  if(length(facs_point) > 0)
  {
    pData(p1_frames)[i,"vb_cd7"] = facs_data[facs_point[1],"vb_cd7"]
    pData(p1_frames)[i,"clone_quadrant"] = facs_data[facs_point[1],"clone_quadrant"]
  }
}

#8) Perform panel 1 gating
gating(all_gating, p1_frames)

#9) Set panel 2 phenotypic data
pData(p2_frames)$name <- row.names(pData(p2_frames))
pData(p2_frames)$tissue <- "none"
pData(p2_frames)$patient_id <- "none"
pData(p2_frames)$vb_cd7 <- "vbeta"
pData(p2_frames)$type <- "none"
for(i in 1:nrow(pData(p2_frames)))
{
  patient_id = strsplit(pData(p2_frames)[i,"name"],"_")[[1]][1]
  pData(p2_frames)[i,"patient_id"] = patient_id
  pData(p2_frames)[i,"tissue"] = strsplit(pData(p2_frames)[i,"name"],"_")[[1]][2]
  
  facs_point = which(facs_data$patientnumber == patient_id)
  if(length(facs_point) > 0)
  {
    pData(p2_frames)[i,"vb_cd7"] = facs_data[facs_point,"vb_cd7"]
  }
}

#10) Perform panel 2 gating
gating(all_gating, p2_frames)

#11) Set panel 3 phenotypic data
pData(p3_frames)$name <- row.names(pData(p3_frames))
pData(p3_frames)$tissue <- "none"
pData(p3_frames)$patient_id <- "none"
pData(p3_frames)$vb_cd7 <- "vbeta"
pData(p3_frames)$type <- "none"
for(i in 1:nrow(pData(p3_frames)))
{
  patient_id = strsplit(pData(p3_frames)[i,"name"],"_")[[1]][1]
  pData(p3_frames)[i,"patient_id"] = patient_id
  pData(p3_frames)[i,"tissue"] = strsplit(pData(p3_frames)[i,"name"],"_")[[1]][2]
  
  facs_point = which(facs_data$patientnumber == patient_id)
  if(length(facs_point) > 0)
  {
    pData(p3_frames)[i,"vb_cd7"] = facs_data[facs_point,"vb_cd7"]
  }
}

#12) Perform panel 3 gating (note different gates)
gating(panel_3_extra_gating, p3_frames)

#13) Set panel 4 phenotypic data
pData(p4_frames)$name <- row.names(pData(p4_frames))
pData(p4_frames)$tissue <- "none"
pData(p4_frames)$patient_id <- "none"
pData(p4_frames)$vb_cd7 <- "vbeta"
pData(p4_frames)$type <- "none"
for(i in 1:nrow(pData(p4_frames)))
{
  patient_id = strsplit(pData(p4_frames)[i,"name"],"_")[[1]][1]
  pData(p4_frames)[i,"patient_id"] = patient_id
  pData(p4_frames)[i,"tissue"] = strsplit(pData(p4_frames)[i,"name"],"_")[[1]][2]
  
  facs_point = which(facs_data$patientnumber == patient_id)
  if(length(facs_point) > 0)
  {
    pData(p4_frames)[i,"vb_cd7"] = facs_data[facs_point,"vb_cd7"]
  }
}

#14) Perform panel 4 gating
gating(all_gating, p4_frames)

#15) Plot gating strategy
plot(p1_gating)

#16) Do a couple of individual plots
png(filename=paste(base_dir, "test1.png", sep="/"), width=1400,height=1000)
print(plotGate(subs_frames, "cd3_pos", gpar=list(nrow=2)))
dev.off()

png(filename=paste(base_dir, "test2.png", sep="/"), width=1400,height=1000)
print(plotGate(p1_frames, c("cd3_pos/FL1+FL2+", "cd3_pos/FL1-FL2-"), gpar=list(nrow=2)))
dev.off()

#17) Retrieve information, second argument is channels you're interested in quantifying
p1_data <- measure_duncan_flow(p1_frames, c("FSC-A","SSC-A","tigit","pd1","pdl2","pdl1"))
p2_data <- measure_duncan_flow(p2_frames, c("FSC-A","SSC-A","gal9","hladr","lag3","tim3"))
p3_data <- measure_duncan_flow(p3_frames, c("FSC-A","SSC-A","il10","cd25","tgfb1","foxp3"))
p4_data <- measure_duncan_flow(p4_frames, c("FSC-A","SSC-A","mhci","fasl","fas"))