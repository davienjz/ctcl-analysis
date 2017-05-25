library(readr)
library(flowCore)
library(flowDensity)
library(flowViz)
library(flowStats)
library(flowClust)
library(parallel)
library(ggplot2)

#Uncomment the transform you want to use
#lgcl <- logicleTransform(transformationId="defaultLogicleTransform", w = 1, t = 1048576,m = 4.5, a = 1)
lgcl <- arcsinhTransform(a=0, b=0.15, c=0)

merge_flow_frames <- function(flow_frames, limit=0)
{
  
  print("Merging...")
  new_flow_frames <- NULL
  for(i in 1:length(flow_frames))
  {
    flow_frame = flow_frames[[i]]
    
    sample_name <- strsplit(keyword(flow_frame)$GUID, "\\.")[[1]][1]

    done = 0
    if(length(new_flow_frames) > 0)
    {
      for(j in 1:length(new_flow_frames))
      {
        new_flow_name <- strsplit(keyword(new_flow_frames[[j]])$GUID, "\\.")[[1]][1]
        
        if(sample_name == new_flow_name)
        {
          exprs(new_flow_frames[[j]]) = rbind(exprs(new_flow_frames[[j]]), exprs(flow_frame))
          done = 1
          print("Merged")
        }
      }
    }
    if(done == 0)
    {
      print("New")
      new_flow_frames <- c(new_flow_frames, flow_frame)
    }
  }
  
  filtered_flow_frames <- NULL
  snames <- c()
  for(i in 1:length(new_flow_frames))
  {
    if(nrow(new_flow_frames[[i]]) > limit)
    {
      filtered_flow_frames <- c(filtered_flow_frames, new_flow_frames[[i]])
      snames <- c(snames, keyword(new_flow_frames[[i]])$GUID)
      print(keyword(new_flow_frames[[i]])$GUID)
    }
  }
  fset <- flowSet(filtered_flow_frames)
  sampleNames(fset) <- snames
  return(fset)
}

load_compensate_transform <- function(folder, cnames, cut_marginals=TRUE, comp_mat = NULL, exclude=c())
{
  set.seed(1)
  #opar <- par()
  #options(warn=-1)
  
  all_LMD_filenames <- list.files(folder, pattern="*.LMD")
  
  f_frames <- NULL

  for(i in 1:length(all_LMD_filenames))
  {
    #load in the FCS
    file = all_LMD_filenames[i]
    print("Loading file...")
    print(file)
    f_set <- try(read.FCS(paste(folder, file, sep=""), dataset=2), silent=TRUE)
    
    if(!inherits(f_set, "try-error"))
    {
      
      fname =  gsub(' ','_',strsplit(file, "\\.")[[1]][1])
      
      if(!fname %in% exclude)
      {
        sample_name = strsplit(file," ")[[1]][1]
        tissue = strsplit(file," ")[[1]][2]
        panel = strsplit(file," ")[[1]][3]
        
        #keyword(f_set)$GUID <- paste(sample_name, tissue, panel, sep="_")
        keyword(f_set)$GUID <- fname
        #sort out the scatter channel names
        
        colnames(f_set)[1:4] <- c("FSC-H", "FSC-A", "SSC-H", "SSC-A")
        parameters(f_set)[[1]][1:4] <- c("FSC-H", "FSC-A", "SSC-H", "SSC-A")
        colnames(f_set)[5:14] <- cnames
        parameters(f_set)[[1]][5:14] <- cnames
        
        #compensate
        print("Compensating...")
        ff_kw <- keyword(f_set)
        cmat <- ff_kw$"$SPILLOVER"
        
        if(!is.null(comp_mat))
        {
          cmat = comp_mat
        }
        
        rownames(cmat) <- colnames(f_set)[5:14]
        colnames(cmat) <- colnames(f_set)[5:14]
        f_set <- compensate(f_set, cmat)
        
        #transform
        print("Transforming...")
        trans <- transformList(cnames, lgcl)
        
        # f_set <- transform(f_set, `FL1-A` = lgcl(`FL1-A`), `FL2-A` = lgcl(`FL2-A`), `FL3-A` = lgcl(`FL3-A`), `FL4-A` = lgcl(`FL4-A`), `FL5-A` = lgcl(`FL5-A`), `FL6-A` = lgcl(`FL6-A`), `FL7-A` = lgcl(`FL7-A`), `FL8-A` = lgcl(`FL8-A`), `FL9-A` = lgcl(`FL9-A`), `FL10-A` = lgcl(`FL10-A`))
        f_set <- transform(f_set, trans)
        
        #set lower boundaries to zero
        parameters(f_set)[[4]] <- rep(0,15)
        
        print("Dealing with marginals...")
        #Remove marginal events from all channels
        boundFilt <- boundaryFilter(filterId = "boundFilt", side="both", x = colnames(f_set)[1:14])
        boundFilt <- filter(f_set, boundFilt)
        
        if(cut_marginals == TRUE)
        {
          f_set <- Subset(f_set, boundFilt) 
        }
        
        f_frames <- c(f_frames, f_set)
      }
    }
    else
    {
      print("Error loading file... moving on")
    }
  }
  samp_names <- c()
  for(i in 1:length(f_frames))
  {
    samp_names <- c(samp_names, keyword(f_frames[[i]])$GUID)
  }
  
  fset <- flowSet(f_frames)
  sampleNames(fset) <- make.unique(samp_names)
  return(fset)
}
  
  
gate_from_file <- function(folder, gating_file, flow_data, names)
{
  options(warn=-1)
  
  #Open gating file
  gating_file_d  <- file(gating_file, open = "r")
  
  #Read file as lines
  gating_lines <- readLines(gating_file_d, warn=FALSE)
  
  flow_list <- list("UNGATED" = flow_data)
  
  returned_data <- data.frame("SAMPLE" = names, UNGATED_CELL_COUNT = rep(NA, length(names)))
  for(i in 1:length(flow_data))
  {
    returned_data$UNGATED_CELL_COUNT[i] = nrow(flow_data[[i]])
  }
  
  #Now do the automated gating
  #for each line (gate)
  for(i in 1:length(gating_lines))
  {
    line = trim(gating_lines[[i]])
    gating_values = strsplit(line, "\t")[[1]]
    
    #get the information out
    population_name <- gating_values[1]
    channel_one <- gating_values[2]
    channel_two <- gating_values[3]
    gate_type <- gating_values[4]
    new_population_name <- gating_values[5]
    
    print(paste("Gating ", new_population_name, " from ", population_name, "..."))
    
    #gate = quadGate
    print("Subsetting...")
    flow_subs = NULL
    
    if(gate_type == "POLY")
    {
      #get coordinates for whichever gate into a matrix
      coordinates <- as.numeric(gating_values[6:length(gating_values)])
      coordinates_mat <- matrix(coordinates,ncol=2,nrow=length(coordinates) / 2)
      colnames(coordinates_mat) <- c(channel_one, channel_two)
      
      gate <- polygonGate(filterId = new_population_name, .gate = coordinates_mat)
      flow_subs <- Subset(flow_list[[population_name]], gate)
    }
    if(gate_type == "QUAD")
    {
      channel_one_lim = as.numeric(gating_values[6])
      channel_two_lim = as.numeric(gating_values[7])
      
      marker_ref = gsub("—", "-",gating_values[8])
      channel_one_ref = substring(marker_ref,1,1)
      channel_two_ref = substring(marker_ref,2,3)
      population = paste(channel_one, channel_one_ref, channel_two, channel_two_ref, sep="")

      gate_values = c(channel_one_lim, channel_two_lim)
      names(gate_values) = c(channel_one, channel_two)
      
      gate <- quadGate(filterId = new_population_name, .gate = gate_values)

      f_res = filter(flow_list[[population_name]], gate)
      flow_split <- split(flow_list[[population_name]], f_res)
      flow_subs <- flow_split[[population]]
    }
    if(gate_type == "1D")
    {
      lower_lim = as.numeric(gating_values[6])
      upper_lim = as.numeric(gating_values[7])
      filter = gsub("—", "-",gating_values[8])
      
      gate_values = list()
      gate_values[[channel_one]] = c(lower_lim, upper_lim)
      gate_values[[channel_two]] = c(-Inf, Inf)
      
      gate <- rectangleGate(filterId="rectGate",.gate = gate_values)
      f_res <- filter(flow_list[[population_name]],gate)
      
      flow_split <- split(flow_list[[population_name]], f_res)

      keep_population = paste("rectGate", filter, sep="")
      flow_subs <- flow_split[[keep_population]]
    }
    
    channel_one_new <- paste("`", channel_one, "`", sep="")
    channel_two_new <- paste("`", channel_two, "`", sep="")
    
    channel_formula = formula(paste(channel_one_new, "~", channel_two_new))

    print("Plotting...")
    png(filename=paste(folder, "/", new_population_name ,".png", sep=""))
    print(xyplot(channel_formula, flow_list[[population_name]], filter = gate))
    dev.off()
    
    flow_list[[new_population_name]] <- flow_subs
    
    new_col_name <- paste(new_population_name, "FROM", population_name, sep="_")
    returned_data[,new_col_name] <- rep("NA", length(names))
    
    rec_col_name <- paste(population_name,channel_one, "EXPRESSION", sep="_")
    if(gate_type == "1D")
    {
      returned_data[,rec_col_name] <- rep("NA", length(names))
    }
    
    for(i in 1:length(flow_subs))
    {
      new_pop_cells = nrow(flow_list[[new_population_name]][[i]])
      old_pop_cells = nrow(flow_list[[population_name]][[i]])
      
      #new_pop_cells = 0
      #old_pop_cells = 0
      if(is.null(new_pop_cells))
      {
        returned_data[,new_col_name][i] = 0
      }
      else if(is.null(old_pop_cells))
      {
        returned_data[,new_col_name][i] = 0
      }
      else
      {
        returned_data[,new_col_name][i] = new_pop_cells / old_pop_cells
      }
      
      if(gate_type == "1D")
      {
        expression_table = exprs(flow_list[[population_name]][[i]])
        recording = mean(expression_table[,channel_one])
        returned_data[,rec_col_name][i] = recording
      }

    }
    
  }
  
  close(gating_file_d)
  
  return(returned_data)
}

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

