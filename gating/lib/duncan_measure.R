measure_duncan_flow <- function(frames, channels_to_measure=c())
{
  panel_measurements <- data.frame(sample_names = sampleNames(frames))
  
  sample_frame <- getData(p1_frames[[1]],"root")
  channel_names <- colnames(exprs(sample_frame))
  p_data <- pData(frames)
  
  for(i in 1:length(frames))
  {
    cur_gates = frames[[i]]
    
    #Get FL1/FL2 populations
    pp_lymph = getData(cur_gates, "pp_lymph")
    np_lymph = getData(cur_gates, "np_lymph")
    pn_lymph = getData(cur_gates, "pn_lymph")
    nn_lymph = getData(cur_gates, "nn_lymph")
    
    #Get cd4/cd8 splits
    pp_cd8 <- getData(cur_gates, "pp_lymph/cd4-cd8+")
    np_cd8 <- getData(cur_gates, "np_lymph/cd4-cd8+")
    pn_cd8 <- getData(cur_gates, "pn_lymph/cd4-cd8+")
    nn_cd8 <- getData(cur_gates, "nn_lymph/cd4-cd8+")
    
    pp_cd4 <- getData(cur_gates, "pp_lymph/cd4+cd8-")
    np_cd4 <- getData(cur_gates, "np_lymph/cd4+cd8-")
    pn_cd4 <- getData(cur_gates, "pn_lymph/cd4+cd8-")
    nn_cd4 <- getData(cur_gates, "nn_lymph/cd4+cd8-")
    
    tum_exprs <- NULL
    til_exprs <- NULL
    if(p_data[i, "vb_cd7"] == "no_clone")
    {
      til_exprs <- rbind(exprs(np_lymph), exprs(pn_lymph), exprs(pp_lymph), exprs(nn_lymph))
      prop_til_cd4 <- (nrow(pn_cd4) + nrow(pp_cd4) + nrow(nn_cd4) + nrow(np_cd4)) / nrow(til_exprs)
      prop_til_cd8 <- (nrow(pn_cd8) + nrow(pp_cd8) + nrow(nn_cd8) + nrow(np_cd8)) / nrow(til_exprs)
      
      panel_measurements[i,"prop_tum_cd4"] = NULL
      panel_measurements[i,"prop_til_cd4"] = prop_til_cd4
      panel_measurements[i,"prop_tum_cd8"] = NULL
      panel_measurements[i,"prop_til_cd8"] = prop_til_cd8  
    }
    if(p_data[i,"vb_cd7"] == "vbeta")
    {
      prop_tum_cd4 <- NULL
      prop_til_cd4 <- NULL
      prop_tum_cd8 <- NULL
      prop_til_cd8 <- NULL
      
      if(p_data[i,"clone_quadrant"] == 1)
      {
        tum_exprs <- exprs(np_lymph)
        til_exprs <- rbind(exprs(pn_lymph), exprs(pp_lymph), exprs(nn_lymph))
        
        prop_tum_cd4 <- nrow(np_cd4) / nrow(tum_exprs)
        prop_til_cd4 <- (nrow(pn_cd4) + nrow(pp_cd4) + nrow(nn_cd4)) / nrow(til_exprs)
        prop_tum_cd8 <- nrow(np_cd8) / nrow(tum_exprs)
        prop_til_cd8 <- (nrow(pn_cd8) + nrow(pp_cd8) + nrow(nn_cd8)) / nrow(til_exprs)
      }
      else if(p_data[i,"clone_quadrant"] == 2)
      {
        tum_exprs <- exprs(nn_lymph)
        til_exprs <- rbind(exprs(np_lymph), exprs(pp_lymph), exprs(pn_lymph))
        
        prop_tum_cd4 <- nrow(nn_cd4) / nrow(tum_exprs)
        prop_til_cd4 <- (nrow(np_cd4) + nrow(pp_cd4) + nrow(pn_cd4)) / nrow(til_exprs)
        prop_tum_cd8 <- nrow(nn_cd8) / nrow(tum_exprs)
        prop_til_cd8 <- (nrow(np_cd8) + nrow(pp_cd8) + nrow(pn_cd8)) / nrow(til_exprs)
      }
      else if(p_data[i,"clone_quadrant"] == 3)
      {
        tum_exprs <- exprs(pn_lymph)
        til_exprs <- rbind(exprs(np_lymph), exprs(pp_lymph), exprs(nn_lymph))
        
        prop_tum_cd4 <- nrow(pn_cd4) / nrow(tum_exprs)
        prop_til_cd4 <- (nrow(np_cd4) + nrow(pp_cd4) + nrow(nn_cd4)) / nrow(til_exprs)
        prop_tum_cd8 <- nrow(pn_cd8) / nrow(tum_exprs)
        prop_til_cd8 <- (nrow(np_cd8) + nrow(pp_cd8) + nrow(nn_cd8)) / nrow(til_exprs)
      }
      else
      {
        tum_exprs <- exprs(pp_lymph)
        til_exprs <- rbind(exprs(np_lymph), exprs(pn_lymph), exprs(nn_lymph))
        
        prop_tum_cd4 <- nrow(pp_cd4) / nrow(tum_exprs)
        prop_til_cd4 <- (nrow(np_cd4) + nrow(pn_cd4) + nrow(nn_cd4)) / nrow(til_exprs)
        prop_tum_cd8 <- nrow(pp_cd8) / nrow(tum_exprs)
        prop_til_cd8 <- (nrow(np_cd8) + nrow(pn_cd8) + nrow(nn_cd8)) / nrow(til_exprs)
      }
    
      panel_measurements[i,"prop_tum_cd4"] = prop_tum_cd4
      panel_measurements[i,"prop_til_cd4"] = prop_til_cd4
      panel_measurements[i,"prop_tum_cd8"] = prop_tum_cd8
      panel_measurements[i,"prop_til_cd8"] = prop_til_cd8  
    }
    else if(p_data[i, "vb_cd7"] == "cd7")
    {
      tum_exprs <- rbind(exprs(nn_lymph), exprs(pn_lymph))
      til_exprs <- rbind(exprs(np_lymph), exprs(pp_lymph))
      
      prop_tum_cd4 <- (nrow(nn_cd4) + nrow(pn_cd4)) / nrow(tum_exprs)
      prop_til_cd4 <- (nrow(np_cd4) + nrow(pp_cd4)) / nrow(til_exprs)
      prop_tum_cd8 <- (nrow(nn_cd8) + nrow(pn_cd8)) / nrow(tum_exprs)
      prop_til_cd8 <- (nrow(np_cd8) + nrow(pp_cd8)) / nrow(til_exprs)
      
      panel_measurements[i,"prop_tum_cd4"] = prop_tum_cd4
      panel_measurements[i,"prop_til_cd4"] = prop_til_cd4
      panel_measurements[i,"prop_tum_cd8"] = prop_tum_cd8
      panel_measurements[i,"prop_til_cd8"] = prop_til_cd8
    }
    
    #get proportion of cd3+ lymphocytes that are tumour
    cd3_lymph = getData(cur_gates, "cd3_lymph")
    prop_lymph_tumour = nrow(tum_exprs) / nrow(cd3_lymph)
    panel_measurements[i,"prop_tum_lymph"] = prop_lymph_tumour

    for(j in 1:length(channels_to_measure))
    {
      channel = channels_to_measure[j]
      
      tum_average = NULL
      til_average = NULL
      if(!is.null(tum_exprs))
      {
        tum_average <- mean(tum_exprs[,channel])
      }
      if(!is.null(til_exprs))
      {
        til_average <- mean(til_exprs[,channel])
      }
      
      panel_measurements[i,paste(channel, "_average_", "tum", sep="")] = tum_average
      panel_measurements[i,paste(channel, "_average_", "til", sep="")] = til_average
    }
  }

  return(panel_measurements)
}

#To do: this 
measure_duncan_panel_3 <- function(frames)
{
  
}