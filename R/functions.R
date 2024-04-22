##### Function to enable printing messages within a pipe
pipe_message = function(.data, status) {message(status); .data}

##### Generate a table of chlorinated paraffin congeners based on input chain length range

CP.build <- function(min_C, max_C){
  CPtable_plusCl <- tibble(Cl = sequence((min_C):(max_C))) %>%
    mutate(C = rep(min_C:max_C, times = min_C:max_C)) %>% 
    mutate(H = (C*2+2-Cl)) %>%
    relocate(Cl, .after = H) %>%
    mutate(Cl = Cl+1) %>%
    mutate(cmpd_id = row_number()) %>%
    gather(key, value, -cmpd_id) %>%
    group_by(cmpd_id) %>%
    summarize (formula_plusCl = paste0(key, value, collapse = ""))
  
  CPtable <- tibble(Cl = sequence((min_C):(max_C))) %>%
    mutate(C = rep(min_C:max_C, times = min_C:max_C)) %>% 
    mutate(H = (C*2+2-Cl)) %>%
    relocate(Cl, .after = H) %>%
    mutate(cmpd_id = row_number()) %>%
    gather(key, value, -cmpd_id) %>%
    group_by(cmpd_id) %>%
    summarize (formula = paste0(key, value, collapse = "")) %>%
    left_join(CPtable_plusCl, by = "cmpd_id")
  
  #Add details
  max_abundance <- NULL
  isomass <- NULL
  cmpd_class <- NULL
  mass_percent_Cl <- NULL
  for(i in 1:nrow(CPtable)){
    isoCPs <- as_tibble(isopattern(iso_list, CPtable$formula_plusCl[i], 0.01))
    max_abundance <- c(max_abundance, max(isoCPs$abundance))
    isomass <- c(isomass, isoCPs$mass[which.max(isoCPs$abundance)])
    
    avgmass <- MolecularWeight(ListFormula(CPtable$formula[i]))
    Clnum <- str_extract(CPtable$formula[i], "Cl[:digit:]+")
    Clmass <- MolecularWeight(ListFormula(Clnum))
    mass_percent_Cl <- c(mass_percent_Cl, (Clmass / avgmass))
    
    if(i%%10==0){print(paste0(i, " out of ", nrow(CPtable)))}
  }
  
  CPtable <- CPtable %>% 
    mutate(mass = isomass, mz = mass+0.00055, iso_fraction = max_abundance, Cnum = as.integer(str_extract(formula_plusCl, "[:digit:]{1,2}"))) %>%
    mutate(cmpd_class = case_when(
      Cnum <= 13 ~ "SCCP",
      between(Cnum, 14, 17) ~ "MCCP",
      Cnum >= 18 ~ "LCCP",
      .default = "CP"
    ), mz_iso1 = mz + 1.997, mz_iso2 = mz + 3.994, mass_percent_Cl = mass_percent_Cl) %>%
    select(-c(Cnum, formula))
    #mz_iso1 and #mz_iso2 will help filter out false positives by checking for the next two 2 isotope peaks
    #This should be robust, as it's uncommon to detect CPs with less than 3 chlorine atoms
  return(CPtable)
}



##### Add specific standards to the table of chlorinated paraffins
std.Append <- function(std_table, formula, std_class){
  
  if(length(formula) != length(std_class)){
    stop("List of formulae must be the same length as list of compound classes")
    }
  
  for(i in 1:length(formula)){
    id <- nrow(std_table) + i
    iso_std <- as_tibble(isopattern(iso_list, formula[i], 0.01))
    mass <- iso_std$mass[which.max(iso_std$abundance)]
    mz <- mass + 0.00055
    iso_fraction <- max(iso_std$abundance)
    cmpd_class <- std_class[i]
    
    halogens <- str_extract_all(formula[i], "(Cl|Br)[:digit:]+")[[1]]
    
    if(any(str_detect(halogens, "Cl"))){
      halo_Cl <- halogens[str_which(halogens, "Cl")]
      Cl_num <- as.integer(str_extract(halo_Cl, "[:digit:]+"))
      if(Cl_num > 2){
        mz_iso1 <- mz + 1.997
        mz_iso2 <- mz_iso1 + 1.997
      }
      
    }else if(any(str_detect(halogens, "Br"))){
      halo_Br <- halogens[str_which(halogens, "Br")]
      Br_num <- as.integer(str_extract(halo_Br, "[:digit:]+"))
      if(Br_num > 2){
        mz_iso1 <- mz + 1.998
        mz_iso2 <- mz_iso1 + 1.998
      }
      
    }else if(all(str_detect(halogens, "Cl|Br"))){
      mz_iso1 <- mz + 1.997
      mz_iso2 <- mz_iso1 + 1.998
    }else{
      mz_iso1 <- mz
      mz_iso2 <- mz
    }
    
    std_table <- add_row(std_table, cmpd_id = id, formula_plusCl = formula[i], mass = mass, mz = mz, iso_fraction = iso_fraction, cmpd_class = cmpd_class, mz_iso1 = mz_iso1, mz_iso2 = mz_iso2, mass_percent_Cl = 1)
  }
  
  return(std_table)
}


##### Match detected peaks to list of CPs

CP.peakmatch <- function(peaklist, CPref, ppm, IC_cutoff){
  ppm <- ppm/(10^6)
  
  CPref <- CPref %>%
    rowwise() %>%
    mutate(mz_lower = (mz-(mz*ppm)), mz_upper = (mz+(mz*ppm))) %>%
    select(-mz) %>%
    ungroup()
  
  peaklist <- as_tibble(peaklist) %>%
    full_join(CPref, join_by(between(mz, mz_lower, mz_upper))) %>%
    filter(!is.na(mz))
  
  keepindex <- NULL
  for(i in 1:nrow(peaklist)){
    if(is.na(peaklist$formula_plusCl[i])){next} #Skip NA values (non-matches)
    
    mziso1 <- peaklist$mz_iso1[i]
    mziso2 <- peaklist$mz_iso2[i]
    #Check for chlorine isotope peaks to reduce false positives
    if(any(between(peaklist$mz, mziso1-(mziso1*ppm), mziso1+(mziso1*ppm)))){
      if(any(between(peaklist$mz, mziso2-(mziso2*ppm), mziso2+(mziso2*ppm)))){
        keepindex <- c(keepindex, i)
      }
    }
  }  
  
  
  peaklist <- peaklist[keepindex,] %>%
    filter(!is.na(formula_plusCl)) %>%
    select(-c(mass, mz_lower, mz_upper, mz_iso1, mz_iso2)) %>%
    group_by(cmpd_id, sample) %>%
    filter(maxo == max(maxo)) %>%
    filter(maxo > IC_cutoff) %>%
    ungroup()

  return(peaklist)
  
}

##### Run this if you get the error "Error in summary.connection(connection) : invalid connection"
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}




##### Peak detection and CP screening

CP.screen <- function(msfiles, CP_RT, CPref, ppm, IC_cutoff){
  
  data <- readMSData(msfiles, mode = "onDisk", msLevel. = 1) %!>% #Imports mzXML (or mzML) files. Similar to 'xcmsSet'. 'onDisk' mode reduces system memory usage by not storing the entire raw file in R, just the spectrum header. This speeds up analysis of large datasets.
    filterRt(rt = CP_RT) %!>% #Optional step to reduce the processing time of findChromPeaks by reducing the retention time range (in seconds) to the region where CPs are expected to elute.
    pipe_message("Data import complete") %!>%
    smooth(method = "SavitzkyGolay", verbose = TRUE, halfWindowSize = 4) %!>% #Smooths data and reduces noise. Note: I have literally no idea what the units of halfWindowSize are, but '4' gave the best IC improvement, so I'm just rolling with it. 
    pipe_message("Smoothing complete") %!>%
    pickPeaks(refineMz = "descendPeak") %!>% #Centroids profile-mode MS data (keeping only the maximum signal for each mass peak).
    pipe_message("Centroiding complete") %!>%
    findChromPeaks(param = CentWaveParam(peakwidth = c(2,10), ppm = ppm, snthresh = 10)) %!>% #Identifies chromatographic peaks (peakwidth is retention time width range in seconds, ppm is mass error, and snthresh is signal-to-noise threshold) CentWaveParam() works only for centroided data (which is why pickPeaks() was used).
      # ^This step may take up to two minutes per sample to run, depending on CPU speed and available RAM
    pipe_message("Chromatographic peak detection and integration complete") %!>%
    groupChromPeaks(param = PeakDensityParam(sampleGroups = msfiles)) %!>% #group chromatographic peaks across samples
    pipe_message("Peak grouping complete") %!>%
    fillChromPeaks() %!>% #Performs peak filling (check for missing/undetected peaks across samples)
    pipe_message("Peak filling complete") %!>%
    adjustRtime(param = PeakGroupsParam(minFraction = 0.5)) %!>% #Performs retention time correction/alignment on grouped peaks
    pipe_message("Retention time correction complete")
    
  peaklist <- data %>%
    chromPeaks() %>%
    CP.peakmatch(CPref, ppm, IC_cutoff) %>% #CPref = table of target CP formulae and mz, ppm = mass error, IC_cutoff = minimum IC to count as a peak detection
    select(c(mz, rt, into, maxo, sample, cmpd_class, formula_plusCl, iso_fraction, mass_percent_Cl))
  
  return(peaklist)
}


##### Filter out peaks which were detected in blanks, based on input S/N ratio (typically 3, 5, or 10)

blank.filter <- function(blk_table, peak_table, SN){
  blk_table <- blk_table %>%
    group_by(formula_plusCl) %>%
    filter(maxo == max(maxo)) %>%
    rowwise() %>%
    mutate(intensity_SN = maxo * SN) %>%
    select(!maxo)
  
  peak_table <- peak_table %>%
    left_join(blk_table, by = "formula_plusCl") %>%
    rowwise() %>%
    mutate(intensity_SN = replace_na(maxo)) %>%
    filter(maxo >= intensity_SN) %>%
    select(!intensity_SN)
  
  return(peak_table)
}
  




## Example XIC plotting code

#mzerror <- 5e-6
#mz_SCCP <- 446.8939
#mzr <- c(mz_SCCP-(mz_SCCP*mzerror), mz_SCCP+(mz_SCCP*mzerror))
#rtr <- c(3.46*60, 3.6*60)
  
#test <- data_centroided %>%
#       filterRt(rtr) %>%
#       filterMz(mzr)
#plot(test, type = "XIC")

#Note: 
#chromPeaks(x) can be used to pull mz, RT, and IC data for chromatographic peaks
#featureValues(x, value, method) can quickly access specific data from chromPeaks, with one column per sample
#plotChromPeakDensity(x, mz, PeakDensityParam(sampleGroups, bw, minFraction)) can be used to test different grouping parameters, especially band width (bw)

