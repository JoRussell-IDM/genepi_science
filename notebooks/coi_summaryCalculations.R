################################################################################
# Purpose: COI summary for EIR sweep
# Author: Jessica Ribado
# Date: 09/2021
################################################################################

################################################################################
# set-up 
################################################################################
# load libraries
for(p in c('data.table', 'dplyr', 'ggplot2','ggpubr')){
  if(!p %in% installed.packages()[,1]){
    install.packages(p)
    library(p, character.only = TRUE)
  }
  library(p, character.only = TRUE)
}

project_dir <-'/mnt/data/malaria/synthetic_genomes/outbreak_eir/root_moi_fixed'

################################################################################
eff_coi_summary <- function(file){
  tmp <- data.table::fread(file) %>% 
    purrr::possibly(rename, otherwise = .)("variants14000" = `interval-genome`) %>%
    dplyr::rename("variants0" = true_coi) 
  tmp <- tidyr::pivot_longer(tmp, cols=starts_with("variants"), names_to = "comparison", values_to = "measured_coi")
  
  tmp_mean <- dplyr::group_by(tmp, year, comparison) %>%
    dplyr::summarise(mean = mean(measured_coi, na.rm=T),
                     median = median(measured_coi, na.rm=T),
                     sd = sd(measured_coi, na.r=T))
  tmp_mono <- dplyr::mutate(tmp, coi_group = ifelse(measured_coi == 1, "mono", 'poly')) %>%
    dplyr::count(year, comparison, coi_group) %>%
    dplyr::group_by(year, comparison) %>%
    dplyr::mutate(percent = n/sum(n)) %>% 
    dplyr::filter(coi_group == "mono") %>%
    dplyr::select(year, comparison, percent)
  tmp_merge <- dplyr::left_join(tmp_mean, tmp_mono)
  tmp_merge$multiplier <- gsub(".*[Hh]abitat_(.*?)/features.*", "\\1", file)
  return(tmp_merge)
}

################################################################################
eff_coi_files  <- list.files(project_dir, pattern = "infIndex-effective-coi-counts.csv" , full.names = T, recursive = T)

eff_summary_df <- dplyr::bind_rows(lapply(eff_coi_files, eff_coi_summary)) %>%
  dplyr::mutate(comparison = gsub("[[:alpha:]]", "", comparison)) %>%
  tidyr::separate(comparison, c("variants", "af", "gtSeed"), sep='_') %>%
  dplyr::mutate(variants = ifelse(variants == 0, "true", ifelse(variants == 14000, "effective", variants))) %>%
  dplyr::rename('percent_mono' = percent)
write.table(eff_summary_df, paste(project_dir, "year_allCOI_summaryStats_addMonoPerc.txt", sep="/"), sep="\t", quote=F, row.names = F)



