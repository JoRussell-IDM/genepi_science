################################################################################
# Purpose: IBx summary functions
# Author: Jessica Ribado
# Date: 10/2021
################################################################################

################################################################################
# set-up 
################################################################################
# load libraries
for(p in c('jsonlite', 'data.table', 'dplyr', 'ggplot2','ggpubr')){
  if(!p %in% installed.packages()[,1]){
    install.packages(p)
    library(p, character.only = TRUE)
  }
  library(p, character.only = TRUE)
}
theme_set(theme_bw())

project_dir <-'/mnt/data/malaria/synthetic_genomes/testing/genotype_subset/output/summaries'
plot_output_dir <- paste(project_dir, "plots", sep="/")
if (!dir.exists(plot_output_dir)){ dir.create(plot_output_dir) }

################################################################################
# functions
################################################################################
# reformatting json
ibx_aggregate <- function(agg_json){
  agg_list <- lapply(names(agg_json), function(i){
    cols = agg_json[[i]]$columns
    agg_json[[i]]$columns <- NULL 
    bind_rows(lapply(as_tibble(agg_json[[i]]), function(j){
      map_df(j, as_tibble, .id="grouping") %>%
        tidyr::separate(grouping, cols, sep="_")}
    ), .id="genotype")
  })
  agg_df <- bind_rows(agg_list)
  return(agg_df)
}

coi_json2df <- function(json){
  full_json <- jsonlite::fromJSON(json)

################################################################################
coi_files <- list.files(project_dir, pattern = "-coiSummary.json" , full.names = T, recursive = T)
coi_list <- lapply(ibx_files, ibx_json2df)
names(ibx_list) <- sapply(ibx_files, function(i) gsub("-ibx.*", "", basename(i)))
