################################################################################
# Purpose: Allele frequency summary functions
# Author: Jessica Ribado
# Date: 01/2022
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

project_dir <-'/mnt/data/malaria/synthetic_genomes/fixed_txn_report_test/output/summaries'
plot_output_dir <- paste(project_dir, "summaries/plots", sep="/")
if (!dir.exists(plot_output_dir)){ dir.create(plot_output_dir) }

################################################################################
# functions
################################################################################
nested_json2df <- function(json){
  cols = ""
  if("columns" %in% names(json)){
    cols = json$columns
    json$columns <- NULL
    tmp_tbl <-  as_tibble(json) %>%
      tibble::rownames_to_column("index") %>%
      tidyr::pivot_longer(-index, names_to = c("grouping"), values_to = c("value")) %>%
      tidyr::separate(grouping, cols, sep="_")
  } else{
    tmp_tbl <-  as_tibble(json) %>%
      tibble::rownames_to_column("index")
  }
  return(tmp_tbl)
}

group_json2df <- function(json){
  tmp_list <- bind_rows(
    lapply(setNames(names(json), names(json)), function(i){
      nested_json2df(json[[i]])   
  }), .id="grouping")
}        

json2groups <- function(json){
  json <- jsonlite::fromJSON(json)
  tmp_group <- lapply(setNames(names(json), names(json)), function(j){
    group_json2df(json[[j]])})
  return(tmp_group)                    
}

# plotting
base_cols <- c("#BF1B0B", "#FFC465", "#66ADE5", "#252A52")
colors <- function(df){
  vars <- sort(as.numeric(unique(df$variants)[!unique(df$variants) %in% "All"]))
  afs <- sort(as.numeric(unique(df$af)[!unique(df$af) %in% c("NA")]))
  shade_cols <- lapply(cumsum(rep(round(1/length(afs), 2), length(afs))*0.5), function(i)
    as.vector(shades::saturation(base_cols[1:length(vars)], delta(-i))))
  names(shade_cols) <- rev(afs)
  shade_df <- reshape2::melt(shade_cols)
  names(shade_df) <- c("color", "af")
  shade_df$var <- vars
  if ("All" %in% unique(df$variants)){
    shade_df <- rbind(shade_df, c("tan4", "NA", "All"))
  }
  cols <- setNames(shade_df$color, paste(shade_df$var, shade_df$af))
}
  
p_af <- function(df){
  p <- df %>% 
    dplyr::mutate(year = as.numeric(year), af = as.numeric(af), variants = as.numeric(variants)) %>%
    dplyr::filter(year > 0, subset_replicate == 0) %>%
    ggplot(aes(x=year, y=value, group = paste(genotype, index, subset_replicate),
               color=as.character(variants))) +
    geom_line(alpha=0.25) +
    scale_color_manual(values = setNames(base_cols[1:length(unique(popDist$variants))], c("24", "100", "All"))) +
    labs(x="Year", y="Alternative allele frequency") +
    facet_grid(variants~af) +
    guides(color=FALSE)
  p
}

################################################################################
het_dir <- '/mnt/data/malaria/synthetic_genomes/fixed_txn_report_test/output/'
het_files <- list.files(het_dir, pattern = "hetDistribution.json" , full.names = T, recursive = T)
het_list <- lapply(het_files, json2groups)
names(het_list) <- gsub("-", "_", basename(het_files))

# plot
af_df <- bind_rows(lapply(het_list, "[[", 2), .id="genotype") %>%
  dplyr::mutate(initial = gsub("[[:alpha:]]", "", genotype)) %>%
  tidyr::separate(initial, c("variants", "af", "seed", "het_prop"), sep="_")
ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_", name, "_perSiteAF.png"), 
       plot = p_af(af_df),
       path = plot_output_dir,
       width = 8, height = 4, units = c("in"), dpi = 300)  
