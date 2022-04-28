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

project_dir <-'/mnt/data/malaria/synthetic_genomes/fixed_txn_report_test/output/summaries'
plot_output_dir <- paste(project_dir, "plots", sep="/")
if (!dir.exists(plot_output_dir)){ dir.create(plot_output_dir) }

################################################################################
# functions
################################################################################
# reformatting json
coi_all <- function(agg_json){
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
  
  fractions <- names(full_json)[!names(full_json) %in% "clones"]
  full_list <- lapply(setNames(fractions, fractions), function(i){
    coi_all(full_json[[i]])
  })
  full_list[['clones']] <- reshape2::melt(full_json[['clones']]) %>%
    dplyr::rename(year=value, hash=L2, genotype=L1)
  return(full_list)
}


df_reframe <- function(df, mapping_file){
  r_df <- df %>%
    dplyr::left_join(., mapping_file) %>%
    dplyr::mutate(genotype = ifelse(grepl("interval", genotype), gsub("interval", "All_NA_NA", genotype), 
                                    gsub("[[:alpha:]]", "", genotype))) %>%
    tidyr::separate(genotype, c("variants", "af", "seed"), sep="_")
  return(r_df)
}


p_base <- function(df, year){
  p <- dplyr::filter(df, is.na(age_bin) & year == !!year) %>%
    ggplot(aes(x=aEIR, y=value, color=af)) +
    geom_point(alpha=0.5) +
    scale_color_manual(values=RColorBrewer::brewer.pal(3, "Set1")) +
    scale_x_log10() +
    labs(x="Average EIR",
         color="Seeding allele\nfrequency")
  return(p)
}

################################################################################
# load and format data
################################################################################
sim_mapping <- fread('/mnt/data/malaria/synthetic_genomes/fixed_txn_report_test/sim_id_mapping.csv')

coi_files <- list.files(project_dir, pattern = "-coiSummary.json" , full.names = T, recursive = T)
coi_list <- lapply(coi_files, coi_json2df)
names(coi_list) <- gsub("-year.*", "", basename(coi_files))

mono_df <- df_reframe(bind_rows(lapply(coi_list, "[[", 1), .id="sim_id"), sim_mapping)
uniq_df <- df_reframe(bind_rows(lapply(coi_list, "[[", 2), .id="sim_id"), sim_mapping)

################################################################################
# plotting
################################################################################
coi_plots <- lapply(setNames(unique(mono_df$year), unique(mono_df$year)), function(y){
  mono_p <- p_base(mono_df, y) + labs(y="Percent monogenomic")
  uniq_p <- p_base(uniq_df, y) + labs(y="Percent unique genomes")
  merged_p <- ggarrange(mono_p, uniq_p,  
                        align = "hv", common.legend = TRUE)
  return(merged_p)
  }
)  

# save plots
plot_dir <- '/mnt/data/malaria/synthetic_genomes/fixed_txn_report_test/obs_layer_plots' 
lapply(names(coi_plots), function(s){
  ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_variants24_year", s, "_percentSummaryStats.svg"), 
        plot = coi_plots[[s]],
        path = plot_dir,
        width = 8, height = 3, units = c("in"), dpi = 300)
  
})
   

