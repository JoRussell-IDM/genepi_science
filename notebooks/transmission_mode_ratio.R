################################################################################
# Purpose: Testing a new metris for identifying transmission mode of polygenomic strains
# Author: Jessica Ribado
# Date: 10/2021
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

################################################################################
# functions
################################################################################
dtk_concID <- function(file){
  print(file)
  awk_string = paste("awk -vFPAT='([^,]*)|(\"[^\"]+\")' -vOFS=, '{print $15,$16}'", file, "| uniq")
  tmp <- data.table::fread(cmd=awk_string)
  names(tmp) <- c("infIndex", "ccID_count")
  tmp <- dplyr::filter(tmp, ccID_count > 0)
  tmp$name <- as.numeric(gsub(".*habitat_(.*?)/genome.*", "\\1", file))
  return(tmp)
}

recursive_concID <- function(file){
  print(file)
  awk_string = paste("awk -vFPAT='([^,]*)|(\"[^\"]+\")' -vOFS=, '{print $1,$2,$3,$4,$7}'", file)
  tmp <- data.table::fread(cmd=awk_string) 
  names(tmp) <- c("pid", "year", "infIndex", "day", "strain_count")
  tmp <- dplyr::filter(tmp, strain_count > 0)
  tmp$name <- as.numeric(gsub(".*habitat_(.*?)/inf.*", "\\1", file))
  return(tmp)
}


################################################################################
# load data
################################################################################
project_dir <-'/mnt/data/malaria/synthetic_genomes/outbreak_eir'

genome_files  <- list.files(project_dir, pattern = "genome-df.csv" , full.names = T, recursive = T)
inf_files <- list.files(project_dir, pattern = "infIndexRecursive-genomes-df.csv", full.names = T, recursive = T)

genome_df <- dplyr::bind_rows(lapply(genome_files, dtk_concID))
inf_df <- dplyr::bind_rows(lapply(inf_files, recursive_concID))
merge_df <- dplyr::inner_join(genome_df, inf_df)
write.table(merge_df, paste(project_dir, "strainCountsAll.txt", sep="/"), sep="\t", quote = F, row.names = F)

merge_df <- dplyr::mutate(merge_df, co_perc = strain_count/ccID_count) 
merge_summary <- dplyr::group_by(merge_df, pid, year, name) %>%
  dplyr::summarise(n = n(),
                   mean = mean(co_perc),
                   median = median(co_perc),
                   sd = sd(co_perc)) %>%
  dplyr::filter(year > 0)

 
coi_ratio <- ggplot(merge_summary, aes(y=mean, x=year, color=name, group=name)) + 
  geom_line() +
  labs(x="Year", y="Mean genomes/concurrent\ninfections", color="Habitat\nmultiplier")
ggsave("coiRatio.png", 
       plot = coi_ratio,
       path = paste(project_dir, "plots", sep="/"),
       width = 8, height = 4, units = c("in"), dpi = 300)


n_infections <- ggplot(merge_summary, aes(x=name, y=n, color=as.character(year))) + 
  geom_point() +
  labs(x="Habitat multiplier", y="Number of infections", color="Year")
ggsave("coiRatio_nInfections.png", 
       plot = n_infections,
       path = paste(project_dir, "plots", sep="/"),
       width = 5, height = 4, units = c("in"), dpi = 300)
  

