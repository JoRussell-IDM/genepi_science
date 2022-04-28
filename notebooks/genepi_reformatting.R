################################################################################
# Purpose: Functions to reformat Snakemake Genepi output for summary statistic plotting. 
# Author: Jessica Ribado
# Date: 02/2022
################################################################################

################################################################################
# shared functions
################################################################################
unnest_json <- function(json){
  #' Converts the genepi summary statistic format to a data frame. Structure example:
  #' "genotype":{
  #'     grouping:{
  #'         column_name: values 
  #'         columns: column_names
  #'     }
  #' }
  json_list <- lapply(names(json), function(i){
    cols = json[[i]]$columns
    json[[i]]$columns <- NULL 
    bind_rows(lapply(as_tibble(json[[i]]), function(j){
      map_df(j, as_tibble, .id="grouping") %>%
        tidyr::separate(grouping, cols, sep="_")}
    ), .id="genotype")
  })
  df <- bind_rows(json_list)
  return(df)
}


load_inf <- function(project_dir){
  inf_files <- list.files(project_dir, full.names = T, recursive = T, pattern="subset-indices.csv")
  inf_list <- lapply(inf_files, data.table::fread)
  names(inf_list) <- basename(gsub("subset-indices.csv", "", inf_files))  
  inf_df <- bind_rows(inf_list, .id="sim_id")
  return(inf_df)
}


df_reframe <- function(df, mapping_file){
  #' Converts column with observational model information to individual columns and merges the relevant EMOD report metadata.  
  r_df <- df %>%
    dplyr::mutate(genotype = ifelse(grepl("interval", genotype), 
        gsub("interval", "All_NA_NA", genotype), 
        gsub("[[:alpha:]]", "", genotype))) %>%
    tidyr::separate(genotype, c("variants", "af", "seed"), sep="_") %>%
    dplyr::left_join(., mapping_file) 
  r_df <- dplyr::mutate(r_df, variants = factor(variants, levels = gtools::mixedsort(unique(r_df$variants))))
  return(r_df)
}

import_reframe <- function(df, mapping_file){
  #' Converts column with observational model information to individual columns, acccount for more than one population
  r_df <- df %>%
    dplyr::mutate(genotype = gsub("pid[0-9]-", "", genotype),
      genotype = ifelse(grepl("interval", genotype), 
                                    gsub("interval", "All_NA_NA_NA", genotype), 
                                    gsub("[[:alpha:]]", "", genotype))) %>%
    tidyr::separate(genotype, c("variants", "af", "seed"), sep="_") %>%
    tidyr::separate(af, c("af_1", "af_2"), sep="-") %>%
    dplyr::left_join(., mapping_file) 
  r_df <- dplyr::mutate(r_df, variants = factor(variants, levels = gtools::mixedsort(unique(r_df$variants))))
  return(r_df)
}


model_rebase <- function(df, feature_column){
  filter_columns <- c("uniq_perc", "mono_perc", "median", "std", "variants", "af", "seed")
  truth <- dplyr::filter(df, variants == "All") %>%
    dplyr::rename('model' = feature_column)
  truth <- dplyr::select(truth, -names(truth)[names(truth) %in% filter_columns])
  observed <- dplyr::filter(df, variants != "All") %>%
    dplyr::inner_join(truth, .) %>%
    dplyr::mutate(model_difference = model - get(feature_column)) %>%
    dplyr::select(-model)
}




################################################################################
# inset json functions
################################################################################
inset_summary <- function(inset_file){
  inset_json <- read_json(inset_file)
  # remove nested list that contains description of column name 
  prune <- rrapply::rrapply(inset_json[['Channels']], condition = function(x, .xparents) !any(.xparents == "Units"), how = "prune")
  # bind all 
  channels <- do.call(cbind, lapply(prune, function(x) as.vector(unlist(x))))
  column_df <- cbind.data.frame(
    day = seq(1, nrow(channels)),
    year = floor(seq(1, nrow(channels))/365.24),
    channels
  )
  names(column_df) <- gsub(" ", "_", names(column_df))
  return(column_df)
}


inset_df <- function(report_path){
  #' Loads inset charts for EMOD output and reformats to data frame
  inset_files <- list.files(report_path, full.names = T, recursive = T, pattern="InsetChart")
  inset_data <- lapply(inset_files, inset_summary)
  names(inset_data) <- sapply(inset_files, function(i) tail(strsplit(dirname(i), "/+")[[1]], 1))
  inset_df <- dplyr::bind_rows(inset_data, .id="sim_id")
  return(inset_df)
} 


################################################################################
# allele frequency and heterozygosity specific functions
################################################################################
obs_group <- c("sampling") 

genome_count <- function(summary_json){
  count_df <- bind_rows(lapply(names(summary_json), function(i){
    nested_json <- summary_json[[i]]
    if(i == "population"){
      cols = list(c("year", "pid"))
    } else{
      cols = nested_json[['columns']]
    }
    nested_json$columns <- NULL
    nested_df <- cbind.data.frame(columns = names(nested_json), 
                                  n_genomes = as.vector(unlist(nested_json))) %>%
      tidyr::separate(columns, unlist(cols), sep="_")
  }))
  return(nested_df)
}


samp_values <- function(nested_json){
  nested_df <- cbind.data.frame(
    infIndex = names(nested_json),
    sample_heterozygosity = do.call(rbind, nested_json))
  return(nested_df)
}


site_values <- function(summary_json){
  count_df <- bind_rows(lapply(names(summary_json), function(i){
    nested_json <- summary_json[[i]]
    if(i == "population"){
      cols = list(c("year", "pid"))
    } else {
      cols = nested_json[['columns']]
    }
    nested_json$columns <- NULL
    nested_df <- data.frame(do.call(rbind, nested_json)) %>%
      tibble::rownames_to_column("columns") %>%
      tidyr::separate(columns, unlist(cols), sep="_")
  }))
  return(count_df)
}


het_json2df <- function(json){
  full_json <- jsonlite::fromJSON(json)
  genotype = basename(json)
  if(grepl("pid", genotype)){
    genotype = gsub("pid[0-9]-", "", genotype)
    name_df <- rbind.data.frame(unlist(
      str_split(ifelse(grepl("interval", genotype), 
                       gsub("interval", "All_NA_NA_NA", genotype), 
                       gsub("[[:alpha:]]", "", gsub(".json", "", genotype))), "_|-")))
    names(name_df) <- c("variants", "af1", "af2", "seed", "het_threshold")
    
  } else{
    name_df <- rbind.data.frame(unlist(
      str_split(ifelse(grepl("interval", genotype), 
                       gsub("interval", "All_NA_NA", genotype), 
                       gsub("[[:alpha:]]", "", gsub(".json", "", genotype))), "_|-")))
    names(name_df) <- c("variants", "af", "seed", "het_threshold")
  }
  
  n_genomes = cbind(name_df, genome_count(full_json[['n_genomes']]))
  allele_freq = cbind(name_df, site_values(full_json[['allele_freq']]))
  site_het = cbind(name_df, site_values(full_json[['site_het']]))
  samp_het = cbind(name_df, samp_values(full_json[['samp_het']]))
  return(list(n_genomes = n_genomes, allele_freq = allele_freq, site_het = site_het, samp_het = samp_het))
}


load_het_files <- function(summary_path){
  het_files <- list.files(summary_path, pattern = "hetDistribution.json" , full.names = T, recursive = T)
  het_list <- lapply(het_files, het_json2df)
  names(het_list) <- basename(gsub("/features.*", "", het_files))
  return(het_list)
} 


################################################################################
# coi specific functions
################################################################################
coi_json2df <- function(json){
  full_json <- jsonlite::fromJSON(json)
  
  coi_json <- full_json[['summary']]
  full_list <- list()
  full_list[['summary']] <- bind_rows(lapply(names(coi_json), function(j){
    bind_rows(lapply(setNames(names(coi_json[[j]]), names(coi_json[[j]])), function(i) {
      data.frame(t(do.call(rbind.data.frame,coi_json[[j]][[i]]))) %>% 
        tibble::rownames_to_column("group")}), .id="genotype") %>%
      tidyr::separate(group, full_json[['columns']][[j]], sep="_") %>%
      dplyr::mutate(year = gsub("X", "", year))
  }))
  
  full_list[['clones']] <- reshape2::melt(full_json[['clones']]) %>%
    dplyr::rename(year=value, hash=L2, genotype=L1)
  return(full_list)
}


load_coi_files <- function(summary_path){
  coi_files <- list.files(summary_path, pattern = "-coiSummary.json" , full.names = T, recursive = T)
  coi_list <- lapply(coi_files, coi_json2df)
  names(coi_list) <- gsub("-year.*", "", basename(coi_files))
  return(coi_list)
}

################################################################################
# ibx specific functions
################################################################################
ibx_inf <- function(inf_json){
  bind_rows(lapply(setNames(names(inf_json), names(inf_json)), function(i) 
    dplyr::bind_rows(map_df(inf_json[[i]], as_tibble, .id="infIndex"))),
    .id="genotype")
}

ibx_json2df <- function(json){
  full_json <- jsonlite::fromJSON(json)
  agg_df <- unnest_json(full_json$aggregate)
  
  inf_df <- "Empty"
  if("per_infIndex" %in% names(full_json)){
    # report_prefix <- gsub("-year.*|-month.*", "", basename(json))
    # inf_meta <- data.table::fread(paste(dirname(json), "..", report_prefix, "infIndexRecursive-genomes-df.csv", sep="/")) %>%
    #   dplyr::select(-recursive_nid) %>%
    #   dplyr::mutate(year = paste0(year, ".0"))
    inf_df <- ibx_inf(full_json$per_infIndex) %>%
      tidyr::separate(infIndex, c("infIndex", "parentInfIndex"), sep="_")  %>%
      dplyr::mutate(infIndex = as.numeric(infIndex),
                    cotxn = ifelse(is.na(parentInfIndex), "Superinfection", "Co-transmission")) #%>%
      #dplyr::left_join(., inf_meta)
  } 
  return(list(agg = agg_df, inf = inf_df))
}  

load_ibx_files <- function(summary_path){
  ibx_files <- list.files(summary_path, pattern = "-ibxSummary.json" , full.names = T, recursive = T)
  ibx_list <- lapply(ibx_files, ibx_json2df)
  names(ibx_list) <- gsub("-year.*", "", basename(ibx_files))
  return(ibx_list)
}




