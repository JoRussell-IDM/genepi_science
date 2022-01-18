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
    bind_rows(lapply(tbl_df(agg_json[[i]]), function(j){
      map_df(j, as_tibble, .id="grouping") %>%
        tidyr::separate(grouping, cols, sep="_")}
    ), .id="genotype")
  })
  agg_df <- bind_rows(agg_list)
  return(agg_df)
}


ibx_inf <- function(inf_json){
  bind_rows(lapply(setNames(names(inf_json), names(inf_json)), function(i) 
    bind_rows(map_df(inf_json[[i]], as_tibble, .id="infIndex"))),
    .id="genotype")
}
  
# load files
ibx_json2df <- function(json){
  full_json <- jsonlite::fromJSON(json)
  agg_df <- ibx_aggregate(full_json$aggregate)
  inf_df <- "Empty"
  if("per_infIndex" %in% names(full_json)){
    inf_df <- ibx_inf(full_json$per_infIndex)
  } 
  return(list(agg = agg_df, inf = inf_df))
}   


ibxDist_json2df <- function(json){
  full_json <- jsonlite::fromJSON(json)
  
  agg_ibx <- full_json$aggregate
  pop_ibx <- lapply(setNames(names(agg_ibx), names(agg_ibx)), function(i){
    columns = agg_ibx[[i]]$columns
    agg_ibx[[i]]$columns <- NULL
    tmp_df <- bind_rows(lapply(setNames(names(agg_ibx[[i]]), names(agg_ibx[[i]])), function(j){
      data.frame(count = unlist(agg_ibx[[i]][[j]])) %>% 
        tibble::rownames_to_column("ibx")
    }), .id="grouping") %>%
      tidyr::separate(grouping, columns, sep="_")
  })
  
  inf_ibx <- full_json$per_infIndex
  ctxm_ibx <- bind_rows(lapply(setNames(names(inf_ibx), names(inf_ibx)), function(j){
    data.frame(count = unlist(inf_ibx[[j]])) %>% 
      tibble::rownames_to_column("ibx")
  }), .id="infIndex")
  return(list(pop_ibx = pop_ibx, inf_ibx = ctxm_ibx))
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

mean_pop_plot <- function(df){
  df <- df %>%
    dplyr::mutate(initial = gsub("[[:alpha:]]", "", genotype)) %>%
    tidyr::separate(initial, c("variants", "af", "seed"), sep="_") %>%
    dplyr::mutate(variants = ifelse(genotype == "interval", "All", variants),
                  af = ifelse(genotype == "interval", "NA", af)) 
    if("All" %in% unique(df$variants)){
      variants = gtools::mixedsort(unique(df$variants))
      var_cols = setNames(c(base_cols[1:length(variants)-1], "tan4"), variants)
    } else {
      variants = sorted(unique(df$variants))
      var_cols = setNames(base_cols[1:length(variants)], variants)
    }
  
   p <- df %>%
     ggplot(aes(x = year, y=mean, ymax = mean + std, ymin = mean - std,
               color = variants, fill = variants, group = paste0(genotype, subset_replicate))) +
    geom_ribbon(alpha = 0.15, color=NA) +
    geom_line(aes(linetype=af)) +
    guides(fill=FALSE) +
    scale_color_manual(values = var_cols) + 
    scale_fill_manual(values = var_cols) +
    labs(x= "Year", 
         y = "Mean population strain relatedness\n(All variants = Root IBD, Other = IBS)",
         color = "Variants",
         linetype = "Seeding\nallele\nfrequency")
 return(p)
}

mean_cotxn_plot <- function(df){
  df <- df %>%
    dplyr::mutate(initial = gsub("[[:alpha:]]", "", genotype)) %>%
    tidyr::separate(initial, c("variants", "af", "seed"), sep="_") %>%
    dplyr::mutate(variants = ifelse(genotype == "interval", "All", variants),
                  af = ifelse(genotype == "interval", "NA", af),
                  cotxn = ifelse(grepl("_", infIndex), "Co-transmission", "Superinfection"),
                  color_var = paste(variants, af),
                  color_var = factor(color_var, levels = gtools::mixedsort(unique(color_var)))) 
    new_colors <- colors(df)
    
    # boxplot of men strain IBx
    p1 <- df %>% 
      ggplot(aes(x = cotxn, y=mean, fill = color_var)) +
      geom_boxplot() +
      scale_fill_manual(values = new_colors) + 
      labs(x= "Transmission mechanism", 
          y = "Mean infection strain relatedness\n(All variants = Root IBD, Other = IBS)",
          fill = "Variants & Seeding AF") 
    # direct line relationship
    p2 <- df %>% 
      tidyr::separate(infIndex, c("infIndex", "parentInfIndex"), sep="_") %>%
      ggplot(aes(x = cotxn, y=mean, color = color_var, group = paste(genotype, infIndex))) +
      geom_line(alpha=0.15) +
      scale_color_manual(values = new_colors) + 
      labs(x= "Transmission mechanism", 
           y = "Mean infection strain relatedness\n(All variants = Root IBD, Other = IBS)",
           color = "Variants & Seeding AF") +
      facet_grid(~color_var)
    
    legend <- get_legend(
      # create some space to the left of the legend
      p1 + guides(color = guide_legend(nrow = 1)) +
        theme(legend.position = "bottom")
    )
    
    p_merge <- cowplot::plot_grid(p1 + theme(legend.position="none"), 
                                  p2 + theme(axis.text = element_blank(), 
                                             axis.title = element_blank(),
                                            strip.background = element_blank(),
                                            strip.text.x = element_blank(),
                                            legend.position="none"),
                                  ncol = 2, axis = "bt", align = "h")
                                  
    p_legend <- cowplot::plot_grid(p_merge, legend,
                                  ncol = 1, 
                                  rel_heights = c(1,0.05))
    return(p_legend)
}

ibx_plots <- function(json_list, name){
  pop_mean <- mean_pop_plot(json_list$agg)
  ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_", name, "_ibxPopulationMean.png"), 
         plot = pop_mean,
         path = plot_output_dir,
         width = 6, height = 3, units = c("in"), dpi = 300)
  
  txn_mean <- mean_cotxn_plot(json_list$inf)
  ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_", name, "_ibxTransmissionMean.png"), 
         plot = txn_mean,
         path = plot_output_dir,
         width = 7, height = 5, units = c("in"), dpi = 300)
}

################################################################################
ibx_files <- list.files(project_dir, pattern = "-ibxSummary.json" , full.names = T, recursive = T)
ibx_list <- lapply(ibx_files, ibx_json2df)
names(ibx_list) <- gsub("-ibx.*", "", basename(ibx_files))

lapply(names(ibx_list), function(x) ibx_plots(ibx_list[[x]], x))

# individual distributions
ibx_outdir <- '/mnt/data/malaria/synthetic_genomes/testing/genotype_subset/output/Habitat_0.572/features/year'
ibxDist_files <- list.files(ibx_outdir, pattern = "-ibxDistribution.json" , full.names = T, recursive = T)
ibxDist_list <- lapply(ibxDist_files, ibxDist_json2df)
names(ibxDist_list) <- gsub("-ibx.*", "", basename(ibxDist_files))

popDist <- bind_rows(sapply(ibxDist_list, "[[", 1), .id="genotype") %>%
  tidyr::separate(genotype, c("genotype", "grouping"), sep="\\.(?=[^\\.]+$)") %>%
  dplyr::mutate(genotype = ifelse(grepl("interval", genotype), gsub("interval", "All_NA_NA", genotype), gsub("[[:alpha:]]", "", genotype))) %>%
  tidyr::separate(genotype, c("variants", "af", "seed"), sep="_") %>%
  dplyr::group_by(variants, af, seed, grouping, year, subset_replicate) %>% 
  mutate(cnt = sum(count)) %>%
  mutate(freq = round(count /cnt, 3))


popDistAll <- popDist %>% 
  dplyr::filter(year != "0.0") %>%
  ggplot(aes(x=as.numeric(ibx), y=freq, group = paste(variants, af, subset_replicate), 
                   color=variants)) +
  #geom_point(alpha=0.25) + 
  geom_line(aes(linetype=af)) +
  scale_color_manual(values=setNames(base_cols[1:length(unique(popDist$variants))], c("24", "100", "All"))) +
  facet_grid(~year) +
  labs(x= "IBx", 
       y = "Percent pairwise comparisons",
       color = "Variants",
       linetype = "Seeding\nallele\nfrequency") 
ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_", name, "_ibxPopulationDist.png"), 
       plot = popDistAll,
       path = plot_output_dir,
       width = 12, height = 4, units = c("in"), dpi = 300)


infDist <- bind_rows(sapply(ibxDist_list, "[", 2), .id="genotype")  %>%
  tidyr::separate(genotype, c("genotype", "grouping"), sep="\\.(?=[^\\.]+$)") %>%
  dplyr::mutate(genotype = ifelse(grepl("interval", genotype), gsub("interval", "All_NA_NA", genotype), gsub("[[:alpha:]]", "", genotype)),
                cotxn = ifelse(grepl("_", infIndex), "Co-transmission", "Superinfection")) %>%
  tidyr::separate(genotype, c("variants", "af", "seed"), sep="_") %>%
  tidyr::separate(infIndex, c("infIndex", "parentInfIndex"), sep="_")  %>%
  dplyr::mutate(infIndex = as.numeric(infIndex))

inf_meta <- fread("/mnt/data/malaria/synthetic_genomes/testing/genotype_subset/output/Habitat_0.572/infIndexRecursive-genomes-df.csv")
infDistAll <- infDist %>% 
  dplyr::left_join(., inf_meta, by="infIndex") %>%
  dplyr::filter(year > 0) %>%
  ggplot(aes(x=as.numeric(ibx), y=count.x, group = paste(variants, af), 
             color=variants)) +
  #geom_point(alpha=0.25) + 
  #geom_line(aes(linetype=af)) +
  geom_bar() + 
  scale_color_manual(values=setNames(base_cols[1:length(unique(popDist$variants))], c("24", "100", "All"))) +
  facet_grid(paste(af, cotxn, sep="\n")~year, scales="free") +
  labs(x= "IBx", 
       y = "Pairwise comparisons",
       color = "Variants",
       linetype = "Seeding\nallele\nfrequency")
ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_", name, "_ibxTransmissionDist.png"), 
       plot = infDistAll,
       path = plot_output_dir,
       width = 12, height = 8, units = c("in"), dpi = 300)

