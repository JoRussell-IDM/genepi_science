################################################################################
# Purpose: Summary statistic plotting functions from Snakemake Genepi output.
# Author: Jessica Ribado
# Date: 02/2022
################################################################################

################################################################################
# set-up 
################################################################################
setwd("/home/jribado/git/genepi_science/notebooks")
      
# load libraries
for(p in c('jsonlite', 'data.table', 'dplyr', 'ggplot2','ggpubr', 'gtools')){
  if(!p %in% installed.packages()[,1]){
    install.packages(p)
    library(p, character.only = TRUE)
  }
  library(p, character.only = TRUE)
}

# load functions
source("genepi_reformatting.R")

# plotting options
options(scipen=10000)
theme_set(theme_bw())
variant_colors <- setNames(c("#B17C51FF", "#A4432DFF", "#274C4FFF", "#C6D4D6FF"),
                           c("24", "100", "1000", "All"))
frequency_colors <- setNames(c("#FA8975FF", "#C86C7CFF", "#8B5975FF", "#4E475FFF"),
                             c("0.05", "0.25", "0.5", "NA"))
transmission_colors <- setNames(#c("#D99755FF", "#933450FF", "#180C18FF"),
                            rev(c("#DE9B71", "#EFBC82", "#FBDFA2", "cornsilk")),
                             c("Very low", "Low", "Moderate", "High"))


# set variables
report_start_year = 25
# set paths
#report_dir <- "/mnt/data/malaria/synthetic_genomes/dtk_downloads/manuscript_10seed_replicates"
report_dir <- "/mnt/data/malaria/synthetic_genomes/dtk_downloads/test_fixed_transmission_reports"
summary_dir <- "/mnt/data/malaria/synthetic_genomes/fixed_txn_report_test/output/summaries"
plot_dir  <- "/mnt/data/malaria/synthetic_genomes/fixed_txn_report_test/obs_layer_plots/large_sweep" 
if (!dir.exists(plot_dir)){ dir.create(plot_dir) }

################################################################################
# functions
################################################################################
n_genome_base <- function(df, x_var="aEIR"){
  dummy_bin <- cbind.data.frame(eir_bin = c("Very low", "Low", "Moderate", "High"))
  bin_count <- dplyr::group_by(df, eir_bin) %>% 
    dplyr::summarise(count=n()) %>%
    dplyr::full_join(dummy_bin, .) %>%
    dplyr::mutate(label = paste0(eir_bin, "\n(N=", count, ")"))
  
  genome_p <- df %>%
    ggplot(aes(x=get(x_var), y = n_genomes)) + 
    scale_x_log10() +
    scale_y_log10() +
    annotate("rect", xmin = c(0, 0.1, 10, 50), xmax = c(0.1, 10, 50, Inf),
             ymin = 0, ymax = Inf,
             alpha = 0.4, fill = transmission_colors) +
    annotate("text", x = c(0.01, 1, 25, 80), y = 5000, label = bin_count$count) +
    geom_point(alpha=0.5, size=1.5, color="grey20") +
    labs(x="Average EIR", y="Total genomes")
  return(genome_p)
}


coi_base <- function(df, x_var="aEIR", y_var){
  p <- df %>%
    ggplot(aes(x=get(x_var), y=get(y_var))) +
    annotate("rect", 
             xmin = c(0, 0.1, 10, 50), xmax = c(0.1, 10, 50, Inf),
             ymin = -Inf, ymax = Inf,
             alpha = 0.4, fill = transmission_colors) +
      geom_point(aes(color=af, shape=variants), alpha=0.5) +
      scale_color_manual(values=frequency_colors) +
      scale_x_log10() +
      labs(x="Average EIR",
           color="Seedingnallele\nfrequency",
           shape="Variants") 
    
    if(!grepl("difference", y_var)){
      p = p + labs(y="Observed proportion") + 
        ylim(0,1)
    } else{
      p = p + labs(y= expression("Modeled"~Delta)) + 
        ylim(-1,1) +
        geom_hline(yintercept = 0) 
    }
    return(p)
}


ibx_base <- function(df, x_var="aEIR", y_var="mean"){
  p <- df %>%
    ggplot(aes(x=get(x_var), y = get(y_var), 
               color=af, shape = variants)) + 
    annotate("rect", xmin = c(0, 0.1, 10, 50), xmax = c(0.1, 10, 50, Inf),
             ymin = -Inf, ymax = Inf,
             alpha = 0.4, fill = transmission_colors) +
    geom_point(alpha=0.5) +
    scale_x_log10() +
    scale_color_manual(values=frequency_colors) +
    labs(x="Average EIR",
         color="Seeding\allele\nfrequency", shape="Variants")
  
  if(y_var != "model_difference"){
    p <- p + ylim(0,1) +
      geom_linerange(aes(ymin = get(y_var) - std, ymax = get(y_var) + std), alpha=0.2) +
      labs(y="Mean population IBx\n(-/+ 1 standard deviation)")
  } else{
    p <- p + ylim(-1,1) +
      geom_hline(yintercept = 0) +
      annotate("text", x = 0, y = 0.05, hjust = 0, vjust = 1, 
               label = "Model and observation parity") +
      labs(y=expression("Modeled"~Delta))
  }
  return(p)
}

ibx_by_inf <- function(df){
  variable <- ifelse(length(names(df)[grepl("model", names(df))]) == 0, "mean", "model_difference")
  p <- df %>% 
    dplyr::mutate(eir_bin = factor(eir_bin, levels = rev(c("Very low", "Low", "Moderate", "High")))) %>%
    ggplot(aes(y=eir_bin, x=get(variable), fill=af, color=af)) +
    geom_density_ridges(alpha=0.75, scale=0.95, quantile_lines = TRUE, quantiles = 2) + 
    facet_grid(variants~cotxn) +
    scale_fill_manual(values=frequency_colors) + 
    scale_color_manual(values=frequency_colors) + 
    labs(y="Transmission intensity bin",
         fill="Seeding\nallele\nfrequency",
         color="Seeding\nallele\nfrequency") 

    if(variable == "mean"){
      p <- p +
        xlim(0,1) +
        geom_vline(xintercept = 0.5, linetype="dashed") + 
        geom_vline(xintercept = 0.75, linetype="dotted") +
        labs(x="Mean infection IBx")
    } else{
      p <- p + xlim(-1,1) +
        geom_vline(xintercept = 0) +
        labs(x=expression("Mean infection IBS modeled"~Delta))
    }  
  return(p)
}


binned_eir_ridgeplots <- function(df, feature_column){
  p <- df %>%
    dplyr::mutate(eir_bin = factor(eir_bin, levels = rev(c("Very low", "Low", "Moderate", "High")))) %>%
    ggplot(aes(x=get(feature_column), y=eir_bin, fill=af, color=af)) +
    geom_density_ridges(alpha=0.75, scale = 0.95, quantile_lines = TRUE, quantiles = 2) +
    facet_grid(variants~.) +
    scale_fill_manual(values=frequency_colors) +
    scale_color_manual(values=frequency_colors, guide = F) +
    labs(x = "Observed proportion", y = "Transmission level bin",
         fill="Seeding\nallele\nfrequency",
         color="Seeding\nallele\nfrequency") 
  return(p)
}


summary_panels_all <- function(year_genomes, year_ibx, year_coi){
  if("model_difference" %in% names(year_ibx)){
    ibx_title <- "Mean population IBS"
    shared_y_lab <- expression("Modeled"~Delta)
    y_var = "model_difference"
  } else{
    ibx_title <- "Mean population IBx"
    shared_y_lab <- "Observed proportion" 
    y_var = "mean"
  }
  gen_p  <- n_genome_base(year_genomes)
  ibx_p  <- ibx_base(year_ibx, y_var = y_var)
  mono_p <- coi_base(year_coi, y_var = names(year_coi)[grepl("mono", names(year_coi))])
  uniq_p <- coi_base(year_coi, y_var = names(year_coi)[grepl("uniq", names(year_coi))])
  
  prop_p <- ggarrange(ibx_p + ggtitle(ibx_title) + 
                        theme(plot.margin = margin(0,0.25,0,0, "cm"),
                              axis.title = element_blank(),
                              plot.title = element_text(hjust = 0.5)),
                      uniq_p + ggtitle("Unique genomes") + 
                        theme(plot.margin = margin(0,0.25,0,0, "cm"),
                              axis.title = element_blank(),
                              axis.text.y = element_blank(),
                              plot.title = element_text(hjust = 0.5)),
                      mono_p + ggtitle("Monogenomic\ninfections") + 
                        theme(plot.margin = margin(0,0.25,0,0, "cm"),
                              axis.title = element_blank(),
                              axis.text.y = element_blank(),
                              plot.title = element_text(hjust = 0.5)),  
                      ncol = 3, nrow=1, widths=c(1.15, 1, 1), align = "h", 
                      common.legend = TRUE,  legend="right") 
  prop_p <- annotate_figure(prop_p, left= text_grob(shared_y_lab, rot=90)) 
  
  merged_p <- ggarrange(gen_p + 
                          ggtitle("Total genomes from\nsampled infections") + 
                          theme(plot.margin = margin(0,0,0,0, "cm"),
                                axis.title = element_blank(),
                                plot.title = element_text(hjust = 0.5)),
                        prop_p,
                        ncol = 2, widths=c(1,3.25), align = "v") 
  merged_p <- annotate_figure(merged_p, bottom = text_grob("Average EIR"))
  return(merged_p)
}
  
eir_ridge_all <- function(year_ibx, year_coi){
  if("model_difference" %in% names(year_ibx)){
    ibx_title <- "Mean population IBS"
    shared_x_lab <- expression("Modeled"~Delta)
    y_var = "model_difference"
  } else{
    ibx_title <- "Mean population IBx"
    shared_x_lab <- "Observed proportion" 
    y_var = "mean"
  }
  
  ibx_p  <- binned_eir_ridgeplots(year_ibx, y_var)
  mono_p <- binned_eir_ridgeplots(year_coi, names(year_coi)[grepl("mono", names(year_coi))])
  uniq_p <- binned_eir_ridgeplots(year_coi, names(year_coi)[grepl("uniq", names(year_coi))])
  
  merged_p <- ggarrange(ibx_p + ggtitle(ibx_title) + 
                        theme(plot.margin = margin(0,0.25,0,0, "cm"),
                              axis.title = element_blank(),
                              plot.title = element_text(hjust = 0.5),
                              strip.text.y = element_blank()),
                      uniq_p + ggtitle("Unique genomes") + 
                        theme(plot.margin = margin(0,0.25,0,0, "cm"),
                              axis.title = element_blank(),
                              axis.text.y = element_blank(),
                              plot.title = element_text(hjust = 0.5),
                              strip.text.y = element_blank()),
                      mono_p + ggtitle("Monogenomic\ninfections") + 
                        theme(plot.margin = margin(0,0,0,0, "cm"),
                              axis.title = element_blank(),
                              axis.text.y = element_blank(),
                              plot.title = element_text(hjust = 0.5)),  
                      ncol = 3, nrow=1, widths=c(1.35, 1, 1), align = "h", 
                      common.legend = TRUE,  legend="right") 
  merged_p <- annotate_figure(merged_p, 
                              bottom = text_grob(shared_x_lab),
                              left = text_grob("Transmission intensity bin", rot = 90))
  return(merged_p)
}


baseline_het <- function(df){
  baseline_df <- df %>% 
    tidyr::pivot_longer(cols = starts_with("X"), names_to = "index", 
                        values_to = "het", values_drop_na = TRUE) %>%
    group_by(eir_bin, variants, af, seed, index) %>% 
    mutate(baseline_diff = het[[1]] - het,
           variants = as.numeric(variants)) %>% ungroup() %>%
    dplyr::group_by(eir_bin, variants, af, seed, year) %>%
    dplyr::summarise(median = median(baseline_diff),
                     ribbon = quantile(baseline_diff, probs = 0.9))
  return(baseline_df)
}


het_plot <- function(df){
  p <- df %>% 
    dplyr::mutate(eir_bin = factor(eir_bin, levels = c("Very low", "Low", "Moderate", "High"))) %>%
    ggplot(aes(x=year, y=median, color = af, fill = af, group = paste(eir_bin, variants, af, seed))) + 
    geom_ribbon(aes(ymin = median - ribbon, ymax = median + ribbon), alpha=0.25, colour = NA) +
    geom_line(size=1.5, alpha=0.8) +
    scale_color_manual(values=frequency_colors) +
    scale_fill_manual(values=frequency_colors) +
    facet_grid(variants~eir_bin) +
    ylim(-1,1) +
    labs(x = "Year", y = "Median heterozygosity\n(90th percentile)",
         fill="Seeding allele\nfrequency", 
         color ="Seeding allele\nfrequency")
  return(p)
}

################################################################################
# load data
################################################################################
emod_inset <- inset_df(report_dir)
inf_df   <- load_inf(paste(summary_dir, "..", sep="/"))
coi_list <- load_coi_files(summary_dir)
het_list <- load_het_files(paste(summary_dir, "..", sep="/"))
ibx_list <- load_ibx_files(summary_dir)

# reformat dataframes
inset_eir <- dplyr::select(emod_inset, sim_id, year, Daily_EIR) %>%
  dplyr::filter(year >= report_start_year - 1) %>%
  dplyr::group_by(sim_id, year) %>%
  dplyr::summarise(aEIR = mean(Daily_EIR)*365.24) %>%
  dplyr::mutate(year = paste0(year-min(year), ".0"),
                eir_bin = ifelse(aEIR < 0.1, "Very low",
                      ifelse(aEIR < 10, "Low",
                                 ifelse(aEIR < 50, "Moderate", "High"))),
                eir_bin = factor(eir_bin, levels = c("Very low", "Low", "Moderate", "High")))

#inset_eir <-  

popIBx_df <- df_reframe(bind_rows(lapply(ibx_list, "[[", 1), .id="sim_id"), inset_eir) %>%
  dplyr::filter(is.na(age_bin))
infIBx_df <- bind_rows(lapply(ibx_list, "[[", 2), .id="sim_id") 
infIBx_df <- df_reframe(inner_join(inf_df, infIBx_df) %>% dplyr::mutate(year = paste0(year, ".0")), inset_eir)
coi_df <- df_reframe(bind_rows(lapply(coi_list, "[[", 1), .id="sim_id"), inset_eir)
clone_df <- df_reframe(bind_rows(lapply(coi_list, "[[", 2), .id="sim_id"), inset_eir)
n_genome_df <- inner_join(
  dplyr::select(bind_rows(lapply(het_list, "[[", 1), .id="sim_id"), sim_id, year, n_genomes) %>% unique(), 
  dplyr::select(inset_eir, sim_id, eir_bin, year, aEIR) %>% unique()
) 


################################################################################
# running and plotting
################################################################################
years <- c("0.0", "1.0", "2.0", "3.0", "4.0", "5.0")   
lapply(setNames(years, years), function(y){
  lapply(c("model_difference", "raw"), function(j){
    year_pop <- dplyr::filter(popIBx_df, year == !!y)
    year_coi <- dplyr::filter(coi_df, year == !!y)
    year_inf <- dplyr::filter(infIBx_df, year == !!y)
    output_suffix <- ".svg"
    
    if(j == "model_difference"){
      year_pop <- model_rebase(dplyr::filter(popIBx_df, year == !!y), "mean")
      year_coi <- inner_join(
        model_rebase(dplyr::filter(year_coi, year == !!y), "mono_perc") %>%
          dplyr::rename("mono_difference"="model_difference") %>%
          dplyr::select(-contains("perc")),
        model_rebase(dplyr::filter(year_coi, year == !!y), "uniq_perc") %>%
          dplyr::rename("uniq_difference"="model_difference") %>%
          dplyr::select(-contains("perc"))
      )
      year_inf <- model_rebase(year_inf, "mean")
      output_suffix <- "ModelRescaled.svg"
    }
    
    summary_p <- summary_panels_all(dplyr::filter(n_genome_df, year == !!y), year_pop, year_coi)
    ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_year", y, "_allSummary", output_suffix),
           plot = summary_p,
           path = plot_dir,
           width = 12, height = 4, units = c("in"), dpi = 300)
    
    inf_p <- ibx_by_inf(year_inf)
    ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_year", y, "_infIBxSummary", output_suffix), 
           plot = inf_p,
           path = plot_dir,
           width = 5, height = 4, units = c("in"), dpi = 300)
    
    eirRidge_p <- eir_ridge_all(year_pop, year_coi)
    ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_year", y, "_eirBinSummary", output_suffix), 
           plot = eirRidge_p,
           path = plot_dir,
           width = 8, height = 4, units = c("in"), dpi = 300)
  })
})  
    
    

################################################################################
# heterozygoisty
site_het <- het_plot(baseline_het(inner_join(bind_rows(lapply(het_list, "[[", 3), .id="sim_id"), 
      dplyr::select(inset_eir, sim_id, eir_bin) %>% unique()))) +
  labs(y="Median site heterozygosity\n(90th percentile)")
ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_siteHeterozygoisty.svg"), 
       plot = site_het,
       path = plot_dir,
       width = 8, height = 4, units = c("in"), dpi = 300)

samp_het <- bind_rows(lapply(het_list, "[[", 4), .id="sim_id") %>%
  dplyr::mutate(infIndex = as.integer(infIndex),
                variants = as.integer(variants)) %>%
  dplyr::inner_join(., inf_df) 


samp_summary <- samp_het %>% 
  dplyr::group_by(sim_id, variants, af, seed, pid, year) %>%
  dplyr::summarise(total_infections = n(),
              potential_het = sum(recursive_count > 1, na.rm = TRUE),
              observed_het = sum(sample_heterozygosity > 0, na.rm = TRUE)) %>%
  dplyr::mutate(heterozygous_percent = potential_het/total_infections,
                percent_observed = observed_het/total_infections,
                percent_difference = heterozygous_percent - percent_observed) %>%
  dplyr::inner_join(., dplyr::select(inset_eir, sim_id, eir_bin) %>% unique())


samp_median <- tidyr::pivot_longer(samp_summary, cols=starts_with("percent"), 
                                   names_to = "category", values_to = "percent") %>%
  dplyr::group_by(category, variants, af, year, eir_bin) %>%
  summarise(sim_median = median(percent),
            sim_sd = sd(percent))

poly_diff <- samp_median %>%
  dplyr::filter(year=="2.0" & variants==100) %>%
  group_by(eir_bin, year, variants, af) %>% 
  mutate(cum_tot = cumsum(sim_median)) %>% 
  #ggplot(aes(x=eir_bin, y=cum_tot, fill=af)) +
  #geom_col(data = . %>% filter(category=="percent_observed"), position = position_dodge(width = 0.9), alpha = 1) +
  #geom_col(data = . %>% filter(category=="percent_difference"), position = position_dodge(width = 0.9), alpha = 0.4) +
  #geom_tile(aes(y=NA_integer_, alpha = category)) + 
  #scale_alpha_manual(values = c(1,0.4)) +
  ggplot(aes(x=af, y=sim_median)) +
  geom_bar(aes(fill=category), stat = "identity", color="white") +
  # geom_errorbar(aes(ymin= cum_tot - sim_sd, 
  #                   ymax= cum_tot + sim_sd),
  #               width=.2,                    # Width of the error bars
  #               position="identity") +
  facet_grid(~eir_bin) + 
  scale_fill_manual(values=c("grey70", "grey10")) +
  labs(x="Seeding allele frequency", y="Median polygenomic percent\n(-/+ 1 SD)") +
  #theme(panel.margin = grid::unit(-1.25, "lines")) +
  guides(fill="none")
  

samp_p <- samp_het %>%
  dplyr::filter(year=="2.0" & sample_heterozygosity > 0) %>%
  dplyr::inner_join(., dplyr::select(inset_eir, sim_id, eir_bin) %>% unique()) %>%
  dplyr::mutate(eir_bin = factor(eir_bin, levels = c("Very low", "Low", "Moderate", "High"))) %>%
  ggplot(aes(x=sample_heterozygosity, y=..scaled.., color=af, group=paste(sim_id, variants, af, seed))) +
  geom_density(alpha=0.8) +
  scale_color_manual(values=frequency_colors) +
  xlim(0,1) +
  labs(x = "Percent heterozygosity\n(Polygenomic samples with at least one heterozygous site only)", y = "Scaled density",
      fill="Seeding allele\nfrequency", 
      color ="Seeding allele\nfrequency") +
  facet_grid(variants~eir_bin)
ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_sampleHeterozygoisty.svg"), 
       plot = samp_p,
       path = plot_dir,
       width = 8, height = 4, units = c("in"), dpi = 300)

lapply(setNames(years, years), function(y){
samp_p_summ <- samp_summary %>%
  dplyr::mutate(variants = as.integer(variants)) %>%
  dplyr::filter(year == !!y) %>%
  #tidyr::pivot_longer(cols = starts_with("percent"), names_to = "index", 
  #                    values_to = "het", values_drop_na = TRUE) %>%
  dplyr::mutate(eir_bin = factor(eir_bin, levels = rev(c("Very low", "Low", "Moderate", "High")))) %>%
  ggplot(aes(x=percent_observed, y=eir_bin, fill=af, color=af)) + 
  geom_density_ridges(alpha=0.75, scale = 0.95, quantile_lines = TRUE, quantiles = 2) +
  facet_grid(variants~.) +
  scale_fill_manual(values=frequency_colors) +
  scale_color_manual(values=frequency_colors) +
  labs(x = "Median samples with\nat least one heterozygous position", y = "Transmission level bin",
       fill="Seeding allele\nfrequency",
       color="Seeding allele\nfrequency" )
  
ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_year", y, "_medianHeterozygousDist.svg"), 
       plot = samp_p_summ,
       path = plot_dir,
       width = 4, height = 4, units = c("in"), dpi = 300)
})

perc <-samp_summary %>% group_by(eir_bin, year, variants, af) %>%
  summarise(median_true = median(percent_polygenomic),
            ribbon_true = quantile(percent_polygenomic, probs = 0.9),
            median_obs = median(percent_observed),
            ribbon_obs = quantile(percent_observed, probs = 0.9)) %>% ungroup() %>%
  dplyr::mutate(eir_bin = factor(eir_bin, levels = rev(c("Very low", "Low", "Moderate", "High"))),
                variants = as.integer(variants)) %>%
  ggplot(aes(x=year, y=median_obs, color=af, fill=af, group=paste(eir_bin, variants, af))) +
  geom_ribbon(aes(ymin = median_true - ribbon_true, ymax = median_true + ribbon_true), alpha=0.15, fill= "#4E475FFF", colour=NA) +
  geom_ribbon(aes(ymin = median_obs - ribbon_obs, ymax = median_obs + ribbon_obs), alpha=0.25, colour = NA) +
  geom_line() +
  geom_line(aes(y = median_true), color = "#4E475FFF") +
  scale_color_manual(values=frequency_colors) +
  scale_fill_manual(values=frequency_colors) +
  labs(x = "Year", y = "Percent polygenomic",
       fill="Seeding allele\nfrequency", 
       color ="Seeding allele\nfrequency") +
  facet_grid(variants~eir_bin) 
ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_percPolyByObservedHet.svg"), 
       plot = perc,
       path = plot_dir,
       width = 8, height = 4, units = c("in"), dpi = 300)  




  

  
