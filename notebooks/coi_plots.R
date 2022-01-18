################################################################################
# Purpose: Comparison plots of Genepi COI and THE REAL McCOIL estimations
# Author: Jessica Ribado
# Date: 05/2021

# Usage example: parallel Rscript ~/git/tskit_test/snakemake_parallelization/scripts/coi_plots.R -i xoutbreak_{}_seed_0/features/real_mccoil/hetPropNone -c xoutbreak_{}_seed_0/features/infIndex-effective-coi-counts.csv -o plots_210810/ ::: 0 1 2 5 10 20
################################################################################

################################################################################
# set-up 
################################################################################
# load libraries
for(p in c('optparse', 'data.table', 'dplyr', 'ggplot2', 'ggridges', 'ggpubr')){
  if(!p %in% installed.packages()[,1]){
    install.packages(p)
    library(p, character.only = TRUE)
  }
  library(p, character.only = TRUE)
}

# plotting options
theme_set(theme_bw())
width=8/2
height=6/1.5
units=c("in") 
dpi=300


################################################################################
# functions
################################################################################
is_all_numeric <- function(x) {
  !any(is.na(suppressWarnings(as.numeric(na.omit(x))))) & is.character(x)
}

read_het_files <- function(het_type){

  if (het_type == "site"){
    het_files <- list.files(paste(opt$input_directory, "..", sep="/"), pattern = "-siteHet.txt", full.names = T)
  } else {
    het_files <- list.files(paste(opt$input_directory, "..", sep="/"), pattern = "-sampleHet.txt", full.names = T)
  }
  
  het_df <- dplyr::bind_rows(lapply(het_files, function(x) {
    tmp <- data.table::fread(x)
    tmp$maf <- sapply(tmp$name, pid_maf_map)
    tmp$name = gsub("_af(.*?)_gt", "_gt", tmp$name)
 
    tmp <- tmp %>%
      dplyr::mutate(name = gsub("[[:alpha:]]", "",  gsub("-", "_", name))) %>%
      tidyr::separate(name, c("variants", "seed", "perc_poly", "subset_replicate", "pid", "year"), sep="_") %>%
      dplyr::mutate(
        perc_poly = ifelse(perc_poly == "", "None", perc_poly)) %>%
      dplyr::mutate_if(is_all_numeric,as.numeric) %>%
      dplyr::mutate_if(is.numeric, round, digits=3)
  })) 
  return(het_df)
}

read_genepi_coi <- function(file){
  tmp <- data.table::fread(file) %>%
    dplyr::rename(perc_poly = poly_fraction) %>%
    dplyr::mutate(name = gsub("\\-perc.*", "", basename(file)),
                  name = gsub("[[:alpha:]]", "", name))
  tmp <- tidyr::separate(tmp, name, c("variants", "maf", "seed"), sep="_|-") %>%
    dplyr::mutate(
      perc_poly = ifelse(perc_poly=="", "None", perc_poly),
      maf = sapply(tmp$name, pid_maf_map)) %>% 
    dplyr::mutate_if(is_all_numeric,as.numeric) %>%
  return(tmp)
}

read_rmc_coi <- function(file){
  awk_string = paste("gawk '$2==\"C\" {print $1,$3,$4,$5}'", file)
  tmp <- data.table::fread(cmd=awk_string) 
  names(tmp) <- c("name", "infIndex", "median_coi", "mean_coi")
  tmp$maf  <- sapply(tmp$name, pid_maf_map)
  tmp$name <- gsub("_af(.*?)_gt", "_gt", tmp$name)
  
  tmp <- dplyr::mutate(tmp, name = gsub("[[:alpha:]]", "",  gsub("-", "_", name)),) %>%
  tidyr::separate(name, c("variants", "seed", "perc_poly", "subset_replicate", "pid", "year", "n_iter", "err_model", "e1", "e2"), sep="_") %>%
  dplyr::mutate(perc_poly = as.character(ifelse(perc_poly=="", "None", perc_poly))) %>%
  return(tmp)
}

pid_maf_map <- function(x){
  pop_mafs <- stringr::str_split(gsub("pid[0-9]-", "", gsub(".*_af|_gt.*", "\\2", x)), pattern="-")[[1]]
  counts <- sapply(1:length(pop_mafs),    function(y){
    stringr::str_count(pattern = paste0("pid", y-1), x)
  })
  id_index <- which.max(counts)
  maf <- pop_mafs[id_index]
  return(maf)
}

year_subset_mccoil <- function(file_list, year){
  year_file <- file_list[grepl(paste0("year", year), file_list)]
  year_df <- dplyr::bind_rows(lapply(year_file, read_rmc_coi)) %>% 
    dplyr::mutate_if(is_all_numeric,as.numeric)
  return(year_df)
}

year_subset_het <- function(het_df, year){
  year_df <- dplyr::filter(het_df, year == !!year)
  return(year_df)
}

combine_coi_plot <- function(coi_merge, file_prefix){
  
  coi_category <- dplyr::mutate(coi_merge, 
    eff_coi_group = ifelse(hash_coi < 3, hash_coi, "> 2"),
    rmc_coi_group = ifelse(mean_coi < 3, mean_coi, "> 2"),   
    coi_group_match = ifelse(eff_coi_group == rmc_coi_group, "Match",
                        ifelse(eff_coi_group > rmc_coi_group, "Overestimate", "Underestimate")))
  coi_parity <- coi_category %>% dplyr::filter(perc_poly == "None") %>% unique() %>%
    dplyr::mutate(eff_coi_group = factor(eff_coi_group, levels=c('1', '2', '> 2'))) %>%
    ggplot(aes(x=eff_coi_group, y=..count.., fill=coi_group_match)) +
    geom_bar(position="identity") +
    facet_grid(maf~variants) +
    scale_fill_manual(name="Estimated parity", values=c('cornsilk2', 'tan3', 'tan')) +
    labs(x="Effective COI\n(Genepi)", y="Sampled infections (max=2000)") +
    ylim(0, 1000) +
    guides(fill = guide_legend(ncol = 1)) +
    theme(legend.position="top")
  ggsave(paste0(file_prefix, "-coiParity_polyNone.png"), path=opt$output_directory, 
         plot=coi_parity, width=width, height=height, units=units, dpi=dpi)
  
  coi_par_percent <- coi_category %>% dplyr::filter(perc_poly == "None") %>% unique() %>%
    dplyr::mutate(eff_coi_group = factor(eff_coi_group, levels=c('1', '2', '> 2')),
                  rmc_coi_group = factor(rmc_coi_group, levels=c('1', '2', '> 2'))) %>%
    dplyr::group_by(eff_coi_group, rmc_coi_group, coi_group_match, maf, variants) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::group_by(eff_coi_group, maf, variants) %>%
    dplyr::mutate(percent = count/sum(count)*100)
  
  coi_circle <- coi_par_percent %>%
    ggplot(aes(x=eff_coi_group, y=rmc_coi_group, fill=coi_group_match, size=percent)) +
    geom_point(position="identity", shape=21) +
    facet_grid(maf~variants) +
    scale_fill_manual(name="Estimated\nparity", values=c('cornsilk2', 'tan3', 'tan')) +
    labs(x="Effective COI category\n(Genepi)", y="Estimated COI category\n(The REAL McCOIL)",
         size="Percent\ninfections\nassigned") +
    guides(fill = guide_legend(ncol = 1),
           size = guide_legend(ncol = 1)) +
    theme(legend.position="top")
  ggsave(paste0(file_prefix, "-coiParity_polyNone_circle.png"), path=opt$output_directory, 
         plot=coi_circle, width=width, height=height+1, units=units, dpi=dpi)
  
  coi_merge <- dplyr::mutate(coi_merge, variants = factor(variants, levels=c("24", "100")))
  coi_plot <- dplyr::filter(coi_merge, perc_poly == "None") %>%
    ggplot(coi_merge, aes(x=hash_coi, y=mean_coi)) +
    #geom_point(aes(color=as.character(maf)), alpha=0.05) +
    geom_abline(slope=1, intercept=0) +
    facet_grid(.~variants) +
    geom_smooth(aes(group=paste(variants, maf, seed, perc_poly, pid, err_model, sep="_"), 
      color=as.character(maf),
      linetype=as.character(perc_poly)), 
      method = "lm", formula = y ~ splines::bs(x, 3), se=TRUE) + 
    scale_color_manual(values=c("#edae49", "#00798c")) +
    labs(x="Effective COI\n(Genepi)", y="Estimated mean COI\n(THE REAL McCOIL)",
         color="Seed MAF",
         shape="Subset\nreplicate", 
         linetype="Percent true\npolygenomic") +
    xlim(0,25) + ylim(0,25) +
    guides(color = guide_legend(ncol = 1), linetype=guide_legend(ncol = 2)) +
    theme(legend.position="top")
  ggsave(paste0(file_prefix, "-meanCoiRMC.png"), path=opt$output_directory, 
        plot=coi_plot, width=width+2, height=height, units=units, dpi=dpi)
  
  coi_plot2 <- dplyr::filter(coi_merge, perc_poly == "None") %>%
    ggplot(aes(x=hash_coi, y=mean_coi)) +
    #geom_point(aes(color=as.character(variants)), alpha=0.05) +
    geom_abline(slope=1, intercept=0) +
    facet_grid(as.numeric(maf)~.) +
    geom_smooth(aes(group=paste(variants, maf, seed, perc_poly, pid, err_model, sep="_"), 
                    color=as.character(variants),
                    linetype=as.character(perc_poly)), 
                method = "lm", formula = y ~ splines::bs(x, 3), se=TRUE) + 
    scale_color_manual(values=c("#2e4057", "#66a182")) +
    labs(x="Effective COI\n(Genepi)", y="Estimated mean COI\n(THE REAL McCOIL)",
         color="Variants",
         shape="Subset\nreplicate", 
         linetype="Percent true\npolygenomic") +
    xlim(0,25) + ylim(0,25) +
    guides(color = guide_legend(ncol = 1), linetype=guide_legend(ncol = 2)) +
    theme(legend.position="top")
  ggsave(paste0(file_prefix, "-meanCoiRMCbyMAF.png"), path=opt$output_directory, 
         plot=coi_plot2, width=width, height=height, units=units, dpi=dpi)
}

het_plot <- function(rmc_year, sample_het, file_prefix){

  het_merge <- dplyr::inner_join(rmc_year, sample_het) %>% 
    dplyr::filter(het_prop == !!het_prop)
  
  p1 <- ggplot(het_merge, aes(x=sample_het, y=mean_coi, color=perc_poly)) + 
    geom_point(alpha=0.25) + 
    stat_smooth(method=lm) +
    labs(x='Heterozygous percent', y="Esitmated COI", color = "Percent\npolygenomic") + 
    facet_grid(maf~variants) 
  ggsave(paste0(file_prefix, "_year", year, "-hetCoiComp.png"), path=opt$output_directory, 
         plot=p1, width=width, height=height, units=units, dpi=dpi)
  
  coi_het <- dplyr::inner_join(het_merge, genepi_coi_df)
  coi_het <- dplyr::mutate(coi_het, 
    eff_coi_group = ifelse(hash_coi < 3, hash_coi, "> 2"),
    rmc_coi_group = ifelse(mean_coi < 3, mean_coi, "> 2"),   
    coi_group_match = ifelse(eff_coi_group == rmc_coi_group, "Match",
        ifelse(eff_coi_group > rmc_coi_group, "Overestimate", "Underestimate")))
  coi2 <- dplyr::filter(coi_het, eff_coi_group == 2 & perc_poly == "None")
  p2 <- ggplot(coi2, aes(x=sample_het, y=perc_poly, fill=coi_group_match)) + 
    stat_density_ridges(position = 'identity', alpha=0.5, quantile_lines = TRUE, quantiles = 2) +
    facet_grid(maf~variants) +
    scale_fill_manual(name="Estimated parity", values=c('cornsilk2', 'tan3')) +
    labs(x="Sample heterozygosity") +
    guides(fill = guide_legend(ncol = 1)) +
    theme(legend.position="top",
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) 
  ggsave(paste0(file_prefix, "_year", year, "-polyPercNone-COI2-sampleHet.png"), path=opt$output_directory, 
         plot=p2, width=width/1.5, height=height, units=units, dpi=dpi)
  
  genome_df <- data.table::fread(paste0(opt$input_directory, "/../../../", "genome-df.csv")) %>%
    dplyr::mutate(infIndex = paste("infIndex", infIndex, sep="_")) 
  
  coi2_gdf <- dplyr::inner_join(coi2, genome_df, by=c("pid", "infIndex")) %>% unique() %>%
    dplyr::mutate(tx_category = ifelse(concInfIds_count == 1 & true_coi > 1, "Co-transmission", "Superinfection"))
  p3 <- ggplot(coi2_gdf, aes(x=sample_het, y=tx_category, fill = coi_group_match)) + 
    geom_density_ridges(alpha=0.5, quantile_lines = TRUE, quantiles = 2) + 
    facet_grid(maf~variants) +
    scale_fill_manual(name="Estimated parity", values=c('cornsilk2', 'tan3')) +
    labs(x="Sample heterozygosity") +
    guides(fill = guide_legend(ncol = 1)) +
    theme(legend.position="top",
          axis.title.y=element_blank())
  ggsave(paste0(file_prefix, "_year", year, "-polyPercNone-COI2-infType.png"), path=opt$output_directory, 
         plot=p3, width=width, height=height, units=units, dpi=dpi)
  
  coi2plus <- dplyr::filter(coi_het, eff_coi_group == "> 2" & perc_poly == "None")
  coi2plus_gdf <- dplyr::inner_join(coi2plus, genome_df, by=c("pid", "infIndex")) %>% unique() %>%
    dplyr::mutate(tx_category = ifelse(concInfIds_count == 1 & true_coi > 1, "Co-transmission", "Superinfection"))
  p4 <- ggplot(coi2plus_gdf, aes(x=sample_het, y=tx_category, fill = coi_group_match)) + 
    geom_density_ridges(alpha=0.5, quantile_lines = TRUE, quantiles = 2) + 
    facet_grid(maf~variants) +
    scale_fill_manual(name="Estimated parity", values=c('cornsilk2', 'tan3')) +
    labs(x="Sample heterozygosity") +
    guides(fill = guide_legend(ncol = 1)) +
    theme(legend.position="top",
          axis.title.y=element_blank())
  ggsave(paste0(file_prefix, "_year", year, "-polyPercNone-COI2plus-infType.png"), path=opt$output_directory, 
         plot=p4, width=width, height=height, units=units, dpi=dpi)
  
}

category_colors <- setNames(c("grey50", ggsci::pal_material("blue")(6)[seq(1,5,3)+1]), c(1, 2, ">2"))
effecive_coi_plots <- function(genepi_coi_df, start, stop, output_file){
  genepi_coi_df <- dplyr::mutate(genepi_coi_df, coi_match = true_coi - hash_coi)
  genepi_long <- tidyr::pivot_longer(genepi_coi_df, ends_with('_coi'), names_to = "coi_category", values_to = "coi") %>%
    dplyr::mutate(coi_group = ifelse(coi > 2, ">2", coi),
                  coi_group = factor(coi_group, levels= c(1, 2, ">2")))
  
  truth_set <- genepi_long %>%
    dplyr::filter(coi_category == "effective_coi") %>% 
    dplyr::select(month, coi_category, coi_group) %>%
    ggplot(aes(x=month, y=..count.., fill=coi_group)) + 
    geom_bar() + 
    scale_fill_manual(name="COI Category", values=category_colors) +
    labs(x="Month", y="Infections")
  
  coi_set <- genepi_long %>%
    filter(coi_category == "hash_coi", month >= start & month <= stop) %>%
    ggplot(aes(x=month, y=..count.., fill=coi_group)) + 
    geom_bar() + 
    scale_fill_manual(name="COI Category", values=category_colors) +
    facet_grid(maf~variants) +
    labs(x="Month", y="Infections")
  
  diff_tmp <- genepi_long %>%
    filter(month >= start & month <= stop) %>%
    dplyr::select(infIndex, month, variants, maf, coi_match) %>% unique() %>%
    dplyr::mutate(coi_match = ifelse(coi_match > 1, "Incorrect", "Correct"))
  diff <- dplyr::group_by_at(diff_tmp, setdiff(names(diff_tmp), "infIndex")) %>% count() %>%
    dplyr::group_by(month, variants, maf) %>%
    dplyr::mutate(percent = round(n/sum(n),2))
  diff_plot <- diff %>%
    ggplot(aes(x=month, y=percent, fill=coi_match)) +
    geom_bar(stat="identity") +
    facet_grid(maf~variants) +
    scale_fill_manual(name="", values=c("#66CC99", "#E69F00")) +
    labs(x="Month", y="COI match proportion") + 
    theme(legend.position="top")
  
  joint_coi <- ggpubr::ggarrange(truth_set + theme(legend.position="top"), 
                                coi_set + theme(legend.position="none"), 
                                diff_plot + theme(legend.position="none"), 
                                ncol = 3, nrow = 1, widths = c(2,1,1))
  ggsave(output_file, path=opt$output_directory, 
         plot=joint_coi, width=width+4*2, height=height, units=units, dpi=dpi)
  
}

sample_het_plots <- function(sample_het, file_prefix){
  
  sample_het <- dplyr::filter(sample_het, het_prop == !!het_prop)
  
  p1 <- ggplot(sample_het, aes(x=sample_het, y=year, group=paste(pid, seed, year, maf, perc_poly), fill=year)) + 
    stat_density_ridges(position = 'identity', alpha=0.5, quantile_lines = TRUE, quantiles = 2) + 
    facet_grid(maf~as.numeric(variants)) +
    labs(x="Sample heterozygosity", y="Year") +
    xlim(c(0,1)) +
    guides(fill=FALSE) +
    coord_flip()
  
  het_counts <- dplyr::group_by(sample_het, variants, seed, perc_poly, subset_replicate, pid, year, het_prop, maf) %>% count()
  p2 <- ggplot(het_counts, aes(x= reorder(variants, desc(variants)), y=n, color=year)) + 
    geom_point() + 
    facet_grid(perc_poly~.) + 
    labs(x="Variants", y="Samples with >= 1\nheterozygous position(s) (max=1000)") + 
    guides(color=FALSE)
  joint_het <- ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1, widths = c(2,1))
  ggsave(paste0(file_prefix, "-sampleHetAllYears.png"), path=opt$output_directory, 
         plot=joint_het, width=width+2, height=height, units=units, dpi=dpi)
}

site_af_plots <- function(file_prefix){
  
  af_files <- list.files(paste(opt$input_directory, "..", sep="/"), pattern = "-siteAF.txt", full.names = T)
  af_df <- dplyr::bind_rows(lapply(af_files, function(x) {
    tmp <- data.table::fread(x) %>%
      dplyr::filter(het_prop == !!het_prop)
    return(tmp) }))
    
  rmc_files <-   af_files <- list.files(opt$input_directory, pattern = "summary.txt", full.names = T)
  rmc_df <- dplyr::bind_rows(lapply(rmc_files, function(file) {
    awk_string = paste("gawk '$2==\"P\" {print $1,$3,$4}'", file)
    tmp <- data.table::fread(cmd=awk_string)
  }))
  names(rmc_df) <- c("name", "site", "median_maf")
  rmc_df <- dplyr::mutate(rmc_df,
    name =  gsub("-nIter.*", "", name),
    site = as.numeric(gsub("site", "", site)))
  
  af_merge <- dplyr::inner_join(af_df, rmc_df)
  af_merge <- dplyr::mutate(af_merge,
    pop_af = sapply(name, pid_maf_map),
    name = gsub("_af(.*?)_gt", "_gt", name)) %>%
    dplyr::mutate(name = gsub("[[:alpha:]]", "",  gsub("-", "_", name))) %>%
    tidyr::separate(name, c("variants", "seed", "perc_poly", "subset_replicate", "pid", "year"), sep="_", remove = F) %>%
    dplyr::mutate(perc_poly = ifelse(perc_poly == "", "None", perc_poly),
                  median_maf = 1 - median_maf)
 
  af_long <- tidyr::pivot_longer(af_merge, ends_with("_maf"), names_to="maf_category", values_to="af") %>% 
   dplyr::mutate(maf_category = factor(maf_category, levels=c("true_maf","adj_maf", "median_maf")))
  
  p1 <- ggplot(af_long, aes(x=af, y=year, fill=maf_category), group=paste(variants, pop_af, perc_poly, het_prop, maf_category)) + 
    stat_density_ridges(position = 'identity', alpha=0.5, quantile_lines = TRUE, quantiles = 2) +
    facet_grid(pop_af~as.numeric(variants)) +
    scale_fill_manual(name="Alelle frequency\nmethod", 
                      limits = c("true_maf", "adj_maf", "median_maf"),
                      labels = c("True", "Strain Collapase", "RMC Estimate"),
                      values = c("deepskyblue4", "lightblue2", "slateblue")) +
    labs(x="Minor allele frequency", y="Year") +
    xlim(0,1) + 
    theme(legend.position="top") +
    coord_flip()
  ggsave(paste0(file_prefix, "-siteHetAllYears.png"), path=opt$output_directory, 
         plot=p1, width=width, height=height, units=units, dpi=dpi)
  
  p2 <- ggplot(dplyr::filter(af_long, perc_poly == "None" & maf_category != "median_maf"),
               aes(x=af, y=year, fill=maf_category)) + 
    stat_density_ridges(position = 'identity', alpha=0.5, quantile_lines = TRUE, quantiles = 2) +
    facet_grid(pop_af~as.numeric(variants)) +
    scale_fill_manual(name="Alelle frequency\nmethod", 
                        limits = c("true_maf", "adj_maf"),
                        labels = c("True", "Strain Collapase"),
                        values = c("deepskyblue4", "lightblue2")) +
    labs(x="Minor allele frequency", y="Year") +
    xlim(0,1) +
    theme(legend.position="top") + 
    coord_flip()
  ggsave(paste0(file_prefix, "-sampleHetAllYears_polyNone_varVmaf.png"), path=opt$output_directory, 
         plot=p2, width=width, height=height, units=units, dpi=dpi)
                 
}


main_run <- function(genepi_df, rmc_files, year){
  
  rmc_year <- year_subset_mccoil(rmc_files, year)
  
  coi_merge <- dplyr::left_join(rmc_year, genepi_df) 
  print("Last lines of plotting file:")
  print(tail(data.frame(coi_merge)))
  
  combine_coi_plot(coi_merge, paste0(file_prefix, "_year", year))
  het_plot(rmc_year, sample_het, file_prefix)
}




################################################################################
# read command line arguments
option_list = list(
  make_option(c("-i", "--input_directory"), type="character", default=getwd(), 
    help="Directory containing THE REAL McCOIL output."),
  make_option(c("-c", "--coi_file"), type="character", 
    help="File with true and effective COI data."),  
  make_option(c("-o", "--output_directory"), type="character",
     help="Output directory for files.")       
)
opt_parser = OptionParser(option_list=option_list)
opt        = parse_args(opt_parser)

het_prop = gsub(".*hetProp(.+)", "\\1", opt$input_directory)
file_prefix = paste(sapply(strsplit(opt$input_directory, "/"), tail, 4)[1,1],
                    paste0("hetCollapse", het_prop), sep="-")
dir.create(file.path(opt$output_directory), showWarnings = FALSE)

################################################################################
# run main commands

# read in coi hash files
if (dir.exists(opt$coi_file)){
  genepi_coi_files <- list.files(paste(opt$input_directory, "..", "..", sep="/"), pattern = "variants(.*)percPoly(.*)-sampleHashCoi.txt", full.names = T)
  print(paste("Merging", length(genepi_coi_files), "Genepi COI hash files."))
  genepi_coi_df <- dplyr::bind_rows(lapply(genepi_coi_files, read_genepi_coi))                   
} else {
  genepi_coi_df <- fread(opt$coi_file) %>% 
    tidyr::pivot_longer(cols = starts_with("variants"),
                        names_to = "name", values_to="hash_coi") %>%
    dplyr::mutate(maf = sapply(name, pid_maf_map),
                  infIndex = paste("infIndex", infIndex, sep="_"),
                  name = gsub("_af(.*?)_gt", "_gt", name)) %>%
  dplyr::mutate(name = gsub("[[:alpha:]]", "",  name)) %>%
  tidyr::separate(name, c("variants", "seed"), sep="_") %>%
  dplyr::rename(effective_coi = `interval-genome`) %>%
  dplyr::mutate_if(is_all_numeric,as.numeric)
}
effecive_coi_plots(genepi_coi_df, 48, 60, paste0(sapply(strsplit(opt$input_directory, "/"), tail, 4)[1,1], "-monthCOI.png"))

# list files in directory
rmc_coi_files <- list.files(opt$input_directory, pattern = "_summary.txt", full.names = T) 

# read in heterozygous files
site_het   <- read_het_files("site")
sample_het <- read_het_files("sample")

# run plots
sample_het_plots(sample_het, file_prefix)
site_af_plots(file_prefix)
# look through the real mccoil files
for (year in unique(genepi_coi_df$year)[unique(genepi_coi_df$year) > 0]){
    main_run(genepi_coi_df, rmc_coi_files, year)
}


###############################################################################
# GR one-off plots
coi_plot2 <- dplyr::filter(coi_merge, perc_poly == "None") %>%
  ggplot(aes(x=hash_coi, y=mean_coi)) +
  #geom_point(aes(color=as.character(variants)), alpha=0.05) +
  geom_abline(slope=1, intercept=0) +
  geom_smooth(aes(group=paste(variants, maf, seed, perc_poly, pid, err_model, sep="_"), 
                  color=as.character(variants),
                  linetype=as.character(maf)), 
              method = "lm", formula = y ~ splines::bs(x, 3), se=TRUE) + 
  scale_color_manual(values=c("#2e4057", "#66a182")) +
  labs(x="Effective COI\n(Genepi)", y="Estimated mean COI\n(THE REAL McCOIL)",
       color="Variants",
       shape="Subset\nreplicate", 
       linetype="Seeding diversity") +
  xlim(0,15) + ylim(0,15) +
  coord_fixed(ratio = 1) +
  guides(color = guide_legend(ncol = 1), linetype=guide_legend(ncol = 1)) +
  theme(legend.position="top")
ggsave(paste0(file_prefix, "-meanCoiCombined_newLimit.png"), path=opt$output_directory, 
       plot=coi_plot2, width=width, height=height, units=units, dpi=dpi)