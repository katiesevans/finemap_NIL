# script for figures for midwest worm meeting poster 2019
library(tidyverse)
library(linkagemapping)
source("~/Dropbox/AndersenLab/LabFolders/Katie/projects/chemos/scripts/NIL_genotype_plots.R")
source("~/Dropbox/AndersenLab/LabFolders/Katie/projects/chemos/scripts/NIL_phenotype_plots.R")


# load RIAIL phenotypes
load("~/Dropbox/AndersenLab/RCode/Linkage mapping/RIAILsMappings/RIAIL_data/allRIAILsregressed.Rda")

# zinc mean.EXT parents
parents <- allRIAILsregressed %>%
    dplyr::filter(condition == "zinc", trait == "mean.EXT", strain %in% c("N2", "CB4856"))

parents %>%
    ggplot(.) +
    aes(x = strain, y = phenotype, fill = strain) +
    geom_jitter(width = 0.1) +
    geom_boxplot(alpha = 0.5) +
    scale_fill_manual(values = c("N2" = "orange", "CB4856" = "blue")) +
    theme_bw(24) +
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.title = element_text(face = "bold", color = "black")) +
    labs(x = "", y = "Animal Size")

ggsave("~/Dropbox/AndersenLab/LabFolders/Katie/projects/applications/midwestworm19/parent_pheno.pdf", height = 5, width = 5)

# riail phenotypes
# only set 2
data("N2xCB4856cross")
set2 <- N2xCB4856cross$pheno %>%
    dplyr::filter(set == 2) %>%
    dplyr::mutate(strain = as.character(strain)) %>%
    dplyr::pull(strain)

riails <- allRIAILsregressed %>%
    dplyr::filter(condition == "zinc", trait == "mean.EXT", strain %in% set2)


allRIAILsregressed %>%
    dplyr::filter(condition == "zinc", trait == "mean.EXT", strain %in% c(set2, "N2", "CB4856")) %>%
    dplyr::group_by(strain) %>%
    dplyr::summarise(meanpheno = mean(phenotype)) %>%
    dplyr::mutate(type = ifelse(strain == "N2", "N2_parent", ifelse(strain == "CB4856", "CB_parent", "RIL"))) %>%
    dplyr::arrange(meanpheno) %>%
    ggplot(.) +
    aes(x = factor(strain, levels = unique(strain)), y = meanpheno, fill = type) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("N2_parent" = "orange", "CB_parent" = "blue", "RIL" = "grey")) +
    theme_bw(24) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(color = "black"),
          axis.title = element_text(face = "bold", color = "black"),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none") +
    labs(x = "Recombinant Lines", y = "Animal Length")

ggsave("~/Dropbox/AndersenLab/LabFolders/Katie/projects/applications/midwestworm19/riailpheno.pdf", height = 5, width = 14)


# load zinc mappings
load("~/Dropbox/AndersenLab/RCode/Linkage mapping/RIAILsMappings/20180829/data/zinc-GWER.chromosomal.annotated.Rda")

# reannotated with proximal CI
load_cross_obj("N2xCB4856cross_full")
map <- annotatedmap %>%
    dplyr::filter(trait == "zinc.mean.EXT") %>%
    dplyr::select(marker:iteration)
cross <- linkagemapping::mergepheno(N2xCB4856cross_full2, riails)
zincmap <- linkagemapping::annotate_lods(map, cross, cutoff = "proximal")

# plot only the two QTL we are interested in
maxlodplot_kt(zincmap %>% dplyr::filter(var_exp > 0.1 | is.na(var_exp))) +
    theme_bw(24) +
    theme(axis.text = element_text(color = "black"),
          axis.title = element_text(face = "bold", color = "black"),
          panel.grid.minor = element_blank()) +
    labs(x = "Genomic Position (Mb)", y = "LOD", title = "")

ggsave("~/Dropbox/AndersenLab/LabFolders/Katie/projects/applications/midwestworm19/lod.pdf", height = 5, width = 14)

### CHRIII ###

# RIAIL split
pxgplot_kt(cross, zincmap %>% dplyr::filter(chr == "III")) +
    theme_bw(24) +
    theme(axis.text = element_text(color = "black"),
          axis.title = element_text(face = "bold", color = "black"),
          legend.position = "none") +
    labs(x = "", y = "Animal Length", title = "")

ggsave("~/Dropbox/AndersenLab/LabFolders/Katie/projects/applications/midwestworm19/splitIII.pdf", height = 5, width = 5)

# NILs
load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/chemos/HTA_sorter/20181112_zincHTA/zincIII_pruned1a.Rda")

regressed <- easysorter::regress(zincIII_pruned1a)

# plot_genopheno(regressed, "zinc", "mean.EXT", "III", conf = c(.116449,.651305))

nilgeno_plot <- nil_plot(unique(regressed$strain), chr = "III", left.bound = 0, background = T)

geno <- nilgeno_plot[[1]] +
    theme_bw(24) +
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black", face = "bold"),
          panel.grid = element_blank()) +
    geom_vline(xintercept = c(.116449,.651305), color = "black")
geno2 <- nilgeno_plot[[4]] +
    theme_bw(24) +
    theme(legend.position = "none",
          axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank())

regressed$strain <- factor(regressed$strain, levels = unique(nilgeno_plot[[2]]$sample), 
                       labels = unique(nilgeno_plot[[2]]$sample), ordered = T)

pheno <- quick_plot_breakup_flip(regressed, "zinc", "mean.EXT") +
    theme_bw(24) +
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.title = element_text(face = "bold"),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(y = "Animal Length") +
    facet_grid(~trait)

cowplot::plot_grid(geno, geno2, pheno, nrow = 1, ncol = 3, rel_widths = c(1, 0.3, 1),  align = "h", axis = "b")

ggsave("~/Dropbox/AndersenLab/LabFolders/Katie/projects/applications/midwestworm19/NIL_III.pdf", height = 10, width = 14)

#### CHRV ####

# RIAIL split
pxgplot_kt(cross, zincmap %>% dplyr::filter(chr == "V")) +
    theme_bw(24) +
    theme(axis.text = element_text(color = "black"),
          axis.title = element_text(face = "bold", color = "black"),
          legend.position = "none") +
    labs(x = "", y = "Animal Length", title = "")

ggsave("~/Dropbox/AndersenLab/LabFolders/Katie/projects/applications/midwestworm19/splitV.pdf", height = 5, width = 5)

# NILs
load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/chemos/HTA_sorter/20190305_mtx_zincV/L1L4_v3assay_pruned.Rda")

# add eca1114 genotype
eca1114 <- nilgeno %>%
    dplyr::filter(sample == "ECA1058") %>%
    dplyr::mutate(sample = "ECA1114")

eca1114[1,]$gt_name <- "CB4856"
eca1114[1,]$gt <- 2

nilgeno <- nilgeno %>%
    dplyr::bind_rows(eca1114)

nilV <- pruned %>%
    dplyr::filter(plate %in% c(3,4), condition %in% c("water", "zinc-250")) %>%
    dplyr::mutate(control = ifelse(condition == "water", NA, control),
                  condition = ifelse(condition == "zinc-250", "zinc", condition))

regressed <- easysorter::regress(nilV)

# plot_genopheno(regressed, "zinc", "mean.EXT", "V", conf = c(8.602178,13.084931))

nilgeno_plot <- nil_plot(unique(regressed$strain), chr = "V", left.bound = 0, background = T)

geno <- nilgeno_plot[[1]] +
    theme_bw(24) +
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black", face = "bold"),
          panel.grid = element_blank()) +
    geom_vline(xintercept = c(8.602178,13.084931), color = "black")
geno2 <- nilgeno_plot[[4]] +
    theme_bw(24) +
    theme(legend.position = "none",
          axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank())

regressed$strain <- factor(regressed$strain, levels = unique(nilgeno_plot[[2]]$sample), 
                           labels = unique(nilgeno_plot[[2]]$sample), ordered = T)

pheno <- quick_plot_breakup_flip(regressed, "zinc", "mean.EXT") +
    theme_bw(24) +
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.title = element_text(face = "bold"),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(y = "Animal Length") +
    facet_grid(~trait)

cowplot::plot_grid(geno, geno2, pheno, nrow = 1, ncol = 3, rel_widths = c(1, 0.3, 1),  align = "h", axis = "b")


ggsave("~/Dropbox/AndersenLab/LabFolders/Katie/projects/applications/midwestworm19/NIL_V.pdf", height = 10, width = 14)
