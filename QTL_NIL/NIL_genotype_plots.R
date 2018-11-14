library(tidyverse)
library(ggplot2)

#load nilgenos
nilgeno1 <- data.table::fread("~/Dropbox/AndersenLab/Reagents/WormReagents/_SEQ/NIL/20170829/vcf/gt_hmm_fill.tsv") %>%
    dplyr::mutate(gt_name = ifelse(gt == 1, "N2", "CB4856")) %>%
    dplyr::select(chrom, start, end, sample, gt, gt_name, supporting_sites, sites, DP, switches)
nilgeno2 <- data.table::fread("~/Dropbox/AndersenLab/Reagents/WormReagents/_SEQ/NIL/20170128/hmm/gt_hmm_fill.tsv") %>%
    dplyr::mutate(gt_name = ifelse(gt == 1, "N2", "CB4856")) %>%
    dplyr::select(chrom, start, end, sample, gt, gt_name, supporting_sites, sites, DP, switches)
nilgeno3 <- data.table::fread("~/Dropbox/AndersenLab/Reagents/WormReagents/_SEQ/NIL/20171212/hmm/gt_hmm_fill.tsv") %>%
    dplyr::select(-rle)
nilgeno4 <- data.table::fread("~/Dropbox/AndersenLab/Reagents/WormReagents/_SEQ/NIL/20180409/hmm/gt_hmm_fill.tsv") %>%
    dplyr::select(-rle)
nilgeno5 <- data.table::fread("~/Dropbox/AndersenLab/Reagents/WormReagents/_SEQ/NIL/20180626/hmm/gt_hmm_fill.tsv") %>%
    dplyr::select(-rle)
nilgeno6 <- data.table::fread("~/Dropbox/AndersenLab/Reagents/WormReagents/_SEQ/NIL/20181016/hmm/gt_hmm_fill.tsv") %>%
    dplyr::select(-rle)
N2_CB_geno <- read.csv("~/Dropbox/AndersenLab/LabFolders/Katie/projects/chemos/HTA_sorter/N2_CB_geno.csv")  %>%
    dplyr::mutate(gt_name = ifelse(gt == 1, "N2", "CB4856")) %>%
    dplyr::select(chrom, start, end, sample, gt, gt_name, supporting_sites, sites, DP, switches)
nilgeno <- rbind(nilgeno1, nilgeno2, nilgeno3, nilgeno4, nilgeno5, nilgeno6, N2_CB_geno)

# df - nil genotype df in segment format
# chr - NIL chromosome
# left.cb - if looking at NIL breakups, this number corresponds to the left flank for CB4856 NILs
# left.n2 - if looking at NIL breakups, this number corresponds to the left flank for N2 NILs
# left.bound - left boundary for region of interest
# right.bound - right boundary for region of interest
# scan.range - cutoff for small genomic regions when identifying left nils
# all.chr - plot all chromosomes or just the one with the NIL
# section - "all" returns all NILs, "N2-NILs" returns only NILs N2 > CB. "CB-NILs" returns only NILs CB > N2
# background - FALSE returns the genotype of just the chromosome or all chromosomes. TRUE returns just the chrom of interest and the "genome" genotype
# ci - default is NA (no lines drawn) otherwise input vector of positions for confidence intervals of QTL

nil_plot <- function(strains, chr, left.cb = 0, left.n2 = 0, left.bound = 1, right.bound = 19e6, scan.range = 2e4, 
                     all.chr=F, section = "all", background = F, ci = 1){
    # # # determine if NILs are CB or N2
    nilsII_sort_type <- nilgeno %>%
        dplyr::filter(sample %in% strains) %>%
        dplyr::group_by(sample)%>%
        dplyr::mutate(size = end - start )%>%
        dplyr::group_by(sample, start)%>%
        dplyr::mutate(gt_ct = sum(size))%>%
        dplyr::ungroup()%>%
        dplyr::group_by(sample, gt)%>%
        dplyr::mutate(major_gt = sum(gt_ct))%>%
        dplyr::ungroup()%>%
        dplyr::arrange(desc(major_gt))%>%
        dplyr::distinct(sample, .keep_all = T)%>%
        dplyr::mutate(nil_type = ifelse(gt == 1, "CB", "N2"))%>%
        dplyr::select(sample, nil_type)

    # # # keep NILs that lost NIL genotype on right side
    nilsII_left <- nilgeno %>%
        dplyr::filter(sample %in% strains) %>%
        dplyr::filter(chrom == chr)%>%
        dplyr::left_join(.,nilsII_sort_type, by = "sample")%>%
        dplyr::filter((start > left.cb - .3e6 & start < left.cb + .3e6 & nil_type == "CB") |
                          (start > left.n2 - .3e6 & start < left.n2 + .3e6 & nil_type == "N2") )%>%
        dplyr::filter(gt == 2 & nil_type == "CB" | gt == 1 & nil_type == "N2")%>%
        dplyr::mutate(side = "LEFT")%>%
        dplyr::mutate(size = end - start )%>%
        dplyr::ungroup()%>%
        dplyr::arrange(desc(size))%>%
        dplyr::distinct(sample, .keep_all = T)%>%
        dplyr::filter(size > scan.range)%>% # # # remove small (likely wrong calls) around interval site
        dplyr::select(sample, side)

    # # # keep NILs that lost NIL genotype on left side
    nilsII_right <- nilgeno %>%
        dplyr::filter(sample %in% strains) %>%
        dplyr::filter(chrom == chr)%>%
        dplyr::left_join(.,nilsII_sort_type, by = "sample")%>%
        dplyr::filter(!sample %in% nilsII_left$sample)%>%
        dplyr::mutate(side = "RIGHT")%>%
        dplyr::distinct(sample,.keep_all = T)%>%
        dplyr::select(sample, side)

    nil_sides <- bind_rows(nilsII_left,nilsII_right)


    nilsII_sort_left <- nilgeno %>%
        dplyr::filter(sample %in% strains) %>%
        dplyr::filter(chrom == chr)%>%
        dplyr::left_join(.,nilsII_sort_type, by = "sample")%>%
        dplyr::left_join(.,nil_sides, by = "sample")%>%
        dplyr::filter(gt == 2 & nil_type == "CB" | gt == 1 & nil_type == "N2")%>%
        dplyr::group_by(sample)%>%
        dplyr::filter(start > left.bound & start < right.bound)%>%
        dplyr::filter(end > left.bound | end < right.bound)%>%
        dplyr::filter(end > left.bound | end < right.bound)%>%
        dplyr::mutate(size = end - start )%>%
        dplyr::group_by(sample, start)%>%
        dplyr::mutate(gt_ct = sum(size))%>%
        dplyr::ungroup()%>%
        dplyr::filter(side == "LEFT")%>%
        dplyr::arrange( desc(gt_ct))%>%
        dplyr::distinct(sample, .keep_all = T)%>%
        dplyr::arrange(nil_type, desc(gt_ct))

    #Here
    nilsII_sort_right <- nilgeno %>%
        dplyr::filter(sample %in% strains) %>%
        dplyr::filter(chrom == chr)%>%
        dplyr::left_join(.,nilsII_sort_type, by = "sample")%>%
        dplyr::left_join(.,nil_sides, by = "sample")%>%
        dplyr::filter(side == "RIGHT")%>%
        dplyr::filter(gt == 2 & nil_type == "CB" | gt == 1 & nil_type == "N2")%>%
        dplyr::group_by(sample)%>%
        dplyr::filter(start > left.bound & start < right.bound)%>%
        dplyr::filter(end > left.bound | end < right.bound)%>%
        dplyr::mutate(size = end - start )%>%
        dplyr::group_by(sample, start)%>%
        dplyr::mutate(gt_ct = sum(size))%>%
        dplyr::ungroup()%>%
        dplyr::arrange( desc(gt_ct))%>%
        dplyr::distinct(sample, .keep_all = T)%>%
        dplyr::ungroup()%>%
        dplyr::arrange( nil_type, gt_ct)

    N2CB <- nilgeno %>%
        dplyr::filter(sample %in% c("N2", "CB4856"), chrom == chr) %>%
        dplyr::mutate(nil_type = "parent", side = NA, size = NA, gt_ct = NA)

    nilsII_sort <- bind_rows(nilsII_sort_right, nilsII_sort_left) %>%
        arrange(desc(nil_type), desc(side), size)

    nilsII_sort <- rbind(nilsII_sort, N2CB)

    if (all.chr == T){

        nilsII <- nilgeno %>%
            dplyr::filter(sample %in% strains) %>%
            dplyr::filter(chrom != "MtDNA")%>%
            dplyr::left_join(.,nilsII_sort_type, by = "sample")

        nilsII$sample <- factor(nilsII$sample, levels = unique(nilsII_sort$sample), labels = unique(nilsII_sort$sample), ordered = T)
        nilsII$gt <- as.character(nilsII$gt)
    } else {
        nilsII <- nilgeno %>%
            dplyr::filter(sample %in% strains) %>%
            dplyr::filter(chrom == chr, chrom != "MtDNA")%>%
            dplyr::left_join(.,nilsII_sort_type, by = "sample")

        nilsII$sample <- factor(nilsII$sample, levels = unique(nilsII_sort$sample), labels = unique(nilsII_sort$sample), ordered = T)
        nilsII$gt <- as.character(nilsII$gt)
    }
    # make fake 'background' genome
    if(background == T) {
        # cannot show all chromosomes and the background
        if(all.chr == T) {
            bgplot <- NULL
        } else {
            bg <- data.frame(chrom = NA, start = NA, end = NA, sample = NA, gt = NA, gt_name = NA, supporting_sites = NA, sites = NA,
                             DP = NA, switches = NA, nil_type = NA)
            for(i in unique(nilsII$sample)) {
                # get the NIL type
                type <- nilsII %>%
                    dplyr::filter(sample == i) %>%
                    dplyr::distinct(nil_type) %>%
                    dplyr::pull(nil_type)
                
                # add a new row with the background genotype in a "new" chromosome
                bg <- rbind(bg, c("Genome", 1, 3e6, i, ifelse(type == "N2", 2, ifelse(type == "CB", 1, NA)), 
                                  ifelse(type == "N2", "CB4856", "N2"), NA, NA, NA, NA, type))
            }
            bg <- bg %>%
                tidyr::drop_na(sample) %>%
                dplyr::mutate(start = as.numeric(start), end = as.numeric(end))
            
            # add "chrom" to chromosome for plotting
            nilsII <- nilsII %>%
                dplyr::mutate(chrom = paste0("Chrom ", chrom))
            
            # background plot
            bgplot <- ggplot(bg)+
                geom_segment(aes(x = start/1e6, y = factor(sample, levels = levels(nilsII$sample)), xend = end/1e6, yend = sample, color = gt, size = 2))+
                scale_color_manual(values=c("1"="orange","2"="blue"))+
                facet_grid(~chrom, scales = "free",  space = "free")+
                theme_bw() +
                theme(axis.text.x = element_blank(),
                      axis.text.y = element_blank(),
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      strip.text = element_text(size = 12, face = "bold", color = "black"),
                      plot.title = element_text(size=24, face="bold"),
                      legend.position = "none",
                      panel.grid.minor = element_blank(),
                      panel.grid.major = element_blank(),
                      axis.ticks = element_blank())
        }
    } else { bgplot <- NULL }
    
    # only plot the N2-NILs
    if(section == "N2-NILs") {
        nl.pl <- ggplot(nilsII %>% dplyr::filter(nil_type == "CB"))+
            geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt, size = 2))+
            facet_grid(~chrom, scales = "free",  space = "free")+
            scale_color_manual(values=c("1"="orange","2"="blue"))+
            theme_bw() +
            theme(axis.text.x = element_text(size=12, face="bold", color="black"),
                  axis.text.y = element_text(size=12, face="bold", color="black"),
                  axis.title.x = element_text(size=14, face="bold", color="black"),
                  axis.title.y = element_text(size=14, face="bold", color="black"),
                  strip.text = element_text(size = 12, face = "bold", color = "black"),
                  plot.title = element_text(size=24, face="bold"),
                  legend.position = "none",
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank())+
            geom_vline(xintercept = ci, color = ifelse(ci != 1, "red", NA)) +
            labs(x = "Genomic Position (Mb)", y = "NIL")
    } else if(section == "CB-NILs") {
        # only plot the CB-NILs
        nl.pl <- ggplot(nilsII %>% dplyr::filter(nil_type == "CB"))+
            geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt, size = 2))+
            facet_grid(~chrom, scales = "free",  space = "free")+
            scale_color_manual(values=c("1"="orange","2"="blue"))+
            theme_bw() +
            theme(axis.text.x = element_text(size=12, face="bold", color="black"),
                  axis.text.y = element_text(size=12, face="bold", color="black"),
                  axis.title.x = element_text(size=14, face="bold", color="black"),
                  axis.title.y = element_text(size=14, face="bold", color="black"),
                  strip.text = element_text(size = 12, face = "bold", color = "black"),
                  plot.title = element_text(size=24, face="bold"),
                  legend.position = "none",
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank())+
            geom_vline(xintercept = ci, color = ifelse(ci != 1, "red", NA)) +
            labs(x = "Genomic Position (Mb)", y = "NIL")
    } else {
        
        # plot all NILs
        nl.pl <- ggplot(nilsII)+
            geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt, size = 2))+
            facet_grid(~chrom, scales = "free",  space = "free")+
            scale_color_manual(values=c("1"="orange","2"="blue"))+
            theme_bw() +
            theme(axis.text.x = element_text(size=12, face="bold", color="black"),
                  axis.text.y = element_text(size=12, face="bold", color="black"),
                  axis.title.x = element_text(size=14, face="bold", color="black"),
                  axis.title.y = element_text(size=14, face="bold", color="black"),
                  strip.text = element_text(size = 12, face = "bold", color = "black"),
                  plot.title = element_text(size=24, face="bold"),
                  legend.position = "none",
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank())+
            geom_vline(xintercept = ci, color = ifelse(ci != 1, "red", NA)) +
            labs(x = "Genomic Position (Mb)", y = "NIL")
    }
    
    # return plots and dataframe
    if(is.null(bgplot)) {
        return(list(nl.pl, nilsII_sort, nilsII))
    } else {
        return(list(cowplot::plot_grid(nl.pl, bgplot, nrow = 1, ncol = 2, align = "h", axis = "b", rel_widths = c(1, 0.3)),
                    nilsII_sort, 
                    nilsII))
    }
}

rm(nilgeno1, nilgeno2, nilgeno3, nilgeno4, nilgeno5, nilgeno6, N2_CB_geno)

# plots genotype on left and phenotype on right for dataframe, condition, trait
# dataframe needs condition and trait columns
plot_genopheno <- function(pheno, cond, trt, chrom, back = F, conf = 1) {
    nilgeno_plot <- nil_plot(unique(pheno$strain), chr = chrom, background = back, ci = conf, left.bound = 0)
    
    pheno$strain <- factor(pheno$strain, levels = unique(nilgeno_plot[[2]]$sample), 
                           labels = unique(nilgeno_plot[[2]]$sample), ordered = T)
    
    # Plot NIL phenotypes in condition-
    plot <- cowplot::plot_grid(nilgeno_plot[[1]], 
                               quick_plot_breakup_flip(pheno, cond, trt) + facet_grid(~trait), nrow = 1, ncol = 2)
    
    return(plot)
    
}
