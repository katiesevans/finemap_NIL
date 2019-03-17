# load genotype data
load("data/nil_genotypes.Rda")

# function for nil genotype plot
nil_plot <- function(strains, chr, all.chr=F, section = "all", background = F, ci = 1){
    
    # add blank genotype for strains that are not sequenced yet
    if(sum(strains %in% unique(nilgeno$sample)) != length(strains)) {
        # make dataframe for each missing strain
        for(i in dplyr::setdiff(strains, unique(nilgeno$sample))) {
            blanks <- nilgeno %>%
                dplyr::filter(sample == "N2") %>%
                dplyr::distinct(chrom, start, .keep_all = T) %>%
                dplyr::mutate(sample = i, gt = 3, gt_name = "unknown")
            nilgeno <- nilgeno %>%
                dplyr::bind_rows(blanks)
        }
    }
    
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
        dplyr::filter((start > .3e6 & start < .3e6 & nil_type == "CB") |
                          (start > .3e6 & start < .3e6 & nil_type == "N2") )%>%
        dplyr::filter(gt == 2 & nil_type == "CB" | gt == 1 & nil_type == "N2")%>%
        dplyr::mutate(side = "LEFT")%>%
        dplyr::mutate(size = end - start )%>%
        dplyr::ungroup()%>%
        dplyr::arrange(desc(size))%>%
        dplyr::distinct(sample, .keep_all = T)%>%
        dplyr::filter(size > 2e4)%>% # # # remove small (likely wrong calls) around interval site
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
        dplyr::group_by(sample)%>%
        dplyr::filter(start > 0 & start < 19e6)%>%
        dplyr::filter(end > 0 | end < 19e6)%>%
        dplyr::filter(end > 0 | end < 19e6)%>%
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
        dplyr::group_by(sample)%>%
        dplyr::filter(start > 0 & start < 19e6)%>%
        dplyr::filter(end > 0 | end < 19e6)%>%
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
        dplyr::mutate(nil_type = "parent", side = NA, size = NA, gt_ct = NA) %>%
        dplyr::distinct(sample, .keep_all = T)
    
    nilsII_sort <- bind_rows(nilsII_sort_right, nilsII_sort_left) %>%
        dplyr::arrange(desc(nil_type), desc(side), size) %>%
        dplyr::bind_rows(N2CB)
    
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
    
    # plot all NILs
    nl.pl <- ggplot(nilsII)+
        geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt, size = 2))+
        facet_grid(~chrom, scales = "free",  space = "free")+
        scale_color_manual(values=c("1"="orange","2"="blue", "3"="grey"))+
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
        labs(x = "Genomic Position (Mb)", y = "")
    
    # return plots and dataframe
    return(list(nl.pl, nilsII_sort, nilsII))
}

# function for nil phenotype plot
quick_plot_breakup_flip <- function(df, cond, pltrt, geno = F, pos = NA, chr = NA) {
    
    # default dataframe if no genotype to plot on pheno plot
    phen_gen <- df %>%
        dplyr::filter(trait == pltrt, condition == cond) %>%
        dplyr::mutate(type = dplyr::case_when(as.character(strain) == "N2" ~ "N2_parent",
                                              as.character(strain) == "CB4856" ~ "CB_parent",
                                              TRUE ~ "NIL"))
    
    # show genotype if geno == TRUE
    if(!is.null(geno) && geno == TRUE) {
        
        df2 <- phen_gen %>%
            dplyr::group_by(strain) %>%
            dplyr::summarize(pheno = max(phenotype))
        
        # if multiple QTL, find the genotype at each QTL
        geno <- NULL
        genos <- NULL
        for(i in pos) {
            genodf <- nilgeno %>%
                dplyr::filter(chrom == chr, 
                              sample %in% df2$strain,
                              start < i*1e6,
                              end > i*1e6) %>%
                dplyr::distinct(sample, start, gt_name) %>%
                dplyr::mutate(geno = dplyr::case_when(gt_name == "N2" ~ "N",
                                                      gt_name == "CB4856" ~ "C",
                                                      TRUE ~ "?")) %>%
                dplyr::select(strain = sample, geno)
            genos <- dplyr::bind_rows(genos, genodf)
        }
        
        # combine genotypes (if more than one QTL)
        genodf <- genos %>%
            dplyr::group_by(strain) %>%
            dplyr::mutate(genotype = paste(geno, collapse = " ")) %>%
            dplyr::select(strain, geno = genotype) %>%
            dplyr::distinct() %>%
            dplyr::left_join(df2)
        
        phen_gen <- phen_gen %>%
            dplyr::left_join(genodf) %>%
            dplyr::mutate(strain = factor(strain, levels = levels(phen_gen$strain)))
        
    } else {
        phen_gen <- phen_gen %>%
            dplyr::mutate(geno = "", pheno = phenotype)
    }
    
    # plot
    phen_gen %>%
        ggplot2::ggplot(.) +
        ggplot2::aes(x = factor(strain),
                     y = phenotype, 
                     fill=factor(type)) +
        ggplot2::geom_jitter(size = 0.5, width = 0.1)+
        ggplot2::geom_text(aes(x = strain, y = pheno, label = geno, vjust = 1.5), size = 3) +
        ggplot2::geom_boxplot(outlier.colour = NA, alpha = 0.7)+
        ggplot2::scale_fill_manual(values = c("N2_parent" = "orange", "CB_parent" = "blue", "NIL" = "gray"))+
        ggplot2::theme_bw()+
        ggplot2::coord_flip()+
        ggplot2::theme(axis.text.x = element_text(size=10, face="bold", color="black"),
                       axis.text.y = element_blank(),
                       axis.title.x = element_text(size=12, face="bold", color="black", vjust=-.3),
                       axis.title.y = element_blank(),
                       strip.text.x = element_text(size=10, face="bold", color="black"),
                       strip.text.y = element_text(size=10, face="bold", color="black"),
                       plot.title = element_text(size=12, face="bold", vjust = 1),
                       legend.position="none",
                       panel.background = element_rect( color="black",size=1.2),
                       strip.background = element_rect(color = "black", size = 1.2),
                       panel.border = element_rect( colour = "black"))+
        ggplot2::labs(y = paste0(cond, ".", pltrt), x = "")
    
}

# function for nil phenotype plot without genotype
quick_plot_breakup_flip2 <- function(df, cond, pltrt) {
    
    # default dataframe if no genotype to plot on pheno plot
    phen_gen <- df %>%
        dplyr::filter(trait == pltrt, condition == cond)
    
    # plot
    phen_gen %>%
        ggplot2::ggplot(.) +
        ggplot2::aes(x = strain,
                     y = phenotype, 
                     fill=factor(strain)) +
        ggplot2::geom_jitter(size = 0.5, width = 0.1)+
        ggplot2::geom_boxplot(outlier.colour = NA, alpha = 0.7)+
        ggplot2::theme_bw()+
        ggplot2::coord_flip()+
        ggplot2::theme(axis.text.x = element_text(size=10, face="bold", color="black"),
                       axis.text.y = element_text(size=10, face="bold", color="black"),
                       axis.title.x = element_text(size=12, face="bold", color="black", vjust=-.3),
                       axis.title.y = element_blank(),
                       strip.text.x = element_text(size=10, face="bold", color="black"),
                       strip.text.y = element_text(size=10, face="bold", color="black"),
                       plot.title = element_text(size=12, face="bold", vjust = 1),
                       legend.position="none",
                       panel.background = element_rect( color="black",size=1.2),
                       strip.background = element_rect(color = "black", size = 1.2),
                       panel.border = element_rect( colour = "black"))+
        ggplot2::labs(y = paste0(cond, ".", pltrt), x = "")
    
}

#Function to get the significance of each strain pair using Tukey HSD
nil_stats <- function(df, cond, trt, pval = 0.05) {
    
    # from quick stats
    stat_df <- df %>%
        dplyr::filter(trait == trt, condition==cond)%>%
        dplyr::select(strain, phenotype)
    aov_res <- aov(stat_df$phenotype ~ stat_df$strain)
    statsdf <- broom::tidy(TukeyHSD(aov_res)) %>%
        dplyr::select(strain_pair = comparison, pvalue = adj.p.value) %>%
        dplyr::mutate(significant = ifelse(pvalue < pval, T, F))
    
    return(statsdf)
}

# returns a list of (1) pca object and (2) pca phenotypes
calc_pc2 <- function(pheno) {
    # calculate PC for linkage
    pc_traits <- pheno %>%
        dplyr::mutate(drugtrait = paste0(condition, ".", trait)) %>%
        dplyr::ungroup()%>%
        dplyr::mutate(well = paste(assay, round, plate, row, col, strain, sep = "-")) %>%
        dplyr::select(well, drugtrait, phenotype)%>%
        unique() %>%
        tidyr::spread(drugtrait, phenotype) %>%
        na.omit()
    
    # keep strains as rownames and remove strain
    row.names(pc_traits) <- pc_traits$well
    pc_traits <- pc_traits %>%
        dplyr::select(-well)
    
    # scale traits
    scales_pc_traits <- as.data.frame(scale(pc_traits))
    
    # create PCs
    pca_obj <- princomp(scales_pc_traits)
    
    # figure out how many PCs to keep that will explain > 90% of the variance
    # pull the total variance explained with each PC and call it "drug.cumsum"
    cumsums <- t(as.matrix(cumsum(pca_obj$sdev^2/sum(pca_obj$sdev^2))))
    cumsums <- as.data.frame(cumsums) %>%
        tidyr::gather(comp, var) %>%
        dplyr::filter(var > 0.9) 
    
    # the first component in the dataframe is the first one that goes over 90%, so we want to keep it and everything else below it
    keep <- as.numeric(stringr::str_split_fixed(cumsums$comp[1], "Comp.", 2)[,2])
    
    # keep only those PCs phenotypes for all RIAIL strains
    colnames(pca_obj$scores) <- paste(pheno$condition[1], colnames(pca_obj$scores), sep = "_")
    
    PCpheno <- data.frame(pca_obj$scores[,1:keep]) %>%
        dplyr::mutate(well = rownames(.)) %>%
        tidyr::separate(well, into = c("assay", "round", "plate", "row", "col", "strain"), by = "-") %>%
        tidyr::gather(trait, phenotype, -c(assay:strain)) %>%
        dplyr::mutate(trait = gsub("Comp.", "PC", trait),
                      phenotype = phenotype) %>%
        tidyr::separate(trait, into = c("condition", "trait"), by = "_") %>%
        dplyr::mutate(phenotype = as.numeric(phenotype))
    
    return(list(pca_obj, PCpheno))
}

# returns a list of (1) pca object and (2) pca phenotypes
predict_pc <- function(pheno, pca_obj, keep = 24) {
    # clean dataframe
    pc_traits <- pheno %>%
        dplyr::mutate(drugtrait = paste0(condition, ".", trait)) %>%
        dplyr::ungroup()%>%
        dplyr::mutate(well = paste(round, assay, plate, row, col, strain, sep = "-")) %>%
        dplyr::select(well, drugtrait, phenotype)%>%
        unique() %>%
        tidyr::spread(drugtrait, phenotype) %>%
        na.omit()
    
    # keep strains as rownames and remove strain
    row.names(pc_traits) <- pc_traits$well
    pc_traits <- pc_traits %>%
        dplyr::select(-well)
    
    # scale traits
    scales_pc_traits <- as.data.frame(scale(pc_traits))
    
    # predict PCA based on pca_obj
    pcapredict <- data.frame(predict(pca_obj, scales_pc_traits))
    colnames(pcapredict) <- paste(pheno$condition[1], colnames(pcapredict), sep = "_")
    
    pcaout <- pcapredict[,1:keep] %>%
        dplyr::mutate(well = rownames(.)) %>%
        tidyr::separate(well, into = c("round", "assay", "plate", "row", "col", "strain"), by = "-") %>%
        tidyr::gather(trait, phenotype, -c(round:strain)) %>%
        dplyr::mutate(trait = gsub("Comp.", "PC", trait),
                      phenotype = phenotype) %>%
        tidyr::separate(trait, into = c("condition", "trait"), by = "_") %>%
        dplyr::mutate(phenotype = as.numeric(phenotype))
    
    return(pcaout)
}
