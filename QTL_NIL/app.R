library(shiny)
library(tidyverse)

# Define UI for application
ui <- fluidPage(
   
   # Application title
   titlePanel("Fine-map QTL with NIL phenotypes"),
   
   # text explanation
   uiOutput("intro"),
   
   sidebarLayout(
       sidebarPanel(
           # Input: Select a file ----
           fileInput("file1", "Choose phenotype file (Rdata)",
                     accept = c(".rda", ".Rda", ".RData")),
           
           # show the options to choose trait and conditions after the user has uploaded a file
           uiOutput("choosetrait"),
           
           # show the options to choose QTL and genotype
           uiOutput("chooseqtl")
           
       ),
       
       mainPanel(
           # Show a plot of the generated distribution
           plotOutput("pheno_plot")
       )

   )
)


# load genotype data
load("data/nil_genotypes.Rda")


# function for nil genotype plot
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
            labs(x = "Genomic Position (Mb)", y = "")
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
            labs(x = "Genomic Position (Mb)", y = "")
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

# function for nil phenotype plot
quick_plot_breakup_flip <- function(df, cond, pltrt, ylab = paste0(cond, ".", pltrt),
                                    textsize = 12, titlesize = 14, pointsize = 0.5, geno = F, pos = NA, chr = NA) {
    phen_gen <- df %>%
        dplyr::filter(trait == pltrt, condition == cond) %>%
        dplyr::mutate(type = ifelse(as.character(strain) == "N2", "N2_parent", ifelse(as.character(strain) == "CB4856", "CB_parent", "NIL")))

    if(geno == F) {
        phen_gen %>%
            ggplot(.) +
            aes(x = factor(strain),
                y = phenotype, 
                fill=factor(type)) +
            geom_jitter(size = pointsize, width = 0.1)+
            geom_boxplot(outlier.colour = NA, alpha = 0.7)+
            scale_fill_manual(values = c("N2_parent" = "orange", "CB_parent" = "blue", "NIL" = "gray"))+
            theme_bw()+
            coord_flip()+
            theme(axis.text.x = element_text(size=textsize, face="bold", color="black"),
                  axis.text.y = element_blank(),
                  axis.title.x = element_text(size=titlesize, face="bold", color="black", vjust=-.3),
                  axis.title.y = element_blank(),
                  strip.text.x = element_text(size=textsize, face="bold", color="black"),
                  strip.text.y = element_text(size=textsize, face="bold", color="black"),
                  plot.title = element_text(size=titlesize, face="bold", vjust = 1),
                  legend.position="none",
                  panel.background = element_rect( color="black",size=1.2),
                  strip.background = element_rect(color = "black", size = 1.2),
                  panel.border = element_rect( colour = "black"))+
            labs(y = ylab)
    } else {
        
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
                                                      TRUE ~ "C")) %>%
                dplyr::select(strain = sample, geno)
            genos <- dplyr::bind_rows(genos, genodf)
        }
        
    # combine genotypes (if more than one QTL)
        geno <- genos %>%
            dplyr::group_by(strain) %>%
            dplyr::mutate(genotype = paste(geno, collapse = "")) %>%
            dplyr::select(strain, geno = genotype) %>%
            dplyr::distinct() %>%
            dplyr::left_join(df2)
        
       phen_gen %>%
            dplyr::left_join(geno) %>%
            dplyr::mutate(strain = factor(strain, levels = levels(phen_gen$strain))) %>%
            ggplot(.) +
            aes(x = factor(strain),
                y = phenotype, 
                fill=factor(type)) +
            geom_jitter(size = pointsize, width = 0.1)+
            geom_text(aes(x = strain, y = pheno, label = geno, vjust = 2)) +
            geom_boxplot(outlier.colour = NA, alpha = 0.7)+
            scale_fill_manual(values = c("N2_parent" = "orange", "CB_parent" = "blue", "NIL" = "gray"))+
            # scale_color_manual(values = c("N" = "orange", "C" = "blue"))+
            theme_bw()+
            coord_flip()+
            theme(axis.text.x = element_text(size=textsize, face="bold", color="black"),
                  axis.text.y = element_blank(),
                  axis.title.x = element_text(size=titlesize, face="bold", color="black", vjust=-.3),
                  axis.title.y = element_blank(),
                  strip.text.x = element_text(size=textsize, face="bold", color="black"),
                  strip.text.y = element_text(size=textsize, face="bold", color="black"),
                  plot.title = element_text(size=titlesize, face="bold", vjust = 1),
                  legend.position="none",
                  panel.background = element_rect( color="black",size=1.2),
                  strip.background = element_rect(color = "black", size = 1.2),
                  panel.border = element_rect( colour = "black"))+
            labs(y = ylab)
    }
    
}


# Define server logic required for application
server <- function(input, output) {
    
    email <- a("Katie", href="mailto:kathrynevans2015@u.northwestern.edu")
    
    # intro text
    output$intro <- renderUI({
        tagList(
            "This shiny app can be used to view NIL phenotypes and genotypes and estimate the location of a QTL. To begin,
            click the button below to open a R data file containing NIL phenotypes. Once the file has loaded, choose a condition,
            trait, and chromosome to plot. To show a veritcal line representing the QTL, check the 'Show QTL?' box and use the slider
            to move the QTL position along the chromosome. Further, click the 'Show genotype?' box to print the genotype of each NIL 
            at the QTL position on the phenotype plot to the right. Please direct all questions or comments to", email
        )
        
    })
    
    # give user options for choosing condition and trait once phenotype dataframe is uploaded
    output$choosetrait <- renderUI({
        # only show if phenodf is not NULL
        req(input$file1)
        
        # load phenotype data
        assign('phenodf', get(load(input$file1$datapath)))
        
        tagList(
            
            # choose condition
            selectInput("condition", "Select condition to plot:", choices = unique(phenodf$condition)),
            
            # choose trait
            selectInput("trait", "Select trait to plot:", choices = unique(phenodf$trait)),
            
            # radio button input for chromosome to plot
            radioButtons("chrom", "Choose chromosome to plot:", 
                         choices = c("I", "II", "III", "IV", "V", "X"),
                         inline = TRUE),
            
            # check box for showing QTL on the phenotye plots
            checkboxInput("showqtl", "Show QTL?"),
            
            uiOutput("number_qtl"),
            
            uiOutput("show_qtl_pos"),
            
            # check box for showing genotypes on the phenotye plots
            checkboxInput("showgeno", "Show genotype?")
            
        )

    })
    
    # how many QTL?
    output$number_qtl <- renderUI({
        # only show slider bar if the show QTL checkbox is checked
        if(input$showqtl == T) {
            
            # how many qtl?
            numericInput("numQTL", "How many QTL in your model?", 1)

        }
        
    })
    
    # show slider bars
    output$show_qtl_pos <- renderUI({
        # only show slider bar if the show QTL checkbox is checked
        if(input$showqtl == T) {
            
            # how many qtl in the model?
            qtls <- input$numQTL
            
            # make slider for each QTL
            lapply(1:qtls, function(i) {
                tagList(
                    sliderInput(glue::glue("qtlpos{i}"), glue::glue("Choose position of QTL {i}"),
                                min = 0, max = 20, value = 10, step = 0.1)
                )
            })
        }
        
    })
    
    
    # plot nil phenotypes given user input data
    output$pheno_plot <- renderPlot({
        
        # only show if phenodf is not NULL
        req(input$file1)
        
        assign('phenodf', get(load(input$file1$datapath)))
        
        # strains from dataframe
        strains <- unique(phenodf$strain)
        
        # plot genotype with input chrom and qtl position as vertical line
        geno <- nil_plot(strains, input$chrom, left.bound = 0)
        
        # show QTL if it is clicked
        if(input$showqtl == T) {
            # how many qtl?
            qtls <- input$numQTL
            
            # get all QTL locations
            vals <- NULL
            for(i in 1:qtls) {
                vals <- c(vals, input[[glue::glue("qtlpos{i}")]])
            }

            genoplot <- geno[[1]] +
                ggplot2::geom_vline(xintercept = vals)
            
        } else {
            genoplot <- geno[[1]]
        }
        
        
        # plot data
        phenodf$strain <- factor(phenodf$strain, levels = unique(geno[[2]]$sample),
                                 labels = unique(geno[[2]]$sample))
        pheno <- quick_plot_breakup_flip(phenodf %>% dplyr::ungroup(), 
                                         input$condition, 
                                         input$trait, 
                                         geno = input$showgeno, 
                                         pos = vals, 
                                         chr = input$chrom) + 
            ggplot2::facet_grid(~trait)
        
        cowplot::plot_grid(genoplot, pheno)
        
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

