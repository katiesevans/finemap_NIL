library(shiny)
library(tidyverse)
library(easysorter)
library(shinythemes)

#########################################
#       Load data & functions           #
#########################################


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


#########################################
#                   UI                  #
#########################################


# Define UI for application
ui <- fluidPage(
    
    theme = shinythemes::shinytheme('yeti'),

   # Application title
   titlePanel("Fine-map QTL with NIL phenotypes"),
   
   # text explanation
   uiOutput("intro"),
   
   # add break
   br(),
   
   # sidebar layout
   sidebarLayout(
       sidebarPanel(
           # use provided sample data
           checkboxInput("sampledata", "Use sample data"),
           
           # Input: Select a file
           fileInput("file1", "Choose phenotype file (Rdata)",
                     accept = c(".rda", ".Rda", ".RData")),
           
           # check if your data is not N2/CB NILs, won't show genotypes
           checkboxInput("n2cb", "My data is not N2/CB NILs"),
           
           # show the options to choose trait and conditions after the user has uploaded a file
           uiOutput("choosecontrol"),
           
           # choose to select certain strains
           uiOutput("choosestrains"),
           
           # show the options to choose QTL and genotype
           uiOutput("chooseqtl")
           
       ),
       
       mainPanel(
           # make tabs for each condition
           tabsetPanel(type = "tabs",
                       tabPanel("Control", uiOutput("control_plot")),
                       tabPanel("Condition", uiOutput("condition_plot")),
                       tabPanel("Regressed", uiOutput("regressed_plot")),
                       tabPanel("NIL Genotypes", uiOutput("nil_genotypes")),
                       tabPanel("Help!", uiOutput("help_page"))),
           
           # # button to output plot as png
           uiOutput("saveButton")
       )

   )
)


#########################################
#               Server                  #
#########################################


# Define server logic required for application
server <- function(input, output) {
    
    # load phenotype data
    loadPhenoData <- reactive({
        if(!is.null(input$sampledata) && input$sampledata == T) {
            assign('phenodf', get(load("data/test_NIL_pheno.Rda")))
        } else {
            req(input$file1)
            assign('phenodf', get(load(input$file1$datapath)))
        }
    })
    
    # intro text
    output$intro <- renderUI({
        p("This shiny app can be used to view NIL phenotypes and genotypes and estimate the location of a QTL. To begin,
        click the 'Browse' button to open a R data file containing pruned (non-regressed) NIL phenotypes or check the
        'Use sample data' checkbox to explore our features. For more instructions, refer to the 'Help!' tab below.")
        
    })
    
    # help page
    output$help_page <- renderUI({
        
        tagList(
            h2("Help Page"),
            p("This shiny app can be used to view NIL phenotypes and genotypes and estimate the location of a QTL."),
            tags$ul(
                tags$li("To begin,
            click the 'Browse' button to open a R data file containing pruned (non-regressed) NIL phenotypes or check the
            'Use sample data' checkbox to explore our features."),
                tags$li(p(em("Please note, if your data contains strains that are not N2/CB4856 NILs, you can click the corresponding
               checkbox to hide irrelavent options related to strain genotypes. However, no errors will result by not checking the box."))),
                tags$li(p("Once the file has loaded, choose a control, condition, and trait to plot.")),
                tags$li(p(em("Pro tip: the application will not let you choose a condition that is already named as your control, so choose your control first!"))),
                tags$li(p("The corresponding phenotype box plot should now show up on your screen, you can use the tabs along the top to 
              rotate between your control phenotype, non-regressed condition phenotype, and your control-regressed condition phenotype.")),
                tags$li(p("The strain genotypes are plotted to the left of the phenotypes. To see a specific chromosome, select a chromosome from the sidebar.")),
                tags$li(p("You may also be interested in only a subset of the strains in your assay. To view this, click the checkbox named 'Show a subset
              of strains? and then unclick the strains you do not wish to show. The application will re-plot the figures and re-run the regression analysis")),
                tags$li(p("To visualize the location of your QTL, you can check the 'Show QTL' box. A vertical line representing the QTL will appear on
              the strain genotype plot to the left. You can add more QTL and use the sliders to move the QTL across the chromosome.")),
                tags$li(p("Further, click the 'Show genotype?' box to print the genotype of each strain at the QTL position on the phenotype plot to the right.")),
                tags$li(p(em("Please direct all questions or comments to", a("Katie", href="mailto:kathrynevans2015@u.northwestern.edu"))))
                )
        )
        
    })
    
    # choose control first
    output$choosecontrol <- renderUI({
        # load phenotype data
        phenodf <- loadPhenoData()
        
        tagList(
            # choose control
            selectInput("control", "Select control:", choices = unique(phenodf$condition)),
            
            # show output for choosetrait
            uiOutput("choosetrait")
        )
        
    })
    
    # give user options for choosing condition and trait once phenotype dataframe is uploaded
    output$choosetrait <- renderUI({
        # load phenotype data
        phenodf <- loadPhenoData()
        
        # all conditions, including control
        condition_choices <- unique(phenodf$condition)
        
        tagList(
            
            # choose condition - don't allow user to choose control
            selectInput("condition", "Select condition:", 
                        choices = condition_choices[condition_choices != input$control]),
            
            # choose trait
            selectInput("trait", "Select trait:", choices = unique(phenodf$trait)),
            
            # choose subset of strains?
            checkboxInput("straininput", "Show a subset of strains?"),
            
            # check box for showing QTL on the phenotye plots
            # only if n2cb is not true
            if(!is.null(input$n2cb) && input$n2cb == FALSE) {
                tagList(
                    # radio button input for chromosome to plot
                    radioButtons("chrom", "Choose chromosome to plot:", 
                                 choices = c("I", "II", "III", "IV", "V", "X"),
                                 inline = TRUE,
                                 selected = "I"),
                    
                    checkboxInput("showqtl", "Show QTL?"),
                    
                    uiOutput("number_qtl"),
                    
                    uiOutput("show_qtl_pos")
                )
            }
        )
    })
    
    # choose to select strains for the output
    output$choosestrains <- renderUI({
        # load phenotype data
        phenodf <- loadPhenoData()
        
        # only show the options if this is checked
        if(!is.null(input$straininput) && input$straininput == T) {
            # which strains to show?
            checkboxGroupInput("whichstrains", 
                               "Which strains to include?", 
                               choices = unique(phenodf$strain), 
                               selected = unique(phenodf$strain))

        }
    })
    
    # how many QTL?
    output$number_qtl <- renderUI({
        # only show slider bar if the show QTL checkbox is checked
        if(!is.null(input$showqtl) && input$showqtl == T) {
            # how many qtl?
            numericInput("numQTL", "How many QTL in your model?", 1)
        }
    })
    
    # make a reactive value to keep track of how many sliders we have
    sliders <- reactiveValues(num = 1)

    # show slider bars
    output$show_qtl_pos <- renderUI({
        # only show slider bar if the show QTL checkbox is checked
        if(!is.null(input$showqtl) && input$showqtl == T) {
            # how many qtl in the model?
            qtls <- input$numQTL
            
            # if(sliders$num > 1) {
            #     
            #     
            # } else {
            #     tagList(
            #         # new slider
            #         sliderInput(glue::glue("qtlpos{sliders$num}"),
            #                     glue::glue("Choose position of QTL {sliders$num}"),
            #                     min = 0,
            #                     max = 20,
            #                     value = 10,
            #                     step = 0.1),
            #         # check box for showing genotypes on the phenotye plots
            #         checkboxInput("showgeno", "Show genotype?")
            #     )
            # }
            
            tagList(

                # make slider for each QTL
                lapply(1:qtls, function(i) {
                    tagList(
                        sliderInput(glue::glue("qtlpos{i}"),
                                    glue::glue("Choose position of QTL {i}"),
                                    min = 0,
                                    max = 20,
                                    value = 10,
                                    step = 0.1)
                    )
                }),
                
                # check box for showing genotypes on the phenotye plots
                checkboxInput("showgeno", "Show genotype?")
            )
        }

    })
    
    # initialize reactive values to tell function to look at control, condition, or regressed
    rv <- reactiveValues(cond = NA)
    
    # function to make geno/pheno plot
    plotInput <- reactive({
        # load phenotype data
        phenodf <- loadPhenoData()
        
        # strains from dataframe
        strains <- unique(phenodf$strain)
        if(input$straininput == T) {
            strains <- input$whichstrains
        }
        
        # regress
        pruned <- phenodf %>%
            dplyr::ungroup() %>%
            dplyr::filter(strain %in% strains)
        regressed <- easysorter::regress(pruned) %>%
            dplyr::mutate(condition = paste0(condition, "-regressed")) %>%
            dplyr::bind_rows(pruned)
        
        # only plot phenotype if n2cb is checked
        if(!is.null(input$n2cb) && input$n2cb == TRUE) {
            phenoplot <- quick_plot_breakup_flip2(regressed,
                                             rv$cond,
                                             input$trait) +
                ggplot2::facet_grid(~trait)
            
        } else {
            # plot genotype with input chrom and qtl position as vertical line
            geno <- nil_plot(strains, input$chrom)
            
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
            regressed$strain <- factor(regressed$strain, levels = unique(geno[[2]]$sample),
                                       labels = unique(geno[[2]]$sample))
            pheno <- quick_plot_breakup_flip(regressed,
                                             rv$cond,
                                             input$trait,
                                             geno = input$showgeno,
                                             pos = vals,
                                             chr = input$chrom) +
                ggplot2::facet_grid(~trait)
            
            phenoplot <- cowplot::plot_grid(genoplot, pheno)
        }
        
        # nil stats
        stat <- nil_stats(regressed, rv$cond, input$trait)

        return(list(phenoplot, stat))
    })

    # plot nil phenotypes given user input data (control)
    output$control_plot <- renderUI({
        rv$cond <- input$control
        vars <- plotInput()
        
        # pheno plot
        output$controlpheno <- renderPlot({
            vars[[1]]
        })
        
        # NIL stats dataframe
        output$controlstat <- renderDataTable({
            vars[[2]]
        })
        
        tagList(
            h3("NIL Phenotype"),
            plotOutput("controlpheno"),
            br(),
            h3("NIL stats"),
            dataTableOutput("controlstat")
        )
        
    })
    
    # plot nil phenotypes given user input data (condition)
    output$condition_plot <- renderUI({
        rv$cond <- input$condition
        vars <- plotInput()
        
        # pheno plot
        output$condpheno <- renderPlot({
            vars[[1]]
        })
        
        # NIL stats dataframe
        output$condstat <- renderDataTable({
            vars[[2]]
        })
        
        tagList(
            h3("NIL phenotype"),
            plotOutput("condpheno"),
            br(),
            h3("NIL stats"),
            dataTableOutput("condstat")
        )
    })
    
    # plot nil phenotypes given user input data (regressed)
    output$regressed_plot <- renderUI({
        rv$cond <- paste0(input$condition, "-regressed")
        vars <- plotInput()
        
        # pheno plot
        output$regpheno <- renderPlot({
            vars[[1]]
        })
        
        # NIL stats dataframe
        output$regstat <- renderDataTable({
            vars[[2]]
        })
        
        tagList(
            h3("NIL phenotype"),
            plotOutput("regpheno"),
            br(),
            h3("NIL stats"),
            dataTableOutput("regstat")
        )
    })
    
    # code to generate nil genotype dataset
    nilgeno_dataset <- reactive({
        # load phenotype data
        phenodf <- loadPhenoData()
        
        # strains from dataframe
        strains <- unique(phenodf$strain)
        if(input$straininput == T) {
            strains <- input$whichstrains
        }
        
        # plot genotypes, by chromosome
        nils <- nil_plot(strains, input$chrom)
        return(nils)
    })
    
    # plot nil genotypes for tab
    output$nil_genotypes <- renderUI({
        # return error message if N2/CB checkbox is clicked
        if(!is.null(input$n2cb) && input$n2cb == TRUE) {
            h4(em("Cannot view genotypes for strains that are not N2/CB NILs at this time."))
        } else {
            # call nilgeno_dataset
            nils <- nilgeno_dataset()
            
            output$nilplot <- renderPlot({
                nils[[1]]
            })
            
            
            # show datatable of genotypes
            output$niltable <- renderDataTable({
                nils[[3]] %>%
                    dplyr::select(chrom, start, end, sample, genotype = gt_name)
            })
            
            tagList(
                h3("NIL genotypes - plot"),
                plotOutput("nilplot"),
                br(),
                h3("NIL genotypes - breakpoints"),
                dataTableOutput("niltable"),
                downloadButton('downloadData', "Download data")
            )
        }

    })
    
    output$saveButton <- renderUI({
        req(input$file1)
        downloadButton('saveImage', 'Save plot')
    })
    
    # output figure as png if button is pressed
    savePlot <- eventReactive(input$saveImage, {
        input$n
    })
    
    output$saveImage <- downloadHandler(
        filename = function() { 'test.png' },
        content = function(file) {
            ggsave(file, plot = plotInput()[[1]], device = "png")
        }
    )
    
    # handle download of dataset for nil geno
    output$downloadData <- downloadHandler(
        filename = "nil_genotypes.csv",
        content = function(file) {
            write.csv(nilgeno_dataset()[[3]], file, row.names = FALSE)
        }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)

