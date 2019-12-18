library(shiny)
library(tidyverse)
library(easysorter)
library(shinythemes)
library(cowplot)

#########################################
#       Load data & functions           #
#########################################

source("functions.R")

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
                       tabPanel("PCA", 
                                #selectInput("pca_trait", "Select PC:", choices = c("PC1", "PC2", "PC3", "PC4", "PC5")),
                                uiOutput("pca_trait"),
                                uiOutput("pca_plot")),
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
    
    # for plotting with PCA step
    pcaInput <- reactive({
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
        regressed <- easysorter::regress(pruned)
        
        # calculate PCA
        pca_regressed <- calc_pc2(regressed)
        regressed <- pca_regressed[[2]]
        
        # show trait loadings of PCA
        loadings <- data.frame(unclass(pca_regressed[[1]]$loadings)) %>%
            dplyr::mutate(trait = rownames(.)) %>% 
            tidyr::gather(comp, cor, -trait) %>%
            dplyr::mutate(pc = as.numeric(stringr::str_split_fixed(comp, "\\.",2)[,2]),
                          comp = as.character(comp)) %>%
            dplyr::arrange(pc, trait)
        
        loadings$comp <- factor(loadings$comp, levels = unique(loadings$comp))
        
        # how to best show the loadings? heatmap?
        loadingplot <- ggplot(loadings) +
            aes(x = comp, y = trait, fill = cor) +
            geom_tile() +
            scale_fill_gradient2() +
            theme(axis.text.x = element_text(angle = 90))
        
        
        # choose trait to plot
        output$pca_trait <- renderUI({
            selectInput("pca_trait", "Select PC", unique(regressed$trait))
        })
        
        # only plot phenotype if n2cb is checked
        if(!is.null(input$n2cb) && input$n2cb == TRUE) {
            phenoplot <- quick_plot_breakup_flip2(regressed,
                                                  rv$cond,
                                                  input$pca_trait) + # HOW TO CHANGE TO MULTIPLE PCs
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
                                             input$pca_trait,
                                             geno = input$showgeno,
                                             pos = vals,
                                             chr = input$chrom) +
                ggplot2::facet_grid(~trait)
            
            phenoplot <- cowplot::plot_grid(genoplot, pheno)
        }
        
        # nil stats
        stat <- nil_stats(regressed, rv$cond, input$pca_trait)
        
        return(list(phenoplot, stat))
    })

    # make PCA tab
    output$pca_plot <- renderUI({
        
        rv$cond <- input$condition
        vars <- pcaInput()
        
        # pheno plot
        output$pcapheno <- renderPlot({
            vars[[1]]
        })
        
        # NIL stats dataframe
        output$pcastat <- renderDataTable({
            vars[[2]]
        })
        
        tagList(
            h3("NIL phenotype"),
            plotOutput("pcapheno"),
            br(),
            h3("NIL stats"),
            dataTableOutput("pcastat")
        )

        
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

