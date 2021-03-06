library(shiny)
library(tidyverse)
library(easysorter)
library(shinythemes)
library(cowplot)
library(DT)
library(plotly)

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
           shiny::checkboxInput("sampledata", "Use sample data"),
           
           # Input: Select a file
           shiny::fileInput("file1", "Choose phenotype file (Rdata)",
                     accept = c(".rda", ".Rda", ".RData", ".csv")),
           
           # check if your data is not N2/CB NILs, won't show genotypes
           shiny::checkboxInput("n2cb", "My data is not N2/CB NILs"),
           
           # show the options to choose trait and conditions after the user has uploaded a file
           shiny::uiOutput("choosecontrol"),
           
           # choose to select certain strains
           shiny::uiOutput("choosestrains"),
           
           # show the options to choose QTL and genotype
           shiny::uiOutput("chooseqtl")
           
       ),
       
       mainPanel(
           # make tabs for each condition
           shiny::tabsetPanel(type = "tabs",
                              shiny::tabPanel("Control", shiny::uiOutput("control_plot")),
                              shiny::tabPanel("Condition", shiny::uiOutput("condition_plot")),
                              shiny::tabPanel("Regressed", shiny::uiOutput("regressed_plot")),
                              shiny::tabPanel("NIL Genotypes", shiny::uiOutput("nil_genotypes")),
                              shiny::tabPanel("Help!", shiny::uiOutput("help_page"))),
           
           # # button to output plot as png
           shiny::uiOutput("saveButton")
       )

   )
)


#########################################
#               Server                  #
#########################################


# Define server logic required for application
server <- function(input, output) {
    
    # load phenotype data
    loadPhenoData <- shiny::reactive({
        if(!is.null(input$sampledata) && input$sampledata == T) {
            assign('phenodf', get(load("data/test_NIL_pheno.Rda")))
        } else {
            req(input$file1)
            # If CSV File
            if(grepl(".csv", input$file1$datapath)) {
                phenodf <- read.csv(input$file1$datapath) %>%
                    dplyr::mutate(condition = as.character(condition))
            } else { # otherwise, RDA file
                assign('phenodf', get(load(input$file1$datapath)))
            }
        }
    })
    
    # intro text
    output$intro <- shiny::renderUI({
        p("This shiny app can be used to view NIL phenotypes and genotypes and estimate the location of a QTL. To begin,
        click the 'Browse' button to open a R data file containing pruned (non-regressed) NIL phenotypes or check the
        'Use sample data' checkbox to explore our features. For more instructions, refer to the 'Help!' tab below.")
        
    })
    
    # help page
    output$help_page <- shiny::renderUI({
        
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
    output$choosecontrol <- shiny::renderUI({
        # load phenotype data
        phenodf <- loadPhenoData()
        
        # choose control: if dmso or water exists, choose this first.
        if(sum(c("DMSO", "water", "None", "Water") %in% unique(phenodf$condition)) > 0) {
            control_choices <- intersect(c("DMSO", "water", "None", "Water"), unique(phenodf$condition))
        } else {
            control_choices <- unique(phenodf$condition)
        }
        
        tagList(
           
            # choose control
            selectInput("control", "Select control:", choices = control_choices),
            
            # show output for choosetrait
            uiOutput("choosetrait")
        )
        
    })
    
    # give user options for choosing condition and trait once phenotype dataframe is uploaded
    output$choosetrait <- shiny::renderUI({
        # load phenotype data
        phenodf <- loadPhenoData()
        
        # all conditions, including control
        condition_choices <- unique(phenodf$condition)
        
        tagList(
            
            # choose condition - don't allow user to choose control
            selectInput("condition", "Select condition:", 
                        choices = condition_choices[condition_choices != input$control]),
            
            # choose experiment if there are multiple...
            selectInput("exp", "Select assay:", choices = unique(phenodf$experiment)),
            
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
    output$choosestrains <- shiny::renderUI({
        # load phenotype data
        phenodf <- loadPhenoData()
        
        # select experiment
        phenodf <- phenodf %>%
            dplyr::filter(experiment == input$exp)
        
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
    output$number_qtl <- shiny::renderUI({
        # only show slider bar if the show QTL checkbox is checked
        if(!is.null(input$showqtl) && input$showqtl == T) {
            # how many qtl?
            numericInput("numQTL", "How many QTL in your model?", 1)
        }
    })
    
    # make a reactive value to keep track of how many sliders we have
    sliders <- shiny::reactiveValues(num = 1)

    # show slider bars
    output$show_qtl_pos <- shiny::renderUI({
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
    rv <- shiny::reactiveValues(cond = NA)
    
    # function to make geno/pheno plot
    # plotInput <- shiny::eventReactive(input$go, {
    plotInput <- shiny::reactive({

        # load phenotype data
        phenodf <- loadPhenoData()
        
        # select experiment
        phenodf <- phenodf %>%
            dplyr::filter(experiment == input$exp)
        
        # strains from dataframe
        strains <- unique(phenodf$strain)
        if(input$straininput == T) {
            strains <- input$whichstrains
        }
        
        # regress
        pruned <- phenodf %>%
            dplyr::ungroup() %>%
            dplyr::filter(strain %in% strains)
        
        # assay regression
        if(length(unique(pruned$assay)) > 1) {
            assreg <- easysorter::regress(pruned, assay = TRUE)
        } else {
            assreg <- pruned
        }
        
        regressed <- easysorter::regress(assreg) %>%
            dplyr::ungroup() %>%
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
            
            phenoplot <- plotly::subplot(plotly::ggplotly(genoplot, tooltip = "text"), 
                                         plotly::ggplotly(pheno +
                                                              aes(text = glue::glue("Strain: {strain}\n Rep: p{plate}_{row}{col} \n Assay: {assay} \n Pheno: {round(phenotype, digits = 3)}")), 
                                         tooltip = "text"))
        }
        
        # nil stats
        stat <- nil_stats(regressed, rv$cond, input$trait)

        return(list(phenoplot, stat, cowplot::plot_grid(genoplot, pheno, nrow = 1)))
    })

    # plot nil phenotypes given user input data (control)
    output$control_plot <- shiny::renderUI({
        rv$cond <- input$control
        
        vars <- plotInput()
        
        # try with plotly
        output$controlpheno <- plotly::renderPlotly({
            vars[[1]]
        })
        
        # NIL stats dataframe
        output$controlstat <- DT::renderDataTable({
            vars[[2]]
        })
        
        tagList(
            h3("NIL Phenotype"),
            plotly::plotlyOutput("controlpheno"),
            br(),
            h3("NIL stats"),
            DT::dataTableOutput("controlstat")
        )
        
    })
    
    # plot nil phenotypes given user input data (condition)
    output$condition_plot <- shiny::renderUI({
        rv$cond <- input$condition
        vars <- plotInput()
        
        # pheno plot
        output$condpheno <- plotly::renderPlotly({
            vars[[1]]
        })
        
        # NIL stats dataframe
        output$condstat <- DT::renderDataTable({
            vars[[2]]
        })
        
        tagList(
            h3("NIL phenotype"),
            plotly::plotlyOutput("condpheno"),
            br(),
            h3("NIL stats"),
            DT::dataTableOutput("condstat")
        )
    })
    
    # plot nil phenotypes given user input data (regressed)
    output$regressed_plot <- shiny::renderUI({
        rv$cond <- paste0(input$condition, "-regressed")
        vars <- plotInput()
        
        # pheno plot
        output$regpheno <- plotly::renderPlotly({
            vars[[1]]
        })
        
        # NIL stats dataframe
        output$regstat <- DT::renderDataTable({
            vars[[2]]
        })
        
        tagList(
            h3("NIL phenotype"),
            plotly::plotlyOutput("regpheno"),
            br(),
            h3("NIL stats"),
            DT::dataTableOutput("regstat")
        )
    })
    
    # code to generate nil genotype dataset
    nilgeno_dataset <- shiny::reactive({
        # load phenotype data
        phenodf <- loadPhenoData()
        
        # strains from dataframe
        strains <- unique(phenodf$strain)
        if(input$straininput == T) {
            strains <- input$whichstrains
        }
        
        # plot genotypes, by chromosome
        nils <- nil_plot(strains, input$chrom, all.chr = T)
        return(nils)
    })
        
    # plot nil genotypes for tab
    output$nil_genotypes <- shiny::renderUI({
        # return error message if N2/CB checkbox is clicked
        if(!is.null(input$n2cb) && input$n2cb == TRUE) {
            h4(em("Cannot view genotypes for strains that are not N2/CB NILs at this time."))
        } else {
            # call nilgeno_dataset
            nils <- nilgeno_dataset()
            
            output$nilplot <- plotly::renderPlotly({
                plotly::ggplotly(nils[[1]], tooltip = "text")
            })
            
            
            # show datatable of genotypes
            output$niltable <- DT::renderDataTable({
                nils[[3]] %>%
                    dplyr::select(chrom, start, end, sample, genotype = gt_name)
            })
            
            tagList(
                h3("NIL genotypes - plot"),
                plotly::plotlyOutput("nilplot"),
                br(),
                h3("NIL genotypes - breakpoints"),
                DT::dataTableOutput("niltable"),
                downloadButton('downloadData', "Download data")
            )
        }

    })
    
    output$saveButton <- shiny::renderUI({
        req(input$file1)
        downloadButton('saveImage', 'Save plot')
    })
    
    # output figure as png if button is pressed
    savePlot <- shiny::eventReactive(input$saveImage, {
        input$n
    })
    
    output$saveImage <- shiny::downloadHandler(
        filename = function() { 'test.png' },
        content = function(file) {
            ggsave(file, plot = plotInput()[[3]], device = "png", height = 5, width = 10)
        }
    )
    
    # handle download of dataset for nil geno
    output$downloadData <- shiny::downloadHandler(
        filename = "nil_genotypes.csv",
        content = function(file) {
            write.csv(nilgeno_dataset()[[3]], file, row.names = FALSE)
        }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)

