library(EnhancedVolcano)
library(tidyverse)
library(DT)
library(shiny)
library(shinyjs)
library(shinythemes) 
library(writexl)
library(plotly)

resNH_shiny <- read.csv("resNH_shiny.csv")
resCH_shiny <- read.csv("resCH_shiny.csv")
resCN_shiny <- read.csv("resCN_shiny.csv")
t_box_shiny <- read.csv("t_box_shiny.csv")
GOBP_CH <- read.csv("GOBP_Genes_CH_SeqCorrected.csv")
GOBP_NH <- read.csv("GOBP_Genes_NH_SeqCorrected.csv")
GOBP_CN <- read.csv("GOBP_Genes_CN_SeqCorrected.csv")
GOBP_CH$Group <- rep("Cirrhosis vs Healthy", length(GOBP_CH$Description))
GOBP_CH$GeneFraction <- sapply(GOBP_CH$GeneRatio, function(x) eval(parse(text=x)))
GOBP_NH$Group <- rep("NAFLD vs Healthy", length(GOBP_NH$Description))
GOBP_NH$GeneFraction <- sapply(GOBP_NH$GeneRatio, function(x) eval(parse(text=x)))
GOBP_CN$Group <- rep("Cirrhosis vs NAFLD", length(GOBP_CN$Description))
GOBP_CN$GeneFraction <- sapply(GOBP_CN$GeneRatio, function(x) eval(parse(text=x)))
t_box_shiny$Group <- factor(t_box_shiny$Group, levels = c("Healthy", "NAFLD", "Cirrhosis"))

ui <- fluidPage(theme = shinytheme("paper"),
            navbarPage(title = "Differential Expression Analysis - PLS",
                       tabPanel(title = "Preface",
                                fluidRow(column(12, wellPanel(tags$h4("Transcriptomic analysis of human liver samples from healthy individuals, NAFLD patients, and cirrhosis patients as presented in:"),
                                                              tags$br(),
                                                              tags$i(tags$h5("A web-based browsable resource of hepatic gene expression in health and liver disease")),
                                                              tags$h6("Josephine Grandt, Christian D. Johansen, Anne Sofie H. Jensen, Mikkel P. Werge, Elias B. Rashu1, Andreas Møller, 
                                                                      Anders E. Junker, Lise Hobolth, Christian Mortensen, Mogens Vyberg, Reza Rafiolsadat Serizawa, Søren Møller, Lise Lotte Gluud, Nicolai J. Wewer Albrechtsen"),
                                                              tags$h6("Paper: X, DOI: X"),
                                                              tags$br(),
                                                              tags$h6("This app contains information on RNA sequencing data from liver samples taken from healthy individuals, NAFLD patients, and cirrhosis patients. The subjects were
                                                              randomized into two subgroups: fasting and postprandial. Individuals in the postpradial subgroup ingested a standardized meal prior to sampling, while the fasting group
                                                              were fasted overnight. When comparing healthy, NAFLD, and cirrhosis, the subgroup difference were accounted for in the bioinformatic analyses in addition to the differences 
                                                              in sequencing runs (Design: ~ Sequencing + Intervention + Groups). For detailed information on study design and methods please refer to the abovementioned publication and the R codes available at", 
                                                                      tags$a("https://github.com/nicwin98/PostprandialLiverStudy"), ".")))),
                                fluidRow(column(12, wellPanel(
                                  tags$h5("The app contains three types of tabs where you can find the following information:"),
                                  tags$ul(
                                    tags$li(tags$strong("Differential Expression Analysis:")),
                                    tags$ul(
                                      tags$li("Results for the following comparisons: NAFLD vs Healthy | Cirrhosis vs Healthy | Cirrhosis vs NAFLD"),
                                      tags$ul(
                                        tags$li("An adjustable volcano plot of all analyzed genes."),
                                        tags$li("A table of all genes displaying the results of the differential expression analysis."),
                                          tags$ul(
                                            tags$li("The information in the table is available for download."))))),
                                  tags$ul(
                                    tags$li(tags$strong("Gene Ontology - Biological Pathways")),
                                    tags$ul(
                                      tags$li("All significant GOBPs are displayed in a table."),
                                      tags$li("Create a dotplot by clicking on your GOBPs of interest."),
                                      tags$ul(
                                        tags$li("The table information and the dotplot are available for download.")))),
                                  tags$ul(
                                    tags$li(tags$strong("Single Gene Expression:")),
                                    tags$ul(
                                      tags$li("Boxplot showing the sample variation for individual genes (for all groups)."),
                                      tags$ul(
                                        tags$li("The boxplot and datapoints are available for download.")))),
                                  tags$br(),
                                  tags$h6("Any inquires related to the application and publication should be directed 
                                                              to the corresponding author, Nicolai J. Wewer Albrechtsen: nicolai.albrechtsen@sund.ku.dk."),
                                  tags$h6("An additional app has been created in relation to the abovementioned publication:", tags$a("https://weweralbrechtsenlab.shinyapps.io/PostprandialLiver/"), 
                                          "."))))),  
                    navbarMenu(title = "Differential Expression Analysis", 
                       tabPanel(title = "NAFLD vs Healthy",
                                fluidRow(column(4,
                                                wellPanel(
                                                  sliderInput(
                                                    inputId = "num1",
                                                    label = "Choose a log2FoldChange cutoff",
                                                    value = 1,
                                                    min = 0.0,
                                                    max = 5.0,
                                                    step = 0.5),
                                                  sliderInput(
                                                    inputId = "xlim1",
                                                    label = "Adjust x-axis",
                                                    value = c(-7,9),
                                                    min = -7,
                                                    max = 9,
                                                    step = 1),
                                                  selectInput(
                                                    inputId = "radio1",
                                                    label = "Choose FDR adjusted p-value:",
                                                    selected = 0.05,
                                                    choices = list("0.05","0.01","0.001","0.0001","0.00001")))),
                                         column(8,plotOutput("volcanoNH"))),
                                fluidRow(column(12, wellPanel(tags$h4("Use the table below to search for genes of interest"),
                                                              tags$h5("You can download the information in the table for selected genes by clicking on your genes of interest and clicking 
                                                                 the download button below the table"),
                                                              tags$h6("The table displays information on:"),
                                                              tags$ul(
                                                                tags$li("ENSEMBL ID: Unique identifier (used to select genes in the Single Gene Expression tab)"),
                                                                tags$li("Gene Symbol"),
                                                                tags$li("Gene Name"),
                                                                tags$li("Base Mean: Mean of normalized expression (not comparable across genes)"),
                                                                tags$li("Log 2 fold changes (green: > 0 | red: < 0)"),
                                                                tags$li("Fold changes (green: > 1 | red: < 1)"),
                                                                tags$li("lfcSE: Standard error of the log 2 fold change"),
                                                                tags$li("Stat: Wald statistic used to calculate the p-values"),
                                                                tags$li("P-value: Not corrected for multiple testing"),
                                                                tags$li("Padj: FDR adjusted p-value (green: padj < 0.05 | black: padj > 0.05)")),
                                                              tags$h6("Sometimes adjusted p-values or p-values and adjusted p-values are missing. This is likely due to:"),
                                                              tags$ul(
                                                                tags$li("Missing padj: The gene is to lowly expressed and thus not tested (see baseMean). 
                                                                   DESeq2's independent filtering function removed it."),
                                                                tags$li("Missing pvalue and padj: The gene has one or more extreme outliers. 
                                                                   As a result the gene is removed from the stastical analysis
                                                                   (see if this is the case by using the Single Gene Expression tab).")),
                                                              tags$h6("Green log 2 fold changes / fold changes means the gene is up-regulated in NAFLD compared
                                                                 to the healthy group.")))),
                                fluidRow(column(12, wellPanel(DTOutput("listNH")))),
                                fluidRow(column(12, wellPanel(actionButton(inputId = "clearRowsDE1",
                                                                           label = "Clear selected rows"),
                                                              downloadButton("printDEselect1", 'Download info on all selected genes'))))),
                  
                  tabPanel(title = "Cirrhosis vs Healthy",
                           fluidRow(column(4,
                                           wellPanel(
                                             sliderInput(
                                               inputId = "num2",
                                               label = "Choose a log2FoldChange cutoff",
                                               value = 1,
                                               min = 0.0,
                                               max = 5.0,
                                               step = 0.5),
                                             sliderInput(
                                               inputId = "xlim2",
                                               label = "Adjust x-axis",
                                               value = c(-10,9),
                                               min = -10,
                                               max = 30,
                                               step = 1),
                                             selectInput(
                                               inputId = "radio2",
                                               label = "Chosse FDR adjusted p-value:",
                                               selected = 0.05,
                                               choices = list("0.05","0.01","0.001","0.0001","0.00001")))),
                                    column(8,plotOutput("volcanoCH"))),
                           fluidRow(column(12, wellPanel(tags$h4("Use the table below to search for genes of interest"),
                                                         tags$h5("You can download the information in the table for selected genes by clicking on your genes of interest and clicking 
                                                                 the download button below the table"),
                                                         tags$h6("The table displays information on:"),
                                                         tags$ul(
                                                           tags$li("ENSEMBL ID: Unique identifier (used to select genes in the Single Gene Expression tab)"),
                                                           tags$li("Gene Symbol"),
                                                           tags$li("Gene Name"),
                                                           tags$li("Base Mean: Mean of normalized expression (not comparable across genes)"),
                                                           tags$li("Log 2 fold changes (green: > 0 | red: < 0)"),
                                                           tags$li("Fold changes (green: > 1 | red: < 1)"),
                                                           tags$li("lfcSE: Standard error of the log 2 fold change"),
                                                           tags$li("Stat: Wald statistic used to calculate the p-values"),
                                                           tags$li("P-value: Not corrected for multiple testing"),
                                                           tags$li("Padj: FDR adjusted p-value (green: padj < 0.05 | black: padj > 0.05)")),
                                                         tags$h6("Sometimes adjusted p-values or p-values and adjusted p-values are missing. This is likely due to:"),
                                                         tags$ul(
                                                           tags$li("Missing padj: The gene is to lowly expressed and thus not tested (see baseMean). 
                                                                   DESeq2's independent filtering function removed it."),
                                                           tags$li("Missing pvalue and padj: The gene has one or more extreme outliers. 
                                                                   As a result the gene is removed from the stastical analysis
                                                                   (see if this is the case by using the Single Gene Expression tab).")),
                                                         tags$h6("Green log 2 fold changes / fold changes means the gene is up-regulated in cirrhosis compared
                                                                 to the healthy group.")))),
                           fluidRow(column(12, wellPanel(DTOutput("listCH")))),
                           fluidRow(column(12, wellPanel(actionButton(inputId = "clearRowsDE2",
                                                                      label = "Clear selected rows"),
                                                         downloadButton("printDEselect2", 'Download info on all selected genes'))))),
                  tabPanel(title = "Cirrhosis vs NAFLD",
                           fluidRow(column(4,
                                           wellPanel(
                                             sliderInput(
                                               inputId = "num3",
                                               label = "Choose a log2FoldChange cutoff",
                                               value = 1,
                                               min = 0.0,
                                               max = 5.0,
                                               step = 0.5),
                                             sliderInput(
                                               inputId = "xlim3",
                                               label = "Adjust x-axis",
                                               value = c(-13,13),
                                               min = -13,
                                               max = 13,
                                               step = 1),
                                             selectInput(
                                               inputId = "radio3",
                                               label = "Chosse FDR adjusted p-value:",
                                               selected = 0.05,
                                               choices = list("0.05","0.01","0.001","0.0001","0.00001")))),
                                    column(8,plotOutput("volcanoCN"))),
                           fluidRow(column(12, wellPanel(tags$h4("Use the table below to search for genes of interest"),
                                                         tags$h5("You can download the information in the table for selected genes by clicking on your genes of interest and clicking 
                                                                 the download button below the table"),
                                                         tags$h6("The table displays information on:"),
                                                         tags$ul(
                                                           tags$li("ENSEMBL ID: Unique identifier (used to select genes in the Single Gene Expression tab)"),
                                                           tags$li("Gene Symbol"),
                                                           tags$li("Gene Name"),
                                                           tags$li("Base Mean: Mean of normalized expression (not comparable across genes)"),
                                                           tags$li("Log 2 fold changes (green: > 0 | red: < 0)"),
                                                           tags$li("Fold changes (green: > 1 | red: < 1)"),
                                                           tags$li("lfcSE: Standard error of the log 2 fold change"),
                                                           tags$li("Stat: Wald statistic used to calculate the p-values"),
                                                           tags$li("P-value: Not corrected for multiple testing"),
                                                           tags$li("Padj: FDR adjusted p-value (green: padj < 0.05 | black: padj > 0.05)")),
                                                         tags$h6("Sometimes adjusted p-values or p-values and adjusted p-values are missing. This is likely due to:"),
                                                         tags$ul(
                                                           tags$li("Missing padj: The gene is to lowly expressed and thus not tested (see baseMean). 
                                                                   DESeq2's independent filtering function removed it."),
                                                           tags$li("Missing pvalue and padj: The gene has one or more extreme outliers. 
                                                                   As a result the gene is removed from the stastical analysis
                                                                   (see if this is the case by using the Single Gene Expression tab).")),
                                                         tags$h6("Green log 2 fold changes / fold changes means the gene is up-regulated in cirrhosis compared
                                                                 to the NAFLD group.")))),
                           fluidRow(column(12, wellPanel(DTOutput("listCN")))),
                           fluidRow(column(12, wellPanel(actionButton(inputId = "clearRowsDE3",
                                                                      label = "Clear selected rows"),
                                                         downloadButton("printDEselect3", 'Download info on all selected genes')))))),
              navbarMenu(title = "Gene Ontology - Biological Pathways",    
                  tabPanel(title = "GOBP - NAFLD vs Healthy",
                           fluidRow(column(12, wellPanel(h5("The table below allows you to click and select the Gene Ontology Biological Pathways (GOBPs) you are interested in. 
                                                            The selected GOBPs will appear in the dotplot below the table."),
                                                         h6("Both the displayed dotplot and the selected table information are available for download. 
                                                            All GOBPs in the table are significantly enriched in the NAFLD group compared to the healthy group")))),
                           fluidRow(column(4,
                                           downloadButton("printGOBPselect_NH", 'Download info on all selected GOBPs'),
                                           downloadButton("printGOBPdotplot_NH", "Download the current dotplot")),
                                    column(4,
                                           numericInput(
                                             inputId = "radioGO_NH1",
                                             label = "Adjust the width of the downloadable dotplot",
                                             value = 8,
                                             min = 1,
                                             max = 15,
                                             width = "200%")),
                                    column(4,
                                           numericInput(
                                             inputId = "radioGO_NH2",
                                             label = "Adjust the height of the downloadable dotplot",
                                             value = 4,
                                             min = 1,
                                             max = 15,
                                             width = "200%"))),
                           fluidRow(column(12, wellPanel(DTOutput("GOBPtable_NH")))),
                           fluidRow(column(12, wellPanel(actionButton(inputId = "clearRowsGOBP_NH",
                                                                      label = "Clear selected rows")))),
                           fluidRow(useShinyjs(),style = "padding-top:20px", column(12, plotOutput("GOBPlot_NH")))),
                  
                  tabPanel(title = "GOBP - Cirrhosis vs Healthy",
                           fluidRow(column(12, wellPanel(h5("The table below allows you to click and select the Gene Ontology Biological Pathways (GOBPs) you are interested in. 
                                                            The selected GOBPs will appear in the dotplot below the table."),
                                                         h6("Both the displayed dotplot and the selected table information are available for download. 
                                                            All GOBPs in the table are significantly enriched in the cirrhosis group compared to the healthy group")))),
                           fluidRow(column(4,
                                           downloadButton("printGOBPselect_CH", 'Download info on all selected GOBPs'),
                                           downloadButton("printGOBPdotplot_CH", "Download the current dotplot")),
                                    column(4,
                                           numericInput(
                                             inputId = "radioGO_CH1",
                                             label = "Adjust the width of the downloadable dotplot",
                                             value = 8,
                                             min = 1,
                                             max = 15,
                                             width = "200%")),
                                    column(4,
                                           numericInput(
                                             inputId = "radioGO_CH2",
                                             label = "Adjust the height of the downloadable dotplot",
                                             value = 4,
                                             min = 1,
                                             max = 15,
                                             width = "200%"))),
                           fluidRow(column(12, wellPanel(DTOutput("GOBPtable_CH")))),
                           fluidRow(column(12, wellPanel(actionButton(inputId = "clearRowsGOBP_CH",
                                                                      label = "Clear selected rows")))),
                           fluidRow(useShinyjs(),style = "padding-top:20px", column(12, plotOutput("GOBPlot_CH")))),
                  
                  tabPanel(title = "GOBP - Cirrhosis vs NAFLD",
                           fluidRow(column(12, wellPanel(h5("The table below allows you to click and select the Gene Ontology Biological Pathways (GOBPs) you are interested in. 
                                                            The selected GOBPs will appear in the dotplot below the table."),
                                                         h6("Both the displayed dotplot and the selected table information are available for download. 
                                                            All GOBPs in the table are significantly enriched in the cirrhosis group compared to the NAFLD group")))),
                           fluidRow(column(4,
                                           downloadButton("printGOBPselect_CN", 'Download info on all selected GOBPs'),
                                           downloadButton("printGOBPdotplot_CN", "Download the current dotplot")),
                                    column(4,
                                           numericInput(
                                             inputId = "radioGO_CN1",
                                             label = "Adjust the width of the downloadable dotplot",
                                             value = 8,
                                             min = 1,
                                             max = 15,
                                             width = "200%")),
                                    column(4,
                                           numericInput(
                                             inputId = "radioGO_CN2",
                                             label = "Adjust the height of the downloadable dotplot",
                                             value = 4,
                                             min = 1,
                                             max = 15,
                                             width = "200%"))),
                           fluidRow(column(12, wellPanel(DTOutput("GOBPtable_CN")))),
                           fluidRow(column(12, wellPanel(actionButton(inputId = "clearRowsGOBP_CN",
                                                                      label = "Clear selected rows")))),
                           fluidRow(useShinyjs(),style = "padding-top:20px", column(12, plotOutput("GOBPlot_CN"))))),
                  
              tabPanel(title = "Single Gene Expression - BoxPlot",
                       fluidRow(column(12, wellPanel(tags$h5("The boxplot below displays expression counts normalized with DESeq2 (plotCount function) across treatment groups"),
                                                     tags$h6("The counts are NOT comparable across genes, as they are not normalized to gene length. The boxplot shows the median, 25", 
                                                             tags$sup("th"),", and 75", tags$sup("th"),"percentiles. Points are displayed as outliers if they are above or below 
                                                                 1.5 times the interquartile range.")))),
                       fluidRow(column(12,
                                       wellPanel(
                                         selectizeInput(
                                           inputId = 'gene',
                                           label = "Select your gene of interest",
                                           choices = NULL),
                                         downloadButton("printboxplot", 'Download plot as a PDF file'),
                                         downloadButton("printdatatable", 'Download expression values from plot')))),
                       fluidRow(column(12, plotlyOutput("expressionboxplot")))),
                  
                  
            )
)

server <- function(input, output, session) {
  
  updateSelectizeInput(session, 'gene', choices = as.vector(unique(t_box_shiny$GeneSymbol)), selected = character(0), server = TRUE)
  data_box <- reactive({filter(t_box_shiny, GeneSymbol == input$gene, .preserve = TRUE)})
  
  reactive_NH <- reactive({resNH_shiny})
  reactive_CH <- reactive({resCH_shiny})
  reactive_CN <- reactive({resCN_shiny})
  
  output$volcanoNH <- renderPlot({
    EnhancedVolcano(reactive_NH(),
                    lab = resNH_shiny$gene_symbol,
                    x = 'log2FoldChange',
                    y = 'padj',
                    title = 'NAFLD vs Healthy',
                    pCutoff = as.numeric(input$radio1),
                    FCcutoff = input$num1,
                    pointSize = 1.0,
                    labSize = 2.0,
                    cutoffLineWidth = 0.5,
                    col=c('darkgray', 'darkgray', 'blue', 'red'),
                    colAlpha = 0.75,
                    legendLabels=c('Not sig.','','Sig. & low log2FC','Sig. & high log2FC'),
                    drawConnectors = TRUE,
                    widthConnectors = 0.5,
                    xlim = as.numeric(input$xlim1),
                    ylim = c(0,7))
  })
  
  output$listNH <- renderDT({datatable(resNH_shiny, 
                      options = list(
                        scrollX = TRUE,
                        pageLength = 10,
                        lengthMenu = c(5,10,25,50,200),
                        filter = "bottom"
                      )) %>%
      formatSignif(., columns = c(9,10), digits = 4) %>%
      formatRound(., columns = c(4:8), digits = 4, interval = 0, mark = "") %>%
      formatStyle('padj', color = styleInterval(0.05,c('green','black'))) %>%
      formatStyle('log2FoldChange', color = styleInterval(0, c('red', 'green'))) %>%
      formatStyle('foldChange', color = styleInterval(1, c('red', 'green')))
  })
  
  DTreset1 <- DT::dataTableProxy("listNH")
  shiny::observeEvent(input$clearRowsDE1, {
    DT::selectRows(DTreset1, NULL)
  })
  
  output$printDEselect1 = downloadHandler('NAFLDvsHealthy_GeneList_Selected.xlsx', content = function(file) 
  {srows_data <- input$listNH_rows_selected
  writexl::write_xlsx(resNH_shiny[srows_data, , drop = FALSE], path = file)
  })
  
  output$volcanoCH <- renderPlot({
    EnhancedVolcano(reactive_CH(),
                    lab = resCH_shiny$gene_symbol,
                    x = 'log2FoldChange',
                    y = 'padj',
                    title = 'Cirrhosis vs Healthy',
                    pCutoff = as.numeric(input$radio2),
                    FCcutoff = input$num2,
                    pointSize = 1.0,
                    labSize = 2.0,
                    cutoffLineWidth = 0.5,
                    col=c('darkgray', 'darkgray', 'blue', 'red'),
                    colAlpha = 0.75,
                    legendLabels=c('Not sig.','','Sig. & low log2FC','Sig. & high log2FC'),
                    drawConnectors = TRUE,
                    widthConnectors = 0.5,
                    xlim = as.numeric(input$xlim2),
                    ylim = c(0,36))
  })
  
  output$listCH <- renderDT({datatable(resCH_shiny, 
                                       options = list(
                                         scrollX = TRUE,
                                         pageLength = 10,
                                         lengthMenu = c(5,10,25,50,200),
                                         filter = "bottom"
                                       )) %>%
      formatSignif(., columns = c(9,10), digits = 4) %>%
      formatRound(., columns = c(4:8), digits = 4, interval = 0, mark = "") %>%
      formatStyle('padj', color = styleInterval(0.05,c('green','black'))) %>%
      formatStyle('log2FoldChange', color = styleInterval(0, c('red', 'green'))) %>%
      formatStyle('foldChange', color = styleInterval(1, c('red', 'green')))
  })
  
  DTreset2 <- DT::dataTableProxy("listCH")
  shiny::observeEvent(input$clearRowsDE2, {
    DT::selectRows(DTreset2, NULL)
  })
  
  output$printDEselect2 = downloadHandler('CirrhosisvsHealthy_GeneList_Selected.xlsx', content = function(file) 
  {srows_data <- input$listCH_rows_selected
  writexl::write_xlsx(resCH_shiny[srows_data, , drop = FALSE], path = file)
  })

  output$volcanoCN <- renderPlot({
    EnhancedVolcano(reactive_CN(),
                    lab = resCN_shiny$gene_symbol,
                    x = 'log2FoldChange',
                    y = 'padj',
                    title = 'Cirrhosis vs NAFLD',
                    pCutoff = as.numeric(input$radio3),
                    FCcutoff = input$num3,
                    pointSize = 1.0,
                    labSize = 2.0,
                    cutoffLineWidth = 0.5,
                    col=c('darkgray', 'darkgray', 'blue', 'red'),
                    colAlpha = 0.75,
                    legendLabels=c('Not sig.','','Sig. & low log2FC','Sig. & high log2FC'),
                    drawConnectors = TRUE,
                    widthConnectors = 0.5,
                    xlim = as.numeric(input$xlim3),
                    ylim = c(0,22.5))
  })
  
  output$listCN <- renderDT({datatable(resCN_shiny, 
                                       options = list(
                                         scrollX = TRUE,
                                         pageLength = 10,
                                         lengthMenu = c(5,10,25,50,200),
                                         filter = "bottom"
                                       )) %>%
      formatSignif(., columns = c(9,10), digits = 4) %>%
      formatRound(., columns = c(4:8), digits = 4, interval = 0, mark = "") %>%
      formatStyle('padj', color = styleInterval(0.05,c('green','black'))) %>%
      formatStyle('log2FoldChange', color = styleInterval(0, c('red', 'green'))) %>%
      formatStyle('foldChange', color = styleInterval(1, c('red', 'green')))
  })
  
  DTreset3 <- DT::dataTableProxy("listCN")
  shiny::observeEvent(input$clearRowsDE3, {
    DT::selectRows(DTreset3, NULL)
  })
  
  output$printDEselect3 = downloadHandler('CirrhosisvsNAFLD_GeneList_Selected.xlsx', content = function(file) 
  {srows_data <- input$listCN_rows_selected
  writexl::write_xlsx(resCN_shiny[srows_data, , drop = FALSE], path = file)
  })
  
  output$expressionboxplot <- renderPlotly({
    validate(
      need(input$gene, "Select a gene above to generate the boxplot.")
    )
    ebplot <- ggplotly(ggplot(data_box(),aes(x=Group, y=count, col = Group)) +
                         geom_boxplot(outlier.shape = NA) +
                         geom_jitter(aes(col = Group), alpha = 0.6, color = "black", width = 0.15) +
                         ylab("Count") +
                         xlab("Group") +
                         ggtitle(paste0(input$gene," human liver RNA expression")) +
                         scale_color_manual(values = c("blue","red","black")) +
                         theme_bw())
    print(ebplot)
  })
 
  
  output$printboxplot <- downloadHandler(
    filename = function() {
      paste0('ExpressionBoxPlot_PLS_' ,input$gene ,"_", Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      ggsave(file, 
             ggplot(data_box(),aes(x=Group, y=count, col = Group)) +
               geom_boxplot(outlier.shape = NA) +
               geom_jitter(aes(col = Group), alpha = 0.6, color = "black", width = 0.15) +
               ylab("Count") +
               xlab("Group") +
               ggtitle(paste0(input$gene," human liver RNA expression")) +
               scale_color_manual(values = c("blue","red","black")) +
               theme_bw(), 
             dpi = 1000, height = 7, width = 7)
    })
  
  output$printdatatable <- downloadHandler(
    filename = function() {paste0('NormCounts_PLS_',input$gene ,"_", Sys.Date(), '.xlsx', sep='')},
    content = function(file) 
    {writexl::write_xlsx(dplyr::filter(data_box()), path = file)
    })

  output$GOBPtable_NH <- renderDT(GOBP_NH,
                               filter = "bottom",
                               options = list(
                                 autoWidth = TRUE,
                                 scrollX = TRUE,
                                 rowCallback = JS(
                                   "function(row, data) {",
                                   "for (i = 1; i < data.length; i++) {",
                                   "if (data[i]<0.01){",
                                   "$('td:eq('+i+')', row).html(data[i].toExponential(4));",
                                   "}",
                                   "}",
                                   "}"),
                                 columnDefs = list(list(width = '250px', targets = c(3)),
                                                   list(width = '100px', targets = c(1)),
                                                   list(width = '55px', targets = c(6,7)),
                                                   list(visible=FALSE, targets=c(9,10))),
                                 pageLength = 10,
                                 lengthMenu = c(5,10,25,50,200))
  )
  
  DTreset_GO_NH <- DT::dataTableProxy("GOBPtable_NH")
  shiny::observeEvent(input$clearRowsGOBP_NH, {
    DT::selectRows(DTreset_GO_NH, NULL)
  })
  
  output$GOBPlot_NH <- renderPlot({
    validate(
      need(input$GOBPtable_NH_rows_selected, "Please click and select your GOBPs of interest to generate the dotplot.")
    )
    srows_data <- input$GOBPtable_NH_rows_selected
    GOBP_NH <- GOBP_NH[srows_data, , drop = FALSE]
    ggplot(GOBP_NH, aes(x = Genes, y = Description)) + 
      geom_point(aes(size = GeneFraction, color = p.adjust), stroke = 2) +
      scale_color_distiller(palette = "Spectral", direction = 1) +
      ylab(NULL) +
      xlab(NULL) +
      facet_wrap(~Group, drop = TRUE) +
      theme_bw(base_size = 14) +
      theme(axis.text = element_text(color = "black", size = 12), 
            legend.key.size = unit(0.7, 'cm'), 
            legend.title = element_text(size=11), 
            legend.text = element_text(size=9)) +
      guides(size = guide_legend(order=1))
  })
  
  output$printGOBPselect_NH = downloadHandler('NAFLDvsHealthy_GOBP_Selected.xlsx', content = function(file) 
  {srows_data <- input$GOBPtable_NH_rows_selected
  writexl::write_xlsx(GOBP_NH[srows_data, , drop = FALSE], path = file)
  })
  
  output$printGOBPdotplot_NH <- downloadHandler(
    filename = function() {
      paste('DotPlot_NAFLDvsHealthy_GOBP', Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      srows_data <- input$GOBPtable_NH_rows_selected
      GOBP_NH <- GOBP_NH[srows_data, , drop = FALSE]
      ggsave(file, 
             ggplot(GOBP_NH, aes(x = Genes, y = Description)) + 
               geom_point(aes(size = GeneFraction, color = p.adjust), stroke = 2) +
               scale_color_distiller(palette = "Spectral", direction = 1) +
               ylab(NULL) +
               xlab(NULL) +
               facet_wrap(~Group, drop = TRUE) +
               theme_bw(base_size = 14) +
               theme(axis.text = element_text(color = "black", size = 12), 
                     legend.key.size = unit(0.7, 'cm'), 
                     legend.title = element_text(size=11), 
                     legend.text = element_text(size=9)) +
               guides(size = guide_legend(order=1)), 
             dpi = 1000, width = as.numeric(input$radioGO_NH1), height = as.numeric(input$radioGO_NH2))
    })
  
  output$GOBPtable_CH <- renderDT(GOBP_CH,
                                  filter = "bottom",
                                  options = list(
                                    autoWidth = TRUE,
                                    scrollX = TRUE,
                                    rowCallback = JS(
                                      "function(row, data) {",
                                      "for (i = 1; i < data.length; i++) {",
                                      "if (data[i]<0.01){",
                                      "$('td:eq('+i+')', row).html(data[i].toExponential(4));",
                                      "}",
                                      "}",
                                      "}"),
                                    columnDefs = list(list(width = '250px', targets = c(3)),
                                                      list(width = '100px', targets = c(1)),
                                                      list(width = '55px', targets = c(6,7)),
                                                      list(visible=FALSE, targets=c(9,10))),
                                    pageLength = 10,
                                    lengthMenu = c(5,10,25,50,200))
  )
  
  DTreset_GO_CH <- DT::dataTableProxy("GOBPtable_CH")
  shiny::observeEvent(input$clearRowsGOBP_CH, {
    DT::selectRows(DTreset_GO_CH, NULL)
  })
  
  output$GOBPlot_CH <- renderPlot({
    validate(
      need(input$GOBPtable_CH_rows_selected, "Please click and select your GOBPs of interest to generate the dotplot.")
    )
    srows_data <- input$GOBPtable_CH_rows_selected
    GOBP_CH <- GOBP_CH[srows_data, , drop = FALSE]
    ggplot(GOBP_CH, aes(x = Genes, y = Description)) + 
      geom_point(aes(size = GeneFraction, color = p.adjust), stroke = 2) +
      scale_color_distiller(palette = "Spectral", direction = 1) +
      ylab(NULL) +
      xlab(NULL) +
      facet_wrap(~Group, drop = TRUE) +
      theme_bw(base_size = 14) +
      theme(axis.text = element_text(color = "black", size = 12), 
            legend.key.size = unit(0.7, 'cm'), 
            legend.title = element_text(size=11), 
            legend.text = element_text(size=9)) +
      guides(size = guide_legend(order=1))
  })
  
  output$printGOBPselect_CH = downloadHandler('Cirrhosis_VS_Healthy_GOBP_Selected.xlsx', content = function(file) 
  {srows_data <- input$GOBPtable_CH_rows_selected
  writexl::write_xlsx(GOBP_CH[srows_data, , drop = FALSE], path = file)
  })
  
  output$printGOBPdotplot_CH <- downloadHandler(
    filename = function() {
      paste('DotPlot_Cirrhosis_VS_Healthy_GOBP', Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      srows_data <- input$GOBPtable_CH_rows_selected
      GOBP_CH <- GOBP_CH[srows_data, , drop = FALSE]
      ggsave(file, 
             ggplot(GOBP_CH, aes(x = Genes, y = Description)) + 
               geom_point(aes(size = GeneFraction, color = p.adjust), stroke = 2) +
               scale_color_distiller(palette = "Spectral", direction = 1) +
               ylab(NULL) +
               xlab(NULL) +
               facet_wrap(~Group, drop = TRUE) +
               theme_bw(base_size = 14) +
               theme(axis.text = element_text(color = "black", size = 12), 
                     legend.key.size = unit(0.7, 'cm'), 
                     legend.title = element_text(size=11), 
                     legend.text = element_text(size=9)) +
               guides(size = guide_legend(order=1)), 
             dpi = 1000, width = as.numeric(input$radioGO_CH1), height = as.numeric(input$radioGO_CH2))
    })
  
  output$GOBPtable_CN <- renderDT(GOBP_CN,
                                  filter = "bottom",
                                  options = list(
                                    autoWidth = TRUE,
                                    scrollX = TRUE,
                                    rowCallback = JS(
                                      "function(row, data) {",
                                      "for (i = 1; i < data.length; i++) {",
                                      "if (data[i]<0.01){",
                                      "$('td:eq('+i+')', row).html(data[i].toExponential(4));",
                                      "}",
                                      "}",
                                      "}"),
                                    columnDefs = list(list(width = '250px', targets = c(3)),
                                                      list(width = '100px', targets = c(1)),
                                                      list(width = '55px', targets = c(6,7)),
                                                      list(visible=FALSE, targets=c(9,10))),
                                    pageLength = 10,
                                    lengthMenu = c(5,10,25,50,200))
  )
  
  DTreset_GO_CN <- DT::dataTableProxy("GOBPtable_CN")
  shiny::observeEvent(input$clearRowsGOBP_CN, {
    DT::selectRows(DTreset_GO_CN, NULL)
  })
  
  output$GOBPlot_CN <- renderPlot({
    validate(
      need(input$GOBPtable_CN_rows_selected, "Please click and select your GOBPs of interest to generate the dotplot.")
    )
    srows_data <- input$GOBPtable_CN_rows_selected
    GOBP_CN <- GOBP_CN[srows_data, , drop = FALSE]
    ggplot(GOBP_CN, aes(x = Genes, y = Description)) + 
      geom_point(aes(size = GeneFraction, color = p.adjust), stroke = 2) +
      scale_color_distiller(palette = "Spectral", direction = 1) +
      ylab(NULL) +
      xlab(NULL) +
      facet_wrap(~Group, drop = TRUE) +
      theme_bw(base_size = 14) +
      theme(axis.text = element_text(color = "black", size = 12), 
            legend.key.size = unit(0.7, 'cm'), 
            legend.title = element_text(size=11), 
            legend.text = element_text(size=9)) +
      guides(size = guide_legend(order=1))
  })
  
  output$printGOBPselect_CN = downloadHandler('Cirrhosis_VS_NAFLD_GOBP_Selected.xlsx', content = function(file) 
  {srows_data <- input$GOBPtable_CN_rows_selected
  writexl::write_xlsx(GOBP_CN[srows_data, , drop = FALSE], path = file)
  })
  
  output$printGOBPdotplot_CN <- downloadHandler(
    filename = function() {
      paste('DotPlot_Cirrhosis_VS_NAFLD_GOBP', Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      srows_data <- input$GOBPtable_CN_rows_selected
      GOBP_CN <- GOBP_CN[srows_data, , drop = FALSE]
      ggsave(file, 
             ggplot(GOBP_CN, aes(x = Genes, y = Description)) + 
               geom_point(aes(size = GeneFraction, color = p.adjust), stroke = 2) +
               scale_color_distiller(palette = "Spectral", direction = 1) +
               ylab(NULL) +
               xlab(NULL) +
               facet_wrap(~Group, drop = TRUE) +
               theme_bw(base_size = 14) +
               theme(axis.text = element_text(color = "black", size = 12), 
                     legend.key.size = unit(0.7, 'cm'), 
                     legend.title = element_text(size=11), 
                     legend.text = element_text(size=9)) +
               guides(size = guide_legend(order=1)), 
             dpi = 1000, width = as.numeric(input$radioGO_CN1), height = as.numeric(input$radioGO_CN2))
    })
  
}

shinyApp(ui, server)