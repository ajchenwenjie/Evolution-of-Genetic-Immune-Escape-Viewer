# Source helper functions
rm(list = ls())

options(warn = -1)

# Define file paths
base_path <- getwd()
data_path <- file.path(base_path, "data")

# Source helper script
source(file.path(base_path, "helpers.R"), local = TRUE)

# UI code
ui <- tagList(
  navbarPage("Evolution of Genetic Immune Escape Viewer",
             tabPanel("Introduction",
                      fluidRow(
                        column(10, offset = 1,
                               h2("Welcome to the Evolution of Genetic Immune Escape Viewer"),
                               p("This Shiny application explores the landscape of immune escape evolution in cancer by integrating data from CRISPR screens, mutation frequencies, and mutation timing across multiple cancer types from the WGS projects."),
                               tags$div(
                                 img(src = "immune_escape_intro.png", height = "400px", style = "display: block; margin-left: auto; margin-right: auto;"),
                                 style = "margin-top: 20px; margin-bottom: 20px;"
                               ),
                               p(tags$b("Identifying Immunomodulatory Genes:")),
                               tags$ul(
                                 tags$li("1. Published CRISPR screens were collected for identifying the regulators of MHC-I expression, response to CD8 T-cell-mediated killing, NK-cell-mediated killing, macrophage-mediated phagocytosis and γδ T-cell-mediated killing."),
                                 tags$li("2. The top 100 positive/negative regulators in each study were respectively unitized for Gene Set Enrichment Analysis."),
                                 tags$li("3. Next, we selected the featured enrichments for further analysis based on the frequencies of studies reporting these enrichments."),
                                 tags$li("4. We then combined the immunomodulatory genes in each frequently enriched pathway as a gene set for that pathway. "),
                               ),
                               p(tags$b("Inferring Mutation Timing:"), "GRITIC is used to estimate when clonal copy number gains happened in a tumor's evolution. (PMID: 38943574)"),
                               tags$ul(
                                 tags$li("1. Posterior gain timing distributions for clonal copy number segments are calculated based on the copy number, tumor purity and the read counts for SNVs in the region of the gain."),
                                 tags$li("2. For each SNV, GRITIC samples when the gain happened, how many copies the SNV has, and then estimates when the SNV occurred."),
                                 tags$li("3. The exact timing of each SNV is sampled within a time window defined by nearby gains. "),
                                 tags$li("4. The timing of SNVs is measured on a “mutation time” scale that goes from 0 (representing conception) to 1 (the end of clonal evolution). Each SNV is sampled 250 times to create a full timing distribution."),
                               ),
                               tags$div(
                                 img(src = "Timing.jpg", height = "200px", width = "750px", style = "display: block; margin-left: auto; margin-right: auto;"),
                                 style = "margin-top: 20px; margin-bottom: 20px;"
                               ),
                               p("Please select the cohort that you are interested in, then you can navigate through the panels to:", style = "margin-top: 20px;"),
                               tags$ul(
                                 tags$li("View the summary of the cohort."),
                                 tags$li("Analyze mutation frequency and timing of immune-related genes."),
                                 tags$li("Upload your own study for immune escape evolution.")
                               ),
                               p("Please use the tabs above to begin exploring."),
                               p(
                                 "For more details, refer to our BioRxiv preprint: https://www.biorxiv.org/content/10.1101/2025.01.17.632799v1",
                                 tags$a(href = "https://github.com/ajchenwenjie/", 
                                        target = "_blank", "GitHub repository"),
                                 "."
                               )
                        )
                      )
             ),
             tabPanel(
               "Immunomodulatory Pathways",
               sidebarLayout(
                 sidebarPanel(
                   width = 2,
                   selectInput(
                     "cell_type",
                     "Immunity Scenarios",
                     choices = names(cell_type_files),
                     selected = "MHC-I Regulators"
                   ),
                   selectInput(
                     "immuno_regulators",
                     "Select Regulator Type",
                     choices = c("Positive" = "pos", "Negative" = "neg"),
                     selected = "pos"
                   ),
                   numericInput("immuno_n", "Number of Top Pathways", value = 10, min = 5, max = 50),
                   helpText("This sets how many top-ranked pathways (based on frequency) will be shown in the plot.")
                 ),
                 mainPanel(
                   width = 10,
                   uiOutput("immuno_pathway_plot_ui"),
                   br(),
                   uiOutput("pathway_interpretation")
                 )
               )
             ),
             tabPanel(
               "Select the WGS Datasetss",
               sidebarLayout(
                 sidebarPanel(
                   width = 2,
                   selectInput("cohort_select", "Select Cohort", choices = c("PCAWG", "POG570")),
                   actionButton("submit_btn", "Submit")
                 ),
                 mainPanel(
                   conditionalPanel(
                     condition = "output.loading == true",
                     h4("Loading data.....", style = "color: #FF5733; text-align:center;")
                   ),
                   tabPanel(
                     "Select the WGS Datasetss",
                     fluidRow(
                       column(width = 9, uiOutput("cohort_summary_left")),
                       column(width = 3, uiOutput("cohort_summary_right")),
                       br(),
                       uiOutput("cohortsummary_interpretation")
                     )
                   )
                 )
               )
             ),
             tabPanel(
               "Immune Escape Evolution",
               sidebarLayout(
                 sidebarPanel(
                   width = 2,
                   # Cell type selector for specific tabs only
                   conditionalPanel(
                     condition = "input.crispr_tab_selected == 'Single Pathway' || input.crispr_tab_selected == 'Multi Pathways' || input.crispr_tab_selected == 'Survival' || input.crispr_tab_selected =='Immune Cell Infiltration'", 
                     uiOutput("cell_type_ui")
                   ),
                   
                   # Timeline specific controls
                   conditionalPanel(
                     condition = "input.crispr_tab_selected == 'Timeline'", 
                     uiOutput("timing_type_ui"), 
                     uiOutput("cell_type_timeline_ui")
                   ),
                   
                   # Single Pathway specific controls
                   conditionalPanel(
                     condition = "input.crispr_tab_selected == 'Single Pathway'", 
                     uiOutput("regulator_selector"), 
                     uiOutput("pathway_selector")
                   ),
                   
                   # Survival specific controls (no cell_type_ui here)
                   conditionalPanel(
                     condition = "input.crispr_tab_selected == 'Survival'", 
                     uiOutput("pathway_selector_surv"), 
                     uiOutput("select_histology_ui")
                   ),
                   
                   # Immune Cell Infiltration specific controls (no cell_type_ui here)
                   conditionalPanel(
                     condition = "input.crispr_tab_selected == 'Immune Cell Infiltration'", 
                     uiOutput("pathway_selector_ciber"), 
                     uiOutput("select_histology_ciber")
                   )
                 ),
                 mainPanel(
                   tabsetPanel(
                     id = "crispr_tab_selected",
                     tabPanel("Single Pathway",
                              h4("Top 10 Cancer Types by Mutation Frequency"),
                              withSpinner(plotOutput("top_cancer_plot", width = "400px", height = "400px")),
                              br(),
                              uiOutput("singlepathway_interpretation")
                     ),
                     tabPanel("Multi Pathways",
                              fluidPage(
                                withSpinner(plotOutput("radarplot", width = "800px", height = "500px")),
                                br(),
                                uiOutput("multipathways_interpretation")
                              )
                     ),
                     tabPanel("Timeline",
                              div(style = "overflow-x: auto; text-align: left;",
                                  withSpinner(plotOutput("timing_plot", width = "2000px", height = "500px"))
                              ),
                              br(),
                              uiOutput("timeline_interpretation")
                     ),
                     tabPanel("Survival",
                              withSpinner(plotOutput("timing_survival_plot", height = "500px", width = "450px")),
                              br(),
                              uiOutput("timingsurvival_interpretation")
                     ),
                     tabPanel("Immune Cell Infiltration",
                              withSpinner(plotOutput("ciber_plot", height = "400px", width = "800px")),
                              br(),
                              uiOutput("cibersort_interpretation")
                     )
                   )
                 )
               )
             ),
             tabPanel(
               "Exploring Your Genes",
               sidebarLayout(
                 sidebarPanel(
                   width = 2,
                   checkboxInput("useGeneList", "Use Custom Gene List", value = FALSE),
                   conditionalPanel(condition = "input.useGeneList == true", textAreaInput("geneList", "Enter Gene Symbols (comma or newline separated)", "")),
                   conditionalPanel(condition = "input.gene_tab_selected == 'Survival'", uiOutput("select_histology_list")),
                   conditionalPanel(condition = "input.gene_tab_selected == 'ICB Response'", selectInput("therapy_select_genes", "Select Therapy:", choices = NULL)),
                   conditionalPanel(condition = "input.gene_tab_selected == 'ICB Survival'", tagList(selectInput("therapy_surv_genes", "Select Therapy:", choices = NULL)))
                 ),
                 mainPanel(
                   tabsetPanel(
                     id = "gene_tab_selected",
                     tabPanel("Mutation Frequency", withSpinner(uiOutput("oncoplot_ui")), br(), uiOutput("mutation_freq_interpretation")),
                     tabPanel("Gene Timing", withSpinner(uiOutput("timing_plot_singlegene_ui")), br(), uiOutput("gene_timing_interpretation")),
                     tabPanel("Pathway Timing", withSpinner(plotOutput("timing_plot_gene", width = "500px", height = "500px")), br(), uiOutput("pathway_timing_gene_interpretation")),
                     tabPanel("Survival", withSpinner(plotOutput("timing_survival_plot_genelist", height = "500px", width = "450px")), br(), uiOutput("timingsurvival_gene_interpretation")),
                     tabPanel("ICB Response",
                              div(
                                style = "margin-top: 0px; padding-top: 0px; display: flex; justify-content: left; align-items: flex-start;",
                                withSpinner(plotOutput("forestplot", width = "450px", height = "600px"))
                              ),
                              br(),
                              uiOutput("ICBresponse_interpretation")
                     ),
                     tabPanel("ICB Survival",
                              withSpinner(uiOutput("ICB_survival_plot")),
                              br(),
                              uiOutput("ICBsurvival_interpretation")
                     )
                   ),
                   br(),
                   uiOutput("tab_description")
                 )
               )
             ),
             tabPanel(
               "Exploring Your Study",
               fluidPage(
                 sidebarLayout(
                   sidebarPanel(
                     width = 2,
                     fileInput("file", "Upload Gene Summary File", accept = ".txt"),
                     helpText("Example data can be chosen from ./data/Example_data."),
                     numericInput("top_n", "Number of Top Genes:", value = 20, min = 10, step = 10),
                     conditionalPanel(condition = "input.study_tab_selected == 'Survival'", uiOutput("select_histology_study")),
                     conditionalPanel(condition = "input.study_tab_selected == 'ICB Response'", selectInput("therapy_select_study", "Select Therapy:", choices = NULL, selected = "All")),
                     conditionalPanel(condition = "input.study_tab_selected == 'ICB Survival'", selectInput("therapy_surv_study", "Select Therapy:", choices = NULL, selected = "All")),
                     conditionalPanel(condition = "input.study_tab_selected == 'ICB Survival'", selectInput("regulator_type", "Select Regulator Type", choices = c("Positive", "Negative"), selected = "Positive"))
                   ),
                   mainPanel(
                     tabsetPanel(
                       id = "study_tab_selected",
                       tabPanel("CRISPR Ranking", withSpinner(plotOutput("crisprPlot", height = "350px", width = "800px")), br(), uiOutput("CRISPRRanking_interpretation")),
                       tabPanel("Mutation Frequency",
                                fluidRow(
                                  column(6, h4("Top Positive Regulators"), withSpinner(uiOutput("oncoplot_pos_ui"))),
                                  column(6, h4("Top Negative Regulators"), withSpinner(uiOutput("oncoplot_neg_ui")))
                                ),
                                br(),
                                uiOutput("mutation_freq_study_interpretation")
                       ),
                       tabPanel("Gene Timing",
                                fluidRow(
                                  column(6, h4("Top Positive Regulators"), withSpinner(uiOutput("timing_plot_singlegene_pos_ui"))),
                                  column(6, h4("Top Negative Regulators"), withSpinner(uiOutput("timing_plot_singlegene_neg_ui")))
                                ),
                                br(),
                                uiOutput("gene_timing_study_interpretation")
                       ),
                       tabPanel("Pathway Timing",
                                fluidRow(
                                  column(6, h4("Top Positive Regulators"), withSpinner(plotOutput("timing_plot_pos", height = "500px"))),
                                  column(6, h4("Top Negative Regulators"), withSpinner(plotOutput("timing_plot_neg", height = "500px")))
                                ),
                                br(),
                                uiOutput("pathway_timing_study_interpretation")
                       ),
                       tabPanel("Survival",
                                fluidPage(
                                  fluidRow(
                                    column(6, h4("Top Positive Regulators"), withSpinner(plotOutput("timing_survival_plot_pos", height = "500px", width = "100%"))),
                                    column(6, h4("Top Negative Regulators"), withSpinner(plotOutput("timing_survival_plot_neg", height = "500px", width = "100%")))
                                  )
                                ),
                                br(),
                                uiOutput("timingsurvival_study_interpretation")
                       ),
                       tabPanel("ICB Response",
                                fluidRow(
                                  column(6, h4("Top Positive Regulators", style = "margin-bottom: 5px;"), withSpinner(plotOutput("forestplot_pos", height = "600px"))),
                                  column(6, h4("Top Negative Regulators", style = "margin-bottom: 5px;"), withSpinner(plotOutput("forestplot_neg", height = "600px")))
                                ),
                                br(),
                                uiOutput("ICBresponse_study_interpretation")
                       ),
                       tabPanel("ICB Survival",
                                fluidPage(
                                  fluidRow(column(12, withSpinner(uiOutput("ICB_survival_plot_study"))))
                                ),
                                br(),
                                uiOutput("ICBsurvival_study_interpretation")
                       )
                     )
                   )
                 )
               )
             ),
             tabPanel("Help",
                      fluidPage(
                        h2("User Guide"),
                        
                        h4("Overview"),
                        p("The Evolution of Genetic Immune Escape Viewer is an interactive Shiny web application designed to explore the landscape of the evolution of genetic immune escape mechanisms across cancer types."),
                        p("It integrates CRISPR screens, mutation frequency, timing information, and immune therapy response from multiple cohorts such as PCAWG, TCGA, and ICB studies."),
                        
                        h4("Installation and Launch"),
                        tags$ol(
                          tags$li("Ensure all required R packages are installed."),
                          tags$li("Set working directory and source helper functions."),
                          tags$li("Run `shinyApp(ui, server)`.")
                        ),
                        
                        h4("Required R Packages"),
                        p("broom, forestplot, shiny, ggplot2, dplyr, tidyr, data.table, readxl, stringr, ggradar, scales, gridExtra, ggsci, patchwork, maftools, ggdist, ggthemes, ggrepel, plotly, textshape, survival, survminer, purrr, shinycssloaders"),
                        
                        h3("Tabs Description"),
                        
                        h4("1. Introduction"),
                        tags$ul(
                          tags$li("Overview of the application"),
                          tags$li("Visual introduction image"),
                          tags$li("Summary of available analysis functionalities")
                        ),
                        
                        h4("2. Immunomodulatory Pathways"),
                        p("Explore top immunomodulatory pathways identified via CRISPR screens across various immune cell scenarios."),
                        tags$b("Controls:"),
                        tags$ul(
                          tags$li("Immunity Scenarios: Choose immune cell category (e.g., MHC-I regulators, T-cell Killing)"),
                          tags$li("Select Regulator Type: Choose between positive or negative regulators"),
                          tags$li("Number of Top Pathways: Set number of top pathways to display")
                        ),
                        tags$b("Output:"),
                        tags$ul(
                          tags$li("Interactive dot plot showing pathway enrichment across studies with tooltip support")
                        ),
                        
                        h4("3. Select the WGS Datasetss"),
                        tags$b("Controls:"),
                        tags$ul(tags$li("Cohort Selection: PCAWG, POG570, TCGA-OV, CRC-Prognosis, TNBC, Glioma")),
                        tags$b("Output:"),
                        tags$ul(
                          tags$li("Bar plot of sample counts by histology"),
                          tags$li("Boxplot of TMB"),
                          tags$li("Pie charts of sample type and gender")
                        ),
                        
                        h4("4. Immune Escape Evolution"),
                        tags$p("Explore mutation frequency, timing, and clinical associations across cancer types."),
                        tags$b("Sub-tabs:"),
                        tags$ul(
                          tags$li("Single Pathway: Top 10 cancer types with most frequent mutations in the selected pathway"),
                          tags$li("Multi Pathways: Radar plots comparing mutation timing"),
                          tags$li("Timeline: Mutation timing landscape across cohorts"),
                          tags$li("Survival: KM plots stratified by mutation timing"),
                          tags$li("Immune Cell Infiltration: Boxplots of CIBERSORT scores grouped by timing")
                        ),
                        
                        h4("5. Exploring Your Genes"),
                        tags$b("Controls:"),
                        tags$ul(
                          tags$li("Use Custom Gene List: Paste gene symbols"),
                          tags$li("Therapy Selection: For ICB-related analyses")
                        ),
                        tags$b("Sub-tabs:"),
                        tags$ul(
                          tags$li("Mutation Frequency: Oncoplots of selected gene list"),
                          tags$li("Gene Timing: Raincloud plots of timing distribution"),
                          tags$li("Pathway Timing: Bar plots of mean mutation timing"),
                          tags$li("Survival: KM plots for timing-based groups"),
                          tags$li("ICB Response: Forest plots by ICB response"),
                          tags$li("ICB Survival: KM plots for ICB survival")
                        ),
                        
                        h4("6. Exploring Your Study"),
                        tags$b("Controls:"),
                        tags$ul(
                          tags$li("Upload .txt File"),
                          tags$li("Top N Genes Slider"),
                          tags$li("Regulator Type Selection"),
                          tags$li("Therapy and Cancer Type Dropdowns")
                        ),
                        tags$b("Sub-tabs:"),
                        tags$ul(
                          tags$li("CRISPR Ranking: Plots of positive/negative CRISPR ranks"),
                          tags$li("Mutation Frequency: Oncoplots for top genes"),
                          tags$li("Gene Timing: Raincloud plots for CRISPR-ranked regulators"),
                          tags$li("Pathway Timing: Bar plots of mutation timing"),
                          tags$li("Survival: KM plots for timing groups"),
                          tags$li("ICB Response: Forest plots for ICB response"),
                          tags$li("ICB Survival: KM plots for ICB survival")
                        ),
                        
                        h3("Additional Tips"),
                        tags$ul(
                          tags$li("Use dropdowns to filter cancer or immune cell type"),
                          tags$li("Click 'Submit' after selecting new cohorts"),
                          tags$li("Hover over plots to view tooltips if enabled")
                        ),
                        
                        h3("Troubleshooting"),
                        tags$ul(
                          tags$li("If plots do not appear, make sure the data has loaded"),
                          tags$li("Use a supported browser (e.g., Chrome or Firefox)"),
                          tags$li("Contact: your.email@domain.com for help")
                        )
                      )
             )
  )
)

# Server logic
server <- function(input, output, session) {
  
  cohort_data <- reactiveVal(NULL)
  loading_flag <- reactiveVal(FALSE)
  
  output$loading <- reactive({
    loading_flag()
  })
  outputOptions(output, "loading", suspendWhenHidden = FALSE)
  
  ###################################################
  ## Page 2
  ###################################################
  ## Page2 -- Pathway Summary
  # Load CRISPR_all data
  filtered_data <- reactive({
    CRISPR_all %>% filter(Celltype == cell_type_files[[input$cell_type]])
  })

  # Calculate the frequencies of studies reporting the pathways
  study_count <- reactive({
    req(input$cell_type, input$immuno_regulators)
    cell_type_code <- cell_type_files[[input$cell_type]]
    reg <- input$immuno_regulators
    
    data_crispr <- data_crispr %>%
      filter(celltype == !!cell_type_code, regulators == !!reg)
    
    length(base::unique(data_crispr$Study))
  })
  
  # Set dynamic height and width for the plots
  output$immuno_pathway_plot_ui <- renderUI({
    req(input$immuno_n)
    
    base_width <- 600      
    extra_width <- 20   
    width_cap <- 1200
    
    base_height <- 450
    extra_height <- 30 
    height_cap <- 1200
    
    study_n <- study_count()
    total_width <- min(base_width + extra_width * study_n, width_cap)
    
    pathway_n <- input$immuno_n
    total_height <- min(base_height + extra_height * pathway_n, height_cap)
    
    plotlyOutput("immuno_pathway_plot",
                 height = paste0(total_height, "px"),
                 width = paste0(total_width, "px"))
  })
  
  # Plot the pathways with high frequencies
  output$immuno_pathway_plot <- renderPlotly({
    req(input$cell_type, input$immuno_regulators, input$immuno_n)
    
    cell_type_code <- cell_type_files[[input$cell_type]]
    reg <- input$immuno_regulators
    num <- input$immuno_n
    textcol <- "black"
    
    # Read data
    data_crispr <- data_crispr %>%
      filter(celltype == !!cell_type_code, regulators == !!reg)
    
    # Prepare plot
    plot_freq <- data_crispr %>%
      arrange(desc(Freq), desc(`-log10(Padj)`)) %>%
      distinct(Description, .keep_all = TRUE) %>%
      dplyr::slice_head(n = num) %>%
      dplyr::select(ID, Description, Freq)
    
    plot_order <- plot_freq %>% arrange(Freq) %>% pull(Description)
    
    plot_data <- inner_join(plot_freq, data_crispr, by = c("ID", "Description")) %>%
      mutate(
        tooltip_text = paste0(
          "Genes: ", geneName, "<br>",
          "Padj: ", format(p.adjust, scientific = TRUE, digits = 3), "<br>",
          "NES: ", round(NES, 1)
        ),
        Description = factor(Description, levels = plot_order)
      )
    
    P <- ggplot(plot_data, aes(x = Study, y = Description, size = `-log10(Padj)`, text = tooltip_text)) +
      geom_point(color = "#4393C3", alpha = 0.8) +
      labs(size = "-Log10(Padj)", title = "Pathway Enrichment", x = "", y = "") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 60, hjust = 1, size = 10, colour = textcol),
        axis.text.y = element_text(size = 10, colour = textcol),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.margin = margin(0.7, 0.4, 0.1, 0.2, "cm")
      )
    
    ggplotly(P, tooltip = "text") 
  })
  
  
  ###################################################
  ## Page 3 -- Select the WGS Datasets
  ###################################################
  # Load results from different cohorts
  observeEvent(input$submit_btn, {
    loading_flag(TRUE)
    
    withProgress(message = "Loading data...", value = 0.1, {
      data <- dynamic_cohort_data(input$cohort_select)
      incProgress(0.4, detail = "Processing data...")
      
      cohort_data(data)
      
      type_list <- data$clinical_data %>% distinct(aliquot_id, histology_abbreviation)
      updateSelectInput(session, "cancer_select", choices = c("All", base::unique(type_list$histology_abbreviation)))
      updateSelectInput(session, "timing_type", choices = base::unique(data$diff_all$histology_abbreviation))
      
      incProgress(0.5, detail = "Finalizing...")
      Sys.sleep(0.5)
    })
    
    loading_flag(FALSE)
  })
  
  # Top-left -- make barplot for sample number for each cancer type
  output$cohort_summary_left <- renderUI({
    tagList(
      plotOutput("histology_bar", height = "300px"),
      plotOutput("total_perMB_log_box", height = "300px")
    )
  })
  
  clinical_data_filtered <- reactive({
    req(cohort_data())
    clinical_data <- cohort_data()$clinical_data
    
    hist_count <- clinical_data %>%
      dplyr::count(histology_abbreviation) %>%
      filter(n >= 5)
    
    clinical_data %>%
      filter(histology_abbreviation %in% hist_count$histology_abbreviation)
  })
  
  output$histology_bar <- renderPlot({
    req(cohort_data())
    clinical_data <- cohort_data()$clinical_data
    
    clinical_data_filtered <- clinical_data_filtered()
    
    hist_order <- clinical_data_filtered %>%
      dplyr::count(histology_abbreviation) %>%
      arrange(desc(n)) %>%
      pull(histology_abbreviation)
    
    ggplot(clinical_data_filtered, aes(x = factor(histology_abbreviation, levels = hist_order), fill = histology_abbreviation)) +
      geom_bar() +
      scale_fill_manual(values = type_colors) +
      theme_minimal(base_size = 14) +
      labs(title = "Sample Frequency", x = "", y = "Sample Count") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none")
  })
  
  # Bottom-left -- make boxplot for the distribution of TMB
  output$total_perMB_log_box <- renderPlot({
    req(cohort_data())
    clinical_data_filtered <- clinical_data_filtered()
    
    hist_order <- clinical_data_filtered %>%
      group_by(histology_abbreviation) %>%
      summarise(median_tmb = median(TMB_log, na.rm = TRUE)) %>%
      arrange(desc(median_tmb)) %>%
      pull(histology_abbreviation)
    
    ggplot(clinical_data_filtered, aes(x = factor(histology_abbreviation, levels = hist_order), y = TMB_log, fill = histology_abbreviation)) +
      geom_boxplot() +
      scale_fill_manual(values = type_colors) +
      theme_minimal(base_size = 14) +
      scale_y_continuous(breaks = log10(c(0.1, 1, 10, 100, 1000, 10000)), 
                         labels = c(0.1, 1, 10, 100, 1000, 10000)) +
      labs(title = "TMB", x = "", y = "TMB") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none")
  })
  
  # Right -- make pieplot for sample type and age
  output$cohort_summary_right <- renderUI({
    tagList(
      plotOutput("gender_pie", height = "300px"),
      plotOutput("primary_pie", height = "300px")
    )
  })
  
  output$primary_pie <- renderPlot({
    req(cohort_data())
    clinical_data_filtered <- clinical_data_filtered()
    
    data <- clinical_data_filtered %>% dplyr::count(Sample_Type)
    ggplot(data, aes(x = "", y = n, fill = Sample_Type)) +
      geom_col() +
      coord_polar("y") +
      theme_void(base_size = 14) +
      labs(title = "Sample Type") + scale_fill_npg()
  })
  
  output$gender_pie <- renderPlot({
    req(cohort_data())
    clinical_data_filtered = clinical_data_filtered()
    
    data <- clinical_data_filtered %>% dplyr::count(Gender)
    ggplot(data, aes(x = "", y = n, fill = Gender)) +
      geom_col() +
      scale_fill_manual(values = Gender_colors) +
      coord_polar("y") +
      theme_void(base_size = 14) +
      labs(title = "Gender")
  })
  
  ###################################################
  ## Page 4 Immune Escape Evolution
  ###################################################
  # Single Pathway -- make bar plot for top 10 cancer types by Freq_mut for selected pathways
  output$cell_type_ui <- renderUI({
    selectInput("cell_type", "Immunity Scenarios", choices = names(cell_type_files))
  })
  
  output$cell_type_timeline_ui <- renderUI({
    selectInput("cell_type_timeline", "Immunity Scenarios", choices = c("All", names(cell_type_files)))
  })
  
  output$regulator_selector <- renderUI({
    selectInput("selected_regulator", "Select Regulator Type",
                choices = c("Positive" = "pos", "Negative" = "neg"),
                selected = "pos")
  })
  
  output$pathway_selector <- renderUI({
    req(input$cell_type, input$selected_regulator)
    
    choices <- CRISPR_all %>%
      filter(Celltype == cell_type_files[[input$cell_type]],
             Regulators == input$selected_regulator) %>%
      pull(pathway_new) %>%
      base::unique() %>%
      sort()
    
    selectInput("selected_pathways", "Select Pathways",
                choices = choices,
                selected = "Antigen Presentation via MHC-I (MHC-I)",
                multiple = FALSE)
  })
  
  output$top_cancer_plot <- renderPlot({
    tryCatch({
    req(input$cell_type, input$selected_pathways, cohort_data())
    
    top_data <- cohort_data()$Freq_all %>% filter(Freq_mut > 4) %>%
      filter(cell_type == cell_type_files[[input$cell_type]], pathway_new %in% input$selected_pathways) %>%
      group_by(histology_abbreviation) %>%
      summarise(
        total_mut = sum(percent_mut, na.rm = TRUE),
        number_mut = sum(Freq_mut, na.rm = TRUE)  # Add raw mutation count
      ) %>%
      arrange(desc(total_mut)) %>%
      dplyr::slice_head(n = 10)
    
    if (nrow(top_data) > 0 ) {
      ggplot(top_data, aes(x = reorder(histology_abbreviation, total_mut), y = total_mut)) +
        geom_col(fill = "steelblue") +
        geom_text(aes(label = number_mut),
                  hjust = -0.1, size = 5) +
        coord_flip() +
        scale_y_continuous(labels = percent_format(accuracy = 1),
                           expand = expansion(mult = c(0, 0.15))) +
        labs(x = "", y = "Mutation Frequency", title = "") +
        theme_minimal(base_size = 14) +
        theme(
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text.x = element_text(hjust = 0.5, size = 14, colour = "black"),
          axis.text.y = element_text(hjust = 1, size = 14, colour = "black")
        )
      
    } else {
      
      ggplot() +
        theme_void() +
        annotate("text", x = 0.5, y = 0.5, label = "No sufficient data", size = 8, hjust = 0.5)
      
    }
  
    }, error = function(e) {
      ggplot() +
        theme_void() +
        annotate("text", x = 0.5, y = 0.5, label = "No sufficient data", size = 8, hjust = 0.5)
    })
  })
  
  # Multi Pathways -- make Radar Plot for the frequencies of mutations in different pathway across cancers
  output$radarplot <- renderPlot({
    req(cohort_data(), input$cell_type)
    
    validate(need(input$cell_type %in% names(cell_type_files), "Loading..."))
    
    cell_type_key <- cell_type_files[[input$cell_type]]
    
    plots <- plot_radar_all(
      Freq_all = cohort_data()$Freq_all,
      cell_type = cell_type_key
    )
    
    gridExtra::grid.arrange(plots$pos, plots$neg, ncol = 2)
  }, res = 100)
  
  # Timeline -- make plots for the temporal orders of immunomodulatory pathways and driver genes in cancers
  output$timing_type_ui <- renderUI({
    selectInput("timing_type", "Select Cancer Type", choices = base::unique(cohort_data()$diff_all$histology_abbreviation), selected = "Breast-AdenoCA")
  })
  
  output$timing_plot <- renderPlot({
    req(input$timing_type, input$cell_type_timeline, cohort_data())
    
    cell_types <- if (input$cell_type_timeline == "All") {
      unname(cell_type_files)
    } else {
      cell_type_files[[input$cell_type_timeline]]  
    }
    
    p <- plot_timing_summary(
      diff_all = cohort_data()$diff_all,
      driver_list = driver_list,
      clinical_data = cohort_data()$clinical_data,
      type = input$timing_type,
      cell_type = cell_types
    )
    
    if (!is.null(p)) {
      return(p)
    } else {
      return(
        ggplot() + theme_void() +
          annotate("text", x = 0.5, y = 0.5, label = "No sufficient data", size = 8, hjust = 0.5)
      )
    }
  })
  
  # Survival -- make survival plot by mutation timing
  output$select_histology_ui <- renderUI({
    req(cohort_data())
    histologies <- sort(base::unique(cohort_data()$diff_all$histology_abbreviation))
    selectInput("histology_type", "Select Cancer Type:", choices = c("All", histologies), selected = "All")
  })
  
  output$pathway_selector_surv <- renderUI({
    req(input$cell_type, input$selected_regulator, cohort_data())
    
    choices <- cohort_data()$diff_all %>%
      filter(celltype == !!cell_type_files[[input$cell_type]]) %>%
      pull(pathway) %>% base::unique()
    
    selectInput("selected_pathways", "Select Pathways",
                choices = choices,
                selected = "Antigen Presentation via MHC-I (MHC-I)",
                multiple = FALSE)
  })
  
  output$timing_survival_plot <- renderPlot({
    tryCatch({
      req(input$histology_type,  input$selected_regulator)
      
      wt_included <- "no"
      histology_type <- input$histology_type
      selected_pathways <- input$selected_pathways
      
      data_survival <- cohort_data()$surv_data
      diff_all <- cohort_data()$diff_all %>%
        mutate(Timing_group = case_when(early_ratio > 0.5 ~ "Mut_Early",
                                        late_ratio > 0.5 ~ "Mut_Late",
                                        TRUE ~ "Undetermined")) %>%
        filter(Timing_group != "Undetermined")
      
      sample_survival <- base::unique(data_survival$aliquot_id)
      diff_all <- diff_all %>% filter(aliquot_id %in% sample_survival)
      
      data_survival_sub <- if (histology_type == "All") data_survival else filter(data_survival, histology_abbreviation == histology_type)
      diff_all_sub <- if (histology_type == "All") diff_all else filter(diff_all, histology_abbreviation == histology_type)
      
      filtered_data <- diff_all_sub %>%
        filter(pathway == selected_pathways) %>%
        left_join(data_survival_sub, .,by = "aliquot_id") %>%
        mutate(
          Timing_group = ifelse(is.na(Timing_group), "WT", Timing_group),
          Timing_group = factor(Timing_group, levels = c("WT", "Mut_Early", "Mut_Late"))
        )
      
      plot_survival_by_timing(
        filtered_data = filtered_data,
        wt_included = "no",
        title_prefix = ""
      )
      
    }, error = function(e) {
      ggplot() +
        theme_void() +
        annotate("text", x = 0.5, y = 0.5, label = "No sufficient data", size = 8, hjust = 0.5)
    })
  }, res = 90)
  
  # Immune cell infiltration -- make boxplot for the composistion of immune cells by mutation timing (Cibersort)
  output$select_histology_ciber <- renderUI({
    req(cohort_data())
    histologies <- sort(base::unique(cohort_data()$diff_all$histology_abbreviation))
    selectInput("histology_ciber", "Select Cancer Type:", choices = c("All", histologies), selected = "All")
  })
  
  output$pathway_selector_ciber <- renderUI({
    req(input$cell_type, input$selected_regulator, cohort_data())
    
    choices <- cohort_data()$diff_all %>%
      filter(celltype == !!cell_type_files[[input$cell_type]]) %>%
      pull(pathway) %>% base::unique()
    
    selectInput("pathway_ciber", "Select Pathways",
                choices = choices,
                selected = "Antigen Presentation via MHC-I (MHC-I)",
                multiple = FALSE)
  })
  
  output$ciber_plot <- renderPlot({
    tryCatch({
      req(input$histology_ciber, input$pathway_ciber, cohort_data())
      
      histology_type <- input$histology_ciber
      selected_pathways <- input$pathway_ciber
      
      ciber_all <- cohort_data()$ciber_all
      diff_all <- cohort_data()$diff_all %>%
        mutate(Timing_group = case_when(early_ratio > 0.5 ~ "Mut_Early", 
                                        late_ratio > 0.5 ~ "Mut_Late",
                                        TRUE ~ "Undetermined")) %>% 
        filter(Timing_group != "Undetermined")
      
      diff_all_sub <- if (histology_type == "All") diff_all else filter(diff_all, histology_abbreviation == histology_type)

      diff_ciber <- inner_join(diff_all_sub, ciber_all, by = "aliquot_id") %>%
        filter(pathway == selected_pathways)
      
      ciber_long <- diff_ciber %>% 
        gather(cell_type, Value_1, -pathway, -histology_abbreviation, -celltype, -regulator, -Timing_group) %>% 
        mutate(Value_1 = as.numeric(Value_1),
               cell_type = as.factor(cell_type)) %>% 
        filter(cell_type %in% cibercom_cols)
      
      validate(need(nrow(ciber_long) > 0, "No sufficient data"))
      
      plot_ciber(ciber_long, histology_type, selected_pathways)
      
    }, error = function(e) {
      ggplot() + theme_void() +
        annotate("text", x = 0.5, y = 0.5, label = "No sufficient data", size = 8, hjust = 0.5)
    })
  })
  
  ###################################################
  ## Page 5
  ###################################################
  # Mutation Frequency -- make the onocoplot for the selected gene list
  maf_obj <- reactive({
    req(cohort_data())
    cohort_data()$maf_data
  })
  
  get_selected_genes <- reactive({
    if (input$useGeneList && nzchar(input$geneList)) {
      unlist(strsplit(input$geneList, "[,\n]")) %>% trimws() %>% base::unique()
    } else {
      c(
        "ACHE", "AMOT", "CDK5R1", "CDK6", "CELSR1", "CNTFR", "CRMP1", "DPYSL2",
        "ETS2", "GLI1", "ADGRG1", "HEY1", "HEY2", "L1CAM", "LDB1", "MYH9", "NF1",
        "NKX6-1", "NRCAM", "NRP1", "TP53")
    }
  })
  
  output$oncoplot_ui <- renderUI({
    req(get_selected_genes())
    
    gene_n <- length(get_selected_genes())
    
    base_height <- 200    
    extra_height <- 20  
    height_cap <- 1200
    
    total_height <- min(base_height + extra_height * gene_n, height_cap)
    
    plotOutput("oncoplot", height = paste0(total_height, "px"), width = "600px")
  })
  
  output$oncoplot <- renderPlot({
  req(maf_obj())

  genes_to_plot <- get_selected_genes()
  req(length(genes_to_plot) > 0)

  oncoplot(
    maf = maf_obj(),
    genes = genes_to_plot,
    bgCol = "#efefef",
    colors = cols,
    removeNonMutated = TRUE
  )
  }, res = 100)
  
  # Mutation Timing -- make the density plot for distribution of mutation timing for the selected gene list
  output$timing_plot_singlegene_ui <- renderUI({
    req(get_selected_genes())
    
    gene_n <- length(get_selected_genes())
    
    base_height <- 150
    extra_height <- 20
    height_cap <- 1500
    
    total_height <- min(base_height + extra_height * gene_n, height_cap)
    
    plotOutput("timing_plot_singlegene", height = paste0(total_height, "px"), width = "400px")
  })
  
  output$timing_plot_singlegene <- renderPlot({
    req(cohort_data())
    selected_genes <- get_selected_genes()
    
    diff_density <- cohort_data()$diff_allgene %>%
      filter(Hugo_Symbol %in% selected_genes) %>%
      mutate(mean_diff = rowMeans(as.matrix(dplyr::select(., all_of(column_names)))))
    
    p <- plot_raincloud(diff_density, 
                        x = "Hugo_Symbol", 
                        y = "mean_diff", 
                        mode = "mean", 
                        "", 
                        width = 4, 
                        height = 5)
    req(p)
    p
  }, res = 70)
  
  # Pathway Timing -- make the density plot for distribution of mean mutation timing for the selected gene list
  output$timing_plot_gene <- renderPlot({
    tryCatch({
      
    req(cohort_data())
    selected_genes <- get_selected_genes()
    p <- plot_timing_bar_by_genes(
      selected_genes = selected_genes,
      diff_data = cohort_data()$diff_allgene,
      type_list = cohort_data()$diff_allgene %>% distinct(sample_id, histology_abbreviation),
      column_names = column_names,
      timing_bin_levels = timing_bin_levels,
      color_vector = color_vector,
      legend = legend
    )
    req(p)
    p
    
    }, error = function(e) {
      ggplot() + theme_void() +
        annotate("text", x = 0.5, y = 0.5, label = "No sufficient data", size = 8, hjust = 0.5)
    })
  })
  
  # Survival -- make survival plot by mutation timing for the selected gene list
  output$select_histology_list <- renderUI({
    req(cohort_data())
    histologies <- sort(base::unique(cohort_data()$diff_all$histology_abbreviation))
    selectInput("histology_type_genelist", "Select Cancer Type:", choices = c("All", histologies), selected = "All")
  })
  
  output$timing_survival_plot_genelist <- renderPlot({
    tryCatch({
      
      req(get_selected_genes(), input$histology_type_genelist, cohort_data())
      
      LOH_type = cohort_data()$clinical_data %>% dplyr::select(aliquot_id, histology_abbreviation)
      
      selected_genes <- get_selected_genes()
      histology_type <- input$histology_type_genelist
      wt_included <- "no"
      
      df_Reg <- cohort_data()$diff_allgene %>%
        filter(Hugo_Symbol %in% selected_genes) %>%
        group_by(sample_id) %>%
        summarise(across(all_of(column_names), ~ mean(., na.rm = TRUE)), .groups = "drop") %>%
        mutate(
          aliquot_id = sample_id,
          early_ratio = rowMeans(across(all_of(column_names)) < 0, na.rm = TRUE),
          late_ratio = rowMeans(across(all_of(column_names)) > 0, na.rm = TRUE),
          Timing_group = case_when(early_ratio > 0.5 ~ "Mut_Early", 
                                   late_ratio > 0.5 ~ "Mut_Late",
                                   TRUE ~ "Undetermined")) %>% 
        filter(Timing_group != "Undetermined") %>% inner_join(., LOH_type, by = "aliquot_id")
      
      data_survival <- cohort_data()$surv_data
      diff_all <- df_Reg %>%
        filter(aliquot_id %in% base::unique(data_survival$aliquot_id))
      
      data_survival_sub <- if (histology_type == "All") data_survival else filter(data_survival, histology_abbreviation == histology_type)
      diff_all_sub <- if (histology_type == "All") diff_all else filter(diff_all, histology_abbreviation == histology_type)
      
      filtered_data <- diff_all_sub %>%
        left_join(data_survival_sub, by = "aliquot_id") %>%
        mutate(
          Timing_group = ifelse(is.na(Timing_group), "WT", Timing_group),
          Timing_group = factor(Timing_group, levels = c("WT", "Mut_Early", "Mut_Late"))
        )
      
      plot_survival_by_timing(
        filtered_data = filtered_data,
        wt_included = "no",
        title_prefix = ""
      )
      
    }, error = function(e) {
      ggplot() + theme_void() +
        annotate("text", x = 0.5, y = 0.5, label = "No sufficient data", size = 8, hjust = 0.5)
    })
  }, res = 90)
  
  # ICB Response - make the forest plots for ICB response for the patients with or withour mutations in the selected gene list
  observe({
    req(ICB_merge)
    therapy_choices <- c("All", base::unique(ICB_merge$Therapy))
    updateSelectInput(session, "therapy_select_genes",
                      choices = therapy_choices,
                      selected = "All")
  })
  
  output$forestplot <- renderPlot({
    tryCatch({
    req(cohort_data(), input$therapy_select_genes)
    selected_genes <- get_selected_genes()
    
    ICB_subset <- if (input$therapy_select_genes == "All") {
      ICB_merge
    } else {
      ICB_merge %>% filter(Therapy == input$therapy_select_genes)
    }
    
    plot_forest_gg(ICB_subset, selected_genes)
    }, error = function(e) {
      ggplot() + theme_void() +
        annotate("text", x = 0.5, y = 0.5, label = "No sufficient data", size = 8, hjust = 0.5)
    })}, res = 100)
  
  # ICB Survival - make the survival for ICB response for the patients with or without mutations in the selected gene list
  observe({
    req(ICB_surv)
    therapy_choices <- c("All", base::unique(ICB_surv$Therapy))
    updateSelectInput(session, "therapy_surv_genes",
                      choices = therapy_choices,
                      selected = "All")
  })
  
  output$ICB_survival_plot <- renderUI({
    tryCatch({
    req(cohort_data(), get_selected_genes(), input$therapy_surv_genes)
    selected_genes <- get_selected_genes()
    
    ICB_survival <- ICB_surv %>%
      filter(!is.na(OS_status)) %>%
      { if (input$therapy_surv_genes == "All") . else filter(., Therapy == input$therapy_surv_genes) }
    
    mutated_samples <- ICB_survival %>%
      filter(Hugo_Symbol %in% selected_genes) %>%
      distinct(SampleName, CohortName) %>%
      mutate(Mutated = 1)
    
    sample_data <- ICB_survival %>%
      distinct(SampleName, CohortName, OS_status, OS_time) %>%
      left_join(mutated_samples, by = c("SampleName", "CohortName")) %>%
      mutate(Mutated = ifelse(is.na(Mutated), 0, 1)) %>%
      filter(!is.na(OS_status), !is.na(OS_time)) %>%
      group_by(CohortName) %>%
      filter(length(base::unique(OS_status)) > 1, length(base::unique(Mutated)) > 1) %>%
      ungroup()
    
    n_plots <- length(base::unique(sample_data$CohortName))
    
    base_height <- 350  
    plots_per_row <- 4 
    n_rows <- ceiling(n_plots / plots_per_row)
    plot_height <- base_height * n_rows
    
    cols_group <- c("WT" = "#4DBBD5FF", "Mutated" = "#E64B35FF")
    legend_labels <- c("WT", "Mutated")
    
    km_plots <- sample_data %>%
      split(.$CohortName) %>%
      purrr::map(~{
        data <- .
        fit <- survfit(Surv(OS_time, OS_status) ~ Mutated, data = data)
        
        max_time <- ceiling(max(data$OS_time, na.rm = TRUE) / 12) * 12
        interval <- ceiling(max_time / 60) * 12
        
        surv_plot <- ggsurvplot(
          fit,
          data = data,
          risk.table = TRUE,
          pval = TRUE,
          conf.int = FALSE,
          legend.title = "Mutation",
          legend.labs = legend_labels,
          palette = cols_group,
          xlab = "",
          break.x.by = interval,
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE
        )
        
        surv_plot$plot <- surv_plot$plot +
          theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)) +
          labs(title = base::unique(data$CohortName))
        
        cowplot::plot_grid(
          surv_plot$plot,
          surv_plot$table,
          align = "v",
          ncol = 1,
          rel_heights = c(4, 1.5)
        )
      })
    
    final_plot <- cowplot::plot_grid(plotlist = km_plots, ncol = 4, align = "v")
    
    renderPlot({
      final_plot
    }, height = plot_height, res = 70)
    
    }, error = function(e) {
      ggplot() + theme_void() +
        annotate("text", x = 0.5, y = 0.5, label = "No sufficient data", size = 8, hjust = 0.5)
    })
  })
  
  ###################################################
  ## Page 6
  ###################################################
  # CRISPR Ranking -- plot the top-rank regulators (positive/negative) by the CRISPR studies
  gstable_data <- reactive({
    req(input$file)
    
    gdata <- read.table(input$file$datapath, header = TRUE, stringsAsFactors = FALSE)
    comm <- intersect(gdata$id, MouseGene)
    req("id" %in% colnames(gdata))
    
    if (length(comm) > 50) {
      gdata <- left_join(gdata, MouseToHuman_gene, by = "id")
    } else {
      gdata$HumanGene <- gdata$id
    }
    
    gdata <- gdata %>% filter(!is.na(HumanGene) & HumanGene != "")
    
    gdata_pos <- gdata %>% dplyr::select(HumanGene, `pos.score`, `pos.rank`) %>%
      group_by(HumanGene) %>%
      dplyr::slice_min(`pos.rank`, with_ties = FALSE) %>%
      ungroup()
    
    gdata_neg <- gdata %>% dplyr::select(HumanGene, `neg.score`, `neg.rank`) %>%
      group_by(HumanGene) %>%
      dplyr::slice_min(`neg.rank`, with_ties = FALSE) %>%
      ungroup()
    
    gdata_rank <- inner_join(gdata_pos, gdata_neg, by = "HumanGene")
    
    return(gdata_rank)
  })
  
  output$crisprPlot <- renderPlot({
    req(input$top_n)
    gstable <- gstable_data()
    p1 <- plot_selection(gstable, "pos.rank", "pos.score", "Top Positive Regulators", input$top_n, color_offset = 0)
    p2 <- plot_selection(gstable, "neg.rank", "neg.score", "Top Negative Regulators", input$top_n, color_offset = 10)
    p1 + p2
  }, res = 80)
  
  # Mutation Frequency -- make the onocoplot for the top-rank regulators (positive/negative) by the CRISPR studies
  output$oncoplot_pos_ui <- renderUI({
    req(input$top_n)
    
    base_height <- 200
    extra_height <- 20
    height_cap <- 1000
    
    height <- min(base_height + extra_height * input$top_n, height_cap)
    
    plotOutput("oncoplot_pos", height = paste0(height, "px"))
  })
  
  output$oncoplot_neg_ui <- renderUI({
    req(input$top_n)
    
    base_height <- 200
    extra_height <- 20
    height_cap <- 1000
    
    height <- min(base_height + extra_height * input$top_n, height_cap)
    
    plotOutput("oncoplot_neg", height = paste0(height, "px"))
  })
  
  output$oncoplot_pos <- renderPlot({
    req(maf_obj())
    pos_top <- gstable_data() %>% arrange(pos.rank) %>% dplyr::slice(1:input$top_n) %>% pull(HumanGene)
    oncoplot(maf = maf_obj(), genes = pos_top, bgCol = "#efefef", colors = cols, removeNonMutated = TRUE)
  }, res = 80)
  
  output$oncoplot_neg <- renderPlot({
    req(maf_obj())
    neg_top <- gstable_data() %>% arrange(neg.rank) %>% dplyr::slice(1:input$top_n) %>% pull(HumanGene)
    oncoplot(maf = maf_obj(), genes = neg_top, bgCol = "#efefef", colors = cols, removeNonMutated = TRUE)
  }, res = 80)
  
  # Gene Timing -- make the density plot for distribution of mutation timing for the top-rank regulators (positive/negative) by the CRISPR studies
  output$timing_plot_singlegene_pos_ui <- renderUI({
    req(input$top_n)
    
    base_height <- 150
    extra_height <- 20
    height_cap <- 1500
    
    height <- min(base_height + extra_height * input$top_n, height_cap)
    
    plotOutput("timing_plot_singlegene_pos", height = paste0(height, "px"), width = "400px")
  })
  
  output$timing_plot_singlegene_neg_ui <- renderUI({
    req(input$top_n)
    
    base_height <- 150
    extra_height <- 20
    height_cap <- 1500
    
    height <- min(base_height + extra_height * input$top_n, height_cap)
    
    plotOutput("timing_plot_singlegene_neg", height = paste0(height, "px"), width = "400px")
  })
  
  output$timing_plot_singlegene_pos <- renderPlot({
    req(cohort_data())
    selected_genes <- gstable_data() %>% arrange(pos.rank) %>% dplyr::slice(1:input$top_n) %>% pull(HumanGene)
    
    diff_density <- cohort_data()$diff_allgene %>%
      filter(Hugo_Symbol %in% selected_genes) %>%
      mutate(mean_diff = rowMeans(as.matrix(dplyr::select(., all_of(column_names)))))
    
    plot_raincloud(
      diff_density,
      x = "Hugo_Symbol",
      y = "mean_diff",
      mode = "mean",
      title = "",
      width = 4,
      height = 5
    )
  }, res = 70)
  
  output$timing_plot_singlegene_neg <- renderPlot({
    req(cohort_data())
    selected_genes <- gstable_data() %>% arrange(neg.rank) %>% dplyr::slice(1:input$top_n) %>% pull(HumanGene)
    
    diff_density <- cohort_data()$diff_allgene %>%
      filter(Hugo_Symbol %in% selected_genes) %>%
      mutate(mean_diff = rowMeans(as.matrix(dplyr::select(., all_of(column_names)))))
    
    plot_raincloud(
      diff_density,
      x = "Hugo_Symbol",
      y = "mean_diff",
      mode = "mean",
      title = "",
      width = 4,
      height = 5
    )
  }, res = 70)
   
  # Pathway Timing -- make the density plot for distribution of mean mutation timing for the top-rank regulators (positive/negative) by the CRISPR studies
  output$timing_plot_pos <- renderPlot({
    tryCatch({
    req(cohort_data())
    selected_genes <- gstable_data() %>% arrange(pos.rank) %>% dplyr::slice(1:input$top_n) %>% pull(HumanGene)
    type_list <- cohort_data()$diff_allgene %>% distinct(sample_id, histology_abbreviation)
    plot_timing_bar(cohort_data()$diff_allgene, selected_genes, type_list)
    
    }, error = function(e) {
      ggplot() + theme_void() +
        annotate("text", x = 0.5, y = 0.5, label = "No sufficient data", size = 8, hjust = 0.5)
    })
  })
  
  output$timing_plot_neg <- renderPlot({
    tryCatch({
    req(cohort_data())
    selected_genes <- gstable_data() %>% arrange(neg.rank) %>% dplyr::slice(1:input$top_n) %>% pull(HumanGene)
    type_list <- cohort_data()$diff_allgene %>% distinct(sample_id, histology_abbreviation)
    plot_timing_bar(cohort_data()$diff_allgene, selected_genes, type_list)
    
    }, error = function(e) {
      ggplot() + theme_void() +
        annotate("text", x = 0.5, y = 0.5, label = "No sufficient data", size = 8, hjust = 0.5)
    })
  })
  
  # Survival -- make survival plot by mutation timing for the top-rank regulators (positive/negative) by the CRISPR studies
  output$select_histology_study <- renderUI({
    req(cohort_data())
    histologies <- cohort_data()$diff_all %>% pull(histology_abbreviation) %>% base::unique(.)
    selectInput("histology_type_study", "Select Cancer Type:", choices = c("All", histologies), selected = "All")
  })
  
  output$timing_survival_plot_pos <- renderPlot({
    tryCatch({
      
      req(cohort_data(), gstable_data(), input$top_n, input$histology_type_study)
      LOH_type = cohort_data()$clinical_data %>% dplyr::select(aliquot_id, histology_abbreviation)
      
      pos_genes <- gstable_data() %>% arrange(pos.rank) %>% slice_head(n = input$top_n) %>% pull(HumanGene)
      
      selected_genes <- pos_genes
      histology_type <- input$histology_type_study
      wt_included <- "no"
      
      # Prepare data
      df_Reg <- cohort_data()$diff_allgene %>%
        filter(Hugo_Symbol %in% selected_genes) %>%
        group_by(sample_id) %>%
        summarise(across(all_of(column_names), ~ mean(., na.rm = TRUE)), .groups = "drop") %>%
        mutate(
          aliquot_id = sample_id,
          early_ratio = rowMeans(across(all_of(column_names)) < 0, na.rm = TRUE),
          late_ratio = rowMeans(across(all_of(column_names)) > 0, na.rm = TRUE),
          Timing_group = case_when(early_ratio > 0.5 ~ "Mut_Early", 
                                   late_ratio > 0.5 ~ "Mut_Late",
                                   TRUE ~ "Undetermined")) %>% 
        filter(Timing_group != "Undetermined")  %>%  inner_join(., LOH_type, by = "aliquot_id")
      
      data_survival <- cohort_data()$surv_data
      diff_all <- df_Reg %>%
        filter(aliquot_id %in% base::unique(data_survival$aliquot_id))
      
      data_survival_sub <- if (histology_type == "All") data_survival else filter(data_survival, histology_abbreviation == histology_type)
      diff_all_sub <- if (histology_type == "All") diff_all else filter(diff_all, histology_abbreviation == histology_type)
      
      filtered_data <- diff_all_sub %>%
        left_join(data_survival_sub, by = "aliquot_id") %>%
        mutate(
          Timing_group = ifelse(is.na(Timing_group), "WT", Timing_group),
          Timing_group = factor(Timing_group, levels = c("WT", "Mut_Early", "Mut_Late"))
        )
      
      plot_survival_by_timing(
        filtered_data = filtered_data,
        wt_included = "no",
        title_prefix = ""
      )
      
    }, error = function(e) {
      ggplot() + theme_void() +
        annotate("text", x = 0.5, y = 0.5, label = "No sufficient data", size = 8, hjust = 0.5)
    })
  }, res = 90)
  
  output$timing_survival_plot_neg <- renderPlot({
    tryCatch({
      
      req(cohort_data(), gstable_data(), input$top_n, input$histology_type_study)
      LOH_type = cohort_data()$clinical_data %>% dplyr::select(aliquot_id, histology_abbreviation)
      
      pos_genes <- gstable_data() %>% arrange(neg.rank) %>% slice_head(n = input$top_n) %>% pull(HumanGene)
      
      selected_genes <- pos_genes
      histology_type <- input$histology_type_study
      wt_included <- "no"
      
      # Prepare data
      df_Reg <- cohort_data()$diff_allgene %>%
        filter(Hugo_Symbol %in% selected_genes) %>%
        group_by(sample_id) %>%
        summarise(across(all_of(column_names), ~ mean(., na.rm = TRUE)), .groups = "drop") %>%
        mutate(
          aliquot_id = sample_id,
          early_ratio = rowMeans(across(all_of(column_names)) < 0, na.rm = TRUE),
          late_ratio = rowMeans(across(all_of(column_names)) > 0, na.rm = TRUE),
          Timing_group = case_when(early_ratio > 0.5 ~ "Mut_Early", 
                                   late_ratio > 0.5 ~ "Mut_Late",
                                   TRUE ~ "Undetermined")) %>% 
        filter(Timing_group != "Undetermined")  %>%  inner_join(., LOH_type, by = "aliquot_id")
      
      data_survival <- cohort_data()$surv_data
      diff_all <- df_Reg %>%
        filter(aliquot_id %in% base::unique(data_survival$aliquot_id))
      
      data_survival_sub <- if (histology_type == "All") data_survival else filter(data_survival, histology_abbreviation == histology_type)
      diff_all_sub <- if (histology_type == "All") diff_all else filter(diff_all, histology_abbreviation == histology_type)
      
      filtered_data <- diff_all_sub %>%
        left_join(data_survival_sub, by = "aliquot_id") %>%
        mutate(
          Timing_group = ifelse(is.na(Timing_group), "WT", Timing_group),
          Timing_group = factor(Timing_group, levels = c("WT", "Mut_Early", "Mut_Late"))
        )
      
      plot_survival_by_timing(
        filtered_data = filtered_data,
        wt_included = "no",
        title_prefix = ""
      )
      
    }, error = function(e) {
      ggplot() + theme_void() +
        annotate("text", x = 0.5, y = 0.5, label = "No sufficient data", size = 8, hjust = 0.5)
    })
  }, res = 90)
  
  # ICB Response -- make the forest plots for ICB response for the patients with or without mutations for the top-rank regulators (positive/negative) by the CRISPR studies
  observe({
    req(ICB_merge)
    therapy_choices <- c("All", base::unique(ICB_merge$Therapy))
    updateSelectInput(session, "therapy_select_study",
                      choices = therapy_choices,
                      selected = "All")
  })
  
  output$forestplot_pos <- renderPlot({
    tryCatch({
    req(cohort_data(), input$therapy_select_study, ICB_merge)
    selected_genes <- gstable_data() %>% arrange(pos.rank) %>% dplyr::slice(1:input$top_n) %>% pull(HumanGene)
    
    ICB_subset <- if (input$therapy_select_study == "All") {
      ICB_merge
    } else {
      ICB_merge %>% filter(Therapy == input$therapy_select_study)
    }
    
    plot_forest_gg(ICB_subset, selected_genes)
    }, error = function(e) {
      ggplot() + theme_void() +
        annotate("text", x = 0.5, y = 0.5, label = "No sufficient data", size = 8, hjust = 0.5)
    })}, res = 100)
  
  output$forestplot_neg <- renderPlot({
    tryCatch({
      req(cohort_data(), input$therapy_select_study, ICB_merge)
      selected_genes <- gstable_data() %>% arrange(neg.rank) %>% dplyr::slice(1:input$top_n) %>% pull(HumanGene)
      
      ICB_subset <- if (input$therapy_select_study == "All") {
        ICB_merge
      } else {
        ICB_merge %>% filter(Therapy == input$therapy_select_study)
      }
      
      plot_forest_gg(ICB_subset, selected_genes)
    }, error = function(e) {
      ggplot() + theme_void() +
        annotate("text", x = 0.5, y = 0.5, label = "No sufficient data", size = 8, hjust = 0.5)
    })}, res = 100)
  
  # ICB Survival - make the survival for ICB response for the patients with or without mutations for the top-rank regulators (positive/negative) by the CRISPR studies
  observe({
    req(ICB_surv)
    therapy_choices <- c("All", base::unique(ICB_surv$Therapy))
    updateSelectInput(session, "therapy_surv_study",
                      choices = therapy_choices,
                      selected = "All")
  })
  
  output$ICB_survival_plot_study <- renderUI({
    req(cohort_data(), input$therapy_surv_study, input$regulator_type)
    
    selected_genes <- if (input$regulator_type == "Positive") {
      gstable_data() %>%
        arrange(pos.rank) %>%
        dplyr::slice(1:input$top_n) %>%
        pull(HumanGene)
    } else {
      gstable_data() %>%
        arrange(neg.rank) %>%
        dplyr::slice(1:input$top_n) %>%
        pull(HumanGene)
    }
    
    # Prepare survival data
    survival_df <- ICB_surv %>%
      filter(!is.na(OS_status)) %>%
      { if (input$therapy_surv_study == "All") . else filter(., Therapy == input$therapy_surv_study) }
    
    mutated_samples <- survival_df %>%
      filter(Hugo_Symbol %in% selected_genes) %>%
      distinct(SampleName, CohortName) %>%
      mutate(Mutated = 1)
    
    sample_data <- survival_df %>%
      distinct(SampleName, CohortName, OS_status, OS_time) %>%
      left_join(mutated_samples, by = c("SampleName", "CohortName")) %>%
      mutate(Mutated = ifelse(is.na(Mutated), 0, 1)) %>%
      group_by(CohortName) %>%
      filter(length(base::unique(OS_status)) > 1, length(base::unique(Mutated)) > 1) %>%
      ungroup()
    
    n_plots <- length(base::unique(sample_data$CohortName))
    plots_per_row <- 4
    n_rows <- ceiling(n_plots / plots_per_row)
    plot_height <- max(400, 350 * n_rows)
    
    cols_group <- c("WT" = "#4DBBD5FF", "Mutated" = "#E64B35FF")
    legend_labels <- c("WT", "Mutated")
    
    km_plots <- sample_data %>%
      split(.$CohortName) %>%
      purrr::map(~{
        data <- .
        data$Mutated <- ifelse(data$Mutated == 1, "Mut", "WT")
        data$Mutated <- factor(data$Mutated, levels = c("WT", "Mut"))
        
        fit <- survfit(Surv(OS_time, OS_status) ~ Mutated, data = data)
        
        max_time <- ceiling(max(data$OS_time, na.rm = TRUE) / 12) * 12
        interval <- ceiling(max_time / 60) * 12
        
        surv_plot <- ggsurvplot(
          fit,
          data = data,
          risk.table = TRUE,
          pval = TRUE,
          conf.int = FALSE,
          legend.title = "Mutation",
          legend.labs = legend_labels,
          palette = cols_group,
          xlab = "",
          break.x.by = interval,
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE
        )
        
        surv_plot$plot <- surv_plot$plot +
          theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)) +
          labs(title = base::unique(data$CohortName))
        
        cowplot::plot_grid(
          surv_plot$plot,
          surv_plot$table,
          align = "v",
          ncol = 1,
          rel_heights = c(4, 1.5)
        )
      })
    
    final_plot <- cowplot::plot_grid(plotlist = km_plots, ncol = 4, align = "v")
    
    renderPlot({
      final_plot
    }, height = plot_height, res = 70)
  })
  
  ## Page 2
  output$pathway_interpretation <- renderUI({
    HTML(
      paste0(
        "<p><strong>Interpretation:</strong> The enrichments are generated based on top 100 positive/negative regulators in each study by MAGeCKFlute. (PMID: 30710114)",
        "<br>The Dot plot shows the most frequently enriched pathways reported by CRISPR screen studies.",
        "<br>X-axis presents each study (Year_Journal_Author_Condition).",
        "<br>Y-axis shows the most frequently reported pathways.",
        "<br>Dot size shows the -log10(adjusted p-value) for the enrichement.",
        "<br>Each dot containes informtion for the enriched genes, adjusted p-value and NES score.",
        "<br>C5BP+KEGG+REACTOME databases are untilized as the reference database. The genesets with 5~500 are included for the analysis</p>"
      )
    )
  })
  
  ## Page 3
  output$cohortsummary_interpretation <- renderUI({
    HTML(
      paste0(
        "<br><strong>PCAWG:</strong> 2658 cancers from Pan-cancer analysis of whole genomes. (PMID: 32025007)",
        "<br><strong>POG570:</strong> 570 advanced and metastatic cancers from patients from the Personalized OncoGenomics (POG) Program. (PMID: 35121966)",
        "<br><strong>CRC_Prognosis:</strong> 1063 colorectal cancers from the Swedish population-based cohort. (PMID: 39112715) Coming Soon...",
        "<br><strong>TCGA-OV:</strong> 314 ovarian cancers from latest TCGC-OV WGS cohort. (PMID: 21720365) Coming Soon...",
        "<br><strong>TNBC:</strong> 254 tripple-negative breast cancers from TNBC WGS cohort. (PMID: 31570822) Coming Soon...</p>"
      )
    )
  })
  
  ## Page 4
  output$singlepathway_interpretation <- renderUI({
    HTML(
      paste0(
        "<p><strong>Interpretation:</strong> The barplot shows the top cancer types with the most frequent mutations in the selected pathway.",
        "<br>The number on the right of each bar represents the number of samples with mutations in that pathway.</p>"
      )
    )
  })
  
  output$multipathways_interpretation <- renderUI({
    HTML(
      paste0(
        "<p><strong>Interpretation:</strong> The radar plot shows the mutation frequency of different immunomodulatory pathway categories across cancer types.</p>"
      )
    )
  })
  
  output$timeline_interpretation <- renderUI({
    HTML(
      paste0(
        "<p><strong>Interpretation:</strong> The plot shows the mutation timing of immunomodulatory pathways and cancer drivers in selected cancer types.</p>"
      )
    )
  })
  
  output$timingsurvival_interpretation <- renderUI({
    HTML(
      paste0(
        "<p><strong>Interpretation:</strong> The Kaplan-Meier (KM) plot shows survival outcomes stratified by mutation timing groups.",
        "<br><strong>Timing Difference</strong> = MeanTiming(Pathway) - BackgroundTiming",
        "<br><strong>Mut_Early:</strong> Samples where more than 50% of the 250 timing difference samplings are less than 0.",
        "<br><strong>Mut_Late:</strong> Samples where more than 50% of the 250 timing difference samplings are greater than 0.",
        "<br>For more details, refer to our preprint (PMID: 39868264).</p>"
      )
    )
  })
  
  output$cibersort_interpretation <- renderUI({
    HTML(
      paste0(
        "<p><strong>Interpretation:</strong> The boxplot shows the composition of immune cell infiltration stratified by mutation timing in different immunomodulatory pathways.",
        "<br><strong>Timing Difference:</strong> = MeanTiming(Pathway) - BackgroundTiming",
        "<br><strong>Mut_Early:</strong> Samples where more than 50% of the 250 timing difference samplings are less than 0.",
        "<br><strong>Mut_Late:</strong> Samples where more than 50% of the 250 timing difference samplings are greater than 0.",
        "<br>Differences between groups were assessed using the Wilcoxon test.",
        "<br>The composition of immune cell infiltration were calculated by CIBERSORT Absolute mode. (PMID: 29344893)",
        "<br>Specifically in PCAWG cohort, only PCAWG-TCGA samples with available RNA-seq data were included.</p>"
      )
    )
  })
  
  output$mutation_freq_interpretation <- renderUI({
    HTML("<p><strong>Interpretation:</strong> The Oncoplot shows the mutation status of the selected genes across samples. </p>")
  })
  
  output$gene_timing_interpretation <- renderUI({
    HTML(
      paste0(
        "<p><strong>Interpretation:</strong> The raincloud plot shows the mutation timing distribution of the selected genes across samples.",
        "<br><strong>Timing Difference:</strong> = MeanTiming(Pathway) - BackgroundTiming",
        "<br>The genes are sorted by the mean timing difference, where negative values indicate earlier mutations and positive values indicate later mutations.</p>"
      )
    )
  })
  
  output$pathway_timing_gene_interpretation <- renderUI({
    HTML(
      paste0(
        "<p><strong>Interpretation:</strong> The barplot shows the average mutation timing of all selected genes across cancer types.",
        "<br><span style='color:green;'><strong>Green</strong></span> indicates samples with mutations likely occurring early (based on the proportion of timing values less than 0).",
        "<br><span style='color:purple;'><strong>Purple</strong></span> indicates samples with mutations likely occurring late (based on the proportion of timing values greater than 0).</p>"
      )
    )
  })
  
  output$timingsurvival_gene_interpretation <- renderUI({
    HTML(
      paste0(
        "<p><strong>Interpretation:</strong> The Kaplan-Meier (KM) plot shows survival outcomes stratified by mutation timing groups.",
        "<br><strong>Timing Difference:</strong> = MeanTiming(Pathway) - BackgroundTiming",
        "<br><strong>Mut_Early:</strong> Samples where more than 50% of the 250 timing difference samplings are less than 0.",
        "<br><strong>Mut_Late:</strong> Samples where more than 50% of the 250 timing difference samplings are greater than 0.",
        "<br>For more details, refer to our preprint (PMID: 39868264).</p>"
      )
    )
  })
  
  output$ICBresponse_interpretation <- renderUI({
    HTML(
      paste0(
        "<p><strong>Interpretation:</strong> The forest plot shows odds ratios (OR) for ICB response in samples with versus without mutations in the selected genes.",
        "<br>OR values are calculated using a generalized linear model (GLM), adjusting for tumor mutation burden (TMB).",
        "<br>For more details, please refer to the summary table. </p>"
      )
    )
  })
  
  output$ICBsurvival_interpretation <- renderUI({
    HTML(
      paste0(
        "<p><strong>Interpretation:</strong> The Kaplan-Meier (KM) plot shows overall survival outcomes for samples with or without mutations in the selected genes.",
        "<br>For more details, please refer to the summary table. </p>"
      )
    )
  })
  
  output$mutation_freq_study_interpretation <- renderUI({
    HTML("<p><strong>Interpretation:</strong> The Oncoplot shows the mutation status of the selected genes across samples. </p>")
  })
  
  output$gene_timing_study_interpretation <- renderUI({
    HTML(
      paste0(
        "<p><strong>Interpretation:</strong> The raincloud plot shows the mutation timing distribution of the selected genes across samples.",
        "<br><strong>Timing Difference:</strong> = MeanTiming(Pathway) - BackgroundTiming",
        "<br>The genes are sorted by the mean timing difference, where negative values indicate earlier mutations and positive values indicate later mutations.</p>"
      )
    )
  })
  
  output$pathway_timing_study_interpretation <- renderUI({
    HTML(
      paste0(
        "<p><strong>Interpretation:</strong> The barplot shows the average mutation timing of all selected genes across cancer types.",
        "<br><span style='color:green;'><strong>Green</strong></span> indicates samples with mutations likely occurring early (based on the proportion of timing values less than 0).",
        "<br><span style='color:purple;'><strong>Purple</strong></span> indicates samples with mutations likely occurring late (based on the proportion of timing values greater than 0).</p>"
      )
    )
  })
  
  output$timingsurvival_study_interpretation <- renderUI({
    HTML(
      paste0(
        "<p><strong>Interpretation:</strong> The Kaplan-Meier (KM) plot shows survival outcomes stratified by mutation timing groups.",
        "<br><strong>Timing Difference:</strong> = MeanTiming(Pathway) - BackgroundTiming",
        "<br><strong>Mut_Early:</strong> Samples where more than 50% of the 250 timing difference samplings are less than 0.",
        "<br><strong>Mut_Late:</strong> Samples where more than 50% of the 250 timing difference samplings are greater than 0.",
        "<br>For more details, refer to our preprint (PMID: 39868264).</p>"
      )
    )
  })
  output$ICBresponse_study_interpretation <- renderUI({
    HTML(
      paste0(
        "<p><strong>Interpretation:</strong> The forest plot shows odds ratios (OR) for ICB response in samples with versus without mutations in the selected genes.",
        "<br>OR values are calculated using a generalized linear model (GLM), adjusting for tumor mutation burden (TMB).",
        "<br>For more details, please refer to the summary table. (link...)</p>"
      )
    )
  })
  
  output$ICBsurvival_study_interpretation <- renderUI({
    HTML(
      paste0(
        "<p><strong>Interpretation:</strong> The Kaplan-Meier (KM) plot shows overall survival outcomes for samples with or without mutations in the selected genes.",
        "<br>For more details, please refer to the summary table. (link...)</p>"
      )
    )
  })
  
  output$CRISPRRanking_interpretation <- renderUI({
    HTML(
      paste0(
        "<p><strong>Interpretation:</strong> The plot shows the top-ranked regulators for the positive and negative selection.</p>"
      )
    )
  })
  
}

# Run app
shinyApp(ui, server)
