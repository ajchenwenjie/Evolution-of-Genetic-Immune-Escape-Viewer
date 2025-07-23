rm(list = ls())

options(warn = -1)

# Source helper functions
rm(list = ls())

options(warn = -1)

# Load libraries
libs <- c("broom", "forestplot", "shiny", "ggplot2", "dplyr", "tidyr", "data.table", 
          "readxl", "stringr", "ggradar", "scales", "gridExtra", "ggsci", "patchwork", 
          "maftools", "ggdist", "ggthemes", "ggrepel", "plotly", "textshape", 
          "survival", "survminer", "purrr", "shinycssloaders")
invisible(lapply(libs, function(pkg) suppressWarnings(library(pkg, character.only = TRUE))))

# Define file paths
base_path <- "./MAF"
data_path <- file.path(base_path, "data")

# Source helper script
source(file.path(base_path, "helpers.R"), local = TRUE)

# Load data
load_data_files <- function() {
  list(
    long_data = fread(file.path(data_path, "CRISPR/GeneMatrix_long.csv")),
    CRISPR_pathways = fread(file.path(data_path, "CRISPR/C5BP_selected_list.csv")),
    CRISPR_all = fread(file.path(data_path, "CRISPR/CRISPR_all.csv")),
    ICB_merge = fread(file.path(data_path, "ICB/ICB_merge.csv")),
    ICB_surv = fread(file.path(data_path, "ICB/ICB_survival.csv")),
    data_crispr = fread(file.path(data_path, "CRISPR/C5BP_all.csv")),
    MouseToHuman_gene = fread(file.path(data_path, "Refence_data/Referecen_MouseToHuman_gene.csv")),
    driver_list = read_xlsx(file.path(data_path, "Refence_data/TableS1_compendium_mutational_drivers_clean.xlsx"),
                            sheet = "TableS1_compendium_mutational_d")
  )
}

data_list <- load_data_files()
list2env(data_list, envir = .GlobalEnv)
MouseGene <- MouseToHuman_gene$id

cell_type_files <- c(
  "MHC-I Regulators" = "MHCregulators",
  "T-cell Killing" = "Tcells",
  "NK-cell Killing" = "NKcells",
  "Macrophage Killing" = "Macrophages",
  "GDT-cell Killing" = "GammadeltaTcells"
)

# UI code
ui <- tagList(
  navbarPage("Evolution of Genetic Immune Escape Viewer",
             tabPanel("Introduction",
                      fluidRow(
                        column(10, offset = 1,
                               h2("Welcome to the Evolution of Genetic Immune Escape Viewer"),
                               p("This Shiny application explores the landscape of immune escape evolution in cancer by integrating data from CRISPR screens, mutation frequencies, and mutation timing across multiple cancer types from the PCAWG project."),
                               tags$div(
                                 img(src = "immune_escape_intro.png", height = "400px", style = "display: block; margin-left: auto; margin-right: auto;"),
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
                                 "For more details, visit the ",
                                 tags$a(href = "https://github.com/YourRepo/immune-escape-viewer", 
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
                   uiOutput("immuno_pathway_plot_ui")
                 )
               )
             ),
             tabPanel("Cohort Summary",
                      sidebarLayout(
                        sidebarPanel(
                          width = 2,
                          selectInput("cohort_select", "Select Cohort", choices = c("PCAWG", "POG570", "TCGA-OV", "CRC-Prognosis", "TNBC", "Glioma")),
                          actionButton("submit_btn", "Submit")
                        ),
                        mainPanel(
                          conditionalPanel(
                            condition = "output.loading == true",
                            h4("Loading data.....", style = "color: #FF5733; text-align:center;")
                          ),
                          tabPanel("Cohort Summary",
                                   fluidRow(
                                     column(
                                       width = 9, 
                                       uiOutput("cohort_summary_left")
                                     ),
                                     column(
                                       width = 3,
                                       uiOutput("cohort_summary_right")
                                     )
                                   )
                          )
                        )
                      ),
             ),
             tabPanel("Immune Escape Evolution",
                      sidebarLayout(
                        sidebarPanel(
                          width = 2,
                          conditionalPanel(
                            condition = "input.crispr_tab_selected == 'Timeline'",
                            uiOutput("timing_type_ui")
                          ),
                          uiOutput("cell_type_ui"),
                          conditionalPanel(
                            condition = "input.crispr_tab_selected == 'Single Pathway'",
                            uiOutput("regulator_selector"),
                            uiOutput("pathway_selector")
                          ),
                          conditionalPanel(
                            condition = "input.crispr_tab_selected == 'Survival'",
                            uiOutput("pathway_selector_surv"),
                            uiOutput("select_histology_ui")
                          ),
                          conditionalPanel(
                            condition = "input.crispr_tab_selected == 'Immune Cell Infiltration'",
                            uiOutput("pathway_selector_ciber"),
                            uiOutput("select_histology_ciber")
                          ),
                        ),
                        mainPanel(
                          tabsetPanel(id = "crispr_tab_selected",
                                      tabPanel("Single Pathway",
                                               h4("Top 10 Cancer Types by Mutation Frequency"),
                                               plotOutput("top_cancer_plot", width = "400px", height = "400px")
                                      ),
                                      tabPanel("Multi Pathways",
                                               mainPanel(plotOutput("radarplot", width = "800px", height = "500px"))
                                      ),
                                      tabPanel("Timeline",
                                               div(style = "overflow-x: auto; text-align: left;",
                                                   plotOutput("timing_plot", width = "2000px", height = "500px")
                                               )
                                      ),
                                      tabPanel("Survival",
                                               withSpinner(plotOutput("timing_survival_plot", height = "500px", width = "450px"))
                                      ),
                                      tabPanel("Immune Cell Infiltration",
                                               plotOutput("ciber_plot", height = "400px", width = "800px")
                                      )
                          )
                        )
                      )
             ),
             tabPanel("Exploring Your Genes",
                      sidebarLayout(
                        sidebarPanel(
                          width = 2,
                          
                          checkboxInput("useGeneList", "Use Custom Gene List", value = FALSE),
                          conditionalPanel(
                            condition = "input.useGeneList == true",
                            textAreaInput("geneList", "Enter Gene Symbols (comma or newline separated)", "")
                          ),
                          conditionalPanel(
                            condition = "input.gene_tab_selected == 'Survival'",
                            uiOutput("select_histology_list")
                          ),
                          conditionalPanel(
                            condition = "input.gene_tab_selected == 'ICB Response'",
                            selectInput("therapy_select_genes", "Select Therapy:", choices = NULL)
                          ),
                          conditionalPanel(
                            condition = "input.gene_tab_selected == 'ICB Survival'",
                            tagList(
                              selectInput("therapy_surv_genes", "Select Therapy:", choices = NULL)
                            )
                          )
                        ),
                        mainPanel(
                          tabsetPanel(id = "gene_tab_selected",
                                      tabPanel("Mutation Frequency", 
                                               uiOutput("oncoplot_ui")
                                      ),
                                      tabPanel("Gene Timing", 
                                               uiOutput("timing_plot_singlegene_ui")
                                      ),
                                      tabPanel("Pathway Timing", 
                                               plotOutput("timing_plot_gene", width = "500px", height = "500px")
                                      ),
                                      tabPanel("Survival",
                                               withSpinner(plotOutput("timing_survival_plot_genelist", height = "500px", width = "450px"))
                                      ),
                                      tabPanel("ICB Response",
                                               div(
                                                 style = "margin-top: 0px; padding-top: 0px; display: flex; justify-content: left; align-items: flex-start;",
                                                 plotOutput("forestplot", width = "600px", height = "500px")
                                               )
                                      ),
                                      tabPanel("ICB Survival",
                                               withSpinner(uiOutput("ICB_survival_plot"))
                                      )
                          ),
                          br(),
                          uiOutput("tab_description")
                        )
                      )
             ),
             tabPanel("Exploring Your Study",
                      fluidPage(
                        sidebarLayout(
                          sidebarPanel(
                            width = 2,
                            fileInput("file", "Upload Gene Summary File", accept = ".txt"),
                            
                            numericInput("top_n", "Number of Top Genes:", value = 20, min = 10, step = 10),
                            helpText("Choose a cancer type to visualize mutation timing."),
                            conditionalPanel(
                              condition = "input.study_tab_selected == 'Survival'",
                              uiOutput("select_histology_study")
                            ),
                            conditionalPanel(
                              condition = "input.study_tab_selected == 'ICB Response'",
                              selectInput("therapy_select_study", "Select Therapy:", choices = NULL, selected = "All")
                            ),
                            conditionalPanel(
                              condition = "input.study_tab_selected == 'ICB Survival'",
                              selectInput("therapy_surv_study", "Select Therapy:", choices = NULL, selected = "All")
                            ),
                            conditionalPanel(
                              condition = "input.study_tab_selected == 'ICB Survival'",
                              selectInput("regulator_type", "Select Regulator Type",
                                          choices = c("Positive", "Negative"),
                                          selected = "Positive")
                            )
                          ),
                          mainPanel(
                            tabsetPanel(
                              id = "study_tab_selected",
                              tabPanel("CRISPR Ranking", plotOutput("crisprPlot", height = "350px", width = "800px")),
                              tabPanel("Mutation Frequency",
                                       fluidRow(
                                         column(6,
                                                h4("Top Positive Regulators"),
                                                uiOutput("oncoplot_pos_ui")  # dynamic height
                                         ),
                                         column(6,
                                                h4("Top Negative Regulators"),
                                                uiOutput("oncoplot_neg_ui")  # dynamic height
                                         )
                                       )
                              ),
                              tabPanel("Gene Timing",
                                       fluidRow(
                                         column(6,
                                                h4("Top Positive Regulators"),
                                                uiOutput("timing_plot_singlegene_pos_ui")  # dynamic height
                                         ),
                                         column(6,
                                                h4("Top Negative Regulators"),
                                                uiOutput("timing_plot_singlegene_neg_ui")  # dynamic height
                                         )
                                       )
                              ),
                              tabPanel("Pathway Timing",
                                       fluidRow(
                                         column(6, h4("Top Positive Regulators"), plotOutput("timing_plot_pos", height = "500px")),
                                         column(6, h4("Top Negative Regulators"), plotOutput("timing_plot_neg", height = "500px"))
                                       )
                              ),
                              tabPanel("Survival",
                                       fluidPage(
                                         fluidRow(
                                           column(
                                             width = 6,
                                             h4("Top Positive Regulators"),
                                             withSpinner(plotOutput("timing_survival_plot_pos", height = "500px", width = "100%"))
                                           ),
                                           column(
                                             width = 6,
                                             h4("Top Negative Regulators"),
                                             withSpinner(plotOutput("timing_survival_plot_neg", height = "500px", width = "100%"))
                                           )
                                         )
                                       )
                              ),
                              tabPanel("ICB Response",
                                       fluidRow(
                                         column(6,
                                                div(style = "margin-top: 0px; padding-top: 0px;",
                                                    h4("Top Positive Regulators", style = "margin-bottom: 5px;"),
                                                    plotOutput("forestplot_pos", height = "500px")
                                                )
                                         ),
                                         column(6,
                                                div(style = "margin-top: 0px; padding-top: 0px;",
                                                    h4("Top Negative Regulators", style = "margin-bottom: 5px;"),
                                                    plotOutput("forestplot_neg", height = "500px")
                                                )
                                         )
                                       )
                              ),
                              tabPanel("ICB Survival",
                                       fluidPage(
                                         fluidRow(
                                           column(12,
                                                  withSpinner(uiOutput("ICB_survival_plot_study"))
                                           )
                                         )
                                       )
                              )
                            )
                          )
                        )
                      )
             ),
             tabPanel("Help",
                      fluidPage(
                        h2("User Guide"),
                        p("This application helps explore immune escape mechanisms across cancer cohorts using multiple data modalities."),
                        
                        fluidRow(
                          column(
                            width = 3,
                            wellPanel(
                              h4("Sections"),
                              tags$ul(style = "list-style-type: none; padding-left: 0;",
                                      tags$li(tags$a(href = "#pathways", "1. Immunomodulatory Pathways")),
                                      tags$li(tags$a(href = "#cohort", "2. Cohort Summary")),
                                      tags$li(tags$a(href = "#genes", "3. Exploring Your Genes")),
                                      tags$li(tags$a(href = "#study", "4. Exploring Your Study")),
                                      tags$li(tags$a(href = "#tips", "4. Additional Tips")),
                                      tags$li(tags$a(href = "#troubleshooting", "5. Troubleshooting"))
                              )
                            )
                          ),
                          
                          column(
                            width = 9,
                            div(id = "pathways",
                                h3("1. Immunomodulatory Pathways"),
                                
                                # Subtitle
                                h4("Introduction"),
                                p("This section describes the key pathways involved in immune modulation, including their roles in tumor-immune interactions."),
                                
                                # Another subtitle
                                h4("Parameters"),
                                tags$ul(
                                  tags$li("Immunity Scenarios: This section provide the different categorites for immunomodualtory pathways in cancers, which were identified CRSIPR screen.
                                          MHC-I regulators:
                                          T-cell Killing:
                                          NK-cell Killing:
                                          gdT-cell Killing: 
                                          Macrophage Killing . "),
                                  tags$li("Select Regulator Type: Click the Submit button selecting direction for the immunomodualtory regulations.
                                          Positive regulation:
                                          Negative regulation: "),
                                  tags$li("Number of Top Pathway: Top freqnuent pathways will be shown by selecting the numbers.
                                          Positive regulation:
                                          Negative regulation: ")
                                ),
                                h4("Results"),
                                tags$ul(
                                  tags$li("Go to the Introduction tab to read about the study background."),
                                  tags$li("Click the Submit button after selecting your cohort to load data.")
                                ),
                            ),
                            tags$hr(),
                            
                            div(id = "cohort",
                                h3("2. Cohort Summary"),
                                tags$ul(
                                  tags$li("Go to the Introduction tab to read about the study background."),
                                  tags$li("Click the Submit button after selecting your cohort to load data.")
                                )
                            ),
                            tags$hr(),
                            
                            div(id = "genes",
                                h3("2. Exploring Your Genes"),
                                tags$ul(
                                  tags$li("Use the checkbox to enter a custom list of gene symbols."),
                                  tags$li("The Study Frequency tab shows CRISPR-based heatmaps."),
                                  tags$li("Mutation Frequency shows gene alterations across samples."),
                                  tags$li("ICB Response provides forest plots showing associations."),
                                  tags$li("Mutation Timing shows when mutations occur during evolution.")
                                )
                            ),
                            tags$hr(),
                            
                            div(id = "study",
                                h3("3. Exploring Your Study"),
                                tags$ul(
                                  tags$li("Upload your gene-level result file in .txt format."),
                                  tags$li("Use the sliders to explore top positive/negative regulators."),
                                  tags$li("Compare mutation patterns and timing between groups.")
                                )
                            ),
                            tags$hr(),
                            
                            div(id = "tips",
                                h3("4. Additional Tips"),
                                tags$ul(
                                  tags$li("Use the dropdowns to filter by cancer type or cell type."),
                                  tags$li("Click 'Submit' after selecting new cohorts to refresh the data."),
                                  tags$li("Hover over plots to see details if tooltips are enabled.")
                                )
                            ),
                            tags$hr(),
                            
                            div(id = "troubleshooting",
                                h3("5. Troubleshooting"),
                                tags$ul(
                                  tags$li("If plots do not appear, make sure the data has finished loading."),
                                  tags$li("Use a supported browser like Chrome or Firefox."),
                                  tags$li("For questions, contact: your.email@domain.com")
                                )
                            )
                          )
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
  filtered_data <- reactive({
    CRISPR_all %>% filter(Celltype == cell_type_files[[input$cell_type]])
  })
  
  study_count <- reactive({
    req(input$cell_type, input$immuno_regulators)
    cell_type_code <- cell_type_files[[input$cell_type]]
    reg <- input$immuno_regulators
    
    data_crispr <- data_crispr %>%
      filter(celltype == !!cell_type_code, regulators == !!reg)
    
    length(unique(data_crispr$Study))
  })
  
  output$immuno_pathway_plot_ui <- renderUI({
    req(input$immuno_n)
    
    base_width <- 600       # minimum width
    extra_width <- 20       # pixels per additional study
    width_cap <- 1200       # maximum width
    
    base_height <- 450      # base height for a small number of pathways
    extra_height <- 30      # height per additional pathway
    height_cap <- 1200      # max height
    
    # Compute width
    study_n <- study_count()
    total_width <- min(base_width + extra_width * study_n, width_cap)
    
    # Compute height
    pathway_n <- input$immuno_n
    total_height <- min(base_height + extra_height * pathway_n, height_cap)
    
    plotlyOutput("immuno_pathway_plot",
                 height = paste0(total_height, "px"),
                 width = paste0(total_width, "px"))
  })
  
  output$immuno_pathway_plot <- renderPlotly({
    req(input$cell_type, input$immuno_regulators, input$immuno_n)
    
    cell_type_code <- cell_type_files[[input$cell_type]]
    reg <- input$immuno_regulators
    num <- input$immuno_n
    textcol <- "black"
    
    # Read data
    data_crispr <- fread("/Users/wchen20/Desktop/Website/Shiny/MAF/data/CRISPR/C5BP_all.csv") %>%
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
  
  observeEvent(input$submit_btn, {
    loading_flag(TRUE)
    
    withProgress(message = "Loading data...", value = 0.1, {
      data <- dynamic_cohort_data(input$cohort_select)
      incProgress(0.4, detail = "Processing data...")
      
      cohort_data(data)
      
      type_list <- data$LOH_MSI %>% distinct(aliquot_id, histology_abbreviation)
      updateSelectInput(session, "cancer_select", choices = c("All", unique(type_list$histology_abbreviation)))
      updateSelectInput(session, "timing_type", choices = unique(data$diff_all$histology_abbreviation))
      
      incProgress(0.5, detail = "Finalizing...")
      Sys.sleep(0.5)
    })
    
    loading_flag(FALSE)
  })
  
  ###################################################
  ## Page 3 -- Cohort Summary
  ###################################################
  output$cohort_summary_left <- renderUI({
    tagList(
      plotOutput("histology_bar", height = "300px"),
      plotOutput("total_perMB_log_box", height = "300px")
    )
  })
  
  output$cohort_summary_right <- renderUI({
    tagList(
      plotOutput("msi_pie", height = "300px"),
      plotOutput("primary_pie", height = "300px")
    )
  })
  
  LOH_MSI_filtered <- reactive({
    req(cohort_data())
    LOH_MSI <- cohort_data()$LOH_MSI
    
    hist_count <- LOH_MSI %>%
      count(histology_abbreviation) %>%
      filter(n >= 5)
    
    LOH_MSI %>%
      filter(histology_abbreviation %in% hist_count$histology_abbreviation)
  })
  
  output$histology_bar <- renderPlot({
    req(cohort_data())
    LOH_MSI <- cohort_data()$LOH_MSI
    
    LOH_MSI_filtered <- LOH_MSI_filtered()
    
    hist_order <- LOH_MSI_filtered %>%
      count(histology_abbreviation) %>%
      arrange(desc(n)) %>%
      pull(histology_abbreviation)
    
    ggplot(LOH_MSI_filtered, aes(x = factor(histology_abbreviation, levels = hist_order), fill = histology_abbreviation)) +
      geom_bar() +
      scale_fill_manual(values = type_colors) +
      theme_minimal(base_size = 14) +
      labs(title = "Sample Frequency", x = "", y = "Sample Count") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none")
  })
  
  output$primary_pie <- renderPlot({
    req(cohort_data())
    LOH_MSI_filtered <- LOH_MSI_filtered()
    
    data <- LOH_MSI_filtered %>% dplyr::count(Sample_Type)
    ggplot(data, aes(x = "", y = n, fill = Sample_Type)) +
      #scale_fill_manual(values = prim_colors) +
      geom_col() +
      coord_polar("y") +
      theme_void(base_size = 14) +
      labs(title = "Sample Type") + scale_fill_npg()
  })
  
  output$msi_pie <- renderPlot({
    req(cohort_data())
    LOH_MSI_filtered = LOH_MSI_filtered()
    
    data <- LOH_MSI_filtered %>% dplyr::count(Gender)
    ggplot(data, aes(x = "", y = n, fill = Gender)) +
      geom_col() +
      scale_fill_manual(values = Gender_colors) +
      coord_polar("y") +
      theme_void(base_size = 14) +
      labs(title = "Gender")
  })
  
  output$total_perMB_log_box <- renderPlot({
    req(cohort_data())
    LOH_MSI_filtered <- LOH_MSI_filtered()
    
    hist_order <- LOH_MSI_filtered %>%
      group_by(histology_abbreviation) %>%
      summarise(median_tmb = median(TMB_log, na.rm = TRUE)) %>%
      arrange(desc(median_tmb)) %>%
      pull(histology_abbreviation)
    
    ggplot(LOH_MSI_filtered, aes(x = factor(histology_abbreviation, levels = hist_order), y = TMB_log, fill = histology_abbreviation)) +
      geom_boxplot() +
      scale_fill_manual(values = type_colors) +
      theme_minimal(base_size = 14) +
      scale_y_continuous(breaks = log10(c(0.1, 1, 10, 100, 1000, 10000)), 
                         labels = c(0.1, 1, 10, 100, 1000, 10000)) +
      labs(title = "TMB", x = "", y = "TMB") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none")
  })
  
  ###################################################
  ## Page 4 
  ###################################################
  ## Mutation Frequency
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
      unique() %>%
      sort()
    
    selectInput("selected_pathways", "Select Pathways",
                choices = choices,
                selected = "Antigen Presentation via MHC-I",
                multiple = FALSE)
  })
  
  # Render top 10 cancer types by Freq_mut for selected pathways
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
  
  ## Page 4 -- Radar Plot
    output$cell_type_ui <- renderUI({
    if (input$crispr_tab_selected == "Timeline") {
      selectInput("cell_type", "Immunity Scenarios", choices = c("All", names(cell_type_files)))
    } else {
      selectInput("cell_type", "Immunity Scenarios", choices = names(cell_type_files))
    }
  })
    
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
  
  ## Page 4 -- Timeline
  output$timing_type_ui <- renderUI({
    selectInput("timing_type", "Select Cancer Type", choices = unique(cohort_data()$diff_all$histology_abbreviation), selected = "Breast-AdenoCA")
  })
  
  output$timing_plot <- renderPlot({
    req(input$timing_type, input$cell_type, cohort_data())
    
    # Determine selected cell types
    cell_types <- if (input$cell_type == "All") {
      unname(cell_type_files)
    } else {
      cell_type_files[[input$cell_type]]  
    }
    
    # Generate plot
    p <- plot_timing_summary(
      diff_all = cohort_data()$diff_all,
      driver_list = driver_list,
      LOH_MSI = cohort_data()$LOH_MSI,
      type = input$timing_type,
      cell_type = cell_types
    )
    
    # Return plot or fallback
    if (!is.null(p)) {
      return(p)
    } else {
      return(
        ggplot() + theme_void() +
          annotate("text", x = 0.5, y = 0.5, label = "No sufficient data", size = 8, hjust = 0.5)
      )
    }
  })
  
  ## Page 4 -- Survival
  output$select_histology_ui <- renderUI({
    req(cohort_data())
    histologies <- sort(unique(cohort_data()$diff_all$histology_abbreviation))
    selectInput("histology_type", "Select Cancer Type:", choices = c("All", histologies), selected = "All")
  })
  
  output$pathway_selector_surv <- renderUI({
    req(input$cell_type, input$selected_regulator, cohort_data())
    
    choices <- cohort_data()$diff_all %>%
      filter(celltype == !!cell_type_files[[input$cell_type]]) %>%
      pull(pathway) %>% unique()
    
    selectInput("selected_pathways", "Select Pathways",
                choices = choices,
                selected = "Antigen Presentation via MHC-I",
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
      
      sample_survival <- unique(data_survival$aliquot_id)
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
  
  ## Cibersort
  output$select_histology_ciber <- renderUI({
    req(cohort_data())
    histologies <- sort(unique(cohort_data()$diff_all$histology_abbreviation))
    selectInput("histology_ciber", "Select Cancer Type:", choices = c("All", histologies), selected = "All")
  })
  
  output$pathway_selector_ciber <- renderUI({
    req(input$cell_type, input$selected_regulator, cohort_data())
    
    choices <- cohort_data()$diff_all %>%
      filter(celltype == !!cell_type_files[[input$cell_type]]) %>%
      pull(pathway) %>% unique()
    
    selectInput("pathway_ciber", "Select Pathways",
                choices = choices,
                selected = "Antigen Presentation via MHC-I",
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
  ## Page 5 -- Mutation Frequency
  maf_obj <- reactive({
    req(cohort_data())
    cohort_data()$maf_data
  })
  
  get_selected_genes <- reactive({
    if (input$useGeneList && nzchar(input$geneList)) {
      unlist(strsplit(input$geneList, "[,\n]")) %>% trimws() %>% unique()
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
  
  ## Page 5 -- Mutation Timing
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
  
  output$select_histology_list <- renderUI({
    req(cohort_data())
    histologies <- sort(unique(cohort_data()$diff_all$histology_abbreviation))
    selectInput("histology_type_genelist", "Select Cancer Type:", choices = c("All", histologies), selected = "All")
  })
  
  output$timing_survival_plot_genelist <- renderPlot({
    tryCatch({
      
      req(get_selected_genes(), input$histology_type_genelist, cohort_data())
      
      LOH_type = cohort_data()$LOH_MSI %>% dplyr::select(aliquot_id, histology_abbreviation)
      
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
        filter(aliquot_id %in% unique(data_survival$aliquot_id))
      
      data_survival_sub <- if (histology_type == "All") data_survival else filter(data_survival, histology_abbreviation == histology_type)
      diff_all_sub <- if (histology_type == "All") diff_all else filter(diff_all, histology_abbreviation == histology_type)
      
      filtered_data <- diff_all_sub %>%
        left_join(data_survival_sub, by = "aliquot_id") %>%
        mutate(
          Timing_group = ifelse(is.na(Timing_group), "WT", Timing_group),
          Timing_group = factor(Timing_group, levels = c("WT", "Mut_Early", "Mut_Late"))
        )
      
      # Call the function here
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
  
  ## Page 5 -- ICB Response
  observe({
    req(ICB_merge)
    therapy_choices <- c("All", unique(ICB_merge$Therapy))
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
  
  ## Page 5 -- ICB Survival
  observe({
    req(ICB_surv)
    therapy_choices <- c("All", unique(ICB_surv$Therapy))
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
      filter(length(unique(OS_status)) > 1, length(unique(Mutated)) > 1) %>%
      ungroup()
    
    n_plots <- length(unique(sample_data$CohortName))
    
    # Calculate dynamic height
    base_height <- 350  # Height per row of plots
    plots_per_row <- 4  # Your ncol = 4
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
          labs(title = unique(data$CohortName))
        
        cowplot::plot_grid(
          surv_plot$plot,
          surv_plot$table,
          align = "v",
          ncol = 1,
          rel_heights = c(4, 1.5)
        )
      })
    
    final_plot <- cowplot::plot_grid(plotlist = km_plots, ncol = 4, align = "v")
    
    # Return the plot with dynamic height
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
  ## Page 6 -- CRISPR Ranking
  gstable_data <- reactive({
    req(input$file)
    
    # Read the uploaded file
    gdata <- read.table(input$file$datapath, header = TRUE, stringsAsFactors = FALSE)
    
    comm <- intersect(gdata$id, MouseGene)
    # Ensure 'id' column exists
    req("id" %in% colnames(gdata))
    
    # Attempt gene symbol conversion if ids look like gene symbols
    if (length(comm) > 50) {
      gdata <- left_join(gdata, MouseToHuman_gene, by = "id")
    } else {
      gdata$HumanGene <- gdata$id
    }
    
    # Filter out invalid gene names
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
  
  # CRISPR ranking plot
  output$crisprPlot <- renderPlot({
    req(input$top_n)
    gstable <- gstable_data()
    p1 <- plot_selection(gstable, "pos.rank", "pos.score", "Top Positive Regulators", input$top_n, color_offset = 0)
    p2 <- plot_selection(gstable, "neg.rank", "neg.score", "Top Negative Regulators", input$top_n, color_offset = 10)
    p1 + p2
  }, res = 80)
  
  # Page 6 -- Mutation frequency plots
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
  
  # Page 6 -- Raincloud plots for single gene timing
  output$timing_plot_singlegene_pos_ui <- renderUI({
    req(input$top_n)
    
    base_height <- 150
    extra_height <- 20  # raincloud plots often need more space
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
  
  # Page 6 -- Mutation Timing bar plots
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
  
  # Page 6 -- Survival by Mutation Timing
  output$select_histology_study <- renderUI({
    req(cohort_data())
    histologies <- cohort_data()$diff_all %>% pull(histology_abbreviation) %>% unique(.)
    selectInput("histology_type_study", "Select Cancer Type:", choices = c("All", histologies), selected = "All")
  })
  
  output$timing_survival_plot_pos <- renderPlot({
    tryCatch({
      
      req(cohort_data(), gstable_data(), input$top_n, input$histology_type_study)
      LOH_type = cohort_data()$LOH_MSI %>% dplyr::select(aliquot_id, histology_abbreviation)
      
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
        filter(aliquot_id %in% unique(data_survival$aliquot_id))
      
      data_survival_sub <- if (histology_type == "All") data_survival else filter(data_survival, histology_abbreviation == histology_type)
      diff_all_sub <- if (histology_type == "All") diff_all else filter(diff_all, histology_abbreviation == histology_type)
      
      filtered_data <- diff_all_sub %>%
        left_join(data_survival_sub, by = "aliquot_id") %>%
        mutate(
          Timing_group = ifelse(is.na(Timing_group), "WT", Timing_group),
          Timing_group = factor(Timing_group, levels = c("WT", "Mut_Early", "Mut_Late"))
        )
      
      # Call the function here
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
      LOH_type = cohort_data()$LOH_MSI %>% dplyr::select(aliquot_id, histology_abbreviation)
      
      neg_genes <- gstable_data() %>% arrange(neg.rank) %>% slice_head(n = input$top_n) %>% pull(HumanGene)
      
      selected_genes <- neg_genes
      histology_type <- input$histology_type_study
      
      # Prepare data
      df_Reg <- cohort_data()$diff_allgene %>%
        filter(Hugo_Symbol %in% selected_genes) %>%
        group_by(sample_id) %>%
        summarise(across(all_of(column_names), ~ mean(., na.rm = TRUE)), .groups = "drop") %>%
        mutate(Timing_group = case_when(early_ratio > 0.5 ~ "Mut_Early", 
                                        late_ratio > 0.5 ~ "Mut_Late",
                                        TRUE ~ "Undetermined")) %>% 
        filter(Timing_group != "Undetermined") %>% inner_join(., LOH_type, by = "aliquot_id")
      
      data_survival <- cohort_data()$surv_data
      diff_all <- df_Reg %>%
        filter(aliquot_id %in% unique(data_survival$aliquot_id))
      
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
  
  # Page 6 -- ICB Response forest plots
  observe({
    req(ICB_merge)
    therapy_choices <- c("All", unique(ICB_merge$Therapy))
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
  
  # KM Curve for Top Positive Regulators
  observe({
    req(ICB_surv)
    therapy_choices <- c("All", unique(ICB_surv$Therapy))
    updateSelectInput(session, "therapy_surv_study",
                      choices = therapy_choices,
                      selected = "All")
  })
  
  output$ICB_survival_plot_study <- renderUI({
    req(cohort_data(), input$therapy_surv_study, input$regulator_type)
    
    # Get top genes based on regulator type
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
      filter(length(unique(OS_status)) > 1, length(unique(Mutated)) > 1) %>%
      ungroup()
    
    # Calculate number of plots and dynamic height
    n_plots <- length(unique(sample_data$CohortName))
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
          labs(title = unique(data$CohortName))
        
        cowplot::plot_grid(
          surv_plot$plot,
          surv_plot$table,
          align = "v",
          ncol = 1,
          rel_heights = c(4, 1.5)
        )
      })
    
    final_plot <- cowplot::plot_grid(plotlist = km_plots, ncol = 4, align = "v")
    
    # Return the plot with dynamic height
    renderPlot({
      final_plot
    }, height = plot_height, res = 70)
  })
  
}

# Run app
shinyApp(ui, server)
