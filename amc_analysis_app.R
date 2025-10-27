library(shiny)
library(rhandsontable)
library(writexl)

source("amc_scripts/f1.R")


ui <- fluidPage(
  theme = shinytheme("cerulean"),

  titlePanel("AMC Data Analysis - MAAP"),

  tabsetPanel(
    id = "steps",

    # Step 1
    tabPanel("Step 1",
             h3("Select AMC Variables"),

             textInput("country_name", "Enter Country Name"),
             numericInput("population", "Enter Population", value = NA, min = 1),
             br(),
             actionButton("reg_1", "Register Country & Population"),
             br(), br(),
             verbatimTextOutput("register_msg"),
             br(),

             rHandsontableOutput("table_1"),
             br(),
             downloadButton("download_1", "Save Data"),
             helpText(paste0("Save in ", amc_updates_dir,"/")),

             br(),br(), br(),
             actionButton("run_script_1", "Initiate Look up of unclear entries"),
             br(),br(), br(),
             verbatimTextOutput("console_1"),

             checkboxInput("completed_1", "I have completed this step"),
             br(),
             # Bottom-left: Previous

             # Bottom-right: Next
             fixedPanel(
               actionButton("next_1", "Next"),
               bottom = 10, right = 10, width = "auto"
             )
    ),

    # Step 2
    tabPanel("Step 2",
             h3("Validate unclear entries"),

             rHandsontableOutput("table_2"),
             br(),
             downloadButton("download_2", "Save Data"),
             helpText(paste0("Save in ", amc_updates_dir,"/")),

             br(),br(),
             actionButton("run_script_2", "Cleanup AMC dataset"),
             br(),br(), br(),
             verbatimTextOutput("console_2"),

             checkboxInput("completed_2", "I have completed this step"),
             br(),br(),
             # Bottom-left: Previous
             fixedPanel(
               actionButton("prev_2", "Previous"),
               bottom = 10, left = 10, width = "auto"
             ),

             # Bottom-right: Next
             fixedPanel(
               actionButton("next_2", "Next"),
               bottom = 10, right = 10, width = "auto"
             )
    ),

    # Step 3
    tabPanel("Step 3",
             h4("Update DDD Information"),

             rHandsontableOutput("table_3"),
             br(),
             downloadButton("download_3", "Save Data"),
             helpText(paste0("Save in ", amc_updates_dir,"/")),
             br(),br(),

             h4("Match strength units"),

             rHandsontableOutput("table_3a"),
             br(),
             downloadButton("download_3a", "Save Data"),
             helpText(paste0("Save in ", amc_updates_dir,"/")),


             br(),br(),
             actionButton("run_script_3", "Subset the ineligible data"),
             br(),br(),
             verbatimTextOutput("console_3"),


             checkboxInput("completed_3", "I have completed this step"),
             br(),br(),
             # Bottom-left: Previous
             fixedPanel(
               actionButton("prev_3", "Previous"),
               bottom = 10, left = 10, width = "auto"
             ),

             # Bottom-right: Next
             fixedPanel(
               actionButton("next_3", "Next"),
               bottom = 10, right = 10, width = "auto"
             )
    ),

    # Step 4
    tabPanel("Step 4",
             h3("Begin analysis"),

             rHandsontableOutput("table_4"),
             br(),
          #   downloadButton("download_4", "Download Modified Data"),
             #br(),br(),
             actionButton("run_script_4", "Begin analysis on the processed files"),
             br(),br(),
            verbatimTextOutput("console_4"),

             checkboxInput("completed_4", "I have completed this step"),
          br(),
          # Bottom-left: Previous
          fixedPanel(
            actionButton("prev_4", "Previous"),
            bottom = 10, left = 10, width = "auto"
          ),

          # Bottom-right: Next
          fixedPanel(
            actionButton("next_4", "Next"),
            bottom = 10, right = 10, width = "auto"
          )
    ),

    # Step 5
    tabPanel("Step 5",
             h3("Update class information"),
             helpText(paste0("Please provide antibiotic classes based on column 1")),

             rHandsontableOutput("table_5"),
             br(),
             downloadButton("download_5", "Save Data"),
             helpText(paste0("Save in ", amc_updates_dir,"/")),

             br(),br(),
             actionButton("run_script_5", "Finalize plots"),
             br(),br(),
             verbatimTextOutput("console_5"),
             br(),

             # Bottom-left: Previous
             fixedPanel(
               actionButton("prev_5", "Previous"),
               bottom = 10, left = 10, width = "auto"
             )

             #,
            # checkboxInput("completed_5", "I have completed this step"),
            # actionButton("next_5", "Next")
    )#,

    # # Step 6
    # tabPanel("Step 6",
    #          h3("Step 6"),
    #          actionButton("run_script_6", "Run Script"),
    #          verbatimTextOutput("console_6"),
    #          rHandsontableOutput("table_6"),
    #          br(),
    #          downloadButton("download_6", "Download Modified Data"),
    #          br(),
    #          actionButton("prev_6", "Previous"),
    #          checkboxInput("completed_6", "I have completed this step"),
    #          actionButton("next_6", "Next")
    # ),
    #
    # # Step 7
    # tabPanel("Step 7",
    #          h3("Step 7"),
    #          actionButton("run_script_7", "Run Script"),
    #          verbatimTextOutput("console_7"),
    #          rHandsontableOutput("table_7"),
    #          br(),
    #          downloadButton("download_7", "Download Modified Data"),
    #          br(),
    #          actionButton("prev_7", "Previous"),
    #          checkboxInput("completed_7", "I have completed this step"),
    #          actionButton("next_7", "Next")
    # ),
    #
    # # Step 8
    # tabPanel("Step 8",
    #          h3("Step 8"),
    #          actionButton("run_script_8", "Run Script"),
    #          verbatimTextOutput("console_8"),
    #          rHandsontableOutput("table_8"),
    #          br(),
    #          downloadButton("download_8", "Download Modified Data"),
    #          br(),
    #          actionButton("prev_8", "Previous"),
    #          checkboxInput("completed_8", "I have completed this step")
    # )
  )
)



server <- function(input, output, session) {
  # Step datasets (initially NULL except step1)
  step1_data <- reactiveVal(empty_amc_df)
  step2_data <- reactiveVal(NULL)
  step3_data <- reactiveVal(NULL)
 #step3_data <- reactiveVal(list( df1=ddd_updates,df2=units_updates))
  step4_data <- reactiveVal(NULL)
  step5_data <- reactiveVal(NULL)
  step6_data <- reactiveVal(NULL)
  step7_data <- reactiveVal(NULL)
  step8_data <- reactiveVal(NULL)

  # Step logs
  step_logs <- lapply(1:8, function(i) reactiveVal(""))

  # ---- Step 1 ----
  output$table_1 <- renderRHandsontable({
    df <- step1_data()
    req(df)
    rhandsontable(df)%>%
      hot_col("Corresponding_variables", type = "dropdown", source = choices1, width = 300) %>%
      hot_col("Required_variables", readOnly = TRUE, width = 300)  # Optional: Make label column readonly
  })

  observeEvent(input$reg_1, {
    req(input$country_name, input$population)

       # Save to global environment
    assign("cntry", input$country_name, envir = .GlobalEnv)
    assign("pop", input$population, envir = .GlobalEnv)

    # Feedback to user
    output$register_msg <- renderText({
      paste0("✅ Registered: ", input$country_name,
             " with population ", input$population)
    })
  })


  observe({ req(input$table_1); step1_data(hot_to_r(input$table_1)) })
  observeEvent(input$run_script_1, {
    df <- step1_data()
    script_file <- "amc_scripts/f2.R"
    if(file.exists(script_file)) {
      msg <- capture.output(tryCatch(source(script_file, local = .GlobalEnv),
                                     error = function(e) cat("Error:", e$message)), type = "output")
      step_logs[[1]](paste(msg, collapse="\n"))
    } else step_logs[[1]]("No script found for Step 1")
    step1_data(df)
  })
  output$console_1 <- renderText({ step_logs[[1]]() })
  output$download_1 <- downloadHandler(filename = "select_amc_variables.xlsx",
                                       content = function(file) writexl::write_xlsx(step1_data(), file))

  # ---- Step 2 lazy-load ----
  observeEvent(input$next_1, {
    if(isTRUE(input$completed_1)) {
      updateTabsetPanel(session, "steps", "Step 2")
      # Load step2 data only now
      step2_data(lookup_df[,-3])
    }
  })
  output$table_2 <- renderRHandsontable({
    df <- step2_data()
    req(df)
    rhandsontable(df) %>%
      hot_col("Verdict", type = "dropdown", source = c('','Correct', 'I have edited'), width = 150) %>%
      hot_col("original_entry", readOnly = TRUE, width = 300)  # Optional: Make label column readonly
  })
  observe({ req(input$table_2); step2_data(hot_to_r(input$table_2)) })
  observeEvent(input$run_script_2, {
    df <- step2_data()
    script_file <- "amc_scripts/f3.R"
    if(file.exists(script_file)) {
      msg <- capture.output(tryCatch(source(script_file, local = .GlobalEnv),
                                     error = function(e) cat("Error:", e$message)), type = "output")
      step_logs[[2]](paste(msg, collapse="\n"))
    } else step_logs[[2]]("No script found for Step 2")
    step2_data(df)
  })
  output$console_2 <- renderText({ step_logs[[2]]() })
  output$download_2 <- downloadHandler(filename = "matching_unclear_antibiotic_entries.xlsx",
                                       content = function(file) writexl::write_xlsx(step2_data(), file))


  # ---- Step 3 lazy-load ----
  # small helper to update one element of the list
  set_step3 <- function(name, value) {
    x <- step3_data()
    x[[name]] <- value
    step3_data(x)
  }

  #to run once the step is opened
  make_step3 <- function() {
       list(
         df1=ddd_updates,df2=units_updates
    )
  }

  observeEvent(input$next_2, {
    if(isTRUE(input$completed_2)) {
      updateTabsetPanel(session, "steps", "Step 3")

      if (is.null(step3_data())) step3_data(make_step3())
    }
  })

  #df1 data

  output$table_3 <- renderRHandsontable({

    x <- step3_data(); req(!is.null(x))
    df <- x$df1;        req(!is.null(df))

    tbl <- rhandsontable(df)

    # guard hot_col calls so there is no error if cols are missing
    if ("atc_level_name" %in% names(df)) {
      tbl <- hot_col(tbl, "atc_level_name", readOnly = TRUE, width = 300)
    }
    tbl

  })


  # capture edits for df1
  observeEvent(input$table_3, ignoreInit = TRUE, {
    set_step3("df1", hot_to_r(input$table_3))
  })

# observe({ req(input$table_3); step3_data(hot_to_r(input$table_3)) })
  output$download_3 <- downloadHandler(
    filename = "DDD_information_updates.xlsx",
    content = function(file) {

    x <- step3_data(); req(!is.null(x)); req(!is.null(x$df1))
    writexl::write_xlsx(x$df1, path = file)
    }
  )

  ##Df2
  output$table_3a <- renderRHandsontable({
    x <- step3_data(); req(!is.null(x))
    df <- x$df2;        req(!is.null(df))

    tbl <- rhandsontable(df)

    # guard hot_col calls so there is no error if cols are missing
    if ("strength_unit_in_dataset" %in% names(df)) {
      tbl <- hot_col(tbl, "strength_unit_in_dataset", readOnly = TRUE, width = 300)
    }

    if ("standardized_units" %in% names(df)) {
      tbl <- hot_col(tbl, "standardized_units", type = "dropdown", source = units_options, width = 200)
    }
    tbl

  })


# capture edits for df1
observeEvent(input$table_3a, ignoreInit = TRUE, {
  set_step3("df2", hot_to_r(input$table_3a))
})

# observe({ req(input$table_3); step3_data(hot_to_r(input$table_3)) })
output$download_3a <- downloadHandler(
  filename = "strength_units_updates.xlsx",
  content = function(file) {

    x <- step3_data(); req(!is.null(x)); req(!is.null(x$df1))
    writexl::write_xlsx(x$df2, path = file)
  }
)




  observeEvent(input$run_script_3, {
    df <- step3_data()
    script_file <- "amc_scripts/f4.R"
    if(file.exists(script_file)) {
      msg <- capture.output(tryCatch(source(script_file, local = .GlobalEnv),
                                     error = function(e) cat("Error:", e$message)), type = "output")
      step_logs[[3]](paste(msg, collapse="\n"))
    } else step_logs[[3]]("No script found for Step 3")
    step3_data(df)
  })
  output$console_3 <- renderText({ step_logs[[3]]() })


  # ---- Step 4 lazy-load ----
  observeEvent(input$next_3, {
    if(isTRUE(input$completed_3)) {
      updateTabsetPanel(session, "steps", "Step 4")
      # Load step4 data only now
      step4_data()
    }
  })
  output$table_4 <- renderRHandsontable({
    df <- step4_data()
    req(df)
    rhandsontable(df)
  })
  observe({ req(input$table_4); step4_data(hot_to_r(input$table_4)) })
  observeEvent(input$run_script_4, {
    df <- step4_data()
    script_file <- "amc_scripts/f5.R"
    if(file.exists(script_file)) {
      msg <- capture.output(tryCatch(source(script_file, local = .GlobalEnv),
                                     error = function(e) cat("Error:", e$message)), type = "output")
      step_logs[[4]](paste(msg, collapse="\n"))
    } else step_logs[[4]]("No script found for Step 4")
    step4_data(df)
  })
  output$console_4 <- renderText({ step_logs[[4]]() })
  # output$download_4 <- downloadHandler(filename = "matching_unclear_antibiotic_entries.xlsx",
  #                                      content = function(file) writexl::write_xlsx(step4_data(), file))

  # ---- Step 5 lazy-load ----
  observeEvent(input$next_4, {
    if(isTRUE(input$completed_4)) {
      updateTabsetPanel(session, "steps", "Step 5")
      # Load step5 data only now
      step5_data(unclassified_abs)
    }
  })
  output$table_5 <- renderRHandsontable({
    df <- step5_data()
    req(df)
    rhandsontable(df) %>%
      hot_col("Class", type = "dropdown", source = c(' ',sort(antibiotic_classes_amc)), width = 200) %>%
      hot_col("Category", type = "dropdown", source = c(' ','Access','Watch', 'Reserve','Uncategorized'), width = 150) %>%
      hot_col("antibiotic_names", readOnly = TRUE, width = 300)  # Optional: Make label column readonly
  })
  observe({ req(input$table_5); step5_data(hot_to_r(input$table_5)) })
  observeEvent(input$run_script_5, {
    df <- step5_data()
    script_file <- "amc_scripts/f6.R"
    if(file.exists(script_file)) {
      msg <- capture.output(tryCatch(source(script_file, local = .GlobalEnv),
                                     error = function(e) cat("Error:", e$message)), type = "output")
      step_logs[[5]](paste(msg, collapse="\n"))
    } else step_logs[[5]]("No script found for Step 5")
    step5_data(df)
  })
  output$console_5 <- renderText({ step_logs[[5]]() })
  output$download_5 <- downloadHandler(filename = "updated_AMC_classes.xlsx",
                                       content = function(file) writexl::write_xlsx(step5_data(), file))


  # Repeat for Step 3 to Step 8: lazy-load dataset only when next button is clicked
  #observeEvent(input$next_2, { if(isTRUE(input$completed_2)) { updateTabsetPanel(session, "steps", "Step 3"); step3_data(head(airquality,5)) }})
 # observeEvent(input$next_3, { if(isTRUE(input$completed_3)) { updateTabsetPanel(session, "steps", "Step 4"); step4_data(head(PlantGrowth,5)) }})


  # Previous buttons
  observeEvent(input$prev_2, { updateTabsetPanel(session, "steps", "Step 1") })
  observeEvent(input$prev_3, { updateTabsetPanel(session, "steps", "Step 2") })
  observeEvent(input$prev_4, { updateTabsetPanel(session, "steps", "Step 3") })
  observeEvent(input$prev_5, { updateTabsetPanel(session, "steps", "Step 4") })
  observeEvent(input$prev_6, { updateTabsetPanel(session, "steps", "Step 5") })
  observeEvent(input$prev_7, { updateTabsetPanel(session, "steps", "Step 6") })
  #observeEvent(input$prev_8, { updateTabsetPanel(session, "steps", "Step 7") })

  # ---- Step 3-8 handsontable and scripts (same as Step 2) ----
  # Can replicate the Step 2 code for each step (table render, editing, run_script, download, console)
}

shinyApp(ui, server)
