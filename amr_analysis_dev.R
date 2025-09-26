# library(shiny)
# library(rhandsontable)
# library(writexl)

source("amr_scripts/f1.R")


ui <- fluidPage(

  theme = shinytheme("cerulean"),

  titlePanel("AMR Data Analysis - MAAP"),

  tabsetPanel(
    id = "steps",

    # Step 1
    tabPanel("Step 1",
             div(
               style = "position: relative; min-height: 600px;",  # container with enough space

               # 🔹 faint logo
               tags$img(
                 src = "logo.jpg",   # put logo.png in your www/ folder
                 style = "
            opacity: 0.07;            /* faint transparency */
            position: absolute;
            top: 3%; left: 45%;     /* adjust placement */
            width: 500px;            /* adjust size */
            z-index: 0;
          "
               ),


               div(
                 style = "position: relative; z-index: 1;",  # ensure content is above the logo

                 h3("Select AMR Variables"),

                 textInput("country_name", "Enter Country Name"),
                 br(),
                 actionButton("reg_1", "Register Country"),
                 br(), br(),
                 verbatimTextOutput("register_msg"),
                 br(),

                 selectInput("os_type", "Select your Operating System:",
                             choices = c("Windows", "Mac", "Linux", "Other")),
                 br(),
                 textOutput("os_msg"),

                 br(),

                 rHandsontableOutput("table_1"),
                 br(),
                 downloadButton("download_1", "Save Data"),
                 helpText(paste0("Save in ", amr_updates_dir,"/")),

                 br(),




                 #actionButton("run_script_1b", "Update Patient location type options"),
                 #br(), br(), br(),
                 #verbatimTextOutput("console_1b"),


                 checkboxInput("completed_1", "I have completed this step"),
                 br(),

                 # Bottom-right: Next
                 fixedPanel(
                   actionButton("next_1", "Next"),
                   bottom = 10, right = 10, width = "auto"
                 )

               )
             )
    ),
    # Step 2
    tabPanel("Step 2",
             checkboxInput("data_format_long", "My dataset have antibiotics in one column"),
             helpText(paste0("Leave blank if antibiotics form multiple columns in the dataset")),

             conditionalPanel(
               condition = "input.data_format_long == true",
               rHandsontableOutput("table_2"),
               downloadButton("download_2", "Save Parameters"),
               helpText(paste0("Save in ", amr_updates_dir,"/"))
             ),
             br(), br(), br(),

             h4("Select Antibiotic columns"),

             rHandsontableOutput("table_2a"),
             br(),
             downloadButton("download_2a", "Save Data"),
             helpText(paste0("Save in ", amr_updates_dir,"/")),
             br(), br(), br(),

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

    #step_3
    tabPanel("Step 3",
             h4("SIR interpretations"),

             actionButton("run_script_3", "SIR Interpretations and Look up tables"),
             br(), br(),
             verbatimTextOutput("console_3"),
             br(),
             #checkboxInput("sir", "SIR interpretations complete"),
             #br(),

             br(),br(), br(),
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

    #Step4
    tabPanel("Step 4",
             #run script 1b
             h4("Update patient location entries and Begin Analysis"),
             rHandsontableOutput("table_4"),
             br(),

             downloadButton("download_4", "Save Data"),
             helpText(paste0("Save in ", amr_updates_dir,"/")),

             actionButton("run_script_4", "Begin analysis"),
             br(), br(), br(),
             verbatimTextOutput("console_4"),

             verbatimTextOutput("console_4"),

             checkboxInput("completed_4", "I have completed this step"),
             br(), br(),
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
    #step_4
    tabPanel("Step 5",
             h4("GLASS specimen and analysis"),

             rHandsontableOutput("table_5"),
             br(),
             downloadButton("download_5", "Save Data"),
             helpText(paste0("Save in ", amr_updates_dir,"/")),

             br(),br(),
             actionButton("run_script_5", "Begin analysis for WHO GLASS combos"),
             br(),br(), br(),
             verbatimTextOutput("console_5"),

             checkboxInput("completed_5", "I have completed this step"),

             br(),br(),
             # Bottom-left: Previous
             fixedPanel(
               actionButton("prev_5", "Previous"),
               bottom = 10, left = 10, width = "auto"
             )
    )
    #
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
  step1_data <- reactiveVal(initial_df)
  step2_data <- reactiveVal(list( df1=long_df_cols,df2=ab_cols_user))
  step3_data <- reactiveVal(NULL)

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
      hot_col("my_dataset", type = "dropdown", source = choices1, width = 300) %>%
      hot_col("man_vars", readOnly = TRUE, width = 300)  # Optional: Make label column readonly
  })

  observeEvent(input$reg_1, {
    req(input$country_name)

    # Save to global environment
    assign("cntry", input$country_name, envir = .GlobalEnv)
    #assign("pop", input$population, envir = .GlobalEnv)

    # Feedback to user
    output$register_msg <- renderText({
      paste0("✅ R Country Successfully Registered: ", input$country_name)
    })
  })

  #operating system for antibiograms
  observeEvent(input$os_type, {
    req(input$os_type)  # ensure user has selected something

    # Save to global environment
    assign("user_os", input$os_type, envir = .GlobalEnv)

    # Feedback to user
    output$os_msg <- renderText({
      paste0("✅ Operating System Selected: ", input$os_type)
    })
  })

  observe({ req(input$table_1); step1_data(hot_to_r(input$table_1)) })

  output$download_1 <- downloadHandler(filename = "select_amr_variables.xlsx",
                                       content = function(file) writexl::write_xlsx(step1_data(), file))



  #---Step 2
  # small helper to update one element of the list
  set_step2 <- function(name, value) {
    x <- step2_data()
    x[[name]] <- value
    step2_data(x)
  }

  observeEvent(input$next_1, {
    if (isTRUE(input$completed_1)) {
      updateTabsetPanel(session, "steps", "Step 2")

    }
  })

  # ---- Table for df1 ----
  observeEvent(input$data_format_long, {
    if (isTRUE(input$data_format_long)) {
      output$table_2 <- renderRHandsontable({
        df <- step2_data()$df1
        req(df)
        tbl <- rhandsontable(df)

        # guard hot_col calls so there is no error if cols are missing
        if ("my_dataset" %in% names(df)) {
          tbl <- hot_col(tbl, "my_dataset", type = "dropdown", source = choices1, width = 200)
        }
        if ("man_vars" %in% names(df)) {
          tbl <- hot_col(tbl, "man_vars", readOnly = TRUE, width = 200)
        }
        tbl
      })
    } else {
      output$table_2 <- renderRHandsontable(NULL)
    }
  })

  # capture edits for df1
  observeEvent(input$table_2, ignoreInit = TRUE, {
    set_step2("df1", hot_to_r(input$table_2))
  })

  output$download_2 <- downloadHandler(
    filename = "long_format_amr_columns.xlsx",
    content = function(file) writexl::write_xlsx(step2_data()$df1, path = file)
  )

  # ---- Table for df2 ----
  output$table_2a <- renderRHandsontable({
    df <- step2_data()$df2
    req(df)
    tbl <- rhandsontable(df)
    if ("my_dataset" %in% names(df)) {
      tbl <- hot_col(tbl, "my_dataset", readOnly = TRUE, width = 200)
    }
    if ("antibiotic_column" %in% names(df)) {
      tbl <- hot_col(tbl, "antibiotic_column", type = "dropdown", source = c(" ", "Yes"), width = 200)
    }
    tbl
  })

  # capture edits for df2
  observeEvent(input$table_2a, ignoreInit = TRUE, {
    set_step2("df2", hot_to_r(input$table_2a))
  })

  output$download_2a <- downloadHandler(
    filename = "present_antibiotic_columns.xlsx",
    content = function(file) writexl::write_xlsx(step2_data()$df2, path = file)
  )



  #---Step 3
  observeEvent(input$next_2, {
    if(isTRUE(input$completed_2)) {
      updateTabsetPanel(session, "steps", "Step 3")
      step3_data()  # Load only now
    }
  })

  observeEvent(input$run_script_3, {
    df <- step3_data()
    script_file <- "amr_scripts/f2.R"
    if(file.exists(script_file)) {
      msg <- capture.output(tryCatch(source(script_file, local = .GlobalEnv),
                                     error = function(e) cat("Error:", e$message)), type = "output")
      step_logs[[3]](paste(msg, collapse="\n"))
    } else step_logs[[3]]("No script found for Step 3")
    step3_data(df)
  })
  output$console_3 <- renderText({ step_logs[[3]]() })




  #---Step 4
  ##the spec data

  observeEvent(input$next_3, {
    if(isTRUE(input$completed_3)) {
      updateTabsetPanel(session, "steps", "Step 4")
      step4_data(loc_options)  # Load only now
    }
  })


  output$table_4 <- renderRHandsontable({
    df <- step4_data()
    req(df)
    rhandsontable(df) %>%
      hot_col("my_dataset", readOnly = TRUE, width = 200) %>%
      hot_col("options", type = "dropdown", source = c(loc_opts,'Not available'), width = 200)
  })


  # Capture user edits from rhandsontable
  observe({
    if (!is.null(input$table_4)) {
      step4_data(hot_to_r(input$table_4))
    }
  })

  output$download_4 <- downloadHandler(filename = "patient_location_type.xlsx",
                                       content = function(file) writexl::write_xlsx(step4_data(), file))


  observeEvent(input$run_script_4, {
    df <- step4_data()
    script_file <- "amr_scripts/f2a.R"
    if(file.exists(script_file)) {
      msg <- capture.output(tryCatch(source(script_file, local = .GlobalEnv),
                                     error = function(e) cat("Error:", e$message)), type = "output")
      step_logs[[4]](paste(msg, collapse="\n"))
    } else step_logs[[4]]("No script found for Step 4")
    step3_data(df)
  })
  output$console_4 <- renderText({ step_logs[[4]]() })





  # # ---- Step 5 lazy-load ----
  observeEvent(input$next_4, {
    if(isTRUE(input$completed_4)) {
      updateTabsetPanel(session, "steps", "Step 5")
      #Load step2 data only now
      step5_data(specimen_check_df)
    }
  })
  output$table_5 <- renderRHandsontable({
    df <- step5_data()
    req(df)
    rhandsontable(df) %>%
      hot_col("matched", type = "dropdown", source = c(sort(unique(glass_opts$specimen)),''), width = 150) %>%
      hot_col("my_dataset", readOnly = TRUE, width = 300)  # Optional: Make label column readonly
  })
  observe({ req(input$table_5); step5_data(hot_to_r(input$table_5)) })
  observeEvent(input$run_script_5, {
    df <- step5_data()
    script_file <- "amr_scripts/f3.R"  #bring in the corrections and do the conversions
    if(file.exists(script_file)) {
      msg <- capture.output(tryCatch(source(script_file, local = .GlobalEnv),
                                     error = function(e) cat("Error:", e$message)), type = "output")
      step_logs[[5]](paste(msg, collapse="\n"))
    } else step_logs[[5]]("No script found for Step 5")
    step5_data(df)
  })
  output$console_5 <- renderText({ step_logs[[5]]() })
  output$download_5 <- downloadHandler(filename = "matching_GLASS_specimen_types.xlsx",
                                       content = function(file) writexl::write_xlsx(step5_data(), file))


  #
  # # ---- Step 3 lazy-load ----
  # observeEvent(input$next_2, {
  #   if(isTRUE(input$completed_2)) {
  #     updateTabsetPanel(session, "steps", "Step 3")
  #     step3_data(ddd_updates)  # Load only now
  #   }
  # })
  # output$table_3 <- renderRHandsontable({
  #   df <- step3_data()
  #   req(df)
  #   rhandsontable(df) %>%
  #     hot_col("ATC level name", readOnly = TRUE, width = 300)
  # })
  # observe({ req(input$table_3); step3_data(hot_to_r(input$table_3)) })
  # observeEvent(input$run_script_3, {
  #   df <- step3_data()
  #   script_file <- "amc_scripts/f4.R"
  #   if(file.exists(script_file)) {
  #     msg <- capture.output(tryCatch(source(script_file, local = .GlobalEnv),
  #                                    error = function(e) cat("Error:", e$message)), type = "output")
  #     step_logs[[3]](paste(msg, collapse="\n"))
  #   } else step_logs[[3]]("No script found for Step 3")
  #   step3_data(df)
  # })
  # output$console_3 <- renderText({ step_logs[[3]]() })
  # output$download_3 <- downloadHandler(filename = "DDD_information_updates.xlsx",
  #                                      content = function(file) writexl::write_xlsx(step3_data(), file))
  #
  # # ---- Step 4 lazy-load ----
  # observeEvent(input$next_3, {
  #   if(isTRUE(input$completed_3)) {
  #     updateTabsetPanel(session, "steps", "Step 4")
  #     # Load step4 data only now
  #     step4_data(ddd_updates)
  #   }
  # })
  # output$table_4 <- renderRHandsontable({
  #   df <- step4_data()
  #   req(df)
  #   rhandsontable(df)
  # })
  # observe({ req(input$table_4); step4_data(hot_to_r(input$table_4)) })
  # observeEvent(input$run_script_4, {
  #   df <- step4_data()
  #   script_file <- "amc_scripts/f5.R"
  #   if(file.exists(script_file)) {
  #     msg <- capture.output(tryCatch(source(script_file, local = .GlobalEnv),
  #                                    error = function(e) cat("Error:", e$message)), type = "output")
  #     step_logs[[4]](paste(msg, collapse="\n"))
  #   } else step_logs[[4]]("No script found for Step 4")
  #   step4_data(df)
  # })
  # output$console_4 <- renderText({ step_logs[[4]]() })
  # # output$download_4 <- downloadHandler(filename = "matching_unclear_antibiotic_entries.xlsx",
  # #                                      content = function(file) writexl::write_xlsx(step4_data(), file))
  #
  # # ---- Step 5 lazy-load ----
  # observeEvent(input$next_4, {
  #   if(isTRUE(input$completed_4)) {
  #     updateTabsetPanel(session, "steps", "Step 5")
  #     # Load step5 data only now
  #     step5_data(unclassified_abs)
  #   }
  # })
  # output$table_5 <- renderRHandsontable({
  #   df <- step5_data()
  #   req(df)
  #   rhandsontable(df) %>%
  #     hot_col("Class", type = "dropdown", source = sort(antibiotic_classes_amc), width = 200) %>%
  #     hot_col("Category", type = "dropdown", source = c(' ','Access','Watch', 'Reserve','Uncategorized'), width = 150) %>%
  #     hot_col("antibiotic_names", readOnly = TRUE, width = 300)  # Optional: Make label column readonly
  # })
  # observe({ req(input$table_5); step5_data(hot_to_r(input$table_5)) })
  # observeEvent(input$run_script_5, {
  #   df <- step5_data()
  #   script_file <- "amc_scripts/f6.R"
  #   if(file.exists(script_file)) {
  #     msg <- capture.output(tryCatch(source(script_file, local = .GlobalEnv),
  #                                    error = function(e) cat("Error:", e$message)), type = "output")
  #     step_logs[[5]](paste(msg, collapse="\n"))
  #   } else step_logs[[5]]("No script found for Step 5")
  #   step5_data(df)
  # })
  # output$console_5 <- renderText({ step_logs[[5]]() })
  # output$download_5 <- downloadHandler(filename = "update_AMC_classes.xlsx",
  #                                      content = function(file) writexl::write_xlsx(step5_data(), file))
  #

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
