#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

# app.R
library(shiny)
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(DT)
library(scales)

theme_set(ggpubr::theme_pubr(base_size = 13))

ui <- fluidPage(
  titlePanel("Bradford: Concentration + Loading (long-format)"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload long-format Excel", accept = c(".xlsx", ".xls")),
      textInput("sheet", "Sheet name", value = "End point"),
      numericInput("linear_cutoff", "Std curve max conc to fit (mg/mL)", value = NA, min = 0, step = 0.1),
      numericInput("cv_flag", "Flag samples with CV% >", value = 30, min = 0, step = 1),
      textInput("loads", "Desired loads (µg), comma-separated", value = "20,25,30"),
      hr(),
      downloadButton("dl_conc", "Download concentrations (CSV)"),
      downloadButton("dl_vol",  "Download loading volumes (CSV)")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Parse & Name",
                 verbatimTextOutput("parse_msg"),
                 DTOutput("map_table")),
        tabPanel("Standard curve",
                 plotOutput("std_plot", height = 440),
                 verbatimTextOutput("std_eq")),
        tabPanel("Concentrations",
                 plotOutput("spread_plot", height = 440),
                 DTOutput("conc_table")),
        tabPanel("Curve + Samples",
                 plotOutput("overlay_plot", height = 440)),
        tabPanel("Loading volumes",
                 DTOutput("vol_table"))
      )
    )
  )
)

server <- function(input, output, session){
  
  loads_vec <- reactive({
    as.numeric(strsplit(gsub("\\s", "", input$loads), ",")[[1]])
  })
  
  # --- READ & PARSE ---
  plate_raw <- reactive({
    req(input$file)
    raw0 <- readxl::read_excel(input$file$datapath, sheet = input$sheet, col_names = FALSE)
    hdr <- which(raw0[[1]] == "Content")[1]
    validate(need(!is.na(hdr), "Couldn't find the 'Content' header. Check sheet name / file."))
    
    dat <- readxl::read_excel(input$file$datapath, sheet = input$sheet, skip = hdr - 1)
    
    nm <- names(dat)
    nm <- str_replace_all(nm, "\\r\\n", "\n")
    nm <- str_replace_all(nm, fixed("Well\nCol"), "col")
    nm <- str_replace_all(nm, fixed("Well\nRow"), "row")
    nm <- str_replace_all(nm, fixed("Raw Data (A-595)"), "Absorbance")
    nm <- str_replace_all(nm, "Standard Concentrations.*", "std.conc")
    names(dat) <- nm
    
    validate(need(all(c("Content","col","row","Absorbance") %in% names(dat)),
                  "Required columns not found: Content, Well Col/Row, Raw Data (A-595)."))
    
    plate <- dat %>%
      select(Content, col, row, std.conc = any_of("std.conc"), Absorbance) %>%
      unite("well_ID", c(row, col), sep = "", remove = TRUE) %>%
      separate(Content, into = c("Type","id"), sep = " ", remove = FALSE, fill = "right")
    
    blank_med <- plate %>% filter(Type == "Blank") %>% pull(Absorbance) %>% median(na.rm = TRUE)
    list(plate = plate, blank_med = blank_med)
  })
  
  output$parse_msg <- renderPrint({
    req(plate_raw())
    cat("Rows:", nrow(plate_raw()$plate), "\n")
    cat("Median blank A595:", round(plate_raw()$blank_med, 3), "\n")
    cat("Samples detected:", length(unique(plate_raw()$plate$Content[plate_raw()$plate$Type=="Sample"])), "\n")
  })
  
  # --- EDITABLE sample naming ---
  map_tbl <- reactiveVal(NULL)
  observeEvent(plate_raw(), {
    plate <- plate_raw()$plate
    smp <- plate %>% filter(Type == "Sample") %>% distinct(Content, id) %>% arrange(Content)
    smp$desc <- ""
    smp$name <- smp$Content
    map_tbl(smp)
  }, ignoreInit = FALSE)
  
  output$map_table <- renderDT({
    req(map_tbl())
    datatable(map_tbl(), editable = list(target = "cell", disable = list(columns = c(0,1))),
              options = list(dom = "tip", pageLength = 10))
  })
  observeEvent(input$map_table_cell_edit, {
    info <- input$map_table_cell_edit
    df <- map_tbl(); df[info$row, info$col + 1] <- info$value; map_tbl(df)
  })
  
  # --- CORRECTED DATA & FIT ---
  corrected <- reactive({
    req(plate_raw())
    pr <- plate_raw()
    plate <- pr$plate
    blank_med <- pr$blank_med
    
    proc <- plate %>% filter(Type != "Blank") %>%
      mutate(cor.Absorbance = Absorbance - blank_med)
    
    std <- proc %>%
      filter(Type == "Standard", !is.na(std.conc)) %>%
      group_by(std.conc) %>%
      summarise(A595_corr_mean = mean(cor.Absorbance, na.rm = TRUE),
                A595_corr_sd   = sd(cor.Absorbance,   na.rm = TRUE),
                n = n(), .groups = "drop")
    
    cutoff <- input$linear_cutoff
    if (is.na(cutoff)) cutoff <- max(std$std.conc, na.rm = TRUE)
    
    fit <- lm(A595_corr_mean ~ std.conc, data = dplyr::filter(std, std.conc <= cutoff))
    list(proc = proc, std = std, fit = fit)
  })
  
  output$std_eq <- renderPrint({
    req(corrected())
    fit <- corrected()$fit
    cat("A = m*C + b\n")
    cat("m =", round(coef(fit)[2], 4), "A per (mg/mL)\n")
    cat("b =", round(coef(fit)[1], 4), "A\n")
    cat("R^2 =", round(summary(fit)$r.squared, 4), "\n")
  })
  
  output$std_plot <- renderPlot({
    req(corrected())
    std <- corrected()$std
    fit <- corrected()$fit
    m <- unname(coef(fit)[2]); b <- unname(coef(fit)[1]); r2 <- summary(fit)$r.squared
    
    p <- ggplot(std, aes(std.conc, A595_corr_mean)) +
      geom_point(size = 3, color = "#2b8cbe") +
      geom_errorbar(aes(ymin = A595_corr_mean - A595_corr_sd,
                        ymax = A595_corr_mean + A595_corr_sd),
                    width = 0.03, color = "#2b8cbe", alpha = 0.6) +
      geom_smooth(method = "lm", se = FALSE, color = "#de2d26", linewidth = 1) +
      labs(title = "Bradford standard curve (averaged)",
           subtitle = paste0("A = ", round(m,3),"·C + ", round(b,3), "   |   R² = ", round(r2,3)),
           x = "BSA concentration (mg/mL)", y = "A595 (blank-corrected)")
    ggpubr::ggpar(p, legend = "none")
  })
  
  # --- CONCENTRATIONS & CVs ---
  conc_tbl <- reactive({
    req(corrected(), map_tbl())
    proc <- corrected()$proc
    fit  <- corrected()$fit
    m <- unname(coef(fit)[2]); b <- unname(coef(fit)[1])
    
    out <- proc %>%
      mutate(Conc_ug_per_uL = ifelse(Type == "Sample", (cor.Absorbance - b) / m, NA_real_)) %>%
      left_join(map_tbl(), by = c("Content","id"))
    
    sum_tbl <- out %>%
      filter(Type == "Sample") %>%
      group_by(Content, name, desc) %>%
      summarise(
        n = n(),
        mean_conc   = mean(Conc_ug_per_uL, na.rm = TRUE),
        sd_conc     = sd(Conc_ug_per_uL, na.rm = TRUE),
        median_conc = median(Conc_ug_per_uL, na.rm = TRUE),
        cv_percent  = ifelse(mean_conc > 0, (sd_conc/mean_conc)*100, NA_real_),
        .groups = "drop"
      ) %>%
      mutate(flag_high_CV = cv_percent > input$cv_flag)
    
    list(per_well = out, summary = sum_tbl, m = m, b = b)
  })
  
  output$spread_plot <- renderPlot({
    req(conc_tbl())
    out <- conc_tbl()$per_well %>% filter(Type == "Sample")
    p <- ggplot(out, aes(x = Conc_ug_per_uL, y = ifelse(is.na(name), Content, name))) +
      geom_boxplot(outlier.shape = NA, fill = "#f0f0f0") +
      geom_jitter(height = 0.15, width = 0, size = 1.8, alpha = 0.85, color = "#4d4d4d") +
      labs(title = "Sample concentration spread (replicates)",
           x = "Concentration (µg/µL)", y = NULL)
    ggpubr::ggpar(p, legend = "none")
  })
  
  output$conc_table <- renderDT({
    req(conc_tbl())
    datatable(conc_tbl()$summary %>%
                mutate(across(c(mean_conc, sd_conc, median_conc, cv_percent), ~round(.x, 3))),
              options = list(pageLength = 10, dom = "tip"))
  })
  
  # --- OVERLAY ON CURVE ---
  output$overlay_plot <- renderPlot({
    req(corrected(), conc_tbl())
    std <- corrected()$std
    fit <- corrected()$fit
    m <- unname(coef(fit)[2]); b <- unname(coef(fit)[1])
    
    samples_plot <- conc_tbl()$summary %>%
      mutate(A595_from_curve = m * median_conc + b,
             label = ifelse(is.na(name) | name=="", Content, name))
    
    p <- ggplot(std, aes(x = std.conc, y = A595_corr_mean)) +
      geom_point(size = 3, color = "#2b8cbe") +
      geom_smooth(method = "lm", se = FALSE, color = "#de2d26", linewidth = 1) +
      geom_point(data = samples_plot,
                 aes(x = median_conc, y = A595_from_curve),
                 shape = 4, size = 3, color = "black", stroke = 1.2) +
      ggrepel::geom_text_repel(
        data = samples_plot,
        aes(x = median_conc, y = A595_from_curve, label = label),
        size = 3, color = "black"
      ) +
      labs(title = "Standard curve with sample medians overlaid",
           x = "Protein concentration (mg/mL = µg/µL)",
           y = "A595 (blank-corrected)")
    ggpubr::ggpar(p, legend = "none")
  })
  
  # --- VOLUMES ---
  vol_tbl <- reactive({
    req(conc_tbl())
    conc_summary <- conc_tbl()$summary %>% mutate(label = ifelse(is.na(name) | name=="", Content, name))
    loads <- loads_vec(); validate(need(all(is.finite(loads)), "Provide numeric loads (comma-separated)."))
    conc_summary %>%
      mutate(Conc_used = median_conc) %>%
      crossing(tibble(Load_ug = loads)) %>%
      mutate(Volume_uL = ifelse(is.finite(Conc_used) & Conc_used > 0, Load_ug / Conc_used, NA_real_)) %>%
      select(label, desc, Conc_used, Load_ug, Volume_uL, cv_percent, flag_high_CV) %>%
      arrange(label, Load_ug)
  })
  
  output$vol_table <- renderDT({
    req(vol_tbl())
    datatable(vol_tbl() %>% mutate(across(c(Conc_used, Volume_uL, cv_percent), ~round(.x, 3))),
              options = list(pageLength = 10, dom = "tip"))
  })
  
  # --- DOWNLOADS ---
  output$dl_conc <- downloadHandler(
    filename = function() "bradford_sample_concentrations.csv",
    content  = function(file) write.csv(conc_tbl()$summary, file, row.names = FALSE)
  )
  output$dl_vol <- downloadHandler(
    filename = function() "bradford_loading_volumes.csv",
    content  = function(file) write.csv(vol_tbl(), file, row.names = FALSE)
  )
}

shinyApp(ui, server)