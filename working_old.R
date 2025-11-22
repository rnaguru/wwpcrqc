# ==========================================================
# app.R — SARS-CoV-2 qPCR Primer/Probe Mutation Heatmap + Mutational Burden + Wastewater Z-score
# ==========================================================

library(shiny)
library(tidyverse)
library(lubridate)
library(httr)
library(jsonlite)
library(DT)
library(future.apply)
library(cowplot)
library(patchwork)
plan(multisession)

# --------------------------
# Primer/probe positions & wildtype bases
# --------------------------
positions_N1 <- c(28287:28306, 28309:28332, 28335:28358)
wt_N1 <- c(
  strsplit("GACCCCAAAATCAGCGAAAT", "")[[1]],
  strsplit("ACCCCGCATTACGTTTGGTGGACC", "")[[1]],
  strsplit("TCTGGTTACTGCCAGTTGAATCTG", "")[[1]]
)
blocks_N1 <- c(rep("Forward", length(28287:28306)),
               rep("Probe", length(28309:28332)),
               rep("Reverse", length(28335:28358)))

positions_N2 <- c(29164:29183, 29188:29210, 29213:29230)
wt_N2 <- c(
  strsplit("TTACAAACATTGGCCGCAAA", "")[[1]],
  strsplit("ACAATTTGCCCCCAGCGCTTCAG", "")[[1]],
  strsplit("GCGCGACATTCCGAAGAA", "")[[1]]
)
blocks_N2 <- c(rep("Forward", length(29164:29183)),
               rep("Probe", length(29188:29210)),
               rep("Reverse", length(29213:29230)))

default_assays <- list(
  N1 = list(positions=positions_N1, wt=wt_N1, blocks=blocks_N1),
  N2 = list(positions=positions_N2, wt=wt_N2, blocks=blocks_N2)
)

# Precompute x_pos for heatmap
precompute_xpos <- function(positions, blocks){
  block_offsets <- c("Forward" = 0, "Probe" = 0.5, "Reverse" = 1)
  xpos <- as.numeric(factor(positions)) + block_offsets[blocks]
  return(xpos)
}
default_assays$N1$x_pos <- precompute_xpos(default_assays$N1$positions, default_assays$N1$blocks)
default_assays$N2$x_pos <- precompute_xpos(default_assays$N2$positions, default_assays$N2$blocks)

# --------------------------
# Fetch mutation counts
# --------------------------
fetch_nonwt_lapis <- function(positions, wt_bases, start_date, end_date) {
  lapis_url <- "https://lapis.cov-spectrum.org/open/v2/sample/aggregated"
  all_dates <- seq.Date(as.Date(start_date), as.Date(end_date), by="week")
  
  df_list <- future_lapply(
    seq_along(positions),
    function(i){
      pos <- positions[i]
      ref <- wt_bases[i]
      variant_query <- paste0(ref, pos)
      
      params <- list(
        variantQuery = variant_query,
        fields = "date",
        dateFrom = start_date,
        dateTo = end_date
      )
      
      resp <- httr::RETRY("GET", lapis_url, query=params, times=2)
      if(resp$status_code != 200) return(tibble(date = all_dates, position = pos, count = 0))
      
      dat <- httr::content(resp, as="parsed", simplifyVector=FALSE)
      if(is.null(dat$data) || length(dat$data) == 0) {
        tibble(date = all_dates, position = pos, count = 0)
      } else {
        df <- tibble(
          date = as.Date(sapply(dat$data, function(x) x$date)),
          count = sapply(dat$data, function(x) x$count)
        )
        tibble(date = all_dates) %>%
          left_join(df, by="date") %>%
          mutate(count = ifelse(is.na(count), 0, count), position = pos)
      }
    },
    future.seed = TRUE
  )
  
  bind_rows(df_list)
}

# Fetch total sequences
fetch_total_sequences <- function(start_date, end_date) {
  lapis_url <- "https://lapis.cov-spectrum.org/open/v2/sample/aggregated"
  all_dates <- seq.Date(as.Date(start_date), as.Date(end_date), by="week")
  
  params <- list(
    fields = "date",
    dateFrom = start_date,
    dateTo = end_date
  )
  
  resp <- httr::RETRY("GET", lapis_url, query=params, times=2)
  if(resp$status_code != 200) return(tibble(date=all_dates, total=0))
  
  dat <- httr::content(resp, as="parsed", simplifyVector=FALSE)
  
  if(is.null(dat$data) || length(dat$data)==0){
    tibble(date=all_dates, total=0)
  } else {
    df <- tibble(
      date = as.Date(sapply(dat$data, function(x) x$date)),
      total = sapply(dat$data, function(x) x$count)
    )
    tibble(date=all_dates) %>%
      left_join(df, by="date") %>%
      mutate(total = ifelse(is.na(total), 0, total))
  }
}

# --------------------------
# Shiny UI
# --------------------------
ui <- fluidPage(
  #add_busy_spinner(spin = "dots", position = "full-page"),
  titlePanel("SARS-CoV-2 qPCR Primer/Probe Mutation Heatmap + Mutational Burden"),
  sidebarLayout(
    sidebarPanel(
      dateRangeInput("daterange","Date range",
                     start=as.Date("2020-01-01"),
                     end=Sys.Date()),
      actionButton("fetch","Fetch Data"),
      verbatimTextOutput("debug"),
      width=3
    ),
    mainPanel(
      plotOutput("heatmap_N1", height="300px"),
      plotOutput("heatmap_N2", height="300px"),
      plotOutput("zscore_plot", height="250px"),
      plotOutput("mutational_burden_plot", height="350px"),
      plotOutput("total_sequences_plot", height="300px"),
      
      tags$div(
        style="margin-top:10px; font-size:14px;",
        HTML("<b>Panel A:</b> N1 heatmap<br>
             <b>Panel B:</b> N2 heatmap<br>
             <b>Z:</b> Wastewater N1/N2 Z-score<br>
             <b>Panel C:</b> Weekly Weighted Mutational Burden + Interpretation<br>
             <b>Panel D:</b> Weekly count of submitted sequences (log10 scale)")
      ),
      
      tabsetPanel(
        tabPanel("Table", DTOutput("table")),
        tabPanel("Download", downloadButton("download","Download CSV"))
      )
    )
  )
)

# --------------------------
# Shiny Server
# --------------------------
server <- function(input, output, session){
  
  data_reactive <- eventReactive(input$fetch,{
    range <- input$daterange
    
    df_N1 <- fetch_nonwt_lapis(
      default_assays$N1$positions,
      default_assays$N1$wt,
      as.character(range[1]),
      as.character(range[2])
    ) %>%
      mutate(block = default_assays$N1$blocks[match(position, default_assays$N1$positions)],
             x_pos = default_assays$N1$x_pos[match(position, default_assays$N1$positions)])
    
    df_N2 <- fetch_nonwt_lapis(
      default_assays$N2$positions,
      default_assays$N2$wt,
      as.character(range[1]),
      as.character(range[2])
    ) %>%
      mutate(block = default_assays$N2$blocks[match(position, default_assays$N2$positions)],
             x_pos = default_assays$N2$x_pos[match(position, default_assays$N2$positions)])
    
    total_seq <- fetch_total_sequences(range[1], range[2])
    
    list(N1 = df_N1, N2 = df_N2, total_seq = total_seq)
  })
  
  output$debug <- renderPrint({
    df <- data_reactive()
    paste("N1 rows:", nrow(df$N1), "N2 rows:", nrow(df$N2))
  })
  
  # Heatmap plotting function
  plot_heatmap_blocks <- function(df, assay_name, total_seq_df){
    req(nrow(df) > 0)
    
    df_summary <- df %>%
      group_by(position, date, block) %>%
      summarise(count = sum(count), .groups="drop") %>%
      left_join(total_seq_df %>% rename(total_sequences=total), by="date") %>%
      mutate(count = count / pmax(total_sequences,1))
    
    x_all <- sort(unique(df_summary$position))
    x_breaks <- x_all[seq(1, length(x_all), by=2)]
    
    all_weeks <- sort(unique(df_summary$date))
    y_breaks <- all_weeks[seq(1, length(all_weeks), by=13)]
    
    ggplot(df_summary, aes(x = position, y = date, fill = count)) +
      geom_tile(color = "grey92") +
      scale_fill_gradient2(
        low = "white",
        mid = "#dda0dd",
        high = "maroon",
        midpoint = 0.5,
        limits = c(0,1),
        na.value = "black"
      ) +
      facet_grid(. ~ block, scales = "free_x", space = "free_x") +
      scale_x_continuous(
        breaks = x_breaks,
        minor_breaks = x_all
      ) +
      scale_y_date(
        breaks = y_breaks,
        labels = function(d) format(d, "%b %Y"),
        expand = c(0,0)
      ) +
      labs(x="nucleotide position in SARS-CoV-2 reference genome",
           y="Submission date",
           fill="Fraction of Sequences",
           title = assay_name) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=0.5),
        axis.text.y = element_text(size=10),
        panel.spacing = unit(0.5, "lines"),
        axis.ticks = element_blank()
      )
  }
  
  # Render heatmaps
  output$heatmap_N1 <- renderPlot({
    p <- plot_heatmap_blocks(
      data_reactive()$N1,
      "N1",
      total_seq_df = data_reactive()$total_seq
    )
    cowplot::ggdraw(p) + cowplot::draw_plot_label("A", x = 0, y = 1, hjust=0, vjust=1, size=14)
  })
  
  output$heatmap_N2 <- renderPlot({
    p <- plot_heatmap_blocks(
      data_reactive()$N2,
      "N2",
      total_seq_df = data_reactive()$total_seq
    )
    cowplot::ggdraw(p) + cowplot::draw_plot_label("B", x = 0, y = 1, hjust=0, vjust=1, size=14)
  })
  
  # --------------------------
  # Wastewater Z-score plot
  # --------------------------
  output$zscore_plot <- renderPlot({
    
    wwopen <- read.csv("https://raw.githubusercontent.com/Big-Life-Lab/PHESD/main/Wastewater/Ottawa/Data/wastewater_virus.csv") %>%
      mutate(Date = as.Date(sampleDate, format= "%Y-%m-%d")) %>%
      mutate(
        N1N2norm = (covN1_nPMMoV_meanNr + covN2_nPMMoV_meanNr)/2,
        differenceratio = (covN1_nPMMoV_meanNr - covN2_nPMMoV_meanNr) /
          (covN1_nPMMoV_meanNr + covN2_nPMMoV_meanNr),
        zscore = (differenceratio - mean(differenceratio, na.rm=TRUE)) / sd(differenceratio, na.rm=TRUE)
      ) %>%
      filter(covN1_nPMMoV_meanNr != 0, covN2_nPMMoV_meanNr != 0)
    
    ggplot(wwopen, aes(x=Date, y=zscore)) +
      geom_line(color="darkred", size=1) +
      geom_hline(yintercept=0, linetype="dashed") +
      labs(
        title="Wastewater N1/N2 Z-score (Ottawa WW data)",
        x=NULL,
        y="Z-score"
      ) +
      theme_minimal(base_size=14)
  })
  
  # --------------------------
  # Mutational burden + interpretation
  # --------------------------
  output$mutational_burden_plot <- renderPlot({
    
    df_all <- bind_rows(
      data_reactive()$N1 %>% mutate(locus="N1"),
      data_reactive()$N2 %>% mutate(locus="N2")
    ) %>%
      left_join(data_reactive()$total_seq %>% rename(total_sequences=total), by="date") %>%
      mutate(burden = count / pmax(total_sequences,1)) %>%
      group_by(date, locus) %>%
      summarise(mean_burden = sum(burden), .groups="drop") %>%
      arrange(date)
    
    # Identify low-sequence weeks
    low_seq_df <- df_all %>%
      group_by(date) %>%
      summarise(total_sequences = max(data_reactive()$total_seq$total[which(data_reactive()$total_seq$date==date)]),
                .groups="drop") %>%
      arrange(date) %>%
      mutate(low_flag = total_sequences <= 50,
             group = cumsum(c(TRUE, diff(low_flag) != 0)))
    
    low_seq_ranges <- low_seq_df %>%
      filter(low_flag) %>%
      group_by(group) %>%
      summarise(xmin = min(date)-3, xmax = max(date)+3, .groups="drop")
    
    waterfall_df <- df_all %>%
      tidyr::pivot_wider(names_from = locus, values_from = mean_burden, values_fill=0) %>%
      tidyr::pivot_longer(cols=c("N1","N2"), names_to="locus", values_to="mean_burden")
    
    # Stacked area plot
    p_burden <- ggplot(waterfall_df, aes(x=date, y=mean_burden, fill=locus)) +
      geom_rect(data=low_seq_ranges, aes(xmin=xmin, xmax=xmax, ymin=0, ymax=Inf),
                inherit.aes=FALSE, fill="darkorange", alpha=0.5) +
      geom_area(position="stack", color="black", linewidth=0.2, alpha=0.8) +
      scale_fill_manual(values=c("N1"="firebrick","N2"="steelblue")) +
      labs(title="Weekly Weighted Mutational Burden (stacked N1+N2)", y="Weighted Mutational Burden", x=NULL, fill="Locus") +
      theme_minimal(base_size=14) + theme(legend.position="top")
    
    cowplot::ggdraw(p_burden) + cowplot::draw_plot_label("C", x=0, y=1, hjust=0, vjust=1, size=14)
  })
  
  # --------------------------
  # Total sequences (Panel D)
  # --------------------------
  output$total_sequences_plot <- renderPlot({
    df <- data_reactive()$total_seq
    req(df)
    
    low_seq_df <- df %>% arrange(date) %>%
      mutate(low_flag = total <= 50, group = cumsum(c(TRUE, diff(low_flag) != 0)))
    
    low_seq_ranges <- low_seq_df %>%
      filter(low_flag) %>%
      group_by(group) %>%
      summarise(xmin=min(date)-3, xmax=max(date)+3, .groups="drop")
    
    y_min <- min(df$total[df$total>0])*0.9
    y_max <- max(df$total)*1.1
    
    p <- ggplot(df, aes(x=date, y=total)) +
      geom_rect(data=low_seq_ranges, aes(xmin=xmin, xmax=xmax, ymin=y_min, ymax=y_max),
                inherit.aes=FALSE, fill="darkorange", alpha=0.5) +
      geom_line(color="darkgreen") +
      geom_point(color="darkgreen") +
      scale_y_log10() +
      labs(title="Weekly count of submitted sequences", x="Week", y="Total Sequences (log10)") +
      theme_minimal()
    
    cowplot::ggdraw(p) + cowplot::draw_plot_label("D", x=0, y=1, hjust=0, vjust=1, size=14)
  })
  
  # --------------------------
  # Table and download
  # --------------------------
  output$table <- renderDT({
    df <- data_reactive()
    df_all <- bind_rows(
      df$N1 %>% mutate(assay="N1"),
      df$N2 %>% mutate(assay="N2")
    ) %>%
      left_join(df$total_seq %>% rename(total_sequences=total), by="date")
    datatable(df_all, options=list(pageLength=20))
  })
  
  output$download <- downloadHandler(
    filename=function(){paste0("mut_counts_", Sys.Date(), ".csv")},
    content=function(file){
      df <- data_reactive()
      df_all <- bind_rows(
        df$N1 %>% mutate(assay="N1"),
        df$N2 %>% mutate(assay="N2")
      ) %>%
        left_join(df$total_seq %>% rename(total_sequences=total), by="date")
      write_csv(df_all, file)
    }
  )
}

# --------------------------
# Run the app
# --------------------------
shinyApp(ui, server)
