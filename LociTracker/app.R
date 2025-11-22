# app.R — SARS-CoV-2 qPCR Primer/Probe Mutation Heatmap + Mutational Burden + Wastewater Z-score
# Prefetched CSV implementation for N1, N2, N200, E_Sarbeco

library(shiny)
library(tidyverse)
library(lubridate)
library(DT)
library(cowplot)
library(patchwork)
library(jsonlite)
library(ggplot2)
library(shinybusy)
library(shinyWidgets)
library(tidyquant)
library(stringr)
library(egg)
library(ggrepel)
library(plotly)
library(ggalt)
library(ggstream)
library(zoo)    # for rollmean

# --------------------------------------------------------------
# Remote prefetched data location (GitHub)
github_base <- "https://raw.githubusercontent.com/rnaguru/wwpcrqc/main/LociTracker_site/data"
manifest_url <- paste0(github_base, "/manifest.json")
manifest <- jsonlite::fromJSON(manifest_url)

load_remote_csv <- function(filename){
  url <- paste0(github_base, "/", filename)
  tryCatch(
    readr::read_csv(url, show_col_types = FALSE),
    error = function(e){
      warning(paste("Could not load:", url))
      tibble(date=as.Date(character()), position=numeric(), count=numeric())
    }
  )
}

prefetched_data <- list(
  N1        = load_remote_csv(manifest$prefetched_N1),
  N2        = load_remote_csv(manifest$prefetched_N2),
  N200      = load_remote_csv(manifest$prefetched_N200),
  E_Sarbeco = load_remote_csv(manifest$prefetched_E_Sarbeco)
)

total_sequences <- load_remote_csv(manifest$total_sequences)

# --------------------------------------------------------------
# Primer/probe positions & wildtype bases + blocks
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

positions_E <- c(26269:26294, 26332:26357, 26360:26381)
blocks_E <- c(rep("Forward", 26), rep("Probe", 26), rep("Reverse", 22)) # dummy

positions_N200 <- c(28840:28861, 28891:28905, 28940:28960)
blocks_N200 <- c(rep("Forward", 22), rep("Probe", 15), rep("Reverse", 21)) # dummy

default_assays <- list(
  N1 = list(positions=positions_N1, wt=wt_N1, blocks=blocks_N1),
  N2 = list(positions=positions_N2, wt=wt_N2, blocks=blocks_N2),
  E_Sarbeco = list(positions=positions_E, blocks=blocks_E),
  N200 = list(positions=positions_N200, blocks=blocks_N200)
)

precompute_xpos <- function(positions, blocks){
  block_offsets <- c("Forward" = 0, "Probe" = 0.5, "Reverse" = 1)
  xpos <- as.numeric(factor(positions)) + block_offsets[blocks]
  return(xpos)
}

default_assays$N1$x_pos <- precompute_xpos(default_assays$N1$positions, default_assays$N1$blocks)
default_assays$N2$x_pos <- precompute_xpos(default_assays$N2$positions, default_assays$N2$blocks)
default_assays$E_Sarbeco$x_pos <- precompute_xpos(default_assays$E_Sarbeco$positions, default_assays$E_Sarbeco$blocks)
default_assays$N200$x_pos <- precompute_xpos(default_assays$N200$positions, default_assays$N200$blocks)

# --------------------------------------------------------------
# Determine global date range
all_dates <- bind_rows(prefetched_data) %>% pull(date)
all_dates <- all_dates[!is.na(all_dates)]
if(length(all_dates)==0){
  global_min_date <- as.Date("2020-01-01")
  global_max_date <- Sys.Date()
} else {
  global_min_date <- min(all_dates)
  global_max_date <- max(all_dates)
}

#---------------------------------------------------------------
# Load and calculate wastewater data
wwopen <- read.csv("https://raw.githubusercontent.com/Big-Life-Lab/PHESD/main/Wastewater/Ottawa/Data/wastewater_virus.csv") %>%
  mutate(Date = as.Date(sampleDate, format= "%Y-%m-%d")) %>%
  mutate(N1N2norm = ((covN1_nPMMoV_meanNr + covN2_nPMMoV_meanNr)/2)) %>%
  mutate(zscorediff = ((covN1_nPMMoV_meanNr - covN2_nPMMoV_meanNr)-mean(covN1_nPMMoV_meanNr - covN2_nPMMoV_meanNr, na.rm=TRUE))/sd(covN1_nPMMoV_meanNr - covN2_nPMMoV_meanNr, na.rm=TRUE)) %>%
  mutate(N1N2 = (covN1_nPMMoV_meanNr - covN2_nPMMoV_meanNr)/(0.5*(covN1_nPMMoV_meanNr + covN2_nPMMoV_meanNr))) %>%
  mutate(N1N2v2 = (covN1_nPMMoV_meanNr - covN2_nPMMoV_meanNr)/(covN1_nPMMoV_meanNr + covN2_nPMMoV_meanNr)) %>%
  mutate(zscorediffdenom = (N1N2-mean(N1N2, na.rm=TRUE))/sd(N1N2, na.rm=TRUE)) %>%
  mutate(zscorediffdenomv2 = (N1N2v2-mean(N1N2v2, na.rm=TRUE))/sd(N1N2v2, na.rm=TRUE)) %>%
  filter(!is.na(covN1_nPMMoV_meanNr), !is.na(covN2_nPMMoV_meanNr), covN1_nPMMoV_meanNr != 0, covN2_nPMMoV_meanNr != 0) %>%
  mutate(differenceratio = (covN1_nPMMoV_meanNr - covN2_nPMMoV_meanNr) / (covN1_nPMMoV_meanNr + covN2_nPMMoV_meanNr)) %>%
  mutate(zscore = (differenceratio - mean(differenceratio, na.rm=TRUE))/sd(differenceratio, na.rm=TRUE)) %>%
  mutate(percentchange = ((covN1_nPMMoV_meanNr - covN2_nPMMoV_meanNr) / covN1_nPMMoV_meanNr)*100)

# create 7d rolling mean column for derivative plot
wwopen <- wwopen %>%
  arrange(Date) %>%
  mutate(roll7d = zoo::rollmean(N1N2norm, 7, fill = NA, align = "right")) %>%
  mutate(dxN1N2 = c(NA, diff(roll7d)))


# --------------------------------------------------------------
# Shiny UI
ui <- fluidPage(
  # add a spinner to show busy while loading
  add_busy_spinner(spin = "dots", position = "full-page"),
  
  titlePanel("SARS-CoV-2 qPCR Primer/Probe mutations in global clinical genomes"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("date_slider", "Date range",
                  min = global_min_date,
                  max = global_max_date,
                  value = c(global_min_date, global_max_date),
                  timeFormat = "%Y-%m-%d"),
      selectInput("panelA_assay", "Panel A assay (heatmap)",
                  choices = names(default_assays), selected = "N1"),
      selectInput("panelB_assay", "Panel B assay (heatmap)",
                  choices = names(default_assays), selected = "N2"),
      selectizeInput("burden_assays",
                     "Assays included in Panel C (mutational burden)",
                     choices = names(default_assays),
                     selected = c("N1","N2"),
                     multiple = TRUE),
      sliderInput("N1N2spanslider", "WW smoothing span (LOESS)", min = 0.1, max = 1, value = 0.5),
      width = 3
    ),
    
    mainPanel(
      tabsetPanel(
        type = "tabs",
        
        # ===============================
        # Clinical genomic surveillance
        # ===============================
        tabPanel(
          "Clinical genomic surveillance",
          
          # Nested subtabs for plots / table / download
          tabsetPanel(
            id = "clinical_subtabs",
            type = "tabs",
            
            tabPanel(
              "Heatmaps",
              plotOutput("heatmap_panelA", height = "300px"),
              plotOutput("heatmap_panelB", height = "300px"),
              plotOutput("mutational_burden_plot", height = "300px"),
              plotOutput("total_sequences_plot", height = "300px")
            ),
            
            
            tabPanel(
              "Table / Download",
              h4("Mutation table"),
              DTOutput("table"),
              br(),
              downloadButton("download", "Download CSV")
            )
          )
        ),
        
        # ===============================
        # Wastewater surveillance
        # ===============================
        tabPanel(
          "Wastewater surveillance",
          plotOutput("ww_combined", height = "900px") 
        )
      )
    )
  )
)


# --------------------------------------------------------------
# Shiny Server
server <- function(input, output, session){

  # --------------------------
  # Reactive filtered data
  data_reactive <- reactive({
    range <- input$date_slider

    df_list <- lapply(prefetched_data, function(df){
      if(nrow(df)==0) return(df)
      df %>% filter(date >= range[1], date <= range[2])
    })

    # Add blocks and x_pos
    df_list <- lapply(names(df_list), function(nm){
      df <- df_list[[nm]]
      if(nrow(df)==0) return(df)
      df$block <- default_assays[[nm]]$blocks[match(df$position, default_assays[[nm]]$positions)]
      df$x_pos <- default_assays[[nm]]$x_pos[match(df$position, default_assays[[nm]]$positions)]
      df
    })
    names(df_list) <- names(prefetched_data)

    list(
      assays = df_list,
      total_seq = total_sequences %>% filter(date >= range[1], date <= range[2])
    )
  })

  # --------------------------
  # Heatmap plotting function
  plot_heatmap_blocks <- function(df, assay_name, total_seq_df){
    req(nrow(df)>0)
    df_summary <- df %>%
      group_by(position, date, block) %>%
      summarise(count=sum(count, na.rm=TRUE), .groups="drop") %>%
      left_join(total_seq_df %>% rename(total_sequences=total), by="date") %>%
      mutate(fraction = count / pmax(total_sequences,1))

    x_all <- sort(unique(df_summary$position))
    x_breaks <- x_all[seq(1,length(x_all), by=2)]
    y_breaks <- sort(unique(df_summary$date))[seq(1,length(unique(df_summary$date)), by=13)]

    ggplot(df_summary, aes(x=position, y=date, fill=fraction)) +
      geom_tile(color="grey92") +
      scale_fill_gradient2(low="white", mid="#dda0dd", high="maroon",
                           midpoint=0.5, limits=c(0,1), na.value="black") +
      facet_grid(.~block, scales="free_x", space="free_x") +
      scale_x_continuous(breaks=x_breaks, minor_breaks=x_all) +
      scale_y_date(breaks=y_breaks, labels=function(d) format(d,"%b %Y"), expand=c(0,0)) +
      labs(x="Nucleotide position", y="Submission date", fill="Fraction of Sequences", title=assay_name) +
      theme_minimal() +
      theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5),
            axis.text.y=element_text(size=10),
            panel.spacing=unit(0.5,"lines"),
            axis.ticks=element_blank())
  }

  # --------------------------
  # Heatmaps
  output$heatmap_panelA <- renderPlot({
    df <- data_reactive()$assays[[input$panelA_assay]]
    req(df)
    cowplot::ggdraw(plot_heatmap_blocks(df, input$panelA_assay, data_reactive()$total_seq)) +
      cowplot::draw_plot_label("A", x=0, y=1, hjust=0, vjust=1, size=14)
  })

  output$heatmap_panelB <- renderPlot({
    df <- data_reactive()$assays[[input$panelB_assay]]
    req(df)
    cowplot::ggdraw(plot_heatmap_blocks(df, input$panelB_assay, data_reactive()$total_seq)) +
      cowplot::draw_plot_label("B", x=0, y=1, hjust=0, vjust=1, size=14)
  })

  # --------------------------
  # Mutational burden (Panel C)
  output$mutational_burden_plot <- renderPlot({

    assays <- data_reactive()$assays
    total_seq <- data_reactive()$total_seq

    selected <- input$burden_assays
    req(length(selected) > 0)

    df_all <- bind_rows(lapply(selected, function(a) {
      if(!is.null(assays[[a]]) && nrow(assays[[a]]) > 0) {
        assays[[a]] %>% mutate(locus = a)
      } else {
        tibble(date=as.Date(character()), position=numeric(), count=numeric(), block=character(), x_pos=numeric(), locus=a)
      }
    }))

    req(nrow(df_all) > 0)

    df_all <- df_all %>%
      left_join(total_seq %>% rename(total_sequences = total), by = "date") %>%
      mutate(burden = count / pmax(total_sequences, 1)) %>%
      group_by(date, locus) %>%
      summarise(mean_burden = sum(burden, na.rm=TRUE), .groups = "drop") %>%
      arrange(date)

    df_all <- df_all %>%
      tidyr::pivot_wider(names_from = locus, values_from = mean_burden, values_fill = 0) %>%
      tidyr::pivot_longer(cols = all_of(selected), names_to = "locus", values_to = "mean_burden")

    low_seq_df <- total_seq %>%
      arrange(date) %>%
      mutate(low_flag = total <= 50,
             group = cumsum(c(TRUE, diff(low_flag) != 0)))

    low_seq_ranges <- low_seq_df %>%
      filter(low_flag) %>%
      group_by(group) %>%
      summarise(xmin = min(date)-3, xmax = max(date)+3, .groups="drop")

    all_colors <- c(
      "N1" = "firebrick",
      "N2" = "steelblue",
      "N200" = "forestgreen",
      "E_Sarbeco" = "goldenrod"
    )
    assay_colors <- all_colors[selected]

    p_burden <- ggplot(df_all, aes(x = date, y = mean_burden, fill = locus)) +
      geom_rect(data = low_seq_ranges, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf),
                inherit.aes = FALSE, fill = "darkorange", alpha = 0.5) +
      geom_area(position = "stack", color = "black", linewidth = 0.2, alpha = 0.8) +
      scale_fill_manual(values = assay_colors) +
      labs(title = "Weekly Average Mutation Burden (stacked by assay)",
           y = "# of mutations across locus", x = NULL, fill = "Assay") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "top")

    cowplot::ggdraw(p_burden) + cowplot::draw_plot_label("C", x = 0, y = 1, hjust = 0, vjust = 1, size = 14)
  })

  # --------------------------
  # Total sequences (Panel D)
  output$total_sequences_plot <- renderPlot({
    df <- data_reactive()$total_seq
    req(df)
    if(all(is.na(df$total)) || nrow(df)==0){
      return(NULL)
    }
    y_min <- min(df$total[df$total>0], na.rm=TRUE)*0.9
    y_max <- max(df$total, na.rm=TRUE)*1.1

    p <- ggplot(df, aes(x=date, y=total)) +
      geom_line(color="darkgreen") +
      geom_point(color="darkgreen") +
      scale_y_log10() +
      labs(title="Weekly count of submitted sequences", x="Week", y="Total Sequences (log10)") +
      theme_minimal()

    cowplot::ggdraw(p) + cowplot::draw_plot_label("D", x=0, y=1, hjust=0, vjust=1, size=14)
  })

  # --------------------------
  # Table + download
  output$table <- renderDT({
    df_all <- bind_rows(
      lapply(names(data_reactive()$assays), function(nm){
        if(nrow(data_reactive()$assays[[nm]])==0) return(tibble())
        data_reactive()$assays[[nm]] %>% mutate(assay=nm)
      })
    ) %>% left_join(data_reactive()$total_seq %>% rename(total_sequences=total), by="date")
    datatable(df_all, options=list(pageLength=20))
  })

  output$download <- downloadHandler(
    filename=function(){paste0("mut_counts_", Sys.Date(), ".csv")},
    content=function(file){
      df_all <- bind_rows(
        lapply(names(data_reactive()$assays), function(nm){
          if(nrow(data_reactive()$assays[[nm]])==0) return(tibble())
          data_reactive()$assays[[nm]] %>% mutate(assay=nm)
        })
      ) %>% left_join(data_reactive()$total_seq %>% rename(total_sequences=total), by="date")
      readr::write_csv(df_all, file)
    })

  # ------------------------------------------------
  # Wastewater surveillance tab — aligned x-axis
  library(cowplot)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  # Reactive x-axis limits
  ww_xaxis <- reactive({
    start <- as.Date(input$date_slider[1])
    end   <- as.Date(input$date_slider[2])
    breaks <- seq(start, end, by="4 weeks")
    list(limits = c(start, end), breaks = breaks)
  })
  
  # Combined plot output
  output$ww_combined <- renderPlot({
    df <- wwopen %>% filter(Date >= input$date_slider[1], Date <= input$date_slider[2])
    req(nrow(df) > 0)
    
    # Top: N1/N2
    p_top <- ggplot(df, aes(x = Date)) +
      geom_line(aes(y = covN1_nPMMoV_meanNr, color = "N1"), size=0.7, alpha=0.8) +
      geom_line(aes(y = covN2_nPMMoV_meanNr, color = "N2"), size=0.7, alpha=0.8) +
      scale_color_manual(values = c(N1="black", N2="gray")) +
      theme_classic() +
      scale_x_date(limits = ww_xaxis()$limits,
                   breaks = ww_xaxis()$breaks,
                   expand = c(0, 0)) +   # remove axis padding
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            legend.position = c(0.05, 0.85),  # x, y in normalized coordinates
            legend.justification = c("left", "top"),
            legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 9),
            legend.key.size = unit(0.8, "lines"),
            ) +
      labs(y = "genome copies per PMMoV copies", color = "Assay") +
      geom_label(
        data = data.frame(
        x = ww_xaxis()$limits[1] + 30, 
        y = max(df$covN1_nPMMoV_meanNr, na.rm = TRUE),          
        label = "OTTAWA WASTEWATER DERIVED METRIC"
      ), 
      aes(x = x, y = y, label = label),
      hjust = 0, 
      size = 4, 
      fill = "lightyellow",
      color = "black",
      label.size = 0
    )
    
    
    # Middle: z-score
    req("zscore" %in% colnames(df))
    p_mid <- ggplot(df, aes(Date, zscore)) +
      geom_smooth(method="loess", se=TRUE, span=input$N1N2spanslider, color="#FC4E07", fill="#FC4E07") +
      geom_point(size=0.5, alpha=0.3) +
      theme_classic() +
      scale_x_date(limits = ww_xaxis()$limits,
                   breaks = ww_xaxis()$breaks,
                   expand = c(0, 0)) +   # remove axis padding
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank()) +
      labs(y="z-score") +
      geom_label(
        data = data.frame(
        x = ww_xaxis()$limits[1] + 30, 
        y = max(df$zscore, na.rm = TRUE),          
        label = "OTTAWA WASTEWATER DERIVED METRIC"
        ), 
       aes(x = x, y = y, label = label),
       hjust = 0, 
       size = 4, 
       fill = "lightyellow",
       color = "black",
       label.size = 0
  )
    
    
    # Bottom: percent change
   # req("percentchange" %in% colnames(df))
  # p_bottom <- ggplot(df, aes(Date, percentchange)) +
  #    geom_smooth(method="loess", se=TRUE, span=input$N1N2spanslider, color="#FC4E07", fill="#FC4E07") +
  #    geom_point() +
  #    theme_classic() +
  #    scale_x_date(limits = ww_xaxis()$limits,
  #                 breaks = ww_xaxis()$breaks,
  #                 date_labels="%b %d",
  #                 expand = c(0, 0)) +   # remove axis padding
  #    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  #   labs(y="% change")
    
    # Bottom: mutational burden
    assays <- data_reactive()$assays
    total_seq <- data_reactive()$total_seq
    
    selected <- input$burden_assays
    req(length(selected) > 0)
    
    df_all <- bind_rows(lapply(selected, function(a) {
      if(!is.null(assays[[a]]) && nrow(assays[[a]]) > 0) {
        assays[[a]] %>% mutate(locus = a)
      } else {
        tibble(date=as.Date(character()), position=numeric(), count=numeric(),
               block=character(), x_pos=numeric(), locus=a)
      }
    }))
    
    req(nrow(df_all) > 0)
    
    df_all <- df_all %>%
      left_join(total_seq %>% rename(total_sequences = total), by = "date") %>%
      mutate(burden = count / pmax(total_sequences, 1)) %>%
      group_by(date, locus) %>%
      summarise(mean_burden = sum(burden, na.rm=TRUE), .groups = "drop") %>%
      arrange(date)
    
    df_all <- df_all %>%
      tidyr::pivot_wider(names_from = locus, values_from = mean_burden, values_fill = 0) %>%
      tidyr::pivot_longer(cols = all_of(selected), names_to = "locus", values_to = "mean_burden")
    
    # low coverage shading (optional)
    low_seq_df <- total_seq %>%
      arrange(date) %>%
      mutate(low_flag = total <= 50,
             group = cumsum(c(TRUE, diff(low_flag) != 0)))
    
    low_seq_ranges <- low_seq_df %>%
      filter(low_flag) %>%
      group_by(group) %>%
      summarise(xmin = min(date)-3, xmax = max(date)+3, .groups="drop")
    
    all_colors <- c(
      "N1" = "firebrick",
      "N2" = "steelblue",
      "N200" = "forestgreen",
      "E_Sarbeco" = "goldenrod"
    )
    assay_colors <- all_colors[selected]
    
    p_bottom <- ggplot(df_all, aes(x = date, y = mean_burden, fill = locus)) +
      # Highlight low coverage periods
      geom_rect(
        data = low_seq_ranges,
        aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf),
        inherit.aes = FALSE,
        fill = "darkorange",
        alpha = 0.5
      ) +
      # Stacked area plot
      geom_area(
        position = "stack",
        color = "black",
        linewidth = 0.2,
        alpha = 0.8
      ) +
      # Assay colors
      scale_fill_manual(values = assay_colors) +
      # Axis and labels
      labs(y = "# of mutations across locus", x = NULL, fill = "Assay") +
      scale_x_date(
        limits = ww_xaxis()$limits,
        breaks = ww_xaxis()$breaks,
        date_labels = "%b %d",
        expand = c(0, 0),  
        sec.axis = dup_axis(    
          breaks = seq(      
            as.Date(format(ww_xaxis()$limits[1], "%Y-01-01")),      
            as.Date(format(ww_xaxis()$limits[2], "%Y-01-01")),      
            by = "1 year"    
            ),    
          labels = function(x) format(x, "%Y"),    
          name = NULL  )
      ) +
      # Theme
      theme_classic(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = c(0.05, 0.60),  # x, y in normalized coordinates
        legend.justification = c("left", "top"),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)
      ) +
      geom_label(
      data = data.frame(
        x = ww_xaxis()$limits[1] + 30, 
        y = max(df_all$mean_burden, na.rm = TRUE),          
        label = "WORLDWIDE CLINICAL DERIVED METRIC"
      ),
      aes(x = x, y = y, label = label),
      hjust = 0, 
      size = 4, 
      fill = "lightyellow",
      color = "black",
      label.size = 0
    )   
    

    
    # Combine plots with labels
    combined <- plot_grid(
      p_top, p_mid, p_bottom,
      ncol = 1, align = "v", axis = "lr",
      rel_heights = c(1, 1, 1),
      labels = c("A", "B", "C"), label_size = 16
    )
    
    # Caption as its own grob
    caption <- ggdraw() +
      draw_label(
        "Figure. SARS-CoV-2 N1 and N2 PCR assay accuracy in Ottawa Wastewater surveillance.\n
    Panel A: PMMoV-normalized concentrations of N1 and N2 gene targets (genome copies per PMMoV copies).\n
    Panel B: Smoothed z-score of N1 vs. N2 signal indicating relative changes over time. \n
             Negative standard deviations indicate underestimated N1 signal relative to N2 \n
    Panel C: Mutational burden across selected loci, stacked by assay, with shaded regions indicating \n
    low sequencing coverage.Genomic seqeunces are aggregated from Genbank via CoVSpectrum.",
        x = 0.05, y = 0.9,  # top-left inside caption box
        hjust = 0, vjust = 1,
        size = 10,
        lineheight = 0.4
      )
    
    # Stack plots and caption
    final_plot <- plot_grid(
      combined, caption,
      ncol = 1,
      rel_heights = c(0.9, 0.1)  # adjust space for caption
    )
    
    final_plot
    
    
    
  })
  
  
}

# --------------------------------------------------------------
shinyApp(ui, server)
