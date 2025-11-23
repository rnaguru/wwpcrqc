# ==========================================================
# Prefetch SARS-CoV-2 qPCR mutation counts + total sequences
# ==========================================================

library(tidyverse)
library(future.apply)
library(httr)
plan(multisession)

# --------------------------
# Define assays (positions + wildtype bases)
# --------------------------
assays <- list(
  N1 = list(
    positions = c(28287:28306, 28309:28332, 28335:28358),
    wt = c(
      strsplit("GACCCCAAAATCAGCGAAAT", "")[[1]],
      strsplit("ACCCCGCATTACGTTTGGTGGACC", "")[[1]],
      strsplit("TCTGGTTACTGCCAGTTGAATCTG", "")[[1]]
    )
  ),
  N2 = list(
    positions = c(29164:29183, 29188:29210, 29213:29230),
    wt = c(
      strsplit("TTACAAACATTGGCCGCAAA", "")[[1]],
      strsplit("ACAATTTGCCCCCAGCGCTTCAG", "")[[1]],
      strsplit("GCGCGACATTCCGAAGAA", "")[[1]]
    )
  ),
  E_Sarbeco = list(
    positions = c(26269:26294, 26332:26357, 26360:26381),
    wt = c(
      strsplit("ACAGGTACGTTAATAGTTAATAGCGT", "")[[1]],
      strsplit("ACACTAGCCATCCTTACTGCGCTTCG", "")[[1]],
      strsplit("ATATTGCAGCAGTACGCACACA", "")[[1]]
    )
  ),
  N200 = list(
    positions = c(28840:28861, 28891:28905, 28940:28960),
    wt = c(
      strsplit("TAGTCGCAACAGTTCAAGAAAT", "")[[1]],
      strsplit("TCCTGCTAGAATGGC", "")[[1]],
      strsplit("CTGGTTCAATCTGTCAAGCAG", "")[[1]]
    )
  )
)

# --------------------------
# Fetch mutation counts from LAPIS
# --------------------------
fetch_nonwt_lapis <- function(positions, wt_bases, start_date, end_date){
  lapis_url <- "https://lapis.cov-spectrum.org/open/v2/sample/aggregated"
  all_dates <- seq.Date(as.Date(start_date), as.Date(end_date), by="week")
  
  df_list <- future_lapply(seq_along(positions), function(i){
    pos <- positions[i]
    ref <- wt_bases[i]
    variant_query <- paste0(ref, pos)
    
    params <- list(
      variantQuery = variant_query,
      fields = "date",
      dateFrom = start_date,
      dateTo = end_date
    )
    
    resp <- httr::RETRY("GET", lapis_url, query = params, times=2)
    if(resp$status_code != 200) return(tibble(date = all_dates, position = pos, count = 0))
    
    dat <- httr::content(resp, as="parsed", simplifyVector = FALSE)
    if(is.null(dat$data) || length(dat$data) == 0){
      tibble(date = all_dates, position = pos, count = 0)
    } else {
      df <- tibble(
        date = as.Date(sapply(dat$data, function(x) x$date)),
        count = sapply(dat$data, function(x) x$count)
      )
      tibble(date = all_dates) %>%
        left_join(df, by="date") %>%
        mutate(count = ifelse(is.na(count), 0, count),
               position = pos)
    }
  }, future.seed = TRUE)
  
  bind_rows(df_list)
}

# --------------------------
# Fetch total sequences (per week)
# --------------------------
fetch_total_sequences <- function(start_date, end_date) {
  lapis_url <- "https://lapis.cov-spectrum.org/open/v2/sample/aggregated"
  all_dates <- seq.Date(as.Date(start_date), as.Date(end_date), by="week")
  
  params <- list(
    fields = "date",
    dateFrom = start_date,
    dateTo = end_date
  )
  
  resp <- httr::RETRY("GET", lapis_url, query = params, times = 2)
  if (resp$status_code != 200) {
    return(tibble(date = all_dates, total = 0))
  }
  
  dat <- httr::content(resp, as = "parsed", simplifyVector = FALSE)
  
  if (is.null(dat$data) || length(dat$data) == 0) {
    tibble(date = all_dates, total = 0)
  } else {
    df <- tibble(
      date = as.Date(sapply(dat$data, function(x) x$date)),
      total = sapply(dat$data, function(x) x$count)
    )
    tibble(date = all_dates) %>%
      left_join(df, by = "date") %>%
      mutate(total = ifelse(is.na(total), 0, total))
  }
}

# --------------------------
# Prefetch and save CSVs
# --------------------------
start_date <- "2020-01-01"
end_date <- Sys.Date()

data_dir <- "LociTracker/data"
if(!dir.exists(data_dir)) dir.create(data_dir)

for(nm in names(assays)){
  df <- fetch_nonwt_lapis(assays[[nm]]$positions, assays[[nm]]$wt, start_date, end_date)
  write_csv(df, file.path(data_dir, paste0("prefetched_", nm, ".csv")))
}

total_seq <- fetch_total_sequences(start_date, end_date)
write_csv(total_seq, file.path(data_dir, "prefetched_total_sequences.csv"))

message("Prefetching complete. CSVs saved in '", data_dir, "'")
