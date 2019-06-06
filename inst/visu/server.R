library(tidyverse)
library(HPVcap)
library(shiny)

options(shiny.maxRequestSize = 1024 ^ 3)

server <- function(input, output, session)
{
  # Read and parse the sam file
  sam <- reactive({
    req(input$sam)

    withProgress({
      incProgress(message = "Reading sam file")
      input$sam$datapath %>%
        HPVcap:::read_sam() -> sam

      incProgress(message = "Parsing CIGARs")
      sam %>%
        mutate(parsed = cigar %>% HPVcap:::parse_cigar()) %>%
        unnest -> sam

      incProgress(message = "Computing nucleotide read depth")
      sam %>%
        HPVcap:::read_depth() -> sam

      incProgress(message = "Downsampling for displaying")
      sam %>%
        HPVcap:::downsample()
    })
  })
}
