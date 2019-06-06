library(tidyverse)
library(HPVcap)
library(shiny)

options(shiny.maxRequestSize = 1024 ^ 3)

UIlabels <- function(variable)
{
  variable %>%
    fct_infreq %>%
    summary -> varsumm

  names(varsumm) %>%
    setNames(str_c(names(varsumm), " (", varsumm, ")"))
}

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

  # Read blat summary
  summ_blat <- reactive({
    req(input$summ)

    input$summ$datapath %>%
      HPVcap:::read_summary()
  })

  # Set UI values
  observe({
    sam() %>%
      group_by(genotype) %>%
      summarise(n = max(n)) -> genotypes

    genotypes$genotype %>%
      setNames(str_c(genotypes$genotype, " (", genotypes$n, ")")) -> genotypes

    updateSelectizeInput(session, "genotype", choices = genotypes) #, selected = genotypes[1])
  })

  observe(updateSelectizeInput(session, "scores", choices = summ_blat()$quality %>% UIlabels))
  observe(updateSelectizeInput(session, "chrs", choices = summ_blat()$chr %>% UIlabels))
  observe(updateSliderInput(session, "nreads", min = summ_blat()$n %>% min, max = summ_blat()$n %>% max))
  observe(updateSliderInput(session, "match", min = summ_blat()$match %>% min, max = summ_blat()$match %>% max))
}
