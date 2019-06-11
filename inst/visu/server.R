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
        HPVcap:::read_depth()
    })
  })

  # Read blat summary
  summ_blat_file <- reactive({
    req(input$summ)

    input$summ$datapath %>%
      HPVcap:::read_summary() %>%
      mutate(quality = quality %>% fct_drop)
  })

  # Set UI values
  observe({
    sam() %>%
      group_by(genotype) %>%
      summarise(n = max(n)) %>%
      mutate(genotype = genotype %>% as.character) -> genotypes

    genotypes$genotype %>%
      setNames(str_c(genotypes$genotype, " (", genotypes$n, ")")) -> genotypes

    updateSelectizeInput(session, "genotype", choices = genotypes, selected = genotypes[1])

  })

  observe({
    updateSelectizeInput(session, "scores", choices = summ_blat_file()$quality %>% UIlabels)
    updateSelectizeInput(session, "chrs", choices = summ_blat_file()$chr %>% UIlabels)
    updateSliderInput(session, "nreads", min = summ_blat_file()$n %>% min, max = summ_blat_file()$n %>% max)
    updateSliderInput(session, "match", min = summ_blat_file()$match %>% min, max = summ_blat_file()$match %>% max)
  })

  depths <- reactive({
    req(input$genotype)

    sam() %>%
      filter(genotype %in% input$genotype)
  })

  summ_blat <- reactive({
    summ_blat_file() %>%
      filter(genotype %in% input$genotype,
             n >= input$nreads,
             match >= input$match) -> summ

    if (! input$chrs %>% is.null)
      summ %>% filter(chr %in% input$chrs) -> summ

    if (! input$scores %>% is.null)
      summ %>% filter(quality %in% input$scores) -> summ

    summ
  })

  output$plot <- renderPlot({
    HPVcap:::ggplot_depth(depths())
  })

  })
}

shinyApp(ui, server)
