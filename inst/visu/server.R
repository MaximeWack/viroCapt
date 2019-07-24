library(tidyverse)
library(HPVcap)
library(shiny)
library(plotly)

options(shiny.maxRequestSize = 1024 ^ 3)
# options(shiny.reactlog = TRUE)

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
  # Read the visualisation object
  visu <- reactive({
    req(input$visu)

    input$visu$datapath %>%
      readRDS
  })

  # Extract the sequencing depth object
  sam <- reactive({
    req(visu())

    visu()$profile
  })

  # Extract the blat summary
  summ_blat_file <- reactive({
    req(visu())

    visu()$summary %>%
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

  downsampled_depths <- reactive({
    sam() %>%
      HPVcap:::downsample(10)
  })

  depths <- reactive({
    req(input$genotype)

    downsampled_depths() %>%
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

  localplot <- reactive({
    HPVcap:::ggplot_depth(depths())
  })

  output$plot <- renderPlotly({
    max(depths()$n) -> max_n

    localplot() +
      geom_vline(data = summ_blat(), aes(xintercept = position, color = chr, alpha = quality)) -> p

    if (summ_blat() %>% filter(feature == "left") %>% nrow > 0)
      p + geom_segment(data = summ_blat() %>% filter(feature == "left"), aes(x = position, xend = position - 100, y = -max_n/20, yend = -max_n/20, color = chr, alpha = quality)) -> p

    if (summ_blat() %>% filter(feature == "right") %>% nrow > 0)
      p + geom_segment(data = summ_blat() %>% filter(feature == "right"), aes(x = position, xend = position + 100, y = -max_n/20, yend = -max_n/20, color = chr, alpha = quality)) -> p

    p +
      scale_alpha_ordinal(range = c(0.5, 1)) +
      theme_classic()
  })

  output$table <- DT::renderDataTable(
  {
    summ_blat() %>%
      mutate_at(vars(feature, chr, quality), factor)
  },
  options = list(dom = "Bfrtip", 
                 buttons = c("copy", "excel"),
                 paging = F,
                 info = F,
                 fixedHeader = T),
  filter = "top")
}

shinyApp(ui, server)
