library(viroCapt)
library(ggplot2)
library(shiny)

options(shiny.maxRequestSize = 1024 ^ 3)

UIlabels <- function(variable)
{
  variable %>%
    forcats::fct_infreq() %>%
    summary -> varsumm

  names(varsumm) %>%
    stats::setNames(paste0(names(varsumm), " (", varsumm, ")"))
}

server <- function(input, output, session)
{
  # Read the visualisation object
  visu <- reactive({
    req(input$visu)

    readRDS(input$visu$datapath)
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
      dplyr::mutate(quality = quality %>% forcats::fct_drop())
  })

  # Set UI values
  observe({
    sam() %>%
      dplyr::group_by(genotype) %>%
      dplyr::summarise(n = max(n)) %>%
      dplyr::mutate(genotype = genotype %>% as.character) -> genotypes

    genotypes$genotype %>%
      stats::setNames(paste0(genotypes$genotype, " (", genotypes$n, ")")) -> genotypes

    updateSelectizeInput(session, "genotype", choices = genotypes, selected = genotypes[1])
  })

  observe({
    updateSelectizeInput(session, "scores", choices = summ_blat_file()$quality %>% UIlabels)
    updateSelectizeInput(session, "chrs", choices = summ_blat_file()$chr %>% UIlabels)
    updateSliderInput(session, "nreads", min = summ_blat_file()$n %>% min, max = summ_blat_file()$n %>% max)
    updateSliderInput(session, "match", min = summ_blat_file()$match %>% min, max = summ_blat_file()$match %>% max)
  })

  downsampled_depths <- reactive({
    sam() %>% viroCapt:::downsample(10)
  })

  depths <- reactive({
    req(input$genotype)

    downsampled_depths() %>%
      dplyr::filter(genotype %in% input$genotype)
  })

  summ_blat <- reactive({
    summ_blat_file() %>%
      dplyr::filter(genotype %in% input$genotype,
                    n >= input$nreads,
                    match >= input$match) -> summ

    if (! input$chrs %>% is.null)
      summ %>% dplyr::filter(chr %in% input$chrs) -> summ

    if (! input$scores %>% is.null)
      summ %>% dplyr::filter(quality %in% input$scores) -> summ

    summ
  })

  localplot <- reactive({
    viroCapt:::ggplot_depth(depths())
  })

  output$plot <- plotly::renderPlotly({
    max(depths()$n) -> max_n

    localplot() +
      geom_vline(data = summ_blat(), aes(xintercept = position, color = chr, alpha = quality)) -> p

    if (summ_blat() %>% dplyr::filter(feature == "left") %>% nrow > 0)
      p + geom_segment(data = summ_blat() %>% dplyr::filter(feature == "left"), aes(x = position, xend = position - 100, y = -max_n/20, yend = -max_n/20, color = chr, alpha = quality)) -> p

    if (summ_blat() %>% dplyr::filter(feature == "right") %>% nrow > 0)
      p + geom_segment(data = summ_blat() %>% dplyr::filter(feature == "right"), aes(x = position, xend = position + 100, y = -max_n/20, yend = -max_n/20, color = chr, alpha = quality)) -> p

    p +
      scale_alpha_ordinal(range = c(0.5, 1)) +
      theme_classic()
  })

  output$table <- DT::renderDataTable(
  {
    summ_blat() %>%
      dplyr::mutate_at(vars(feature, chr, quality), factor)
  },
  options = list(dom = "Bfrtip",
                 buttons = c("copy", "excel"),
                 paging = T,
                 info = F,
                 fixedHeader = T),
  filter = "top",
  extensions = "Buttons")
}
