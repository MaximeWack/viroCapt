library(shiny)

ui <- bootstrapPage(
  sidebarLayout(
    sidebarPanel(
      titlePanel("Viral integration pipeline results"),
      inputPanel(fileInput("visu", "Visualisation file", accept = "*.rds")),
      inputPanel(selectizeInput("genotype",
                                "Genotypes",
                                "",
                                multiple = T),
                 selectizeInput("scores",
                                "Scores",
                                "",
                                multiple = T),
                 selectizeInput("chrs",
                                "Chromosomes",
                                "",
                                multiple = T),
                 sliderInput("nreads", "Minimum number of reads", 0, 0, 0, step = 1, round = T),
                 sliderInput("match", "Minimum match score", 0, 0, 0, step = 1, round = T)),
      width = 3),
    mainPanel(plotly::plotlyOutput("plot"),
              DT::dataTableOutput("table"),
              width = 9)))
