library(tidyverse)
library(HPVcap)
library(shiny)
library(plotly)


ui <- bootstrapPage(
        sidebarPanel(
          titlePanel("HPV integration pipeline results"),
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
                     sliderInput("match", "Minimum match score", 0, 0, 0, step = 1, round = T))),
        mainPanel(plotlyOutput("plot"),
                  DT::dataTableOutput("table")))
