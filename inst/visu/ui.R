library(tidyverse)
library(HPVcap)
library(shiny)

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
                     sliderInput("nreads", "Nombre de reads minimum", 0, 0, 0, step = 1, round = T),
                     sliderInput("match", "Score minimal de matching", 0, 0, 0, step = 1, round = T))),
        mainPanel(plotlyOutput("plot"),
                  DT::dataTableOutput("table")))
