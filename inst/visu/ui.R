library(tidyverse)
library(HPVcap)
library(shiny)

ui <- bootstrapPage(
        sidebarPanel(
          titlePanel("HPV integration pipeline results"),
          inputPanel(fileInput("sam", "SAM file", accept = "*.sam"),
                     fileInput("summ", "Blat summary file", accept = "*_summary.tsv")),
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
        mainPanel(plotOutput("plot"),
                  DT::dataTableOutput("table")))
