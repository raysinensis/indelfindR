library(shiny)
library(shinyjs)
source("helpers.R")
shinyUI(fluidPage(
  useShinyjs(),
  tags$style(appCSS),
  titlePanel("indelFindR"),
  sidebarLayout(
    sidebarPanel(
    tabsetPanel(id='navbar',
      tabPanel(title='INPUT',value='INPUT',
	  fluidRow(
        column(12,
          fileInput("file", label = h5("select file(s)"), multiple = TRUE),
          textInput("gRNA", label = h5("gRNA target sequence:"), 
            value = "CCAAATAGGGTCATCTTCTA"))
	  ),
	  fluidRow(
	    column(12, 
          textAreaInput("seq", label = h5("genomic sequence around cut site:"), rows = 10, resize = "none",
            value = "GGACAAGAAGACCGTTGACTGTGTTTAAGAACAGCATTTACAAGGAAATAAGCAACGCTAATGTTAAGAAATTAATATATTACAATATAAGTGGTTCTCCTTTTTTATTTTTCCTTTCCCAAATAGGGTCATCTTCTATGGCATCTGCATGTGGTGGAAGTTTGGCATTAATGGATGCAGGTAAAGATTATGTCTCAAAACTGTTAGTATTTTACATATATATTCATTAGTACAGAATATAAAATTTGTTTAGATTGGATTTAATAAAAGTCAAAACCGTGGGCTAGATACAATATAATTACATCAGAAAGATATTAGCTGAAGTTCTTATTTTTAGCTAAAAGACATTGGAAATTTATAATATACTGAAAATTATAAACACTTGCCTTGGTGACCAAATGCACTTTCTCACC"))
	  ),
	  fluidRow(
	    column(12, 
   		  withBusyIndicatorUI(actionButton("button", label = "Submit", class = "btn-primary", style="color: #fff; background-color: #337ab7; border-color: #2e6da4")))
      )
	  
	),
      tabPanel(title='Settings',value='Settings',
        fluidRow(
	    column(12,
	    textInput("threshold", label=h5("set threshold (0.05 works well), or 0 for autothreshold"),value = "0"),
	    textInput("range", label=h5("set indel search range"),value = "35"),
	    textInput("length", label=h5("set search length"),value = "20")))	
	))),
    mainPanel(
	  div(style="display: inline-block;vertical-align:top; width: 60px;",tags$button(id = 'close', type = "button", class = "btn action-button", onclick = "setTimeout(function(){window.close();},23);", "Close"),
          p(" ")),
	  div(style="display: inline-block;vertical-align:top; width: 150px;",actionButton("How", "How it works")),
	  p("indelFindR will search for potential indels (or WT) from sequencing results"),
	  br(),
	  p("Please select AB1 file(s)"),
	  p("Please also input the gRNA target sequence, and the desired genomic area ~400bp"),
	  br(),
	  h2("Analysis Results"),
	  h5(verbatimTextOutput("text1"))
    )
  )
))
