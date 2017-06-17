library(shiny)
source("indel.R")
shinyServer(function(input, output) {
  gRNA <- eventReactive(input$button, {
    input$gRNA
    })
  seq <- eventReactive(input$button, {
    input$seq
    })
  filenames <- eventReactive(input$button, {
    input$file$name
    })
  filenumbers <- eventReactive(input$button, {
    length(input$file$name)
    })
  filelist <- eventReactive(input$button, {
    input$file$datapath
    })
  found <- eventReactive(input$button, {
    indel(filelist(),filenames(),gRNA(),seq())
	})
  output$text1<-renderText({found()
    })
})