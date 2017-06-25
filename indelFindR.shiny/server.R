library(shiny)
source("indel.R")
shinyServer(function(input, output) {
  observeEvent(input$How, {
      showModal(modalDialog(
        title = "How it works:",
        "The simple script tries to locate the first half of the gRNA sequence, and reads 100 bp downstream to assess frameshift. Indels up to 35 bp are compared to sequencing signal peaks at each call position ",
        easyClose = TRUE,
        footer = NULL
      ))
    })
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
  observe({
    if (input$close > 0) stopApp()                             # stop shiny
    })
})