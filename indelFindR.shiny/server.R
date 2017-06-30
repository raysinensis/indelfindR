library(shiny)
source("indel.R")
source("helpers.R")
shinyServer(function(input, output, session) {
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
  sthreshold <- eventReactive(input$button, {
    input$threshold
    })
  srange <- eventReactive(input$button, {
    input$range
    })
  slength <- eventReactive(input$button, {
    input$length
    })
  found <- eventReactive(input$button, {
    withBusyIndicatorServer("button", {indel(filelist(),filenames(),gRNA(),seq(),as.numeric(sthreshold()),as.numeric(srange()),as.numeric(slength()))})
	})
  output$text1<-renderText({found()
    })
  observe({
    if (input$close > 0) stopApp()                             # stop shiny
    })
  session$onSessionEnded(stopApp)
})
