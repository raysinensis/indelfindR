library(shiny)
library(ggplot2)
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
  observeEvent(input$default, {
      updateTextInput(session, "threshold", value = "0")
      updateTextInput(session, "range", value = "35")
      updateTextInput(session, "length", value = "20")
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
  output$text1<-renderText({found()[[1]]
    })
  output$table1<-renderTable({found()[[3]]},colnames=FALSE,na=''
    )
  fig <- eventReactive(input$button, {
    df=as.data.frame(as.matrix(found()[[2]]))
    out=ggplot(data=df,aes(x=row.names(df), y=freq))+
    geom_bar(stat="identity", width = .7)+scale_fill_grey() +
    theme_bw()+theme(axis.title.x=element_blank(),text=element_text(family="serif",size=15))
    out
    })
  output$plot1<-renderPlot({
	#barplot(t(as.matrix(found()[[2]])),main="genotyping results", xlab="mutation type",col=c("darkblue"))
	fig()
    })
  output$Download <- downloadHandler(
    filename = "seq.png",
    content = function(file) {
        ggsave(file, plot = fig(), device = "png")
    })
  observe({
    if (input$close > 0) stopApp()                             # stop shiny
    })
  session$onSessionEnded(stopApp)
})
