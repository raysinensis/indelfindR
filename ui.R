library(shiny)
shinyUI(fluidPage(
  titlePanel("indelFindR"),
  sidebarLayout(
    sidebarPanel(
      p("INPUT"),
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
   		  actionButton("button", label = "Submit", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))
      )
	  
	),
    mainPanel(
	  p("indelFindR will search for potential indels (or WT) from sequencing results"),
	  br(),
	  p("Please select ABI file(s)"),
	  p("Please also input the gRNA target sequence, and desired genomic area ~500bp"),
	  br(),
	  h2("Analysis Results"),
	  h5(verbatimTextOutput("text1"))
    )
  )
))