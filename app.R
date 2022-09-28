library(shiny)
library(DT)
load("AllAppData.R")

data_descr = p("codes ",tags$b(".fxf")," and ",tags$b(".dtc")," refers to egg study and high resolution time course study respectively.")

fxf_ref= tags$li("Mother-Specific Signature in the Maternal Transcriptome Composition of Mature, Unfertilized Zebrafish Eggs",
            a(href="http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0147151","publication" ))
dtc_ref= tags$li("Transcriptome dynamics in early zebrafish embryogenesis determined by high-resolution time course analysis",a(href="http://","publication"),"(to be added later)")


ui <- shinyUI(fluidPage(
     tags$style(type="text/css",".datatables { font-size: 9px; } p { margin: 3px; }"),
     titlePanel("Zebrafish embryogenesis gene-expression browser"),
     fluidRow(column(12,p("Cite one of the following papers if the data is used."),tags$ul(fxf_ref,dtc_ref),data_descr)),
     fluidRow( column(12,
             actionButton("tabReset","Reset gene table"),"use after selection by correlation",
             downloadButton('downloadData','Get gene table'))),
     fluidRow(
          column(3,plotOutput('x2',width="100%", height=360 )), 
          column(5,plotOutput('x3',width="100%",height=360)),
          column(4,tableOutput("geninfo"))),
     fluidRow(
             column(4, actionButton("fxfSim","100 most correlating genes"), 
             p("(calculation takes about 10 seconds)")) ,
             column(4, actionButton("dtcSim","100 most correlating genes"), 
             p("(calculation takes about 20 seconds)")) 
             ),
     fluidRow(column(12,h4("Use searchbox to find your gene of interest and select it in the table to show graph",
                            style="color: red; text-align: right"),
                        DT::dataTableOutput('x1')))
             
))


colour <- rainbow(8)[c(2, 3, 5, 7, 8)]

getCorr.fxf <- function(n) { apply(fxf.data,1,function(x) {cor(x,as.numeric(fxf.data[n,]))} ) }
getCorr.dtc <- function(n) { apply(dtc.data,1,function(x) {cor(x,as.numeric(dtc.data[n,]))} ) }

plot_fxf <- function(ydat=rep(1,24),title="",note="") {
   plot(x=1:24, y=ydat,  ylim=c(6, 18), col=colour[design$Mother], pch=19, cex=1, 
        main=title, xlab="mother", ylab="Signal intensity (log2)", labels=F, bty="l", xaxt='n', yaxt='n')
   if (note != "") { text(12,12,note) } 
}
dtc.colors=unique(dtc.design$SpawnCol)
plot_dtc <- function(ydat=rep(0,180),title="",note="") {

   plot(ydat,pch=20,main=title,col=dtc.design$SpawnCol,ylab="Signal intensity (log2)",xlab="Embryos in developmental order",ylim=c(8,18),xlim=c(0,180))
   if (note != "") { text(90,12,note) } 
   legend("topleft",legend=1:9,col=dtc.colors,pch=20,cex=0.7,title="spawn",horiz=T)
}

server <- shinyServer(function(input, output,session) {
   curTable <- Annotation; 

   observeEvent(input$tabReset,{ 
    curTable <<- Annotation; 
    output$x1 <- renderDataTable(curTable,server=TRUE,selection='single')  
   })  

   observeEvent(input$fxfSim,{ 
    sel.row =  input$x1_row_last_clicked;
    gen = rownames(curTable)[sel.row]
    if  (!is.null(sel.row) & (gen %in% rownames(fxf.data)) ) {
      corv = sort(getCorr.fxf(gen),decreasing=T)[1:100]
      nti <- names(corv)
      sel <- !nti %in%rownames(Annotation)
      nti[sel]  <- paste(nti[sel],".fxf",sep="");
      newTable <- Annotation[nti,]; 
      colname=paste("correlation_to",gen,sep="_")
      newTable[,colname] <- corv
      curTable <<- newTable
      output$x1 <- renderDataTable(curTable,server=TRUE,selection='single')  
    } 
   })  

   observeEvent(input$dtcSim,{ 
    sel.row =  input$x1_row_last_clicked;
    gen = rownames(curTable)[sel.row]
    if  (!is.null(sel.row) & (gen %in% rownames(dtc.data)) ) {
      corv = sort(getCorr.dtc(gen),decreasing=T)[1:100]
      nti <- names(corv)
      sel <- !nti %in%rownames(Annotation)
      nti[sel]  <- paste(nti[sel],".dtc",sep="");
      newTable <- Annotation[nti,]; 
      colname=paste("correlation_to",gen,sep="_")
      newTable[,colname] <- corv
      curTable <<- newTable
      output$x1 <- renderDataTable(curTable,server=TRUE,selection='single')  
    } 
   })  

  output$downloadData <- downloadHandler(
    filename = "Dr-gene-table.tsv",
    content = function(file) { write.table(curTable, file,row.names=FALSE,quote=FALSE,sep='\t') }
  )

  output$x2 <- renderPlot({ 
      sel.row = input$x1_row_last_clicked
      if (!is.null(sel.row)) {
          sel.gen = rownames(curTable)[sel.row]
          s.sel.gen = gsub(".fxf$","",sel.gen)
          d.sel.gen = paste(gsub(".dtc","",sel.gen),".fxf",sep="")
          if (sel.gen %in% rownames(fxf.data))  
               {  plot_fxf(ydat=fxf.data[sel.gen,],title=paste(s.sel.gen,"in Egg")) } 
          else { 
                  if  (d.sel.gen %in% rownames(fxf.data)) { 
	          plot_fxf(note="different transcript expressed in egg")
                  } else {
          
                  plot_fxf(note=paste("gene",sel.gen,"not expressed in egg")); 
	 	  }
          }
      }
      else {
          plot_fxf(note="no gene selected")
      }
      axis(side = 1, at= c((seq(1,20,by=5)+2), 22.5), las=1,lwd=1,lwd.ticks=1,cex.axis=1,labels=1:5)
      axis(side = 2, las=1, cex.axis=1, lwd=1, lwd.ticks=1)
  });

  output$x3 <- renderPlot({
      sel.row = input$x1_row_last_clicked
      if (!is.null(sel.row)) {
        sel.gen = rownames(curTable)[sel.row]
        d.sel.gen = gsub(".dtc$","",sel.gen)
        e.sel.gen = paste(gsub(".fxf","",sel.gen),".dtc",sep="")

        if (sel.gen %in% rownames(dtc.data)) { 
         	plot_dtc(ydat=as.numeric(dtc.data[sel.gen,]),title=paste(d.sel.gen,"in time course")) }
        else {
          if (e.sel.gen %in% rownames(dtc.data)) { 
	          plot_dtc(note="different transcript expressed in time course")
          } else {
	          plot_dtc(note=paste("gene ",sel.gen,"not expressed in time course")) 
          }
        }
      } else {
          plot_dtc(note="no gene selected")
      }
   
  })

  output$x1 <- renderDataTable(curTable,server=TRUE,selection='single')
  output$geninfo <- renderTable({ 
      if (!is.null(input$x1_row_last_clicked)) {
            t(Annotation[input$x1_row_last_clicked,])
      } 
 },rownames=TRUE ) 

})

shinyApp(ui = ui, server = server)