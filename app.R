#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(BiocManager)
options(repos = BiocManager::repositories())
library(devtools)
#install_github("sheffield-bioinformatics-core/sbcMisc")
library(ggbeeswarm)
library(tidyverse)
library(RColorBrewer)

anno <- readRDS("anno.rds")
bcode_info <- readRDS("barcodes.rds")


# Define UI for application that draws a histogram
ui <- navbarPage("Explore the single-cell data...",
  
   # Sidebar with a slider input for number of bins 
      tabPanel("Compare Features",
      sidebarLayout(         
      sidebarPanel(
        img(src="logo-sm.png"), br(),a("sbc.shef.ac.uk",href="http://sbc.shef.ac.uk"),
         selectInput("sel_genes","Select your genes..",choices = anno$mgi_symbol,selected = "",selectize = TRUE,multiple = TRUE),
         checkboxInput("exclude_NA","Exclude Missing Values",value=TRUE),
         checkboxInput("flip_axis","Rotate Axis?",value=FALSE),
         checkboxGroupInput(inputId = "plot_type",label="Type of Plot", choices=c("Boxplot","Points","Violin"),selected="Boxplot"),
        selectInput(inputId = "selected_clus",label="Select Clusters to Compare",0:17, multiple = TRUE,selectize = TRUE,selected=c(0,1,3:5,7:15,17)),
        selectInput("umap_type","Color UMAP by...",choices=c("expression","cluster"),selected = "expresesion"),
        checkboxInput("include_type","Breakdown into Cell Type", value=TRUE),
        selectInput("cell_types", "Choose Cell Types...",choices = unique(bcode_info$Cell.Group),selected = unique(bcode_info$Cell.Group), selectize = TRUE,multiple = TRUE),
        h2("Output options"),
        textInput("profileBasename", label = "What to call the output files",value="boxplot"),
        radioButtons("profilePlotFormat", "File format for plots", choices=c("pdf","png"), selected="pdf"),
        helpText("PDF can be imported into Illustrator (or similar) for editing. PNG plots are suitable for presentation"),
        helpText("PDF dimensions are measured in inches, and PNG dimensions are measured in pixels"),
        textInput("profileWidth", label="Width of plot ",value = 1200),
        textInput("profileHeight", label="Height of plot ",value=600),
        downloadButton("dl_boxplot","Export current boxplot(s)...."),
        downloadButton("dl_umap","Export current UMAP(s)....")
        # selectInput("boxplotTheme", "Pick a plot style",choices=c("bw","grey","classic","minimal","light","Wall Street Journal","Economist", "Excel","solarized","stata","calc","dark","fivethirtyeight","tufte"),selected="bw"),
        # selectInput("pal_name", "Choose a colour palette", choices = c("UOS","Set1","Set2","Set3","Pastel1","Pastel2","Paired","Dark2","Accent"),selected="UOS")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        tabsetPanel(
         tabPanel("Boxplots", plotOutput("gene_boxplot")),
         tabPanel("UMAP",plotOutput("umap_plot")),
         tabPanel("About",h3("App developed by Sheffield Bioinformatics Core"),br(),
                  img(src="logo-sm.png"), br(),a("sbc.shef.ac.uk",href="http://sbc.shef.ac.uk"),br(),
                 helpText("For queries about the app and its usage, please contact Sheffield Bioinformatics Core"),br(),
                 helpText("For queries about the data and its interpretation, please contact the study authors"),
                 h2("Chen, P., Qin, L., Li, G. et al. Endothelial TGF-β signalling drives vascular inflammation and atherosclerosis. Nat Metab 1, 912–926 (2019) doi:10.1038/s42255-019-0102-3"))
        )
      )
      )
      )
   )

# Define server logic required to draw a histogram
server <- function(input, output,session) {
  library(tidyverse)   

  counts <- readRDS("counts.rds")
  bcode_info <- readRDS("barcodes.rds")
  anno <- readRDS("anno.rds")

  observeEvent(input$profilePlotFormat,{
    
    if(input$profilePlotFormat == "pdf"){
      updateTextInput(session, "profileWidth",value=12)
      updateTextInput(session, "profileHeight",value=6.3)
    } else{
      updateTextInput(session, "profileWidth",value=1200)
      updateTextInput(session, "profileHeight",value=600)
      
    }
    
  })


  output$gene_boxplot <- renderImage({
    
    p1 <- prepareBoxplot()
    width  <- session$clientData$output_gene_boxplot_width
    height <- session$clientData$output_gene_boxplot_height
    
    nrows <- length(input$sel_genes)
    
    ### stop an error occurring if nothing is selected
    
    if(is.na(nrows) | nrows == 0) nrows <- 1
    height <- height * nrows
    
    # For high-res displays, this will be greater than 1
    pixelratio <- session$clientData$pixelratio
    
    # A temp file to save the output.
    outfile <- tempfile(fileext='.png')
    
    # Generate the image file
    png(outfile, width=width*pixelratio, height=height*pixelratio,
        res=72*pixelratio)
    
    print(p1)
    dev.off()
    list(src = outfile,
         width = width,
         height = height,
         alt = "This is alternate text")
  }
  )
 
  output$dl_boxplot <- downloadHandler(
    filename = function(){
      paste0(input$profileBasename,".",input$profilePlotFormat)
    },
    content = function(file) {
      
      if(input$profilePlotFormat == "pdf") pdf(file, width=as.numeric(input$profileWidth), height=as.numeric(input$profileHeight))
      else png(file, width=as.numeric(input$profileWidth),height=as.numeric(input$profileHeight))
      
      
      
      p1 <- prepareBoxplot()
      
      print(p1)
      dev.off()
    }
    
    
  )
  
  
  
  
#   output$gene_boxplot <- renderPlot({

  
    prepareBoxplot <- reactive({
  
     p <- ggplot()

#  sel_genes <- ifelse(is.null(input$gene_tbl_rows_selected),"ACTB", 
 # dplyr::slice(anno, input$gene_tbl_rows_selected) %>% pull(hgnc_symbol) %>%  as.character)

  if(!is.null(input$sel_genes)){
    sel_genes <- filter(anno, mgi_symbol %in% input$sel_genes) %>% pull(ensembl_gene_id)
     

  
  data <-  data.frame(counts) %>% tibble::rownames_to_column("ensembl_gene_id") %>% 
    filter(ensembl_gene_id %in% sel_genes) %>% 
    tidyr::gather("Cell.ID","Expression",-ensembl_gene_id) %>% 
    left_join(bcode_info) %>% left_join(anno) %>% 
    filter(Cell.ClusterLabel %in% as.numeric(input$selected_clus)) %>% 
    filter(Cell.Group %in% input$cell_types)

      p <- data %>% 
      ggplot(aes(x=Cell.ClusterLabel,y=Expression,fill=factor(Cell.ClusterLabel),group=Cell.ClusterLabel)) + 
        scale_fill_manual("Cluster",values= c("0"=rgb(3,4,4,maxColorValue = 255),#0
                                    "1"=rgb(242,235,15, maxColorValue = 255),#1
                                    "2"=rgb(89,202,233,maxColorValue = 255),#2
                                    "3"=rgb(186,80,153,maxColorValue = 255),#3
                                    "4"=rgb(253,40,58,maxColorValue = 255),#4
                                    "5"=rgb(7,123,55,maxColorValue = 255),#5
                                    "6"=rgb(26,87,131,maxColorValue = 255),#6
                                    "7"=rgb(154,22,73,maxColorValue = 255),#7
                                    "8"=rgb(249,209,221,maxColorValue = 255),#8
                                    "9"=rgb(113,58,15,maxColorValue = 255),#9
                                    "10"=rgb(45,52,160,maxColorValue = 255),#10
                                    "11"=rgb(121,203,127,maxColorValue = 255),#11
                                    "12"=rgb(173,135,68,maxColorValue = 255),#12
                                    "13"=rgb(2,60,43,maxColorValue = 255),#13
                                    "14"=rgb(129, 147, 209, maxColorValue = 255),#14
                                    "15"=rgb(128, 98,106, maxColorValue = 255),#15
                                    "16"=rgb(64,19,16,maxColorValue = 255),#16
                                    "17"=rgb(109,137,134,maxColorValue = 255))#17
                          )
    
      if(input$include_type) {
        p <- p + facet_wrap(mgi_symbol~Cell.Group,ncol=length(input$cell_types)) + ylim(-5,5)
      } else p <- p + facet_wrap(~mgi_symbol,ncol=1) + ylim(-5,5)
  
  
  if("Boxplot" %in% input$plot_type) p <- p + geom_boxplot(alpha=0.6)
  if("Violin" %in% input$plot_type) p  <- p + geom_violin(alpha=0.6)
  if("Beeswarm" %in% input$plot_type) p  <- p + geom_beeswarm(alpha=0.6,col="black")
  if("Points" %in% input$plot_type) p  <- p + geom_jitter(width=0.1,col="black")
  
  if(input$flip_axis) p <- p + coord_flip()
      p <- p + theme_bw()
      p <- p + theme(axis.title.x = element_text(size=20),
                     axis.text.y = element_text(size=20),
                     axis.text.x = element_text(size=20, face="italic",angle = 90),
                     axis.title = element_text(size=20),
                     strip.text.x = element_text(size=20, face="bold.italic"))
              
  }
      
  p
  
  
   })
   
   
    
    output$umap_plot <- renderImage({
      
      p1 <- prepareUmap()
      width  <- session$clientData$output_gene_boxplot_width
      height <- session$clientData$output_gene_boxplot_height
      
      nrows <- ceiling(length(input$sel_genes)/2)
      
      ### stop an error occurring if nothing is selected
      
      if(is.na(nrows) | nrows == 0) nrows <- 1
      height <- height * nrows
      
      # For high-res displays, this will be greater than 1
      pixelratio <- session$clientData$pixelratio
      
      # A temp file to save the output.
      outfile <- tempfile(fileext='.png')
      
      # Generate the image file
      png(outfile, width=width*pixelratio, height=height*pixelratio,
          res=72*pixelratio)
      
      print(p1)
      dev.off()
      list(src = outfile,
           width = width,
           height = height,
           alt = "This is alternate text")
    }
    )
    
    
  prepareUmap <- reactive({  
    
##   output$umap_plot <- renderPlot({
     
     p <- ggplot()
     
     #  sel_genes <- ifelse(is.null(input$gene_tbl_rows_selected),"ACTB", 
     # dplyr::slice(anno, input$gene_tbl_rows_selected) %>% pull(hgnc_symbol) %>%  as.character)
     if(!is.null(input$sel_genes)){
       sel_genes <- filter(anno, mgi_symbol %in% input$sel_genes) %>% pull(ensembl_gene_id)
     
     
     data <-  data.frame(counts) %>% tibble::rownames_to_column("ensembl_gene_id") %>% 
       filter(ensembl_gene_id %in% sel_genes) %>% 
       tidyr::gather("Cell.ID","Expression",-ensembl_gene_id) %>% 
       left_join(bcode_info) %>% left_join(anno) %>% 
       filter(Cell.ClusterLabel %in% as.numeric(input$selected_clus)) %>% 
       filter(Cell.Group %in% input$cell_types)
     
     if(input$umap_type == "expression"){
     p <- data %>% 
       ggplot(aes(x=UMAP.1,y=UMAP.2,col=Expression)) + geom_point(alpha=0.4) + scale_color_distiller(palette = "RdBu")
     p <- p + facet_wrap(~mgi_symbol)
     }     
     else {
       p <- data %>% 
         ggplot(aes(x=UMAP.1,y=UMAP.2,col=factor(Cell.ClusterLabel))) + geom_point() + scale_color_manual("Cluster",values= c("0"=rgb(3,4,4,maxColorValue = 255),#0
                                                                                                                            "1"=rgb(242,235,15, maxColorValue = 255),#1,
                                                                                                                            "2"=rgb(89,202,233,maxColorValue = 255),#2
                                                                                                                            "3"=rgb(186,80,153,maxColorValue = 255),#3
                                                                                                                            "4"=rgb(253,40,58,maxColorValue = 255),#4
                                                                                                                            "5"=rgb(7,123,55,maxColorValue = 255),#5
                                                                                                                            "6"=rgb(26,87,131,maxColorValue = 255),#6
                                                                                                                            "7"=rgb(154,22,73,maxColorValue = 255),#7
                                                                                                                            "8"=rgb(249,209,221,maxColorValue = 255),#8
                                                                                                                            "9"=rgb(113,58,15,maxColorValue = 255),#9
                                                                                                                            "10"=rgb(45,52,160,maxColorValue = 255),#10
                                                                                                                            "11"=rgb(121,203,127,maxColorValue = 255),#11
                                                                                                                            "12"=rgb(173,135,68,maxColorValue = 255),#12
                                                                                                                            "13"=rgb(2,60,43,maxColorValue = 255),#13
                                                                                                                            "14"=rgb(129, 147, 209, maxColorValue = 255),#14
                                                                                                                            "15"=rgb(128, 98,106, maxColorValue = 255),#15
                                                                                                                            "16"=rgb(64,19,16,maxColorValue = 255),#16
                                                                                                                            "17"=rgb(109,137,134,maxColorValue = 255))) + scale_size_manual(3)#17
                                                                                                                            
     }


     p <- p + theme_bw()
     p <- p + theme(axis.title.x = element_text(size=20),
               axis.text.y = element_text(size=20),
               axis.text.x = element_text(size=20, face="italic"),
               axis.title = element_text(size=20),
               strip.text.x = element_text(size=20, face="bold.italic"),
               legend.text = element_text(size=20))
     }
     
     p
   })
   
  output$dl_umap <- downloadHandler(
    filename = function(){
      paste0(input$profileBasename,".",input$profilePlotFormat)
    },
    content = function(file) {
      
      if(input$profilePlotFormat == "pdf") pdf(file, width=as.numeric(input$profileWidth), height=as.numeric(input$profileHeight))
      else png(file, width=as.numeric(input$profileWidth),height=as.numeric(input$profileHeight))
      
      
      
      p1 <- prepareUmap()
      
      print(p1)
      dev.off()
    }
    
    
  )
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)

