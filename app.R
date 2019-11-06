#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
#library(BiocManager)
#options(repos = BiocManager::repositories())
library(devtools)
#install_github("sheffield-bioinformatics-core/sbcMisc")
library(ggbeeswarm)
library(ggthemes)
library(tidyverse)
library(RColorBrewer)
if(!require(ggthemr)) devtools::install_github('cttobin/ggthemr');library(ggthemr)

uos_pal <- function(){
  uos_pal <- c("Process_Cyan"=rgb(0,159,218,maxColorValue = 255),
               "Pantone_274"=rgb(31,20,93,maxColorValue = 255),
               "Process_Yellow"=rgb(249,227,0,maxColorValue = 255),
               
               "Pantone_347"=rgb(0,155,72,maxColorValue = 255),
               "Pantone_382"=rgb(190,214,0,maxColorValue = 255),
               "Process_Magenta"=rgb(209,0,116,maxColorValue = 255),
               "Pantone_Orange_021"=rgb(255,88,0,maxColorValue = 255),
               
               "Pantone_512"=rgb(119,33,111,maxColorValue = 255),
               "Pantone_485"=rgb(213,43,30,maxColorValue = 255),
               "Pantone_Black"=rgb(30,30,30,maxColorValue = 255),
               "Pantone_161"=rgb(98,60,27,maxColorValue = 255),
               
               "Pantone_7501"=rgb(219,206,172,maxColorValue = 255),
               "Pantone_343"=rgb(3,86,66,maxColorValue = 255),
               "Pantone_322"=rgb(0,116,122,maxColorValue = 255),
               "Pantone_202"=rgb(130,36,51,maxColorValue = 255)
  )
  
  uos_pal
}

uos_colours <- as.character(uos_pal())
# you have to add a colour at the start of your palette for outlining boxes, we'll use a grey:
uos_colours <- c("#555555", uos_colours)
# remove previous effects:

ggthemr_reset()
# Define colours for your figures with define_palette
uos <- define_palette(
  swatch = uos_colours, # colours for plotting points and bars
  gradient = c(lower = uos_colours[1L], upper = uos_colours[2L]), #upper and lower colours for continuous colours
  background = "white" #defining a grey-ish background 
)
# set the theme for your figures:
#ggthemr(uos)

# Define UI for application that draws a histogram
ui <- navbarPage("Explore the single-cell data...",
  
   # Sidebar with a slider input for number of bins 
      tabPanel("Compare Features",
      sidebarLayout(         
      sidebarPanel(
        img(src="logo-sm.png"), br(),a("sbc.shef.ac.uk",href="http://sbc.shef.ac.uk"),
          DT::dataTableOutput("gene_tbl"),
         checkboxInput("exclude_NA","Exclude Missing Values",value=TRUE),
         checkboxInput("flip_axis","Rotate Axis?",value=FALSE),
         checkboxGroupInput(inputId = "plot_type",label="Type of Plot", choices=c("Boxplot","Points","Violin","Beeswarm"),selected="Boxplot"),
        selectInput(inputId = "selected_clus",label="Select Clusters to Compare",0:17, multiple = TRUE,selectize = TRUE,selected=c(0,1,3:5,7:15,17)),
        checkboxInput("include_type","Include Cell Type", value=TRUE)
        # selectInput("boxplotTheme", "Pick a plot style",choices=c("bw","grey","classic","minimal","light","Wall Street Journal","Economist", "Excel","solarized","stata","calc","dark","fivethirtyeight","tufte"),selected="bw"),
        # selectInput("pal_name", "Choose a colour palette", choices = c("UOS","Set1","Set2","Set3","Pastel1","Pastel2","Paired","Dark2","Accent"),selected="UOS")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("gene_boxplot"),
         plotOutput("umap_plot")
      )
      )
      )
   )

# Define server logic required to draw a histogram
server <- function(input, output) {
  library(tidyverse)   

  counts <- readRDS("counts.rds")
  bcode_info <- readRDS("barcodes.rds")
  anno <- readRDS("anno.rds")

  

  output$gene_tbl <- DT::renderDataTable(anno,options=list(pageLength=5))
  
   output$gene_boxplot <- renderPlot({


#  sel_genes <- ifelse(is.null(input$gene_tbl_rows_selected),"ACTB", 
 # dplyr::slice(anno, input$gene_tbl_rows_selected) %>% pull(hgnc_symbol) %>%  as.character)
  sel_genes <- "ENSMUSG00000000001"
  
  if(!is.null(input$gene_tbl_rows_selected)) sel_genes <- dplyr::slice(anno, input$gene_tbl_rows_selected) %>% pull(ensembl_gene_id)
     

  
  data <-  data.frame(counts) %>% tibble::rownames_to_column("ensembl_gene_id") %>% 
    filter(ensembl_gene_id %in% sel_genes) %>% 
    tidyr::gather("FeatureID","Expression",-ensembl_gene_id) %>% 
    left_join(bcode_info) %>% left_join(anno) %>% 
    filter(VAE_ClusterLabel %in% as.numeric(input$selected_clus))

      p <- data %>% 
      ggplot(aes(x=VAE_ClusterLabel,y=Expression,fill=factor(VAE_ClusterLabel),group=VAE_ClusterLabel)) + 
        scale_fill_manual(values= c("0"=rgb(3,4,4,maxColorValue = 255),#0
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
        p <- p + facet_wrap(mgi_symbol~MouseGenotype,ncol=4) + ylim(-5,5)
      } else p <- p + facet_wrap(~mgi_symbol) + ylim(-5,5)
  
  
  if("Boxplot" %in% input$plot_type) p <- p + geom_boxplot(alpha=0.6)
  if("Violin" %in% input$plot_type) p  <- p + geom_violin(alpha=0.6)
  if("Beeswarm" %in% input$plot_type) p  <- p + geom_beeswarm(alpha=0.6,col="black")
  if("Points" %in% input$plot_type) p  <- p + geom_jitter(width=0.1,col="black")
  
  if(input$flip_axis) p <- p + coord_flip()
              

      
  p
  
  
   })
   
   
   output$umap_plot <- renderPlot({
     
     
     #  sel_genes <- ifelse(is.null(input$gene_tbl_rows_selected),"ACTB", 
     # dplyr::slice(anno, input$gene_tbl_rows_selected) %>% pull(hgnc_symbol) %>%  as.character)
     sel_genes <- "ENSMUSG00000000001"
     
     if(!is.null(input$gene_tbl_rows_selected)) sel_genes <- dplyr::slice(anno, input$gene_tbl_rows_selected) %>% pull(ensembl_gene_id)
     
     
     
     data <-  data.frame(counts) %>% tibble::rownames_to_column("ensembl_gene_id") %>% 
       filter(ensembl_gene_id %in% sel_genes) %>% 
       tidyr::gather("FeatureID","Expression",-ensembl_gene_id) %>% 
       left_join(bcode_info) %>% left_join(anno)
     
     
     p <- data %>% 
       ggplot(aes(x=UMAP_1,y=UMAP_2,col=Expression)) + geom_point(alpha=0.4) + scale_color_distiller(palette = "RdBu")
     
     p <- p + facet_wrap(~mgi_symbol)
     p
   })
   
}

# Run the application 
shinyApp(ui = ui, server = server)

