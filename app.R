#This code is for gene-level analysis (LTG preferred)
library(tidyr)
library(ggplot2)
library(shiny)
library(pheatmap)
library(factoextra)
library(FactoMineR)
library(EnhancedVolcano)
library(stringr)
library(tibble)
library(dplyr)

data_url <- "https://raw.githubusercontent.com/arhaeusler/LTG_shiny/main/LTG_TEA_KCL_datafile.csv"
anno_url <- "https://raw.githubusercontent.com/arhaeusler/LTG_shiny/main/LTG_TEA_KCL_annotation.csv"

# Set file size (for 50MB)
options(shiny.maxRequestSize = 85*1024^2)

# Define UI
ui <- fluidPage(
  titlePanel("i3 Neuronal Activity-Dependent Transcriptomic Landscape"),
  sidebarLayout(
    sidebarPanel(
      uiOutput("treatmentSelect"),
      uiOutput("timeSelect"),
      uiOutput("groupSelect"),
      checkboxInput("selectAllGenes", "Select All Genes", FALSE),
      textInput("genes", "Enter Gene Symbol (comma-separated)", ""),
      actionButton("update", "Update"),
      br(),
      br(),
      #downloadButton("downloadPlots", "Downlaod All Plots"),   ## UNDER CONSTRUCTION ##
      br(),br(),
      actionButton("githubButton", label = "Go to Github Repository")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Heatmaps",
                 fluidRow(
                   column(12, plotOutput("heatmapPlot_TEA")),
                   column(12, plotOutput("heatmapPlot_KCl"))
                 )
        ),
        tabPanel("TEA Volcano",
                 fluidRow(
                   column(10, plotOutput("volcanoPlot_C9_UT_TEA")),
                   column(10, plotOutput("volcanoPlot_C9_TTX_TEA"))
                 )
        ),
        tabPanel("KCL Volcano",
                 fluidRow(
                   column(10, plotOutput("volcanoPlot_KCL_C902")),
                   column(10, plotOutput("volcanoPlot_KCL_C906")),
                   column(10, plotOutput("volcanoPlot_KCL_WT02")),
                   column(10, plotOutput("volcanoPlot_KCL_WT06"))
                 )
        ),
        tabPanel("Expression Plots",
                 fluidRow(
                   column(6, plotOutput("expressionPlot_UT")),
                   column(6, plotOutput("expressionPlot_TTX")),
                   column(6, plotOutput("expressionPlot_TEA"))
                 )
        ),
        tabPanel("KCL Expression over Time",
                 fluidRow(
                   column(12, plotOutput("expressionPlot_KCL"))
                 )
        ),
        tabPanel("TEA Pseudotime",
                 fluidRow(
                   column(12, plotOutput("pseudotime_TEA"))
                 )
        ),
        tabPanel("PCA",
                 fluidRow(
                   column(6, plotOutput("pcaPlot_TTX")),
                   column(6, plotOutput("pcaPlot_UT")),
                   column(6, plotOutput("pcaPlot_TEA"))
                 )
        )
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  observeEvent(input$githubButton, {
    browseURL("https://github.com/arhaeusler/LTG_shiny")
  })
  ########## UNDER CONSTRUCTION #####
  # output$downloadPlots <- downloadHandler(
  #   filename = function() {
  #     paste("all_plots", Sys.Date(), ".pdf", sep = "")
  #   },
  #   content = function(file) {
  #     # Open a PDF device to save the plots
  #     pdf(file, width = 8, height = 6)
  #     # Create and save each plot
  #     print(heatmapPlot_TEA)
  #     print(heatmapPlot_KCl)
  #     # Close the PDF device
  #     dev.off()
  #   }
  # )
  ####################################
  data <- reactive({
    read.csv(data_url)
  })
  
  annotation_data <- reactive({
    read.csv(anno_url)
  })
  
  observe({
    df <- data()
    updateSelectInput(session, "variable", choices = colnames(df))
  })
  output$variableSelect <- renderUI({
    df <- data()
    selectInput("variable", "Select Variable to Compare",
                choices = colnames(df))
  })
  
  gene_list <- reactive({
    if (input$selectAllGenes) {
      unique(data()$Gene)              
    } else {
      req(input$genes)
      strsplit(input$genes, ",")[[1]]
    }
  })
  
  unfiltered_data_tea <- reactive({  #for pca plots
    df <- data()
    unfiltered_df <- df %>%
      dplyr::select(-Length,-GeneID,-(starts_with("KCL"))) %>%
      group_by(Gene) %>%
      summarise(across(everything(), median, na.rm = TRUE))
    unfiltered_df
  })
  
  filtered_data <- reactive({
    genes <- gene_list()
    df <- data()
    #new (get median if there are multiple matches to a gene symbol)
    filtered_df <- df %>%
      dplyr::select(-Length, -GeneID, -(starts_with("KCL"))) %>%
      filter(Gene %in% genes) %>%
      group_by(Gene) %>%
      summarize(across(everything(), median, na.rm = TRUE))
    filtered_df
  })
  
  
  kcl_filtered <- reactive({  #For KCL expression 
    genes <- gene_list()
    df <- data()
    kcl_filtered_df <- df %>%
      dplyr::select(Gene, starts_with("KCL")) %>%
      dplyr::filter(Gene %in% genes) %>%
      mutate(across(ends_with("rlog") & where(is.character), as.numeric)) %>%
      group_by(Gene) %>%
      summarize(across(everything(), median, na.rm = TRUE))
    kcl_filtered_df
  })
  
  
  #start volcano plots
  #C9: UT vs TEA
  output$volcanoPlot_C9_UT_TEA <- renderPlot({
    df <- filtered_data()
    #top_genes <- df[order(df$pvalue), "Gene"]
    #top_genes <- ifelse(length(top_genes) > 20, top_genes[1:20], top_genes)
    
    df_fc <- df %>%
      dplyr::select(Gene,ends_with("_log2FoldChange")) %>%
      pivot_longer(cols= ends_with("_log2FoldChange"), names_to = "sample", values_to = "log2FoldChange")
    df_fc$sample = gsub("_log2FoldChange$", "", df_fc$sample)
    df_padj <- df %>%
      dplyr::select(Gene,ends_with("_padj")) %>%
      pivot_longer(cols= ends_with("_padj"), names_to = "sample", values_to = "pvalue")
    df_padj$sample = gsub("_padj$", "", df_padj$sample)
    
    df_volcano <- left_join(df_fc,df_padj, by = c("sample", "Gene"))
    df_volcano[is.na(df_volcano)] = 0
    
    df_C9_UT_TEA <- df_volcano %>%
      filter(str_detect(sample, "C9_UT_TEA"))
    
    EnhancedVolcano(
      df_C9_UT_TEA,
      lab = df_C9_UT_TEA$Gene,
      x = 'log2FoldChange',
      y = 'pvalue',
      selectLab = df_C9_UT_TEA$Gene,
      xlab = bquote(~Log[2]~ 'fold change'),
      pCutoff = 0.05,
      FCcutoff = 2.0,
      pointSize = 4.0,
      labSize = 4,
      labCol = 'black',
      labFace = 'bold',
      boxedLabels = TRUE,
      colAlpha = 4/5,
      legendPosition = 'right',
      legendLabSize = 14,
      legendIconSize = 4.0,
      drawConnectors = TRUE,
      widthConnectors = 1.0,
      colConnectors = 'black') +
      ggtitle("C9: UT vs TEA") +
      theme(plot.title = element_text(hjust = 0.5))
  })
  #end C9:UT vs TEA
  
  #C9:TTX vs TEA
  output$volcanoPlot_C9_TTX_TEA <- renderPlot({
    df <- filtered_data()
    #top_genes <- df[order(df$pvalue), "Gene"]
    #top_genes <- ifelse(length(top_genes) > 20, top_genes[1:20], top_genes)
    
    df_fc <- df %>%
      dplyr::select(Gene,ends_with("_log2FoldChange")) %>%
      pivot_longer(cols= ends_with("_log2FoldChange"), names_to = "sample", values_to = "log2FoldChange")
    df_fc$sample = gsub("_log2FoldChange$", "", df_fc$sample)
    df_padj <- df %>%
      dplyr::select(Gene,ends_with("_padj")) %>%
      pivot_longer(cols= ends_with("_padj"), names_to = "sample", values_to = "pvalue")
    df_padj$sample = gsub("_padj$", "", df_padj$sample)
    
    df_volcano <- left_join(df_fc,df_padj, by = c("sample", "Gene"))
    df_volcano[is.na(df_volcano)] = 0
    
    df_C9_TTX_TEA <- df_volcano %>%
      filter(str_detect(sample, "C9_TTX_TEA"))
    
    EnhancedVolcano(
      df_C9_TTX_TEA,
      lab = df_C9_TTX_TEA$Gene,
      x = 'log2FoldChange',
      y = 'pvalue',
      selectLab = df_C9_TTX_TEA$Gene,
      xlab = bquote(~Log[2]~ 'fold change'),
      pCutoff = 0.05,
      FCcutoff = 2.0,
      pointSize = 4.0,
      labSize = 4,
      labCol = 'black',
      labFace = 'bold',
      boxedLabels = TRUE,
      colAlpha = 4/5,
      legendPosition = 'right',
      legendLabSize = 14,
      legendIconSize = 4.0,
      drawConnectors = TRUE,
      widthConnectors = 1.0,
      colConnectors = 'black') +
      ggtitle("C9: TTX vs TEA") +
      theme(plot.title = element_text(hjust = 0.5))
  })
  # end TEA volcano plots
  
  ## KCl volcano plots
  output$volcanoPlot_KCL_C902 <- renderPlot({
    df <- kcl_filtered()
    #top_genes <- df[order(df$pvalue), "Gene"]
    #top_genes <- ifelse(length(top_genes) > 20, top_genes[1:20], top_genes)
    
    df_fc <- df %>%
      dplyr::select(Gene,ends_with("_log2FoldChange")) %>%
      pivot_longer(cols= ends_with("_log2FoldChange"), names_to = "sample", values_to = "log2FoldChange")
    df_fc$sample = gsub("_log2FoldChange$", "", df_fc$sample)
    df_padj <- df %>%
      dplyr::select(Gene,ends_with("_padj")) %>%
      pivot_longer(cols= ends_with("_padj"), names_to = "sample", values_to = "pvalue")
    df_padj$sample = gsub("_padj$", "", df_padj$sample)
    
    df_volcano <- left_join(df_fc,df_padj, by = c("sample", "Gene"))
    df_volcano[is.na(df_volcano)] = 0
    
    df_KCL_C902 <- df_volcano %>%
      filter(str_detect(sample, "KCL_C90_C92"))
    
    EnhancedVolcano(
      df_KCL_C902,
      lab = df_KCL_C902$Gene,
      x = 'log2FoldChange',
      y = 'pvalue',
      selectLab = df_KCL_C902$Gene,
      xlab = bquote(~Log[2]~ 'fold change'),
      pCutoff = 0.05,
      FCcutoff = 2.0,
      pointSize = 4.0,
      labSize = 4,
      labCol = 'black',
      labFace = 'bold',
      boxedLabels = TRUE,
      colAlpha = 4/5,
      legendPosition = 'right',
      legendLabSize = 14,
      legendIconSize = 4.0,
      drawConnectors = TRUE,
      widthConnectors = 1.0,
      colConnectors = 'black') +
      ggtitle("KCl: C90 vs C92") +
      theme(plot.title = element_text(hjust = 0.5))
  })
  
  output$volcanoPlot_KCL_C906 <- renderPlot({
    df <- kcl_filtered()
    #top_genes <- df[order(df$pvalue), "Gene"]
    #top_genes <- ifelse(length(top_genes) > 20, top_genes[1:20], top_genes)
    
    df_fc <- df %>%
      dplyr::select(Gene,ends_with("_log2FoldChange")) %>%
      pivot_longer(cols= ends_with("_log2FoldChange"), names_to = "sample", values_to = "log2FoldChange")
    df_fc$sample = gsub("_log2FoldChange$", "", df_fc$sample)
    df_padj <- df %>%
      dplyr::select(Gene,ends_with("_padj")) %>%
      pivot_longer(cols= ends_with("_padj"), names_to = "sample", values_to = "pvalue")
    df_padj$sample = gsub("_padj$", "", df_padj$sample)
    
    df_volcano <- left_join(df_fc,df_padj, by = c("sample", "Gene"))
    df_volcano[is.na(df_volcano)] = 0
    
    df_KCL_C906 <- df_volcano %>%
      filter(str_detect(sample, "KCL_C90_C96"))
    
    EnhancedVolcano(
      df_KCL_C906,
      lab = df_KCL_C906$Gene,
      x = 'log2FoldChange',
      y = 'pvalue',
      selectLab = df_KCL_C906$Gene,
      xlab = bquote(~Log[2]~ 'fold change'),
      pCutoff = 0.05,
      FCcutoff = 2.0,
      pointSize = 4.0,
      labSize = 4,
      labCol = 'black',
      labFace = 'bold',
      boxedLabels = TRUE,
      colAlpha = 4/5,
      legendPosition = 'right',
      legendLabSize = 14,
      legendIconSize = 4.0,
      drawConnectors = TRUE,
      widthConnectors = 1.0,
      colConnectors = 'black') +
      ggtitle("KCl: C90 vs C96") +
      theme(plot.title = element_text(hjust = 0.5))
  })
  
  output$volcanoPlot_KCL_WT02 <- renderPlot({
    df <- kcl_filtered()
    #top_genes <- df[order(df$pvalue), "Gene"]
    #top_genes <- ifelse(length(top_genes) > 20, top_genes[1:20], top_genes)
    
    df_fc <- df %>%
      dplyr::select(Gene,ends_with("_log2FoldChange")) %>%
      pivot_longer(cols= ends_with("_log2FoldChange"), names_to = "sample", values_to = "log2FoldChange")
    df_fc$sample = gsub("_log2FoldChange$", "", df_fc$sample)
    df_padj <- df %>%
      dplyr::select(Gene,ends_with("_padj")) %>%
      pivot_longer(cols= ends_with("_padj"), names_to = "sample", values_to = "pvalue")
    df_padj$sample = gsub("_padj$", "", df_padj$sample)
    
    df_volcano <- left_join(df_fc,df_padj, by = c("sample", "Gene"))
    df_volcano[is.na(df_volcano)] = 0
    
    df_KCL_WT02 <- df_volcano %>%
      filter(str_detect(sample, "KCL_WT0_WT2"))
    
    EnhancedVolcano(
      df_KCL_WT02,
      lab = df_KCL_WT02$Gene,
      x = 'log2FoldChange',
      y = 'pvalue',
      selectLab = df_KCL_WT02$Gene,
      xlab = bquote(~Log[2]~ 'fold change'),
      pCutoff = 0.05,
      FCcutoff = 2.0,
      pointSize = 4.0,
      labSize = 4,
      labCol = 'black',
      labFace = 'bold',
      boxedLabels = TRUE,
      colAlpha = 4/5,
      legendPosition = 'right',
      legendLabSize = 14,
      legendIconSize = 4.0,
      drawConnectors = TRUE,
      widthConnectors = 1.0,
      colConnectors = 'black') +
      ggtitle("KCl: WT0 vs WT2") +
      theme(plot.title = element_text(hjust = 0.5))
  })
  
  output$volcanoPlot_KCL_WT06 <- renderPlot({
    df <- kcl_filtered()
    #top_genes <- df[order(df$pvalue), "Gene"]
    #top_genes <- ifelse(length(top_genes) > 20, top_genes[1:20], top_genes)
    
    df_fc <- df %>%
      dplyr::select(Gene,ends_with("_log2FoldChange")) %>%
      pivot_longer(cols= ends_with("_log2FoldChange"), names_to = "sample", values_to = "log2FoldChange")
    df_fc$sample = gsub("_log2FoldChange$", "", df_fc$sample)
    df_padj <- df %>%
      dplyr::select(Gene,ends_with("_padj")) %>%
      pivot_longer(cols= ends_with("_padj"), names_to = "sample", values_to = "pvalue")
    df_padj$sample = gsub("_padj$", "", df_padj$sample)
    
    df_volcano <- left_join(df_fc,df_padj, by = c("sample", "Gene"))
    df_volcano[is.na(df_volcano)] = 0
    
    df_KCL_WT06 <- df_volcano %>%
      filter(str_detect(sample, "KCL_WT0_WT6"))
    
    EnhancedVolcano(
      df_KCL_WT06,
      lab = df_KCL_WT06$Gene,
      x = 'log2FoldChange',
      y = 'pvalue',
      selectLab = df_KCL_WT06$Gene,
      xlab = bquote(~Log[2]~ 'fold change'),
      pCutoff = 0.05,
      FCcutoff = 2.0,
      pointSize = 4.0,
      labSize = 4,
      labCol = 'black',
      labFace = 'bold',
      boxedLabels = TRUE,
      colAlpha = 4/5,
      legendPosition = 'right',
      legendLabSize = 14,
      legendIconSize = 4.0,
      drawConnectors = TRUE,
      widthConnectors = 1.0,
      colConnectors = 'black') +
      ggtitle("KCl: WT0 vs WT6") +
      theme(plot.title = element_text(hjust = 0.5))
  })
  ## end KCL volcano plots
  
  ##heatmaps
  output$heatmapPlot_TEA <- renderPlot({
    df <- filtered_data()
    df_heatmap <- df # Adjusting df for heatmap of rlog values
    df_heatmap = df_heatmap[,grepl("_rlog",colnames(df_heatmap))] #<--- THIS IS THE LINE I SHOULD HAVE HAD PREVIOUSLY
    colnames(df_heatmap) = gsub("_rlog$", "", colnames(df_heatmap))
    rownames(df_heatmap) = df$Gene
    pheatmap(as.matrix(df_heatmap), main="TEA Heatmap (rlog)")
  })
  
  output$heatmapPlot_KCl <- renderPlot({
    df <- kcl_filtered()
    df_heatmap <- df # Adjusting df for heatmap of rlog values
    df_heatmap = df_heatmap[,grepl("_rlog",colnames(df_heatmap))] 
    colnames(df_heatmap) = gsub("_rlog$", "", colnames(df_heatmap))
    rownames(df_heatmap) = df$Gene
    pheatmap(as.matrix(df_heatmap), main="KCl Heatmap (rlog)")
  })
  
  
  ## start expression plots ##NOTE: use boxplot w/points, not barplots. label C9 w/shades of purple. WT w/shades of grey
  output$expressionPlot_UT <- renderPlot({
    df <- filtered_data()
    df_exp <- df %>% 
      dplyr::select(Gene, ends_with("_rlog")) %>%
      pivot_longer(cols = ends_with("_rlog"), names_to = "sample", values_to = "expression")
    df_exp$sample = gsub("_rlog$", "", df_exp$sample)
    
    anno_exp <- annotation_data()
    anno_exp <- anno_exp %>%
      dplyr::select(SampleID, Group) %>%
      dplyr::rename(sample = SampleID)
    df_exp <- df_exp %>%
      left_join(., anno_exp, by = 'sample')
    df_exp_ut <- df_exp %>%
      filter(str_detect(sample, "UT"))
    
    group_colors <- c("C9_UT" = "#9370DB",  # Medium purple
                      "WT_UT" = "#A9A9A9")  # Dark grey
    
    point_colors <- c("C9_UT" = "#6A0DAD",  # Dark purple
                      "WT_UT" = "#2F4F4F")  # Dark slate grey
    
    ggplot(df_exp_ut, aes(x=Gene, y=expression, fill=Group)) +
      geom_boxplot(outlier.shape = NA, position = position_dodge(1)) +  # Boxplot without displaying outliers
      geom_point(aes(color = Group),position = position_jitterdodge(), alpha=0.7) +
      theme_minimal() +
      ggtitle("Expression Plot UT") +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_fill_manual(values = group_colors) +  # Use custom colors for boxes
      scale_color_manual(values = point_colors)   # Use custom colors for points
  })
  
  output$expressionPlot_TTX <- renderPlot({
    df <- filtered_data()
    df_exp <- df %>% # Adjusting df for expression plot.  
      dplyr::select(Gene, ends_with("_rlog")) %>%
      pivot_longer(cols = ends_with("_rlog"), names_to = "sample", values_to = "expression")
    df_exp$sample = gsub("_rlog$", "", df_exp$sample)
    
    anno_exp <- annotation_data()
    anno_exp <- anno_exp %>%
      dplyr::select(SampleID, Group) %>%
      dplyr::rename(sample = SampleID)
    df_exp <- df_exp %>%
      left_join(., anno_exp, by = 'sample')
    df_exp_ttx <- df_exp %>%
      filter(str_detect(sample, "TTX"))
    
    group_colors <- c("C9_TTX" = "#9370DB",  # Medium purple
                      "WT_TTX" = "#A9A9A9")  # Dark grey
    
    point_colors <- c("C9_TTX" = "#6A0DAD",  # Dark purple
                      "WT_TTX" = "#2F4F4F")  # Dark slate grey
    
    ggplot(df_exp_ttx, aes(x=Gene, y=expression, fill=Group)) +
      geom_boxplot(outlier.shape = NA, position = position_dodge(1)) +  # Boxplot without displaying outliers
      geom_point(aes(color = Group),position = position_jitterdodge(), alpha=0.7) +
      theme_minimal() +
      ggtitle("Expression Plot TTX") +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_fill_manual(values = group_colors) +  # Use custom colors for boxes
      scale_color_manual(values = point_colors)   # Use custom colors for points
  })
  output$expressionPlot_TEA <- renderPlot({
    df <- filtered_data()
    df_exp <- df %>% 
      dplyr::select(Gene, ends_with("_rlog")) %>%
      pivot_longer(cols = ends_with("_rlog"), names_to = "sample", values_to = "expression")
    df_exp$sample = gsub("_rlog$", "", df_exp$sample)
    
    anno_exp <- annotation_data()
    anno_exp <- anno_exp %>%
      dplyr::select(SampleID, Group) %>%
      dplyr::rename(sample = SampleID)
    df_exp <- df_exp %>%
      left_join(., anno_exp, by = 'sample')
    df_exp_tea <- df_exp %>%
      filter(str_detect(sample, "TEA"))
    
    group_colors <- c("C9_TEA" = "#9370DB",  # Medium purple
                      "WT_TEA" = "#A9A9A9")  # Dark grey
    
    point_colors <- c("C9_TEA" = "#6A0DAD",  # Dark purple
                      "WT_TEA" = "#2F4F4F")  # Dark slate grey
    
    ggplot(df_exp_tea, aes(x=Gene, y=expression, fill=Group)) +
      geom_boxplot(outlier.shape = NA, position = position_dodge(1)) +  # Boxplot without displaying outliers
      geom_point(aes(color = Group),position = position_jitterdodge(), alpha=0.7) +
      theme_minimal() +
      ggtitle("Expression Plot TEA") +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_fill_manual(values = group_colors) +  # Use custom colors for boxes
      scale_color_manual(values = point_colors)   # Use custom colors for points
  })
  ## end expression plots
  
  ## KCL Expression over time
  #use rlog, not masigpro data
  
  output$expressionPlot_KCL <- renderPlot({
    df <- kcl_filtered()
    df_kcl_plot <- df %>%
      dplyr::select(Gene, ends_with("rlog")) %>%
      pivot_longer(cols = ends_with("_rlog"), names_to = "sample", values_to = "expression")
    df_kcl_plot$sample = gsub("_rlog$", "", df_kcl_plot$sample)
    df_kcl_plot$sample = gsub("KCL_", "", df_kcl_plot$sample)
    
    anno_kcl <- annotation_data()
    anno_kcl <- anno_kcl %>%
      dplyr::filter(Time != "NA") %>%
      dplyr::rename(sample = SampleID)
    anno_kcl$sample = gsub("KCL_", "", anno_kcl$sample)
    df_kcl_plot <- df_kcl_plot %>%
      left_join(., anno_kcl, by = 'sample') %>%
      dplyr::select(-Group,-Condition) %>%
      group_by(Gene, Time, Line) %>% 
      summarize(across(expression, median, na.rm = TRUE))
    
    
    #boxplot with connectors. rlog expression over time. C9/WT with subcolors matching other plots.
    custom_colors <- c("C9_Line1" = "#6A0DAD",  # Dark purple
                       "C9_Line2" = "#9370DB",  # Medium purple
                       "C9_Line3" = "#b27faa",  # Light purple
                       "WT_Line1" = "#2F4F4F",  # Dark slate grey
                       "WT_Line2" = "#A9A9A9",  # Dark grey
                       "WT_Line3" = "#D3D3D3")  # Light grey
    
    #boxplot with connectors. color by Line. x = Time, y = 'expression', and fill = Line.
    expressionPlot_KCL <- ggplot(df_kcl_plot, aes(x=as.factor(Time), y=expression, color=Line, group=Line)) + 
      geom_line(size = 1) + 
      geom_point(size = 3) + 
      scale_color_manual(values = custom_colors) +
      facet_wrap(~Gene) +
      ggtitle("KCl expression over time")
    expressionPlot_KCL
  })
  
  ## TEA pseudotime plot
  output$pseudotime_TEA <- renderPlot({
    df <- filtered_data()
    pseudotime_df <- df %>%
      dplyr::select(Gene, ends_with("rlog")) %>%
      pivot_longer(cols = ends_with("_rlog"), names_to = "sample", values_to = "expression")
    pseudotime_df$sample = gsub("_rlog$", "", pseudotime_df$sample)
    
    anno_tea <- annotation_data()
    anno_tea <- anno_tea %>%
      dplyr::rename(sample = SampleID)
    pseudotime_df <- pseudotime_df %>%
      left_join(., anno_tea, by = 'sample')
    
    ### Make a pseudotime (TEA -> UT -> TTX)
    pseudotime_df$Group = gsub("C9_", "", pseudotime_df$Group)
    pseudotime_df$Group = gsub("WT_", "", pseudotime_df$Group)
    pseudotime_df$pseudotime <- factor(pseudotime_df$Group, levels = c("TEA","UT","TTX"))
    pseudotime_df <- pseudotime_df %>%
      dplyr::select(-Condition,-Group,-Time) %>%
      group_by(Gene, pseudotime, Line) %>% #need to change to numeric? 
      summarize(across(expression, median, na.rm = TRUE))
    
    #boxplot with connectors. rlog expression over time. C9/WT with subcolors matching other plots.
    custom_colors <- c("C9_Line1" = "#6A0DAD",  # Dark purple
                       "C9_Line2" = "#9370DB",  # Medium purple
                       "C9_Line3" = "#b27faa",  # Light purple
                       "WT_Line1" = "#2F4F4F",  # Dark slate grey
                       "WT_Line2" = "#A9A9A9",  # Dark grey
                       "WT_Line3" = "#D3D3D3")  # Light grey
    
    #boxplot with connectors. color by Line. x = Time, y = 'expression', and fill = Line.
    ggplot(pseudotime_df, aes(x=pseudotime, y=expression, color=Line, group = Line)) + #,group=Line
      geom_line(size = 1) + 
      geom_point(size = 3) + 
      scale_color_manual(values = custom_colors) +
      facet_wrap(~Gene) +
      ggtitle("TEA expression over pseudotime")
    
  })
  ## End KCl expression over time
  
  ## PCA   
  #TTX
  output$pcaPlot_TTX <- renderPlot({ 
    df <- unfiltered_data_tea()
    df_pca_ttx <- df %>%
      dplyr::select(Gene, ends_with("TTX_rlog"))
    colnames(df_pca_ttx) = gsub("_rlog$", "", colnames(df_pca_ttx))
    df_pca_ttx2 <- df_pca_ttx %>%
      dplyr::filter(Gene != "NA") %>%
      column_to_rownames("Gene") %>%
      t() %>%
      as.data.frame()
    
    anno <- annotation_data()
    pca_anno_ttx <- anno %>%
      dplyr::filter(str_detect(SampleID, "TTX"))
    
    group_ID <- as.factor(pca_anno_ttx$Line)
    group_ID2 <- as.factor(pca_anno_ttx$Condition)
    #colors purple for C9 and grey for WT
    custom_colors <- c("C9_Line1" = "#6A0DAD",  # Dark purple
                       "C9_Line2" = "#9370DB",  # Medium purple
                       "C9_Line3" = "#b27faa", #light purple
                       "WT_Line1" = "#2F4F4F",  # Dark slate grey
                       "WT_Line2" = "#A9A9A9",  # Dark grey
                       "WT_Line3" = "#D3D3D3")  # Light grey
    
    combined_colors <- c("C9" = "#9370DB",  # Medium purple
                         "WT" = "#A9A9A9")  # Dark grey
    
    res.pca <- PCA(df_pca_ttx2, graph = FALSE) 
    p <- fviz_pca_ind(res.pca,
                      geom.ind = "point",
                      fill.ind = group_ID, col.ind = "black",
                      mean.point = FALSE,
                      pointshape = 21, pointsize = "cos2",
                      palette = custom_colors,
                      addEllipses = TRUE,
                      col.var = "contrib",
                      repel = TRUE,
                      select.var = list(cos2 = 10), legend.title = list(fill = "Group", color="Contrib"),
                      habillage = group_ID2,
                      palette.ellipse = combined_colors) +
      ggtitle("TTX silenced") +
      theme(plot.title = element_text(hjust = 0.5))
    
    p2 <- p + scale_color_manual(values = combined_colors, name = "Condition")
    p2
    
  })
  
  #TEA
  output$pcaPlot_TEA <- renderPlot({ 
    df <- unfiltered_data_tea()
    df_pca_tea <- df %>%
      dplyr::select(Gene, ends_with("TEA_rlog"))
    colnames(df_pca_tea) = gsub("_rlog$", "", colnames(df_pca_tea))
    df_pca_tea2 <- df_pca_tea %>%
      dplyr::filter(Gene != "NA") %>%
      column_to_rownames("Gene") %>%
      t() %>%
      as.data.frame()
    
    anno <- annotation_data()
    pca_anno_tea <- anno %>%
      dplyr::filter(str_detect(SampleID, "TEA"))
    
    custom_colors <- c("C9_Line1" = "#6A0DAD",  # Dark purple
                       "C9_Line2" = "#9370DB",  # Medium purple
                       "C9_Line3" = "#b27faa", #light purple
                       "WT_Line1" = "#2F4F4F",  # Dark slate grey
                       "WT_Line2" = "#A9A9A9",  # Dark grey
                       "WT_Line3" = "#D3D3D3")  # Light grey
    
    combined_colors <- c("C9" = "#9370DB",  # Medium purple
                         "WT" = "#A9A9A9")  # Dark grey
    
    group_ID <- as.factor(pca_anno_tea$Line)
    group_ID2 <- as.factor(pca_anno_tea$Condition)
    res.pca <- PCA(df_pca_tea2, graph = FALSE) 
    p <- fviz_pca_ind(res.pca,
                      geom.ind = "point",
                      fill.ind = group_ID, col.ind = "black",
                      mean.point = FALSE,
                      pointshape = 21, pointsize = "cos2",
                      palette = custom_colors,
                      addEllipses = TRUE,
                      col.var = "contrib",
                      repel = TRUE,
                      select.var = list(cos2 = 10), legend.title = list(fill = "Group", color="Contrib"),
                      habillage = group_ID2,
                      palette.ellipse = combined_colors) +
      ggtitle("TEA activation") +
      theme(plot.title = element_text(hjust = 0.5))
    
    p2 <- p + scale_color_manual(values = combined_colors, name = "Condition")
    p2
  })
  
  
  #UT
  output$pcaPlot_UT <- renderPlot({
    df <- unfiltered_data_tea()
    df_pca_ut <- df %>%
      dplyr::select(Gene, ends_with("UT_rlog"))
    colnames(df_pca_ut) = gsub("_rlog$", "", colnames(df_pca_ut))
    df_pca_ut2 <- df_pca_ut %>%
      dplyr::filter(Gene != "NA") %>%
      column_to_rownames("Gene") %>%
      t() %>%
      as.data.frame()
    
    anno <- annotation_data()
    pca_anno_ut <- anno %>%
      dplyr::filter(str_detect(SampleID, "UT"))
    
    custom_colors <- c("C9_Line1" = "#6A0DAD",  # Dark purple
                       "C9_Line2" = "#9370DB",  # Medium purple
                       "C9_Line3" = "#b27faa", #light purple
                       "WT_Line1" = "#2F4F4F",  # Dark slate grey
                       "WT_Line2" = "#A9A9A9",  # Dark grey
                       "WT_Line3" = "#D3D3D3")  # Light grey
    
    combined_colors <- c("C9" = "#9370DB",  # Medium purple
                         "WT" = "#A9A9A9")  # Dark grey
    
    group_ID <- as.factor(pca_anno_ut$Line)
    group_ID2 <- as.factor(pca_anno_ut$Condition)
    res.pca <- PCA(df_pca_ut2, graph = FALSE) 
    
    p <- fviz_pca_ind(res.pca,
                      geom.ind = "point",
                      fill.ind = group_ID, col.ind = "black",
                      mean.point = FALSE,
                      pointshape = 21, pointsize = "cos2",
                      palette = custom_colors,
                      addEllipses = TRUE,
                      col.var = "contrib",
                      repel = TRUE,
                      select.var = list(cos2 = 10), legend.title = list(fill = "Group", color="Contrib"),
                      habillage = group_ID2,
                      palette.ellipse = combined_colors) +
      ggtitle("Spontaneous firing") +
      theme(plot.title = element_text(hjust = 0.5))
    
    p2 <- p + scale_color_manual(values = combined_colors, name = "Condition")
    p2
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
