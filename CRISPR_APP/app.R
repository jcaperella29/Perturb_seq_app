library(shiny)
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)
library(reshape2)
library(mixtools)
library(SeuratData)
library(SeuratWrappers)
library(presto)
library(plotly)
library(DT)
library(markdown)

options(shiny.maxRequestSize = 1024*1024^2) # 1 G

# Helper: Custom plotting theme
custom_theme <- theme(
  plot.title = element_text(size=16, hjust = 0.5), 
  legend.key.size = unit(0.7, "cm"), 
  legend.text = element_text(size = 14))

     
ui <- fluidPage(
  titlePanel("JCAP CRISPR Mixscape Pipeline"),
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("countsfile", "Upload counts matrix (csv)", accept = ".csv"),
      fileInput("metadatafile", "Upload metadata (csv)", accept = ".csv"),
      actionButton("run_basic", "Run Normalization & UMAP"),
      sliderInput("num_neighbors", "Number of Neighbors (for Mixscape/UMAP):",
                  min = 2, max = 50, value = 20, step = 1),
      actionButton("show_summary", "Show Summary Stats"),
      actionButton("run_mixscape", "Run Mixscape Analysis"),
      downloadButton("downloadSummary", "Download Summary Stats"),
      
      downloadButton("downloadKOgenes", "Download KO Gene Table")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("How To", includeMarkdown("How_To.txt")),
        
        tabPanel("UMAP Plots", plotOutput("umapPlots")),
        tabPanel("Summary table", DT::dataTableOutput("summaryTable")),
        tabPanel("KO Genes", DT::dataTableOutput("koGenes")),
        tabPanel("KO % Barplots", plotlyOutput("koBarPlot")),
        tabPanel("Perturbation Score", plotlyOutput("perturbPlot")),
        tabPanel("Posterior Violin", plotlyOutput("vlnPlot")),
        tabPanel("Heatmap", plotlyOutput("mixscapeHeatmap"))
      )
    )
  )
)


server <- function(input, output, session) {
  vals <- reactiveValues(seurat = NULL, p1 = NULL, p2 = NULL, p3 = NULL)
  
  # Step 1: Run normalization & UMAP
  observeEvent(input$run_basic, {
    req(input$countsfile)
    req(input$metadatafile)
    showModal(modalDialog("Running normalization and UMAP. Please wait...", footer=NULL))
    
    counts <- read.csv(input$countsfile$datapath, row.names=1)
    metadata <- read.csv(input$metadatafile$datapath, row.names=1)
    shared <- intersect(colnames(counts), rownames(metadata))
    counts <- counts[, shared]
    metadata <- metadata[shared, ]
    
    so <- CreateSeuratObject(counts = counts, meta.data = metadata)
    DefaultAssay(so) <- 'RNA'
    so <- NormalizeData(so) %>% FindVariableFeatures() %>% ScaleData()
    n_cells <- ncol(so)
    n_genes <- nrow(so)
    npcs <- min(40, n_cells - 1, n_genes - 1)
    if (npcs < 2) {
      showNotification("Not enough genes/cells for PCA! Please upload a bigger subset.", type = "error")
      return(NULL)  # or stop processing this event
    }
  
    # ...after creating Seurat object and running PCA:
    n_cells <- ncol(so)
    n_genes <- nrow(so)
    npcs <- min(40, n_cells - 1, n_genes - 1)
    if (npcs < 2) {
      showNotification("Not enough genes/cells for PCA! Please upload a bigger subset.", type = "error")
      removeModal()
      return(NULL)
    }
    so <- RunPCA(so, npcs = npcs)
    
    # Now, RunUMAP with correct dimensions and neighbors:
    umap_dims <- 1:npcs
    umap_neighbors <- min(input$num_neighbors, n_cells - 1)
    if (umap_neighbors < 2) umap_neighbors <- 2
    
    so <- RunUMAP(so, dims = umap_dims, n.neighbors = umap_neighbors)
    
    
 
    
    # Save to reactive values
    vals$seurat <- so
    vals$p1 <- DimPlot(so, group.by = 'replicate', pt.size = 0.2, reduction = "umap", repel = TRUE)
    vals$p2 <- DimPlot(so, group.by = 'Phase', pt.size = 0.2, reduction = "umap", repel = TRUE)
    vals$p3 <- DimPlot(so, group.by = 'crispr', pt.size = 0.2, reduction = "umap", repel = TRUE)
    
    # Dynamically update neighbors max/min after loading data
    n_cells <- ncol(so)
    updateNumericInput(session, "num_neighbors",
                       max = max(2, n_cells - 1),
                       value = min(input$num_neighbors, n_cells - 1))
    
    removeModal()
    showNotification("Basic normalization and UMAP complete!", type = "message")
  })
  observeEvent(input$run_mixscape, {
    req(vals$seurat)
    showModal(modalDialog("Running Mixscape analysis. Please wait...", footer = NULL))
    so <- vals$seurat
    
    # Check required columns
    req("gene" %in% colnames(so@meta.data))
    req("replicate" %in% colnames(so@meta.data))
    req("NT" %in% so@meta.data$gene)  # Make sure NTs exist
    
    min_group_cells <- 3  # Change as needed
    
    # Iteratively filter small groups
    repeat {
      group_counts <- table(so@meta.data$gene, so@meta.data$replicate)
      too_small <- group_counts < min_group_cells & group_counts > 0
      if (!any(too_small)) break
      drop_cells <- rownames(so@meta.data)[
        mapply(function(g, r) group_counts[g, r] < min_group_cells,
               so@meta.data$gene, so@meta.data$replicate)
      ]
      if (length(drop_cells) == 0) break
      so <- subset(so, cells = setdiff(rownames(so@meta.data), drop_cells))
    }
    
    # Remove replicates missing any group
    group_counts <- table(so@meta.data$gene, so@meta.data$replicate)
    valid_reps <- colnames(group_counts)[colSums(group_counts > 0) == length(unique(so@meta.data$gene))]
    so <- subset(so, cells = rownames(so@meta.data)[so@meta.data$replicate %in% valid_reps])
    
    cat("---- Group sizes AFTER filtering ----\n")
    print(table(so@meta.data$gene, so@meta.data$replicate))
    cat("-------------------------------------\n")
    
    # Check for enough groups and replicates
    if (length(unique(so@meta.data$gene)) < 2 || length(unique(so@meta.data$replicate)) < 2) {
      removeModal()
      showNotification(
        sprintf("Not enough groups or replicates after filtering (min %d cells/group required).", min_group_cells),
        type = "error"
      )
      return(NULL)
    }
    
    n_cells <- ncol(so)
    n_genes <- nrow(so)
    npcs <- min(40, n_cells - 1, n_genes - 1)
    if (npcs < 2) {
      removeModal()
      showNotification("Not enough genes/cells for PCA!", type = "error")
      return(NULL)
    }
    ndims <- npcs
    
    group_counts <- table(so@meta.data$gene, so@meta.data$replicate)
    min_group_size <- min(group_counts[group_counts > 0])
    safe_neighbors <- max(2, min(if(!is.null(input$num_neighbors)) input$num_neighbors else 7, min_group_size - 1))
    
    # Run CalcPerturbSig
    success <- TRUE
    tryCatch({
      so <- CalcPerturbSig(
        object = so,
        assay = "RNA",
        slot = "data",
        gd.class = "gene",
        nt.cell.class = "NT",
        reduction = "pca",
        ndims = ndims,
        num.neighbors = safe_neighbors,
        split.by = "replicate",
        new.assay.name = "PRTB"
      )
    }, error = function(e) {
      print(e)
      showNotification(paste("Error in CalcPerturbSig:", e$message), type = "error")
      success <<- FALSE
    })
    
    if (!success || !("PRTB" %in% Assays(so))) {
      removeModal()
      showNotification("CalcPerturbSig failed or did not add 'PRTB' assay—cannot run Mixscape.", type = "error")
      return(NULL)
    }
    
    # Run Mixscape
    tryCatch({
      so <- RunMixscape(
        object = so,
        assay = "PRTB",
        slot = "scale.data",
        labels = "gene",
        nt.class.name = "NT",
        min.de.genes = 5,
        iter.num = 10,
        de.assay = "RNA",
        verbose = FALSE,
        prtb.type = "KO"
      )
    }, error = function(e) {
      print(e)
      showNotification(paste("Error in RunMixscape:", e$message), type = "error")
      removeModal()
      return(NULL)
    })
    
    print(table(so@meta.data$mixscape_class))
    
    mc <- so@meta.data$mixscape_class
    if (all(c("IFNGR2 KO", "NT") %in% mc)) {
      # Differential gene expression: KO vs NT
      de_tab <- tryCatch({
        FindMarkers(so, ident.1 = "IFNGR2 KO", ident.2 = "NT", group.by = "mixscape_class", assay = "RNA")
      }, error=function(e) data.frame(message = "DGE failed"))
      vals$koGenes <- de_tab
    } else {
      vals$koGenes <- data.frame(message = "KO or NT class not present")
    }
    
    # Build KO % Barplot Data Frame (vals$df3)
    df <- as.data.frame(table(
      gene = so@meta.data$gene,
      replicate = so@meta.data$replicate,
      mixscape_class = so@meta.data$mixscape_class
    ))
    
    # Calculate percentages *within* each (gene, replicate) group
    df <- df %>%
      dplyr::group_by(gene, replicate) %>%
      dplyr::mutate(percent = Freq / sum(Freq)) %>%
      dplyr::ungroup()
    
    # For barplot code, use Var1 to match your plotting code
    vals$df3 <- df %>%
      dplyr::rename(Var1 = mixscape_class, value = percent)
    
    # --- THIS IS CRUCIAL! ---
    vals$seurat <- so  # update the seurat object with all new results
    removeModal()
  })
  

  library(plotly)
  output$umapPlots <- renderPlot({
    req(vals$p1, vals$p2, vals$p3)
    (vals$p1 / vals$p2 + plot_layout(guides = 'auto')) | vals$p3
  })
  

  
  output$koBarPlot <- renderPlotly({
    req(vals$df3)
    p1 <- ggplot(vals$df3, aes(x = replicate, y = value*100, fill = Var1)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      theme_classic() +
      scale_fill_manual(values = c("grey49", "grey79", "coral1")) +
      ylab("% of cells") +
      xlab("Replicate") +
      theme(axis.text.x = element_text(size = 18, hjust = 1),
            axis.text.y = element_text(size = 18),
            axis.title = element_text(size = 16),
            strip.text = element_text(size = 16, face = "bold")) +
      facet_wrap(vars(gene), ncol = 5, scales = "free") +
      labs(fill = "mixscape class") +
      theme(legend.title = element_text(size = 14),
            legend.text = element_text(size = 12))
    ggplotly(p1)
  })
  
  output$perturbPlot <- renderPlotly({
    req(vals$seurat)
    tryCatch({
      p <- PlotPerturbScore(
        object = vals$seurat, 
        target.gene.ident = "IFNGR2", 
        mixscape.class = "mixscape_class", 
        col = "coral2"
      ) + labs(fill = "mixscape class")
      ggplotly(p)
    }, error = function(e) {
      ggplotly(ggplot() + ggtitle("Error: IFNGR2 or mixscape_class column missing"))
    })
  })
  
  output$vlnPlot <- renderPlotly({
    req(vals$seurat)
    tryCatch({
      p <- VlnPlot(vals$seurat, "mixscape_class_p_ko", idents = c("NT", "IFNGR2 KO", "IFNGR2 NP")) +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
              axis.text = element_text(size = 16), plot.title = element_text(size = 20)) +
        NoLegend() + ggtitle("mixscape posterior probabilities")
      ggplotly(p)
    }, error=function(e) {
      ggplotly(ggplot() + ggtitle("Error: Violin plot variables missing"))
    })
  })
  
  output$mixscapeHeatmap <- renderPlotly({
    req(vals$seurat)
    tryCatch({
      classes <- unique(vals$seurat@meta.data$mixscape_class)
      if (!all(c("NT", "IFNGR2 KO") %in% classes)) {
        stop("Both 'NT' and 'IFNGR2 KO' classes are required for the heatmap.")
      }
      p <- MixscapeHeatmap(
        object = vals$seurat,
        ident.1 = "NT",
        ident.2 = "IFNGR2 KO",
        balanced = FALSE,
        assay = "RNA",
        max.genes = 20,
        angle = 0,
        group.by = "mixscape_class",
        max.cells.group = 300,
        size = 6.5
      ) + NoLegend() + theme(axis.text.y = element_text(size = 16))
      ggplotly(p)
    }, error = function(e) {
      ggplotly(ggplot() + ggtitle(paste("Error: MixscapeHeatmap failed\n", e$message)))
    })
  })
  
  output$koGenes <- DT::renderDataTable({
    req(vals$koGenes)
    DT::datatable(vals$koGenes, options = list(pageLength = 10))
  })

  output$downloadKOgenes <- downloadHandler(
    filename = function() {
      paste0("KO_genes_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(vals$koGenes, file, row.names = TRUE)
    }
  )

  observeEvent(input$show_summary, {
    req(vals$seurat)
    
    showNotification("Generating summary statistics...", type = "message", duration = 2)
    
    meta <- vals$seurat@meta.data
    
    # Compute stats
    n_cells <- nrow(meta)
    n_genes <- nrow(vals$seurat)
    n_ko <- sum(meta$mixscape_class == "IFNGR2 KO")
    n_nt <- sum(meta$mixscape_class == "NT")
    n_np <- sum(meta$mixscape_class == "IFNGR2 NP")
    n_reps <- length(unique(meta$replicate))
    n_guides <- length(unique(meta$guide_ID))
    n_hvgs <- length(VariableFeatures(vals$seurat))
    
    stats <- data.frame(
      Metric = c("Total Cells", "Total Genes", 
                 "KO Cells", "NT Cells", "NP Cells",
                 "Replicates", "sgRNAs", "Variable Genes"),
      Value = c(n_cells, n_genes, n_ko, n_nt, n_np, n_reps, n_guides, n_hvgs)
    )
    
    vals$summaryStats <- stats
    
    showNotification("Summary stats table generated!", type = "message", duration = 2)
  })
  
  output$summaryTable <- DT::renderDataTable({
    req(vals$summaryStats)   # Or however you compute your stats
    DT::datatable(vals$summaryStats, options = list(pageLength = 10))
  })
  
  output$downloadSummary <- downloadHandler(
    filename = function() {
      paste0("summary_stats_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(vals$summaryStats)
      write.csv(vals$summaryStats, file, row.names = FALSE)
    }
  )
  output$downloadSummary <- downloadHandler(
    filename = function() {
      paste0("summary_stats_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(vals$summaryStats)
      write.csv(vals$summaryStats, file, row.names = FALSE)
    }
  )
  

}
  

shinyApp(ui, server)
