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
library(enrichR)       # for Human (Enrichr)
library(gprofiler2)    # for Mouse/Fly/Zebrafish (g:Profiler)


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
      selectInput("ko_label", "KO group to analyze", choices = NULL),
      numericInput("min_de_genes", "min.de.genes", value = 3, min = 1, max = 50, step = 1),
      numericInput("iter_num",     "iter.num",     value = 20, min = 1, max = 200, step = 1),
      actionButton("run_mixscape", "Run Mixscape Analysis"),
      # ---- Pathway Enrichment ----
      hr(),
      h4("Pathway Enrichment"),
      selectInput(
        "species", "Species",
        choices = c(
          "Human" = "hsapiens",
          "Mouse" = "mmusculus",
          "Fly (D. melanogaster)" = "dmelanogaster",
          "Zebrafish (D. rerio)" = "drerio"
        ),
        selected = "hsapiens"
      ),
      helpText("Auto-routing: Human → Enrichr; Mouse/Fly/Zebrafish → g:Profiler"),
      actionButton("enrich_all",  "Enrich ALL DE genes"),
      actionButton("enrich_up",   "Enrich UP genes"),
      actionButton("enrich_down", "Enrich DOWN genes"),
      # Downloads
      downloadButton("downloadSummary", "Download Summary Stats"),
      downloadButton("downloadKOgenes", "Download KO Gene Table"),
      
      downloadButton("dl_enrich_all",  "Download: All (CSV)"),
      downloadButton("dl_enrich_up",   "Download: Up (CSV)"),
      downloadButton("dl_enrich_down", "Download: Down (CSV)")
      
    ),
    mainPanel(
      tabsetPanel(
        id = "tabs",
        tabPanel(
          "How To",
          if (file.exists("How_To.txt")) includeMarkdown("How_To.txt") else div("How_To.txt not found.")
        ),
        tabPanel("UMAP Plots", plotOutput("umapPlots")),
        tabPanel("Summary table", DT::dataTableOutput("summaryTable")),
        tabPanel("KO Genes", DT::dataTableOutput("koGenes")),
        tabPanel("KO % Barplots", plotlyOutput("koBarPlot")),
        tabPanel("Perturbation (Safe Posterior)", plotlyOutput("perturbPlot")),
        tabPanel("Posterior Violin", plotlyOutput("vlnPlot")),
        tabPanel("Heatmap (Safe DE)", plotlyOutput("mixscapeHeatmap")),  # <-- close this tabPanel and add comma
        tabPanel(
          "Pathway Enrichment",
          tabsetPanel(
            id = "enrich_tabs",
            tabPanel("All DE (table)",  DT::dataTableOutput("enrich_all_tbl")),
            tabPanel("All DE (barplot)", plotlyOutput("enrich_all_bar")),
            tabPanel("Up (table)",       DT::dataTableOutput("enrich_up_tbl")),
            tabPanel("Up (barplot)",     plotlyOutput("enrich_up_bar")),
            tabPanel("Down (table)",     DT::dataTableOutput("enrich_down_tbl")),
            tabPanel("Down (barplot)",   plotlyOutput("enrich_down_bar"))
          )
        )
      )
    )))
    

    
    



server <- function(input, output, session) {
  vals <- reactiveValues(seurat = NULL, p1 = NULL, p2 = NULL, p3 = NULL)
  `%||%` <- function(a, b) {
    if (is.null(a)) return(b)
    if (length(a) == 0) return(b)
    if (length(a) == 1 && is.na(a)) return(b)
    a
  }
  vals$enrich_all  <- NULL
  vals$enrich_up   <- NULL
  vals$enrich_down <- NULL
  
  
  
  # ---------- Enrichment helpers ----------
  # Map organism code -> which backend to use
  .which_backend <- function(org) {
    if (tolower(org) == "hsapiens") "enrichr" else "gprof"
  }
  
  # g:Profiler (returns standardized table or NULL)
  .run_gprof <- function(genes, org, sources = c("GO:BP","GO:MF","GO:CC","REAC","KEGG")) {
    if (length(genes) < 3) return(NULL)
    gp <- tryCatch({
      gprofiler2::gost(
        query = unique(genes),
        organism = org,
        sources = sources,
        correction_method = "fdr",
        evcodes = FALSE
      )
    }, error = function(e) NULL)
    if (is.null(gp) || is.null(gp$result) || nrow(gp$result) == 0) return(NULL)
    res <- gp$result
    # standardize to Term / Adjusted.P.value columns
    adj <- if ("p_adjusted" %in% names(res)) res$p_adjusted else p.adjust(res$p_value, "fdr")
    out <- data.frame(
      Term = paste0(res$term_name, " (", res$source, ")"),
      Overlap = paste0(res$intersection_size, "/", res$term_size),
      P.value = res$p_value,
      Adjusted.P.value = adj,
      stringsAsFactors = FALSE, check.names = FALSE
    )
    out[order(out$Adjusted.P.value, out$P.value), , drop = FALSE]
  }
  
  # Enrichr (returns standardized table or NULL)
  # Uses a compact, broadly-available set of libs (works best for Human)
  .enrichr_libs <- c("GO_Biological_Process_2021",
                     "GO_Molecular_Function_2021",
                     "GO_Cellular_Component_2021",
                     "Reactome_2022",
                     "KEGG_2021_Human")
  
  .run_enrichr <- function(genes, libs = .enrichr_libs) {
    if (length(genes) < 3) return(NULL)
    out <- tryCatch(enrichR::enrichr(unique(genes), libs), error = function(e) NULL)
    if (is.null(out)) return(NULL)
    # bind and standardize
    tbl <- do.call(rbind, lapply(names(out), function(lib) {
      df <- out[[lib]]
      if (is.null(df) || nrow(df) == 0) return(NULL)
      df$source <- lib
      df
    }))
    if (is.null(tbl) || nrow(tbl) == 0) return(NULL)
    std <- tbl[, c("Term","Overlap","P.value","Adjusted.P.value","source"), drop = FALSE]
    std$Term <- paste0(std$Term, " (", std$source, ")")
    std <- std[order(std$Adjusted.P.value, std$P.value), , drop = FALSE]
    rownames(std) <- NULL
    std[, c("Term","Overlap","P.value","Adjusted.P.value"), drop = FALSE]
  }
  
  # Unified entry point based on species
  .do_enrichment <- function(gene_symbols, org_code) {
    bk <- .which_backend(org_code)
    if (bk == "enrichr") .run_enrichr(gene_symbols) else .run_gprof(gene_symbols, org_code)
  }
  
  # Pick DE gene symbols from your existing vals$koGenes
  # Accept both Seurat v4/v5 names: avg_log2FC vs avg_logFC; p_val_adj vs p_val_adj
  .pick_de_genes <- function(tab, mode = c("all","up","down")) {
    mode <- match.arg(mode)
    if (is.null(tab) || !nrow(tab)) return(character(0))
    fc_col <- if ("avg_log2FC" %in% names(tab)) "avg_log2FC" else if ("avg_logFC" %in% names(tab)) "avg_logFC" else NA
    padj_col <- if ("p_val_adj" %in% names(tab)) "p_val_adj" else if ("p_adj" %in% names(tab)) "p_adj" else NA
    if (is.na(fc_col) || is.na(padj_col)) return(character(0))
    keep <- tab[is.finite(tab[[padj_col]]) & tab[[padj_col]] < 0.05, , drop = FALSE]
    if (nrow(keep) == 0) return(character(0))
    if (mode == "up")   keep <- keep[keep[[fc_col]] >  0, , drop = FALSE]
    if (mode == "down") keep <- keep[keep[[fc_col]] <  0, , drop = FALSE]
    rownames(keep)
  }
  
  # Simple barplot (top 10 by FDR)
  .enrich_bar <- function(df, title_txt) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df <- df[is.finite(df$Adjusted.P.value) & df$Adjusted.P.value > 0, , drop = FALSE]
    if (!nrow(df)) return(NULL)
    df$mlog10 <- -log10(df$Adjusted.P.value)
    df <- head(df[order(df$Adjusted.P.value, -df$mlog10), , drop = FALSE], 10)
    plotly::plot_ly(df, x = ~mlog10, y = ~reorder(Term, mlog10),
                    type = "bar", orientation = "h") %>%
      plotly::layout(title = title_txt,
                     xaxis = list(title = "-log10(FDR)"),
                     yaxis = list(title = "", automargin = TRUE))
  }
  
  
  # Is this assay a Seurat v5 layered assay?
  .is_assay5 <- function(so, assay) inherits(so[[assay]], "Assay5")
  
  # Ensure an assay has a counts layer; for "derived" assays (e.g., PRTB),
  # it can safely be an all-zero sparse matrix with correct dimnames.
  .ensure_counts_layer <- function(so, assay = DefaultAssay(so)) {
    if (!assay %in% Assays(so)) return(so)
    
    has_counts <- tryCatch({
      lyr <- GetAssayData(so, assay = assay, layer = "counts")
      !is.null(lyr) && length(lyr) > 0
    }, error = function(...) FALSE)
    
    if (!has_counts) {
      # Try reconstruct from data (log1p), else make zeros
      dat <- tryCatch(GetAssayData(so, assay = assay, layer = "data"), error = function(...) NULL)
      if (is.null(dat) || nrow(dat) == 0)
        dat <- tryCatch(GetAssayData(so, assay = assay, slot = "data"), error = function(...) NULL)
      
      if (!is.null(dat) && nrow(dat) > 0) {
        cnt <- Matrix::Matrix(round(pmax(exp(dat) - 1, 0)), sparse = TRUE)
        dimnames(cnt) <- list(rownames(dat), colnames(dat))
      } else {
        rn <- rownames(so[[assay]])
        cn <- colnames(so)
        cnt <- Matrix::Matrix(0, nrow = length(rn), ncol = length(cn), sparse = TRUE,
                              dimnames = list(rn, cn))
      }
      
      if (.is_assay5(so, assay)) {
        so <- SetAssayData(so, assay = assay, layer = "counts", new.data = cnt)
      } else {
        so[[assay]]@counts <- cnt  # v4 fallback
      }
    }
    
    # make sure a "data" layer exists too (some v5 paths expect it)
    has_data <- tryCatch({
      lyr <- GetAssayData(so, assay = assay, layer = "data"); !is.null(lyr) && length(lyr) > 0
    }, error = function(...) FALSE)
    if (!has_data) {
      dat <- tryCatch(GetAssayData(so, assay = assay, slot = "data"), error = function(...) NULL)
      if (!is.null(dat) && length(dat) > 0 && .is_assay5(so, assay)) {
        so <- SetAssayData(so, assay = assay, layer = "data", new.data = dat)
      }
    }
    so
  }
  .harden_prtb <- function(so) {
    if (!"PRTB" %in% Assays(so)) return(so)
    
    # Prefer v5 layered reads; fall back to slot for v4
    get_layer <- function(a, layer) {
      tryCatch(GetAssayData(so, assay = a, layer = layer), error = function(...) NULL)
    }
    get_slot <- function(a, slot) {
      tryCatch(GetAssayData(so, assay = a, slot = slot), error = function(...) NULL)
    }
    
    dat <- get_layer("PRTB", "data")
    if (is.null(dat) || nrow(dat) == 0) dat <- get_slot("PRTB", "data")
    
    # Row/col names to use for layer writes
    rn <- if (!is.null(dat) && nrow(dat) > 0) rownames(dat) else rownames(so[["PRTB"]])
    cn <- colnames(so)
    
    # 1) Give PRTB a harmless "counts" layer (zeros) to satisfy layer-seeking helpers
    if (!is.null(rn) && length(rn) > 0 && length(cn) > 0) {
      zeros <- Matrix::Matrix(0, nrow = length(rn), ncol = length(cn), sparse = TRUE,
                              dimnames = list(rn, cn))
      so <- tryCatch(
        SetAssayData(so, assay = "PRTB", layer = "counts", new.data = zeros),
        error = function(e) { so[["PRTB"]]@counts <- zeros; so } # v4 fallback
      )
    }
    
    # 2) If scale.data exists with different features, clear it before rescaling
    sd <- tryCatch(so[["PRTB"]]@scale.data, error = function(...) NULL)
    if (!is.null(sd) && length(sd)) {
      sd_rn <- rownames(sd)
      # If we can't compare safely, just clear to be safe
      if (is.null(dat) || is.null(sd_rn) || !identical(sd_rn, rownames(dat))) {
        so[["PRTB"]]@scale.data <- matrix(numeric(0), nrow = 0, ncol = 0)
      }
    }
    
    # 3) Ensure a proper 'data' layer exists for v5 callers
    if (!is.null(dat) && nrow(dat) > 0) {
      so <- tryCatch(
        SetAssayData(so, assay = "PRTB", layer = "data", new.data = dat),
        error = function(e) so
      )
    }
    
    # 4) Scale PRTB cleanly
    DefaultAssay(so) <- "PRTB"
    so <- ScaleData(so, assay = "PRTB", verbose = FALSE)
    DefaultAssay(so) <- "RNA" # restore
    
    so
  }
  
  
  
  
  # Harden an assay: ensure counts & data exist as BOTH layer (v5) and slot (v4)
  .ensure_counts_layer <- function(so, assay = DefaultAssay(so)) {
    if (!assay %in% Assays(so)) return(so)
    
    # ---- ensure counts ----
    has_counts_layer <- tryCatch({
      lyr <- GetAssayData(so, assay = assay, layer = "counts")
      !is.null(lyr) && length(lyr) > 0
    }, error = function(...) FALSE)
    
    if (!has_counts_layer) {
      # prefer data (log1p) to reconstruct counts
      dat <- tryCatch(GetAssayData(so, assay = assay, layer = "data"), error = function(...) NULL)
      if (is.null(dat) || nrow(dat) == 0)
        dat <- tryCatch(GetAssayData(so, assay = assay, slot = "data"), error = function(...) NULL)
      
      cnt <- NULL
      if (!is.null(dat) && nrow(dat) > 0) {
        cnt <- Matrix::Matrix(round(pmax(exp(dat) - 1, 0)), sparse = TRUE)
        dimnames(cnt) <- list(rownames(dat), colnames(dat))
      } else {
        rn <- tryCatch(rownames(so[[assay]]), error = function(...) NULL)
        cn <- tryCatch(colnames(so), error = function(...) NULL)
        if (!is.null(rn) && !is.null(cn)) {
          cnt <- Matrix::Matrix(0, nrow = length(rn), ncol = length(cn), sparse = TRUE,
                                dimnames = list(rn, cn))
        }
      }
      
      if (!is.null(cnt)) {
        so <- tryCatch(
          SetAssayData(so, assay = assay, layer = "counts", new.data = cnt),
          error = function(e) { so[[assay]]@counts <- cnt; so }
        )
      }
    }
    
    # ---- ensure data ----
    has_data_layer <- tryCatch({
      lyr <- GetAssayData(so, assay = assay, layer = "data")
      !is.null(lyr) && length(lyr) > 0
    }, error = function(...) FALSE)
    
    if (!has_data_layer) {
      dat <- tryCatch(GetAssayData(so, assay = assay, slot = "data"), error = function(...) NULL)
      if (!is.null(dat) && length(dat) > 0) {
        so <- tryCatch(
          SetAssayData(so, assay = assay, layer = "data", new.data = dat),
          error = function(e) so
        )
      }
    }
    
    so
  }
  
  # Quick debug line to know what's missing when it fails
  .debug_layers <- function(so, assays = c("RNA", "PRTB")) {
    for (a in assays[assays %in% Assays(so)]) {
      has_cnt <- tryCatch({
        lyr <- GetAssayData(so, assay = a, layer = "counts")
        !is.null(lyr) && length(lyr) > 0
      }, error = function(...) FALSE)
      has_dat <- tryCatch({
        lyr <- GetAssayData(so, assay = a, layer = "data")
        !is.null(lyr) && length(lyr) > 0
      }, error = function(...) FALSE)
      message(sprintf("[assay %s] counts_layer=%s, data_layer=%s, class=%s",
                      a, has_cnt, has_dat, class(so[[a]])[1]))
    }
  }
  
  
  
  # ---- Robustly ensure RNA has BOTH a v5 layer and a v4 slot for counts ----
  .ensure_counts_layer <- function(so, assay = "RNA") {
    if (!assay %in% Assays(so)) return(so)
    
    # Try v5 "layer" read first
    has_counts_layer <- FALSE
    has_counts_layer <- tryCatch({
      lyr <- GetAssayData(so, assay = assay, layer = "counts")
      !is.null(lyr) && length(lyr) > 0
    }, error = function(...) FALSE)
    
    if (!has_counts_layer) {
      # Try to derive counts from "data" (log1p-normalized by default)
      dat <- tryCatch(GetAssayData(so, assay = assay, layer = "data"), error = function(...) NULL)
      if (is.null(dat) || nrow(dat) == 0)
        dat <- tryCatch(GetAssayData(so, assay = assay, slot = "data"), error = function(...) NULL)
      
      if (!is.null(dat) && nrow(dat) > 0) {
        # inverse of log1p: exp(dat) - 1
        cnt <- Matrix::Matrix(round(pmax(exp(dat) - 1, 0)), sparse = TRUE)
        dimnames(cnt) <- list(rownames(dat), colnames(dat))
      } else {
        # last resort: make an all-zero counts matrix with correct dimnames
        rn <- tryCatch(rownames(so[[assay]]), error = function(...) NULL)
        cn <- tryCatch(colnames(so), error = function(...) NULL)
        if (is.null(rn) || is.null(cn)) return(so)
        cnt <- Matrix::Matrix(0, nrow = length(rn), ncol = length(cn), sparse = TRUE,
                              dimnames = list(rn, cn))
      }
      
      # Write to v5 "layer" if available; mirror to v4 slot as fallback
      so <- tryCatch(
        SetAssayData(so, assay = assay, layer = "counts", new.data = cnt),
        error = function(e) {
          # v4 fallback
          so[[assay]]@counts <- cnt
          so
        }
      )
    }
    
    # Also make sure there is a "data" layer/slot (some functions expect it)
    has_data_layer <- tryCatch({
      lyr <- GetAssayData(so, assay = assay, layer = "data"); !is.null(lyr) && length(lyr) > 0
    }, error = function(...) FALSE)
    if (!has_data_layer) {
      dat <- tryCatch(GetAssayData(so, assay = assay, slot = "data"), error = function(...) NULL)
      if (!is.null(dat) && length(dat) > 0) {
        so <- tryCatch(
          SetAssayData(so, assay = assay, layer = "data", new.data = dat),
          error = function(e) so
        )
      }
    }
    
    so
  }
  
  # safe ggplotly wrapper that prints a friendly message when we don't have a plot
  as_plotly_or_msg <- function(p, msg) {
    if (inherits(p, "ggplot")) return(plotly::ggplotly(p))
    plotly::ggplotly(ggplot() + ggtitle(msg))
  }
  
  ko_label_choices <- function(so) {
    if (is.null(so)) return(character(0))
    
    # Prefer whatever Mixscape actually produced
    class_cols <- intersect(c("mixscape_class_active","mixscape_class","mixscape_class_global"),
                            colnames(so@meta.data))
    ko_from_class <- character(0)
    if (length(class_cols)) {
      cls <- as.character(so@meta.data[[class_cols[1]]])
      cls <- cls[!is.na(cls)]
      if (length(cls)) {
        ko_from_class <- unique(c(
          sort(unique(cls[grepl(" KO$", cls, ignore.case = TRUE)])),
          if (any(toupper(cls) == "KO")) "KO" else NULL
        ))
      }
    }
    if (length(ko_from_class)) return(ko_from_class)
    
    # Fallback: propose "<GENE> KO" for each gene label present (excluding NT/controls)
    if ("gene" %in% colnames(so@meta.data)) {
      gs <- setdiff(sort(unique(toupper(so@meta.data$gene))), c("NT","CTRL","CONTROL"))
      if (length(gs)) return(paste0(gs, " KO"))
    }
    character(0)
  }
  
  observeEvent(vals$seurat, {
    choices <- ko_label_choices(vals$seurat)
    default_sel <- if ("IFNGR2 KO" %in% choices) "IFNGR2 KO" else (choices[1] %||% NULL)
    updateSelectInput(session, "ko_label", choices = choices, selected = default_sel)
  }, ignoreInit = FALSE)
  
 
  
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
    so <- CreateSeuratObject(counts = counts, meta.data = metadata)
    DefaultAssay(so) <- "RNA"
    
    # NEW: make sure counts & data layers exist
    so <- .ensure_counts_layer(so, "RNA")
    
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
    updateSliderInput(session, "num_neighbors",
                       max = max(2, n_cells - 1),
                       value = min(input$num_neighbors, n_cells - 1))
    
    removeModal()
    showNotification("Basic normalization and UMAP complete!", type = "message")
  })
  observeEvent(input$run_mixscape, {
    req(vals$seurat)
    showModal(modalDialog("Running Mixscape analysis...", footer = NULL))
    
    so <- vals$seurat
    DefaultAssay(so) <- "RNA"
    
    ## -- 0) Basic metadata checks/cleanup
    req(all(c("gene","replicate") %in% colnames(so@meta.data)))
    so$gene      <- toupper(as.character(so$gene))
    so$replicate <- as.character(so$replicate)
    
    keep <- rownames(so@meta.data)[nzchar(so$gene) & nzchar(so$replicate)]
    if (length(keep) < 50) { removeModal(); showNotification("Too few valid cells after metadata cleanup.", type="error"); return(NULL) }
    so <- subset(so, cells = keep)
    
    if (!any(so$gene == "NT")) { removeModal(); showNotification("No 'NT' cells present in 'gene'.", type="error"); return(NULL) }
    if (length(unique(so$gene)) < 2) { removeModal(); showNotification("Need ≥2 perturbation labels in 'gene'.", type="error"); return(NULL) }
    
    ## -- 1) Filter tiny (gene,replicate) groups
    min_group_cells <- 3
    repeat {
      gc <- as.matrix(table(so$gene, so$replicate))
      tiny <- gc < min_group_cells & gc > 0
      if (!any(tiny)) break
      drop_cells <- rownames(so@meta.data)[
        mapply(function(g,r) gc[g, r, drop=FALSE] < min_group_cells, so$gene, so$replicate)
      ]
      if (!length(drop_cells)) break
      so <- subset(so, cells = setdiff(Cells(so), drop_cells))
    }
    
    # keep only replicates that still have NT
    gc <- as.matrix(table(so$gene, so$replicate))
    reps_with_nt <- colnames(gc)[gc["NT", ] > 0]
    if (!length(reps_with_nt)) { removeModal(); showNotification("No replicate contains NT after filtering.", type="error"); return(NULL) }
    so <- subset(so, cells = rownames(so@meta.data)[so$replicate %in% reps_with_nt])
    ## -- 2) Ensure RNA has counts+data; redo Normalize/Scale/PCA (v5-safe)
    # Get normalized 'data' if present; otherwise NormalizeData will create it.
    dat <- tryCatch(GetAssayData(so, assay="RNA", layer="data"), error=function(e) NULL)
    if (is.null(dat) || nrow(dat) == 0)
      dat <- tryCatch(GetAssayData(so, assay="RNA", slot="data"), error=function(e) NULL)
    
    if (is.null(dat) || nrow(dat) == 0) {
      so <- NormalizeData(so) %>% FindVariableFeatures() %>% ScaleData()
      dat <- GetAssayData(so, assay="RNA", slot="data")
    }
    
    # If counts are missing, write them to the **existing** assay (don’t replace it)
    raw <- Matrix::Matrix(round(pmax(exp(dat) - 1, 0)), sparse = TRUE)
    dimnames(raw) <- list(rownames(dat), colnames(dat))
    if (.is_assay5(so, "RNA")) {
      so <- SetAssayData(so, assay = "RNA", layer = "counts", new.data = raw)
    } else {
      so[["RNA"]]@counts <- raw  # v4 fallback
    }
    DefaultAssay(so) <- "RNA"
    
    # Now guarantee v5 layers/slots present
    so <- .ensure_counts_layer(so, "RNA")
    
    so <- NormalizeData(so) %>% FindVariableFeatures()
    so <- ScaleData(so, features = VariableFeatures(so), verbose = FALSE)
    

    
    # PCA
    n_cells <- ncol(so); n_genes <- nrow(so)
    npcs <- max(2, min(40, n_cells - 1, n_genes - 1))
    so <- RunPCA(so, npcs = npcs)
    
    emb <- Embeddings(so, "pca")
    if (is.null(emb) || nrow(emb) == 0 || !identical(rownames(emb), colnames(so))) {
      removeModal(); showNotification("PCA embeddings missing/misaligned after rebuild.", type="error"); return(NULL)
    }
    ndims <- max(2, min(npcs, ncol(emb)))
    
    ## -- 3) Pick safe neighbors
    gc <- as.matrix(table(so$gene, so$replicate))
    min_group_size_per_rep <- apply(gc, 2, function(col) suppressWarnings(min(col[col > 0])))
    global_min_group_size <- suppressWarnings(min(min_group_size_per_rep[is.finite(min_group_size_per_rep)]))
    if (!is.finite(global_min_group_size)) { removeModal(); showNotification("No non-empty (gene,replicate) groups remain.", type="error"); return(NULL) }
    nt_n <- sum(so$gene == "NT")
    safe_neighbors <- max(2, min(if (!is.null(input$num_neighbors)) input$num_neighbors else 7,
                                 global_min_group_size - 1, nt_n - 1))
    
    
    
    ## -- 4) CalcPerturbSig  (ensure counts exist *right before* this call)
    so <- .ensure_counts_layer(so, "RNA")
    Idents(so) <- factor(so$gene)
    so <- .ensure_counts_layer(so, "RNA")
    .debug_layers(so, c("RNA"))
    
    ok <- TRUE
    so <- tryCatch({
      Seurat::CalcPerturbSig(
        object         = so,
        assay          = "RNA",
        slot           = "data",
        gd.class       = "gene",
        nt.cell.class  = "NT",
        reduction      = "pca",
        ndims          = ndims,
        num.neighbors  = safe_neighbors,
        new.assay.name = "PRTB",
        split.by       = "replicate",
        verbose        = TRUE
      )
    }, error = function(e) { message("CalcPerturbSig failed: ", e$message); ok <<- FALSE; return(so) })
    
    if (!ok || !("PRTB" %in% Assays(so))) {
      removeModal()
      showNotification("CalcPerturbSig did not create PRTB (check inputs).", type="error", duration=NULL)
      return(NULL)
    }
    
    
    if ("PRTB" %in% Assays(so)) {
      # PRTB typically lacks "counts"; create a harmless zero layer to satisfy helpers
      prtb_dat <- tryCatch(GetAssayData(so, assay = "PRTB", layer = "data"), error=function(...) NULL)
      if (is.null(prtb_dat) || nrow(prtb_dat) == 0)
        prtb_dat <- tryCatch(GetAssayData(so, assay = "PRTB", slot = "data"), error=function(...) NULL)
      
      ## -- AFTER CalcPerturbSig created PRTB and BEFORE RunMixscape:
      if ("PRTB" %in% Assays(so)) {
        so <- .harden_prtb(so)  # <<< ensures counts layer + aligned features + fresh scale.data
      }
      
      
      rn <- rownames(so[["PRTB"]]); cn <- colnames(so)
      zeros <- Matrix::Matrix(0, nrow = length(rn), ncol = length(cn), sparse = TRUE,
                              dimnames = list(rn, cn))
      if (.is_assay5(so, "PRTB")) {
        so <- SetAssayData(so, assay = "PRTB", layer = "counts", new.data = zeros)
      } else {
        so[["PRTB"]]@counts <- zeros
      }
      
      # Also ensure "data" layer mirrors whatever exists in v4 slot (for v5 callers)
      if (.is_assay5(so, "PRTB") && !is.null(prtb_dat) && nrow(prtb_dat) > 0) {
        so <- SetAssayData(so, assay = "PRTB", layer = "data", new.data = prtb_dat)
      }
    }
    
    
    ## -- 5) RunMixscape
    DefaultAssay(so) <- "PRTB"
    so <- ScaleData(so, assay = "PRTB", verbose = FALSE)
    DefaultAssay(so) <- "RNA"
    
    so <- tryCatch({
      Seurat::RunMixscape(
        object        = so,
        assay         = "PRTB",
        slot          = "scale.data",
        labels        = "gene",
        nt.class.name = "NT",
        min.de.genes  = input$min_de_genes %||% 3,
        iter.num      = input$iter_num %||% 20,
        de.assay      = "RNA",
        verbose       = FALSE,
        prtb.type     = "KO"
      )
    }, error = function(e) {
      print(e); showNotification(paste("Error in RunMixscape:", e$message), type="error"); removeModal(); return(NULL)
    })
    if (is.null(so)) return(NULL)
    
    ## -- 6) Choose active class & posterior columns
    class_candidates <- intersect(c("mixscape_class","mixscape_class_global"), colnames(so@meta.data))
    if (!length(class_candidates)) { removeModal(); showNotification("Mixscape finished but no class column was added.", type="error"); return(NULL) }
    
    cc_div <- sapply(class_candidates, function(x) length(na.omit(unique(so@meta.data[[x]]))))
    class_col <- class_candidates[which.max(cc_div)]
    if (cc_div[which.max(cc_div)] < 2 && length(class_candidates) > 1) {
      class_col <- setdiff(class_candidates, class_col)[1]
    }
    so@meta.data$mixscape_class_active <- so@meta.data[[class_col]]
    
    post_candidates <- intersect(c("mixscape_class_p_ko","mixscape_class_global_p_ko"), colnames(so@meta.data))
    pcol <- if (class_col == "mixscape_class_global" && "mixscape_class_global_p_ko" %in% post_candidates)
      "mixscape_class_global_p_ko" else if (length(post_candidates)) post_candidates[1] else NA_character_
    if (!is.na(pcol)) so@meta.data$mixscape_p_active <- so@meta.data[[pcol]]
    
    ## -- 7) Populate KO selector
    mc <- so@meta.data$mixscape_class_active
    has_NT <- any(mc == "NT", na.rm = TRUE)
    ko_labels <- sort(unique(mc[!is.na(mc) & (mc == "KO" | grepl(" KO$", mc))]))
    try({
      if (length(ko_labels)) {
        updateSelectInput(session, "ko_label",
                          choices = ko_labels,
                          selected = if ("IFNGR2 KO" %in% ko_labels) "IFNGR2 KO" else ko_labels[1])
      } else {
        updateSelectInput(session, "ko_label", choices = character(0), selected = NULL)
      }
    }, silent = TRUE)
    
    ## -- 8) Downstream tables (ensure counts before DE, just in case)
    df <- as.data.frame(table(
      gene = so@meta.data$gene,
      replicate = so@meta.data$replicate,
      mixscape_class = so@meta.data$mixscape_class_active
    )) %>%
      dplyr::group_by(gene, replicate) %>%
      dplyr::mutate(percent = Freq / sum(Freq)) %>%
      dplyr::ungroup() %>%
      dplyr::rename(Var1 = mixscape_class, value = percent)
    vals$df3 <- df
    
    if (has_NT && length(ko_labels) > 0) {
      ko_pick <- ko_labels[1]
      so <- .ensure_counts_layer(so, "RNA")  # extra safety for DE
      so <- .ensure_counts_layer(so, "RNA")
      
      vals$koGenes <- tryCatch({
        FindMarkers(so, ident.1 = ko_pick, ident.2 = "NT",
                    group.by = "mixscape_class_active", assay = "RNA")
      }, error = function(e) data.frame(message = paste("DGE failed:", e$message)))
    } else {
      vals$koGenes <- data.frame(message = "No KO class detected. Try increasing iter.num or lowering min.de.genes, or check guide quality.")
    }
    
    vals$seurat <- so
    removeModal()
    showNotification(paste0("Mixscape complete — using class column: ", class_col), type="message")
  })
  
  

  library(plotly)
  output$umapPlots <- renderPlot({
    req(vals$p1, vals$p2, vals$p3)
    (vals$p1 / vals$p2 + plot_layout(guides = 'auto')) | vals$p3
  })
  

  output$koBarPlot <- renderPlotly({
    req(vals$df3)
    df <- vals$df3
    req(nrow(df) > 0)
    
    # Make sure Var1 is a factor and build a palette that fits exactly
    lvl <- sort(unique(df$Var1))
    pal <- setNames(scales::hue_pal()(length(lvl)), lvl)
    df$Var1 <- factor(df$Var1, levels = lvl)
    
    p <- ggplot(df, aes(x = replicate, y = value * 100, fill = Var1)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      theme_classic() +
      scale_fill_manual(values = pal) +      # << dynamic, no more “insufficient values”
      ylab("% of cells") +
      xlab("Replicate") +
      facet_wrap(vars(gene), ncol = 5, scales = "free_y") +
      labs(fill = "mixscape class") +
      theme(
        axis.text.x = element_text(size = 14, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title  = element_text(size = 16),
        strip.text  = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12)
      )
    
    ggplotly(p)
  })

  output$perturbPlot <- renderPlotly({
    req(vals$seurat)
    so <- vals$seurat
    
    validate(need("mixscape_class_active" %in% colnames(so@meta.data),
                  "Perturbation: run Mixscape first (no mixscape_class_active)."))
    
    # Pick a posterior column that exists in metadata
    pcol <- if ("mixscape_p_active" %in% colnames(so@meta.data)) {
      "mixscape_p_active"
    } else {
      cand <- intersect(c("mixscape_class_p_ko", "mixscape_class_global_p_ko"),
                        colnames(so@meta.data))
      if (length(cand)) cand[1] else NA_character_
    }
    validate(need(!is.na(pcol), "Perturbation: no posterior probability column found."))
    
    # KO label selection (accepts "KO" or "<GENE> KO"); fallback to any KO present
    cls <- as.character(so@meta.data$mixscape_class_active)
    ko_choices <- sort(unique(cls[!is.na(cls) & (cls == "KO" | grepl(" KO$", cls))]))
    validate(need(length(ko_choices) > 0, "Perturbation: no KO classes detected."))
    
    ko <- input$ko_label
    if (is.null(ko) || !(ko %in% ko_choices)) ko <- ko_choices[1]
    
    Idents(so) <- factor(so@meta.data$mixscape_class_active)
    keep_id <- intersect(levels(Idents(so)), c("NT", ko))
    validate(need(length(keep_id) == 2, "Perturbation: need both NT and selected KO."))
    
    p <- VlnPlot(
      so, features = pcol, idents = keep_id, group.by = "mixscape_class_active", pt.size = 0
    ) + ggtitle(paste0("Posterior KO vs NT (", ko, ")")) + NoLegend()
    
    ggplotly(p)
  })
  
   
  output$vlnPlot <- renderPlotly({
    req(vals$seurat)
    so <- vals$seurat
    
    # Ensure RNA counts/data present; some geom/stat helpers poke the assay
    so <- .ensure_counts_layer(so, "RNA")
    
    # pick the best available posterior column
    pcol <- if ("mixscape_p_active" %in% colnames(so@meta.data)) {
      "mixscape_p_active"
    } else {
      cand <- intersect(c("mixscape_class_p_ko","mixscape_class_global_p_ko"),
                        colnames(so@meta.data))
      if (length(cand)) cand[1] else NULL
    }
    req(!is.null(pcol))
    
    gby <- if ("mixscape_class_active" %in% colnames(so@meta.data))
      "mixscape_class_active" else "gene"
    
    p <- VlnPlot(
      so, features = pcol, group.by = gby, pt.size = 0
    ) +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
            axis.text   = element_text(size = 14),
            plot.title  = element_text(size = 18)) +
      ggtitle("Mixscape posterior") + NoLegend()
    
    ggplotly(p)
  })
  
  output$perturbPlot <- renderPlotly({
    req(vals$seurat)
    so <- vals$seurat
    
    validate(need("mixscape_class_active" %in% colnames(so@meta.data),
                  "Perturbation: run Mixscape first (no mixscape_class_active)."))
    
    # Pick a posterior column that exists in metadata
    pcol <- if ("mixscape_p_active" %in% colnames(so@meta.data)) {
      "mixscape_p_active"
    } else {
      cand <- intersect(c("mixscape_class_p_ko", "mixscape_class_global_p_ko"),
                        colnames(so@meta.data))
      if (length(cand)) cand[1] else NA_character_
    }
    validate(need(!is.na(pcol), "Perturbation: no posterior probability column found."))
    
    # KO label selection (accepts "KO" or "<GENE> KO"); fallback to any KO present
    cls <- as.character(so@meta.data$mixscape_class_active)
    ko_choices <- sort(unique(cls[!is.na(cls) & (cls == "KO" | grepl(" KO$", cls))]))
    validate(need(length(ko_choices) > 0, "Perturbation: no KO classes detected."))
    
    ko <- input$ko_label
    if (is.null(ko) || !(ko %in% ko_choices)) ko <- ko_choices[1]
    
    Idents(so) <- factor(so@meta.data$mixscape_class_active)
    keep_id <- intersect(levels(Idents(so)), c("NT", ko))
    validate(need(length(keep_id) == 2, "Perturbation: need both NT and selected KO."))
    
    p <- VlnPlot(
      so, features = pcol, idents = keep_id, group.by = "mixscape_class_active", pt.size = 0
    ) + ggtitle(paste0("Posterior KO vs NT (", ko, ")")) + NoLegend()
    
    ggplotly(p)
  })
  
  
  output$mixscapeHeatmap <- renderPlotly({
    req(vals$seurat)
    so <- vals$seurat
    
    validate(need("mixscape_class_active" %in% colnames(so@meta.data),
                  "Heatmap: run Mixscape first (no mixscape_class_active)."))
    
    cls <- as.character(so@meta.data$mixscape_class_active)
    validate(need(any(cls == "NT", na.rm = TRUE),
                  "Heatmap: no NT cells present in mixscape_class_active."))
    
    # KO label selection
    ko_choices <- sort(unique(cls[!is.na(cls) & (cls == "KO" | grepl(" KO$", cls))]))
    validate(need(length(ko_choices) > 0, "Heatmap: no KO classes detected."))
    ko <- input$ko_label
    if (is.null(ko) || !(ko %in% ko_choices)) ko <- ko_choices[1]
    
    # Use existing DE if available; otherwise compute a small DE table quickly
    de_tab <- NULL
    if (!is.null(vals$koGenes) && is.data.frame(vals$koGenes) && !"message" %in% names(vals$koGenes)) {
      de_tab <- vals$koGenes
    } else {
      de_tab <- suppressWarnings(tryCatch(
        FindMarkers(so, ident.1 = ko, ident.2 = "NT",
                    group.by = "mixscape_class_active", assay = "RNA"),
        error = function(e) NULL))
    }
    validate(need(!is.null(de_tab) && nrow(de_tab) > 0,
                  "Heatmap: no DE genes available (KO vs NT)."))
    
    # Take top 30 by adjusted p-value (then p-value) for signal-dense heatmap
    ord <- order(with(de_tab, if ("p_val_adj" %in% names(de_tab)) p_val_adj else p_val),
                 de_tab$p_val %||% Inf, na.last = NA)
    top_genes <- rownames(de_tab)[head(ord, 30)]
    validate(need(length(top_genes) >= 2, "Heatmap: not enough significant genes."))
    
    DefaultAssay(so) <- "RNA"
    p <- suppressWarnings(tryCatch(
      DoHeatmap(so, features = top_genes, group.by = "mixscape_class_active") + NoLegend(),
      error = function(e) NULL))
    validate(need(!is.null(p), "Heatmap: DoHeatmap failed."))
    
    ggplotly(p)
  })
  
    
  output$koGenes <- DT::renderDataTable({
    req(vals$seurat)
    so  <- vals$seurat
    cls <- so@meta.data$mixscape_class_active
    if (is.null(cls)) {
      return(DT::datatable(data.frame(message = "Mixscape class column not found")))
    }
    
    # any KO labels present?
    ko_choices <- sort(unique(cls[!is.na(cls) & (cls == "KO" | grepl(" KO$", cls))]))
    if (!"NT" %in% cls || length(ko_choices) == 0) {
      return(DT::datatable(data.frame(
        message = "No KO class detected. Try increasing iter.num or lowering min.de.genes, or check guide quality."
      )))
    }
    
    # use selection if available; otherwise first KO present
    ko <- input$ko_label
    if (is.null(ko) || !(ko %in% ko_choices)) ko <- ko_choices[1]
    
    de_tab <- tryCatch({
      FindMarkers(
        so, ident.1 = ko, ident.2 = "NT",
        group.by = "mixscape_class_active", assay = "RNA"
      )
    }, error = function(e) data.frame(message = paste("DGE failed:", e$message)))
    
    DT::datatable(de_tab, options = list(pageLength = 10))
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
  
  # ---- ALL ----
  observeEvent(input$enrich_all, {
    req(vals$koGenes)
    showNotification("Running enrichment (ALL)...", id = "enrich_all_busy", type = "message", duration = NULL)
    on.exit(removeNotification("enrich_all_busy"), add = TRUE)
    
    genes <- .pick_de_genes(vals$koGenes, "all")
    if (!length(genes)) {
      showNotification("No DE genes (ALL) pass FDR < 0.05.", type = "warning")
      return(NULL)
    }
    
    tbl <- .do_enrichment(genes, input$species)
    vals$enrich_all <- tbl
    
    if (is.null(tbl)) {
      showNotification("No terms returned (ALL).", type = "warning")
    } else {
      showNotification("Enrichment (ALL) complete.", type = "message")
    }
  })
  
  # ---- UP ----
  observeEvent(input$enrich_up, {
    req(vals$koGenes)
    showNotification("Running enrichment (UP)...", id = "enrich_up_busy", type = "message", duration = NULL)
    on.exit(removeNotification("enrich_up_busy"), add = TRUE)
    
    genes <- .pick_de_genes(vals$koGenes, "up")
    if (!length(genes)) {
      showNotification("No UP genes pass FDR < 0.05.", type = "warning")
      return(NULL)
    }
    
    tbl <- .do_enrichment(genes, input$species)
    vals$enrich_up <- tbl
    
    if (is.null(tbl)) {
      showNotification("No terms returned (UP).", type = "warning")
    } else {
      showNotification("Enrichment (UP) complete.", type = "message")
    }
  })
  
  # ---- DOWN ----
  observeEvent(input$enrich_down, {
    req(vals$koGenes)
    showNotification("Running enrichment (DOWN)...", id = "enrich_down_busy", type = "message", duration = NULL)
    on.exit(removeNotification("enrich_down_busy"), add = TRUE)
    
    genes <- .pick_de_genes(vals$koGenes, "down")
    if (!length(genes)) {
      showNotification("No DOWN genes pass FDR < 0.05.", type = "warning")
      return(NULL)
    }
    
    tbl <- .do_enrichment(genes, input$species)
    vals$enrich_down <- tbl
    
    if (is.null(tbl)) {
      showNotification("No terms returned (DOWN).", type = "warning")
    } else {
      showNotification("Enrichment (DOWN) complete.", type = "message")
    }
  })
  
  # Tables
  output$enrich_all_tbl  <- DT::renderDataTable({ req(vals$enrich_all);  DT::datatable(vals$enrich_all,  options = list(pageLength = 12)) })
  output$enrich_up_tbl   <- DT::renderDataTable({ req(vals$enrich_up);   DT::datatable(vals$enrich_up,   options = list(pageLength = 12)) })
  output$enrich_down_tbl <- DT::renderDataTable({ req(vals$enrich_down); DT::datatable(vals$enrich_down, options = list(pageLength = 12)) })
  
  # Barplots
  output$enrich_all_bar  <- plotly::renderPlotly({ req(vals$enrich_all);  .enrich_bar(vals$enrich_all,  "Top Enriched (ALL)") })
  output$enrich_up_bar   <- plotly::renderPlotly({ req(vals$enrich_up);   .enrich_bar(vals$enrich_up,   "Top Enriched (UP)") })
  output$enrich_down_bar <- plotly::renderPlotly({ req(vals$enrich_down); .enrich_bar(vals$enrich_down, "Top Enriched (DOWN)") })
  
  # ---- Enrichment downloads ----
  output$dl_enrich_all <- downloadHandler(
    filename = function() sprintf("enrichment_all_%s.csv", Sys.Date()),
    content  = function(file) {
      if (is.null(vals$enrich_all) || nrow(vals$enrich_all) == 0) {
        write.csv(data.frame(message = "No enrichment results for ALL"), file, row.names = FALSE)
      } else {
        write.csv(vals$enrich_all, file, row.names = FALSE)
      }
    }
  )
  
  output$dl_enrich_up <- downloadHandler(
    filename = function() sprintf("enrichment_up_%s.csv", Sys.Date()),
    content  = function(file) {
      if (is.null(vals$enrich_up) || nrow(vals$enrich_up) == 0) {
        write.csv(data.frame(message = "No enrichment results for UP"), file, row.names = FALSE)
      } else {
        write.csv(vals$enrich_up, file, row.names = FALSE)
      }
    }
  )
  
  output$dl_enrich_down <- downloadHandler(
    filename = function() sprintf("enrichment_down_%s.csv", Sys.Date()),
    content  = function(file) {
      if (is.null(vals$enrich_down) || nrow(vals$enrich_down) == 0) {
        write.csv(data.frame(message = "No enrichment results for DOWN"), file, row.names = FALSE)
      } else {
        write.csv(vals$enrich_down, file, row.names = FALSE)
      }
    }
  )
  
  
  
}

  

shinyApp(ui, server)
