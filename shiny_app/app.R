library(shiny)
library(ggplot2)
library(DESeq2)
library(dplyr)
library(reshape2)
library(bslib)
library(DT)
library(tibble)
library(RColorBrewer)

# Define tissues
tissue_list <- c("brain", "heart", "muscle", "spleen")

# Load data
norm_counts <- readRDS("data/normalized_counts.rds")
coldata <- readRDS("data/colData.rds")

# Load full dataset DE results
res_age_df <- readRDS("data/res_age.rds") %>% as.data.frame() %>% tibble::rownames_to_column("gene")
res_sex_df <- readRDS("data/res_sex.rds") %>% as.data.frame() %>% tibble::rownames_to_column("gene")
res_age_female_df <- readRDS("data/res_age_female.rds") %>% as.data.frame() %>% tibble::rownames_to_column("gene")
res_age_male_df <- readRDS("data/res_age_male.rds") %>% as.data.frame() %>% tibble::rownames_to_column("gene")

# UI
ui <- fluidPage(
  theme = bs_theme(bootswatch = "flatly", version = 5),
  titlePanel("Killifish RNA-seq: Gene Expression Explorer"),
  sidebarLayout(
    sidebarPanel(
      helpText("Select gene(s) below to explore DE results and expression patterns."),
      selectizeInput("genes", "Select gene(s):",
                  choices = NULL,
                  #selected = c("mtor", "hdac1"),
                  multiple = TRUE,
                  options = list(placeholder = "Type to search genes...")),
      helpText(
        "Need help identifying your gene of interest?",
        a("View genome annotation at NCBI",
          href = "https://www.ncbi.nlm.nih.gov/datasets/gene/GCF_043380555.1/",
          target = "_blank"))
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("DE Table",
                 downloadButton("download_de", "Download CSV"),
                 br(),
                 HTML("<p><i>Note:</i> <span style='color:red;'>Red</span> log2FC values indicate upregulation, 
       <span style='color:blue;'>Blue</span> indicate downregulation. 
       Genes with <code>padj &lt; 0.05</code> are marked with <b>*</b> in the Significance column.</p>"),
                 DTOutput("de_stats")
        ),
        tabPanel("Expression by Age",
                 downloadButton("download_age_plot", "Download Plot"),
                 plotOutput("age_plot", height = "700px")
        ),
        tabPanel("Expression by Sex",
                 downloadButton("download_sex_plot", "Download Plot"),
                 plotOutput("sex_plot", height = "700px")
        ),
        tabPanel("Expression by Sex & Age",
                 downloadButton("download_combined_plot", "Download Plot"),
                 plotOutput("combined_plot", height = "700px")
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  # Server-side gene list for selectize
  observe({
    updateSelectizeInput(
      session,
      "genes",
      choices = rownames(norm_counts),
      server = TRUE
    )
  })
  
  output$theme_ui <- renderUI({
    bs_themer()
  })
  
  # Show DE results
  output$de_stats <- renderPrint({
    req(input$genes)
    
    gene_info_age <- dplyr::filter(res_age_df, gene %in% input$genes)
    gene_info_sex <- dplyr::filter(res_sex_df, gene %in% input$genes)
    
    for (g in input$genes) {
      cat("Gene:", g, "\n")
      age_row <- gene_info_age %>% filter(gene == g)
      sex_row <- gene_info_sex %>% filter(gene == g)
      
      cat("  log2 Fold Change (Age):", round(age_row$log2FoldChange, 2), "\n")
      cat("  Adjusted p-value (Age):", signif(age_row$padj, 3), "\n")
      cat("  log2 Fold Change (Sex):", round(sex_row$log2FoldChange, 2), "\n")
      cat("  Adjusted p-value (Sex):", signif(sex_row$padj, 3), "\n\n")
    }
  })

  # Load all DE results (full + tissue-specific)
  all_de_results <- reactive({
    req(input$genes)
    results_combined <- list()
    
    for (t in tissue_list) {
      age_path <- file.path("data", paste0("res_age_", t, ".rds"))
      sex_path <- file.path("data", paste0("res_sex_", t, ".rds"))
      
      if (file.exists(age_path) && file.exists(sex_path)) {
        age_df <- readRDS(age_path) %>% as.data.frame() %>%
          tibble::rownames_to_column("gene") %>%
          filter(gene %in% input$genes) %>%
          mutate(Tissue = t, Contrast = "Age",
                 log2FC = log2FoldChange, padj = padj) %>%
          dplyr::select(Tissue, Contrast, gene, log2FC, padj)
        
        sex_df <- readRDS(sex_path) %>% as.data.frame() %>%
          tibble::rownames_to_column("gene") %>%
          filter(gene %in% input$genes) %>%
          mutate(Tissue = t, Contrast = "Sex",
                 log2FC = log2FoldChange, padj = padj) %>%
          dplyr::select(Tissue, Contrast, gene, log2FC, padj)
        
        results_combined[[t]] <- bind_rows(age_df, sex_df)
      }
    }
    
    
    # Full dataset DE results
    full_age <- res_age_df %>%
      filter(gene %in% input$genes) %>%
      mutate(Tissue = "FullData", Contrast = "Age",
             log2FC = log2FoldChange, padj = padj) %>%
      dplyr::select(Tissue, Contrast, gene, log2FC, padj)
    
    full_sex <- res_sex_df %>%
      filter(gene %in% input$genes) %>%
      mutate(Tissue = "FullData", Contrast = "Sex",
             log2FC = log2FoldChange, padj = padj) %>%
      dplyr::select(Tissue, Contrast, gene, log2FC, padj)
    
    female_age_df <- res_age_female_df %>%
      filter(gene %in% input$genes) %>%
      mutate(Tissue = "FullData", Contrast = "Age (females)", log2FC = log2FoldChange) %>%
      select(Tissue, Contrast, gene, log2FC, padj)
    
    male_age_df <- res_age_male_df %>%
      filter(gene %in% input$genes) %>%
      mutate(Tissue = "FullData", Contrast = "Age (males)", log2FC = log2FoldChange) %>%
      select(Tissue, Contrast, gene, log2FC, padj)
    
    # Combine everything
    final_df <- bind_rows(full_age, full_sex, female_age_df, male_age_df, bind_rows(results_combined)) %>%
      arrange(gene, Tissue, Contrast) %>%
      mutate(Significance = ifelse(padj < 0.05, "*", ""))
    
    return(final_df)
  })
  
  
  # DE summary table
  output$de_stats <- DT::renderDataTable({
    df <- all_de_results()
    
    DT::datatable(
      df,
      rownames = FALSE,
      escape = FALSE,
      options = list(
        pageLength = 20,
        autoWidth = TRUE,
        order = list(list(2, 'asc'))
      ),
      caption = htmltools::tags$caption(
        style = 'caption-side: top; text-align: left; font-weight: bold;',
        "Differential Expression Summary (Full Dataset + Tissue-Specific)"
      )
    ) %>%
      DT::formatRound(columns = c("log2FC"), digits = 3) %>%
      DT::formatSignif(columns = "padj", digits = 3) %>%
      DT::formatStyle(
        "log2FC",
        color = DT::styleInterval(0, c("blue", "red")),
        fontWeight = "bold"
      ) %>%
      DT::formatStyle(
        "Tissue",
        backgroundColor = DT::styleEqual(
          unique(df$Tissue),
          RColorBrewer::brewer.pal(n = length(unique(df$Tissue)), "Pastel1")
        )
      )
  })
  
  # Download DE table
  output$download_de <- downloadHandler(
    filename = function() {
      paste0("DE_summary_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(all_de_results(), file, row.names = FALSE)
    }
  )
  
  # Expression plot by age
  output$age_plot <- renderPlot({
    req(input$genes)
    
    sub_counts <- norm_counts[rownames(norm_counts) %in% input$genes, , drop = FALSE]
    sub_counts$Gene <- rownames(sub_counts)
    long_counts <- reshape2::melt(sub_counts, id.vars = "Gene", variable.name = "Sample", value.name = "Expression")
    
    # Merge with coldata
    df_coldata <- coldata %>% mutate(Sample = rownames(coldata))
    merged_df <- merge(long_counts, df_coldata, by = "Sample")
    
    # Plot
    ggplot(merged_df, aes(x = Age, y = Expression, fill = Gene)) +
      geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.8) +
      geom_jitter(aes(fill = Gene), color = "black", shape = 21, position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.8), size = 1.2, alpha = 0.6) +
      facet_wrap(~Tissue, scales = "free_y") +
      scale_fill_brewer(palette = "Spectral") + 
      scale_color_brewer(palette = "Spectral") +
      labs(title = "Expression by Age", x = "Age", y = "Normalized Expression") + 
      theme_minimal() +
      theme(
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 13)  # for facet labels
      )
  })

  output$sex_plot <- renderPlot({
    req(input$genes)
    
    sub_counts <- norm_counts[rownames(norm_counts) %in% input$genes, , drop = FALSE]
    sub_counts$Gene <- rownames(sub_counts)
    
    long_counts <- reshape2::melt(sub_counts, id.vars = "Gene", variable.name = "Sample", value.name = "Expression")
    df_coldata <- coldata %>% mutate(Sample = rownames(coldata))
    merged_df <- merge(long_counts, df_coldata, by = "Sample")
    
    ggplot(merged_df, aes(x = Sex, y = Expression, fill = Gene)) +
      geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.8) +
      geom_jitter(aes(fill = Gene), color = "black", shape = 21, position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.8), size = 1.2, alpha = 0.6) +
      facet_wrap(~Tissue, scales = "free_y") +
      scale_fill_brewer(palette = "Spectral") + 
      scale_color_brewer(palette = "Spectral") +
      labs(title = "Expression by Sex", x = "Sex", y = "Normalized Expression") + 
      theme_minimal() +
      theme(
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 13)  # for facet labels
      )
  })
  
  output$combined_plot <- renderPlot({
    req(input$genes)
    
    # Extract and melt normalized counts
    norm_counts <- norm_counts[rownames(norm_counts) %in% input$genes, , drop = FALSE]
    norm_counts$Gene <- rownames(norm_counts)
    
    long_counts <- reshape2::melt(norm_counts, id.vars = "Gene",
                                  variable.name = "Sample", value.name = "Expression")
    
    df_coldata <- coldata %>% mutate(Sample = rownames(coldata))
    long_counts <- merge(long_counts, df_coldata, by = "Sample")
    
    # Create Label: Sex_Age (e.g. female_6w) and factor to control order
    long_counts$Label <- paste(long_counts$Sex, long_counts$Age, sep = "_")
    long_counts$Label <- factor(long_counts$Label,
                                levels = c("female_6w", "female_16w", "male_6w", "male_16w"))
    
    ggplot(long_counts, aes(x = Label, y = Expression, fill = Gene)) +
      geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, alpha = 0.8) +
      geom_jitter(aes(fill = Gene), color = "black", shape = 21,
                  position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
                  size = 1.5, alpha = 0.7) +
      facet_wrap(~Tissue, scales = "free_y") +
      scale_fill_brewer(palette = "Spectral") +
      labs(title = "Expression by Sex & Age", x = "Sex + Age", y = "Normalized Expression") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(size = 12, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 13)
      )
  })
  
  # Save Age Plot as PDF
  output$download_age_plot <- downloadHandler(
    filename = function() { paste0("expression_age_", Sys.Date(), ".pdf") },
    content = function(file) {
      pdf(file, width = 10, height = 7)
      selected_genes <- input$genes
      req(selected_genes)
      
      sub_counts <- norm_counts[rownames(norm_counts) %in% selected_genes, , drop = FALSE]
      sub_counts$Gene <- rownames(sub_counts)
      long_counts <- melt(sub_counts, id.vars = "Gene", variable.name = "Sample", value.name = "Expression")
      df_coldata <- coldata %>% mutate(Sample = rownames(coldata))
      merged_df <- merge(long_counts, df_coldata, by = "Sample")
      
      p <- ggplot(merged_df, aes(x = Age, y = Expression, fill = Gene)) +
        geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, alpha = 0.8) +
        geom_jitter(aes(fill = Gene), color = "black", shape = 21,
                    position = position_jitterdodge(0.05, 0.8), size = 1.2, alpha = 0.6) +
        facet_wrap(~Tissue, scales = "free_y") +
        scale_fill_brewer(palette = "Spectral") +
        labs(title = "Expression by Age", x = "Age", y = "Normalized Expression") +
        theme_minimal()
      print(p)
      dev.off()
    }
  )
  
  # Save Sex Plot as PDF
  output$download_sex_plot <- downloadHandler(
    filename = function() { paste0("expression_sex_", Sys.Date(), ".pdf") },
    content = function(file) {
      pdf(file, width = 10, height = 7)
      selected_genes <- input$genes
      req(selected_genes)
      
      sub_counts <- norm_counts[rownames(norm_counts) %in% selected_genes, , drop = FALSE]
      sub_counts$Gene <- rownames(sub_counts)
      long_counts <- melt(sub_counts, id.vars = "Gene", variable.name = "Sample", value.name = "Expression")
      df_coldata <- coldata %>% mutate(Sample = rownames(coldata))
      merged_df <- merge(long_counts, df_coldata, by = "Sample")
      
      p <- ggplot(merged_df, aes(x = Sex, y = Expression, fill = Gene)) +
        geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, alpha = 0.8) +
        geom_jitter(aes(fill = Gene), color = "black", shape = 21,
                    position = position_jitterdodge(0.05, 0.8), size = 1.2, alpha = 0.6) +
        facet_wrap(~Tissue, scales = "free_y") +
        scale_fill_brewer(palette = "Spectral") +
        labs(title = "Expression by Sex", x = "Sex", y = "Normalized Expression") +
        theme_minimal()
      print(p)
      dev.off()
    }
  )
  # Save Combined Plot (Sex + Age) as PDF
  output$download_combined_plot <- downloadHandler(
    filename = function() { paste0("expression_sex_age_", Sys.Date(), ".pdf") },
    content = function(file) {
      pdf(file, width = 10, height = 7)
      selected_genes <- input$genes
      req(selected_genes)
      
      sub_counts <- norm_counts[rownames(norm_counts) %in% selected_genes, , drop = FALSE]
      sub_counts$Gene <- rownames(sub_counts)
      long_counts <- reshape2::melt(sub_counts, id.vars = "Gene",
                                    variable.name = "Sample", value.name = "Expression")
      
      df_coldata <- coldata %>% mutate(Sample = rownames(coldata))
      merged_df <- merge(long_counts, df_coldata, by = "Sample")
      
      merged_df$Label <- paste(merged_df$Sex, merged_df$Age, sep = "_")
      merged_df$Label <- factor(merged_df$Label,
                                levels = c("female_6w", "female_16w", "male_6w", "male_16w"))
      
      p <- ggplot(merged_df, aes(x = Label, y = Expression, fill = Gene)) +
        geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, alpha = 0.8) +
        geom_jitter(aes(fill = Gene), color = "black", shape = 21,
                    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
                    size = 1.5, alpha = 0.7) +
        facet_wrap(~Tissue, scales = "free_y") +
        scale_fill_brewer(palette = "Spectral") +
        labs(title = "Expression by Sex & Age", x = "Sex + Age", y = "Normalized Expression") +
        theme_minimal()
      
      print(p)
      dev.off()
    }
  )
}

# Launch the app
shinyApp(ui = ui, server = server)