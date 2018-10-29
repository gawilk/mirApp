# Interactive Shiny App!
# Plots miRNA-specific manhattan plots and SNPs modulating miRNA-gene interactions
# Plots ellipse clouds & linear regression for each genotype
# Obfuscates original data so sample genotypes cannot be re-engineered
# Genotype stats are in matrix form for more efficient data storage and memory loading

# change to repo first if necessary, eg: setwd("path/to/mirApp")
load("data/SNPdb.RData")
library(shiny)
library(ggplot2)

#==============================================================================
# UX
#==============================================================================

ui <- shinyUI(fluidPage(
  titlePanel("miRNA manhattan plotter"),
  sidebarLayout(sidebarPanel(
    helpText("Create miRNA-specific manhattan 
             plots for a tumor type."), 
    selectInput("organ", 
                label = "Choose a cancer type to display",
                choices = list("Breast", "Liver",
                               "Lung", "Prostate"),
                selected = "Breast"),
    textInput("mir", label = "Choose a miRNA", 
    	value = "type in miRNA...(hsa-mir-XXX or hsa-let-XXX)"), 
    helpText("Data may take a while to load...")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Plot", 
          textOutput("text1"),
          textOutput("text2"),
          textOutput("text3"),
          plotOutput("plot1", click = "plot_click"), 
          verbatimTextOutput("info"), 
          plotOutput("plot2")
        ),
        tabPanel("About", 
          h4("mirApp Info"), 
          p("This interactive Shiny app plots SNPs modulating miRNA-gene 
            interactions throughout the genome in tumors. When a user inputs 
            a cancer type and miRNA of choice, the app plots:"),
          p("a) all miRNA-gene 
            interactions modulated by SNPs in the genome in a miRNA-specific 
            manhattan plot. The miRNA-specific manhattan plot appears in 
            the upper panel of the plotting output. All observations are
            FDR-corrected, colored by chromosome, and plotted by genomic location
            of the SNP."),
          p("b) a specific miRNA-gene interaction modulated 
            by a SNP when the user hovers over and clicks on an observation. This figure 
            appears in the lower panel of the plotting output. Each genotype group is plotted
            as an ellipse with a superimposed linear regression line. Genotype sample sizes 
            are shown in the legend."),
          h4("Ellipse Information"), 
          p("Genotype groups are visualized as ellipses, rather than individual data points, 
            to protect the privacy of cancer patients in the study. The ellipses encapsulate 
            the entire portion of the genotype group. In this way, no single tumor sample can 
            be reverse-engineered and reidentified, since revealing genotype information
            can compromise someone's privacy."),
          h4("Notes"), 
          p("***miRNA names must be consistent with TCGA convention. That is, they are
            either of the form \"hsa-mir-XXX\" or \"hsa-let-XXX\", where XXX is replaced with
            whichever miRNA suffix of your choosing."))
      )
    )
)
))

#==============================================================================
# SERVER CODE
#==============================================================================

server <- function(input, output) {
  
  output$text1 <- renderText({ 
    paste("You have selected", input$mir, "in", input$organ, "cancer")
  })
  
  selectedTrios <- reactive({
    switch(input$organ, 
           "Breast" = {load("data/BRCAregQTLPCA.RData", envir = .GlobalEnv)
             BRCAregQTLPCA},
           "Liver" = {load("data/LIHCregQTLPCA.RData", envir = .GlobalEnv)
             LIHCregQTLPCA},
           "Lung" = {load("data/LUSCregQTLPCA.RData", envir = .GlobalEnv)
             LUSCregQTLPCA},
           "Prostate" = {load("data/PRADregQTLPCA.RData", envir = .GlobalEnv)
             PRADregQTLPCA})
  })
  
  selectedStats <- reactive({
    switch(input$organ, 
           "Breast" = {load("data/BRCAstatsMat.RData", envir = .GlobalEnv)
             BRCAstatsMat},
           "Liver" = {load("data/LIHCstatsMat.RData", envir = .GlobalEnv)
             LIHCstatsMat},
           "Lung" = {load("data/LUSCstatsMat.RData", envir = .GlobalEnv)
             LUSCstatsMat},
           "Prostate" = {load("data/PRADstatsMat.RData", envir = .GlobalEnv)
             PRADstatsMat})
  })
  
  getCHRcoords <- function(df) {
    # recomputes observations in absolute coordinates in genome 
    # 
    # Args: 
    #   df: df, look for position (pos) of mutation
    # Returns: 
    #   absolute coords of chromosomes
    chrs <- paste0("chr", c(as.character(seq(1, 22)), "X", "Y"))
    DFsplit <- split(df, as.factor(df$chrom))
    sizes <- sapply(chrs, function(X) max(DFsplit[[X]]$pos))
    sizes <- sizes[is.finite(sizes)]
    sums <- cumsum(as.numeric(sizes))
    starts <- c(0, sums[(1:length(sums)) - 1])
    chrmap <- data.frame(chr = paste0("chr", c(as.character(seq(1, 22)), "X")), 
                         start = starts)
    rownames(chrmap) <- chrmap$chr
    return(chrmap)
  }
  
  chrData <- reactive({
    getCHRcoords(selectedTrios())
  })
  
  subMIR <- reactive({
    inmir <- gsub("miR", "mir", input$mir)
    DF <- subset(selectedTrios(), miRNA == inmir)
    DF$logp <- -log10(DF$pfdrONE)
    DF
  })
  
  subTrio <- reactive({
    subset(subMIR(), logp > 0.25)
  })
  
  output$text2 <- renderText({ 
    paste0("Plotting ", nrow(subTrio()), " out of ", nrow(subMIR()), " trios ",
           "(those above ", 0.25, " on a -log10 scale)")
  })
  
  output$text3 <- renderText({ 
    paste("Click on an individual observation to view SNP effect on miRNA-gene interaction")
  })
  
  output$plot1 <- renderPlot({
    gg <- ggplot(subTrio())
    gg <- gg + geom_point(aes(x = abspos, y = logp, color = chrom), 
                          size = 1.5, alpha = 0.5)
    gg <- gg + scale_x_continuous(labels = rownames(chrData()), breaks = chrData()$start)
    gg <- gg + theme(panel.background = element_rect(fill = "white", colour = "white"), legend.position = "none", 
                     axis.text.x = element_text(size = 16, angle = 90, hjust = 1), 
                     axis.text.y = element_text(size = 16), 
                     axis.title = element_text(size = 16))
    gg <- gg + xlab("") + ylab(expression(-log[10](p[FDR])))
    gg
  })
  
  clickDF <- reactive({
    nearPoints(subTrio(), input$plot_click, 
               threshold = 20, maxpoints = 1,
               addDist = FALSE)
  })
  
  output$info <- renderPrint({
    clickDF <- clickDF()[, c("miRNA", "gene", "snp", "rsID", "chrom", "pfdrONE", "maf")]
    rownames(clickDF) <- NULL
    colnames(clickDF) <- c("miRNA", "gene", "snp", "rsID", "chrom", "pFDR", "MAF")
    clickDF <- as.data.frame(clickDF, stringsAsFactors=FALSE)
    clickDF
  })
  
  reshapeIntoList <- function(statsVector) {
    #converts trio statistics vector into list
    #vector includes ellipse stats and genotype lm regressions
    genoList <- split(statsVector, ceiling(seq_along(statsVector)/10))
    names(genoList) <- c("0","1","2")
    lapply(genoList, function(X) genoStats(X))
  }
  
  genoStats <- function(vect) {
    # called into reshapeIntoList
    statsList <- list()
    cent <- vect[1:2]
    lab <- c("mir", "gene")
    names(cent) <- lab
    statsList[["center"]] <- cent
    statsList[["radius"]] <- vect[3]
    statsList[["shape"]] <- matrix(vect[4:7], 2, 2, dimnames=list(lab, lab))
    statsList[["coeffs"]] <- vect[8:9]
    statsList[["count"]] <- vect[10]
    statsList
  }
  
  plotTrio <- function(miRNA, gene, mut, rsID, 
                             data, db, segments=51) {
    # plots ellipses and lms for each trio
    #
    # Args:
    #  miRNA: miRNA ID string 
    #  gene: gene name string
    #  mut: SNP Affy ID
    #  rsID: SNP rsID
    #  data: cancer genotype stats (huge list including ellipse shapes)
    #  db: SNP database to get SNP genotypes
    #  segments: num of points to calc ellipse
    require(ggplot2)
    angles <- (0:segments) * 2 * pi/segments
    unit.circle <- cbind(cos(angles), sin(angles))
    ggList <- vector("list", length=3)
    trio <- paste(miRNA, gene, mut, sep="+")
    statsList <- reshapeIntoList(data[trio, ])
    snpinfo <- db[mut, ]
    genotypes <- c(paste0(rep(snpinfo$allele_a, 2), collapse = ""), 
                   paste0(snpinfo$allele_a, snpinfo$allele_b), 
                   paste0(rep(snpinfo$allele_b, 2), collapse = ""))
    names(statsList) <- genotypes
    ellipses <- lapply(statsList, function(X) {
      chol_decomp <- chol(X$shape)
      ellipse <- t(X$center + X$radius * t(unit.circle %*% chol_decomp))
      ellipse <- as.data.frame(ellipse)
    })
    lms <- lapply(statsList, function(X) X$coeffs)
    genoCounts <- sapply(statsList, function(X) X$count)
    labs <- mapply(function(X, Y) {
      paste0(X, " (", Y, ")")
    }, names(statsList), genoCounts)
    names(labs) <- c("0","1","2")
    g <- ggplot() + geom_polygon(aes(x=mir, y=gene, fill="0"), 
                                 data=ellipses[[1]], alpha=0.2) + 
      geom_polygon(aes(x=mir, y=gene, fill="1"), 
                   data=ellipses[[2]], alpha=0.2) + 
      geom_polygon(aes(x=mir, y=gene, fill="2"), 
                   data=ellipses[[3]], alpha=0.2) + 
      scale_fill_manual(name=as.character(rsID), 
                        values = c("0"="red", "1"="blue", "2"="forestgreen"), 
                        labels = labs) +
      geom_abline(intercept=lms[[1]][1], slope=lms[[1]][2], col="red") + 
      geom_abline(intercept=lms[[2]][1], slope=lms[[2]][2], col="blue") + 
      geom_abline(intercept=lms[[3]][1], slope=lms[[3]][2], col="forestgreen") +
      xlab(miRNA) + ylab(gene) +
      theme(panel.background = element_rect(fill = "white", colour = "white"), 
            axis.title.x = element_text(size = 16), 
            axis.title.y = element_text(size = 16), 
            title = element_text(size = 14),
            axis.line.x = element_line(size = 0.5),
            axis.line.y = element_line(size = 0.5))
    return(g)
  }
  
  output$plot2 <- renderPlot({
    plotTrio(miRNA = clickDF()$miRNA, gene = clickDF()$gene, 
                   mut = clickDF()$snp, rsID = clickDF()$rsID, 
                   data = selectedStats(), db = SNPdb)
  })
  
}

#==============================================================================
# CALL
#==============================================================================

shinyApp(ui = ui, server = server)