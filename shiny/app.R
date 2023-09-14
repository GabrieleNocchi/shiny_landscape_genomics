library(shiny)
library(lfmm)
library(data.table)
library(tools)
library(scales)
library(LEA)
library(ggplot2)
library(caret)
library(gtools)

ui <- fluidPage(
 titlePanel("LFMM analysis"),
 sidebarLayout(
  position = "left",


  sidebarPanel(
   fileInput(
    inputId = "ped_file",
    label = ".ped file",
    multiple = FALSE,
    accept = NULL,
    width = NULL,
    buttonLabel = "Browse...",
    placeholder = "No file selected",
    capture = NULL
   ),


   fileInput(
    inputId = "map_file",
    label = ".map file",
    multiple = FALSE,
    accept = NULL,
    width = NULL,
    buttonLabel = "Browse...",
    placeholder = "No file selected",
    capture = NULL
   ),


   fileInput(
    inputId = "coo_file",
    label = "Samples coordinates (ID,LON,LAT) (.txt)",
    multiple = FALSE,
    accept = NULL,
    width = NULL,
    buttonLabel = "Browse...",
    placeholder = "No file selected",
    capture = NULL
   ),


   textInput(
    inputId = "var_cor_cutoff",
    label = "Correlation cutoff for environmental variables",
   ),

   textInput(
    inputId = "user_k",
    label = "Number of latent factors (K)",
   ),


   textInput(
    inputId = "FDR_cutoff",
    label = "FDR cutoff",
   )


  ),
  mainPanel(
   plotOutput(outputId = "snmf_plot"),
   downloadButton(outputId = "downloadSNMF", label = "Download SNMF Plot"),
   plotOutput(outputId = "correlation"),
   downloadButton(outputId = "downloadcorr", label = "Download Correlation Plot"),
   plotOutput(outputId = "manh_plot", height = 800),


   downloadButton(outputId = "downloadMANH", label = "Download Manhattan Plot"),
   downloadButton(outputId = "downloadTable", label = "Download Results") # Add this line

  )


 )
)





server <- function(input, output) {

 lfmm.res <- reactiveVal(NULL)
 obj_snmf <- reactiveVal(NULL)
 saved_manh <- reactiveVal(NULL)
 corr_matrix_saved <- reactiveVal(NULL)

 options(shiny.maxRequestSize = 3000*1024^2)

 output$downloadSNMF <- downloadHandler(
  filename = function() {
   paste("snmf_plot_", Sys.Date(), ".png", sep = "")
  },
  content = function(file) {
   png(file)
   plot(obj_snmf(), pch = 16, col = "blue")  # Access obj.snmf using obj.snmf()
   dev.off()
  },
  contentType = "image/png"
 )


 output$snmf_plot <- renderPlot({

  my_prefix <- basename(input$ped_file$name)
  my_prefix <- tools::file_path_sans_ext(my_prefix)

  ped <- input$ped_file$datapath
  ped <- data.table::fread(ped, sep = "\t")
  map <- input$map_file$datapath
  map <- data.table::fread(map, sep = "\t")

  # Build and execute the system command
  system(paste(
   "./plink.exe ",
   " --file ", my_prefix,
   " --recodeA --allow-extra-chr --out ", my_prefix,
   sep = "")
  )

  Y <- data.table::fread(paste(my_prefix, ".raw", sep = ""), na.strings = "NA", header = T)
  Y <- Y[, -c(1:6)]
  fwrite(Y, paste(my_prefix, ".lfmm", sep = ""), col.names = F, row.names = F, sep = "\t", na = "9")
  Y <- LEA::lfmm2geno(paste(my_prefix, ".lfmm",
  sep = ""))

  # sNMF (sparse nonnegative matrix factorization)
  obj.snmf <- LEA::snmf(Y, K = 1:10, entropy = T, ploidy = 2, project = "new")
  obj_snmf(obj.snmf)
  # let’inspect the values of the cross-entropy criterion for each K:
  plot(obj.snmf, pch = 16, col = "blue")
 })

 output$downloadMANH <- downloadHandler(

  filename = function() {
   paste("manh_plot_", Sys.Date(), ".png", sep = "")
  },
  content = function(file) {
   req(lfmm.res())
   for_plot <- length(unique(sort(lfmm.res()$ENV)))
   for_plot <- round(ceiling(for_plot/3))
   req(saved_manh())
   png(file)
   par(mfrow = c(for_plot, 4), oma = c(1, 1, 1, 1), mar = c(4, 6, 3, 1), mgp = c(3, 1, 0))
   for (i in 4:ncol(saved_manh())) {
    manh_i <- saved_manh()[, c(1:3, i)]
    SNPs <- lfmm.res()$SNP[lfmm.res()$ENV == colnames(manh_i)[4]]
    colnames(manh_i)[4] <- "P"
    qqman::manhattan(
     manh_i, genomewideline = F, suggestiveline = F,
     main = paste("LFMM - Env. variable ", colnames(saved_manh())[i], sep = ""),
     highlight = SNPs,
     col = c(alpha("gray70", 0.7), alpha("gray30", 0.7)),
     xlab = "Chromosome", cex.lab = 1.5, cex.axis = 1
    )

   }
   dev.off()
  },
  contentType = "image/png"
 )

 output$downloadcorr <- downloadHandler(
  filename = function() {
   paste("corr_plot_", Sys.Date(), ".png", sep = "")
  },
  content = function(file) {
   png(file)
   corrplot::corrplot(
    corr_matrix_saved(),
    order = "original",
    type = "upper", diag = T,
    tl.cex = 0.4,
    tl.col = c(rep("red", 11),
    rep("blue",8)),
    tl.srt=45, addCoef.col = "darkgray", addCoefasPercent = T
   )
   dev.off()
  },
  contentType = "image/png"
 )

 output$correlation <- renderPlot({
  req(corr_matrix_saved())

  corrplot::corrplot(
   corr_matrix_saved(),
   order = "original",
   type = "upper", diag = T,
   tl.cex = 0.4,
   tl.col = c(rep("red", 11),
   rep("blue",8)),
   tl.srt=45, addCoef.col = "darkgray", addCoefasPercent = T
  )

 })


 output$manh_plot <- renderPlot({

  my_prefix <- basename(input$ped_file$name)
  my_prefix <- tools::file_path_sans_ext(my_prefix)
  ped <- input$ped_file$datapath
  ped <- data.table::fread(ped, sep = "\t")
  map <- input$map_file$datapath
  map <- data.table::fread(map, sep = "\t")


  # Build and execute the system command
  system(paste(
   "./plink.exe ",
   " --file ", my_prefix,
   " --recodeA --allow-extra-chr --out ", my_prefix,
   sep = "")
  )

  Y <- data.table::fread(paste(my_prefix, ".lfmm",
  sep = ""), header = F)

  #     env_file <- input$env_file$datapath
  coo <- input$coo_file$datapath
  coo <- read.table(coo, sep = "\t", h = T, stringsAsFactors = F)

  sample_size <- length(coo$ID)
  # Environmental matrix
  env <- matrix(rep(NA, sample_size * 19), sample_size)
  colnames(env) <- paste("Bio", 1:19, sep = "")
  rownames(env) <- coo$ID


  # Let's fill the matrix
  for (i in 1:nrow(env)) {
   bio <- raster::getData(
    "worldclim", download = T,
    lon = coo$LON[i], lat = coo$LAT[i],
    var = "bio", res = 0.5
   )
   env[i, 1:19] <- raster::extract(
    x = bio, y = coo[i, c("LON", "LAT")]
   )
  }

  # Calculate correlation matrix
  corr_matrix <- cor(env)
  corr_matrix_saved(corr_matrix)
  # Set a correlation threshold (e.g., 0.8)
  threshold <- as.numeric(input$var_cor_cutoff)

  # Find highly correlated columns
  highly_correlated <- findCorrelation(corr_matrix, cutoff = threshold)

  # Remove highly correlated columns from the dataframe
  env <- env[,-highly_correlated]

  for_par <- ncol(env)
  for_par <- round(ceiling(for_par/4))

  #       X <- read.table(env_file, h = TRUE, stringsAsFactors = FALSE)
  X <- cbind(coo, env)
  X <- as.matrix(X[, c(4:ncol(X))])


  # You need to define 'K' before running this code
  K <- as.numeric(input$user_k)

  mod.lfmm <- lfmm::lfmm_ridge(Y = Y, X = X, K = K)
  pv <- lfmm::lfmm_test(Y = Y, X = X, lfmm = mod.lfmm, calibrate = "gif")
  pvalues <- pv$calibrated.pvalue

  my_qvalue <- function(x) {
   q <- qvalue::qvalue(x)
   q <- q$qvalues
   return(q)
  }

  qvalues <- apply(pvalues, 2, my_qvalue)


  qvalues <- as.data.frame(qvalues)
  mylabs <- unique(mixedsort(map$V1)) # we need this to preserve the order, otherwise it will order alphabetically
  map$V1 <- as.numeric(factor(map$V1, levels = mylabs))
  qvalues$SNP <- as.character(map$V2)
  qvalues$CHR <- map$V1
  qvalues$BP <- map$V4
  qcut <- as.numeric(input$FDR_cutoff)
  source("fn-landgen1.R", echo = F, keep.source = TRUE)
  lfmm_res <- lfmm_qvalcut(qvalues = qvalues, cutoff = qcut)
  lfmm.res(lfmm_res)


  pvalues <- as.data.frame(pvalues)
  pvalues$SNP <- as.character(map$V2)
  pvalues$CHR <- map$V1
  pvalues$BP <- map$V4

  a <- ncol(X) + 1
  b <- a + 2
  manh <- pvalues[, c(a:b, 1:ncol(X))]
  saved_manh(manh)
  #### PLOT SINGLE ####
  par(mfrow = c(for_par, 4), oma = c(1, 1, 1, 1), mar = c(4, 6, 3, 1), mgp = c(3, 1, 0))
  for (i in 4:ncol(manh)) {
   manh_i <- manh[, c(1:3, i)]
   SNPs <- lfmm_res$SNP[lfmm_res$ENV == colnames(manh_i)[4]]
   colnames(manh_i)[4] <- "P"
   qqman::manhattan(
    manh_i, genomewideline = F, suggestiveline = F,
    main = paste("LFMM - Env. variable ", colnames(manh)[i], sep = ""),
    highlight = SNPs,
    col = c(alpha("gray70", 0.7), alpha("gray30", 0.7)),
    xlab = "Chromosome", cex.lab = 1.5, cex.axis = 1
   )
   box()
  }

 })



 output$downloadTable <- downloadHandler(
  filename = function() {
   paste("lfmm_results_", Sys.Date(), ".csv", sep = "")
  },
  content = function(file) {
   write.csv(lfmm.res(), file, row.names = F, quote = F)  # Assuming lfmm.res() returns a data frame
  },
  contentType = "text/csv"
 )

}

shinyApp(ui = ui, server = server)
