library(shiny)
library(lfmm)
library(data.table)
library(tools)
library(scales)
library(LEA)

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



      textInput(
        inputId = "user_k",
        label = "Number of latent factors (K)",
      ),


      fileInput(
        inputId = "env_file",
        label = "Environmental data (.txt)",
        multiple = FALSE,
        accept = NULL,
        width = NULL,
        buttonLabel = "Browse...",
        placeholder = "No file selected",
        capture = NULL
      )

    ),
     mainPanel(
        plotOutput(outputId = "snmf_plot"),plotOutput(outputId = "manh_plot",  height=800),tableOutput(outputId = "lfmm_res_table")
    )
  )
)





server <- function(input, output) {

 lfmm.res <- reactiveVal(NULL)
 options(shiny.maxRequestSize = 3000*1024^2)

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

      # let’inspect the values of the cross-entropy criterion for each K:
      plot(obj.snmf, pch = 16, col = "blue")

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

       env_file <- input$env_file$datapath

       X <- read.table(env_file, h = TRUE, stringsAsFactors = FALSE)
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
 mylabs <- unique(map$V1) # we need this to preserve the order, otherwise it will order alphabetically
 map$V1 <- as.numeric(factor(map$V1, levels = mylabs))
 qvalues$SNP <- as.character(map$V2)
 qvalues$CHR <- map$V1
 qvalues$BP <- map$V4
 qcut <- 0.01
 source("fn-landgen1.R", echo = F, keep.source = TRUE)
 lfmm_res <- lfmm_qvalcut(qvalues = qvalues, cutoff = qcut)
print(lfmm_res)
lfmm.res(lfmm_res)


 pvalues <- as.data.frame(pvalues)
 pvalues$SNP <- as.character(map$V2)
 pvalues$CHR <- map$V1
 pvalues$BP <- map$V4

a <- ncol(X) + 1
b <- a + 2
manh <- pvalues[, c(a:b, 1:ncol(X))]

#### PLOT SINGLE ####
par(mfrow = c(4, 4), oma = c(1, 1, 1, 1), mar = c(4, 6, 3, 1), mgp = c(3, 1, 0))
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

output$lfmm_res_table <- renderTable({
   req(lfmm.res())  # Require that lfmm.res is not NULL
   return(lfmm.res())
})

}

shinyApp(ui = ui, server = server)
