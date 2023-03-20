library(shiny)
library(DT)

ui <- fluidPage(
  titlePanel("DINGO single SNP analysis"),
  align = "center",
  fluidRow(
    column(
      width = 10,
      offset = 1,
      sidebarLayout(
        sidebarPanel(
          width = 12,
          h3("Input: Fetal GWAS and Maternal GWAS"),
          fluidRow(
            column(6, offset = 3, align = "left", "Bivariate LDSC intercept:"),
            column(6, offset = 3, align = "left", numericInput("ldsc_int", NULL, value = 0.1, step = 0.01)),
          ),
          fluidRow(
            column(6, offset = 3, align = "left", "Beta from fetal GWAS:"),
            column(6, offset = 3, align = "left", numericInput("beta_fetal_1", NULL, value = 0.1, step = 0.01)),
          ),
          fluidRow(
            column(6, offset = 3, align = "left", "Standard error from fetal GWAS:"),
            column(6, offset = 3, align = "left", numericInput("se_fetal_1", NULL, value = 0.1, step = 0.01))
          ),
          fluidRow(
            column(6, offset = 3, align = "left", "Beta from maternal GWAS:"),
            column(6, offset = 3, align = "left", numericInput("beta_maternal_1", NULL, value = 0.1, step = 0.01)),
            
          ),
          fluidRow(
            column(6, offset = 3, align = "left", "Standard error from maternal GWAS:"),
            column(6, offset = 3, align = "left", numericInput("se_maternal_1", NULL, value = 0.1, step = 0.01))
          ),
          h3("Output: WLM estiamtes, DINGO Meta analysis, and DINGO 2DF test"),
          DT::dataTableOutput("output1"),

        ),
        mainPanel()
      )
    )
  )
)
server <- function(input, output) {
  output$output1 <- DT::renderDataTable({
    #Two degree of freedom DINGO test
    #Compute unbiased estimates of maternal and fetal effect from the unconditional regression coefficients 
    beta_fetal_adj_1 <- ((4/3)*input$beta_fetal_1) - ((2/3)*input$beta_maternal_1) #Unbiased fetal effect
    beta_maternal_adj_1 <- ((4/3)*input$beta_maternal_1) - ((2/3)*input$beta_fetal_1) #Unbiased maternal effect
    
    #Estimate variance of unbiased fetal and maternal effects and their covariance
    fetal_var_adj_1 <- ((16/9)*(input$se_fetal_1)^2) + ((4/9)*(input$se_maternal_1)^2) - ((16/9)*input$ldsc_int*input$se_fetal_1*input$se_maternal_1) #Variance of unbiased fetal effect
    maternal_var_adj_1 <- ((16/9)*(input$se_maternal_1)^2) + ((4/9)*(input$se_fetal_1)^2) - ((16/9)*input$ldsc_int*input$se_fetal_1*input$se_maternal_1) #Variance of unbiased maternal effect
    covar <- ((20/9)*input$ldsc_int*input$se_fetal_1*input$se_maternal_1) - ((8/9)*(input$se_fetal_1)^2)-((8/9)*(input$se_maternal_1)^2) #Covariance between unbiased maternal and fetal effects
    
    #Two degree of freedom DINGO test
    effects <- matrix(c(beta_fetal_adj_1, beta_maternal_adj_1), nrow=1, ncol=2, byrow=TRUE)
    sigma <- matrix(c(fetal_var_adj_1, covar, covar, maternal_var_adj_1), nrow=2, ncol=2, byrow=TRUE)
    chisq_2df <- effects%*%solve(sigma) %*%t(effects)
    results_2df <- pchisq(chisq_2df, df = 2, ncp = 0, lower.tail = FALSE, log.p = FALSE)
    
    #One degree of freedom conditional tests of association
    results_fetal <- pchisq(q = (beta_fetal_adj_1^2/fetal_var_adj_1), df = 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)		#1 df conditional fetal test
    results_maternal <- pchisq(q = (beta_maternal_adj_1^2/maternal_var_adj_1), df = 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)		#1 df conditional maternal test
    
    #Perform meta-analysis of fetal effect
    b_f1 <- input$beta_fetal_1
    b_f2 <- 2*input$beta_maternal_1 #The estimate of the fetal effect estimated from the maternal meta-analysis is two times bm
    var_bm <- input$se_maternal_1^2
    
    #The variance of the fetal effect estimated from the maternal meta-analysis is four times the variance of bm
    var_b_f2 <- 4*var_bm
    var_b_f1 <- input$se_fetal_1^2 #var_bf
    
    w1 <- var_b_f1 #var_bf
    w2 <- var_b_f2 #4*var_bm
    
    #Perform inverse variance weighted meta-analysis of fetal effect estimates
    beta_meta_fetal <- (b_f1/w1 + b_f2/w2)/(1/w1 + 1/w2)
    cov_b_f1_b_f2 <- input$ldsc_int*sqrt(var_b_f1)*sqrt(var_b_f2)
    
    #Calculate variance of beta_meta
    var_beta_meta_fetal <- 1/((1/w1)+(1/w2)) + 2*cov_b_f1_b_f2*((1/w1)/(1/w1+1/w2))*((1/w2)/(1/w1+1/w2))
    
    chisq_fetal <- (beta_meta_fetal)^2/var_beta_meta_fetal
    
    results_meta_fetal <- pchisq(chisq_fetal, df = 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
    
    #Perform meta-analysis of maternal effect
    b_m1 <- input$beta_maternal_1
    b_m2 <- 2*input$beta_fetal_1 #The estimate of the maternal effect estimated from the fetal meta-analysis is two times bf
    var_bf <- input$se_fetal_1^2
    
    #The variance of the maternal effect estimated from the fetal meta-analysis is four times the variance of bf
    var_b_m2 <- 4*var_bf
    var_b_m1 <- input$se_maternal_1^2
    
    w1 <- var_b_m1
    w2 <- var_b_m2
    
    #Perform inverse variance weighted meta-analysis of maternal effect estimates
    beta_meta_maternal <- (b_m1/w1 + b_m2/w2)/(1/w1 + 1/w2)
    
    cov_b_m1_b_m2 <- input$ldsc_int*sqrt(var_b_m1)*sqrt(var_b_m2)
    
    #Calculate variance of beta_meta
    var_beta_meta_maternal <- 1/((1/w1)+(1/w2)) + 2*cov_b_f1_b_f2*((1/w1)/(1/w1+1/w2))*((1/w2)/(1/w1+1/w2))
    
    
    chisq_maternal <- (beta_meta_maternal)^2/var_beta_meta_maternal
    
    results_meta_maternal <- pchisq(chisq_maternal, df = 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
    
    data.frame(
      "Calculated Values" = c("Fetal beta WLM", "Fetal standard error WLM", "Fetal p-value WLM",
                              "Maternal beta WLM", "Maternal standard error WLM", "Maternal p-value WLM",
                              #"Fetal beta meta", "Fetal standard error meta", "Fetal p-value meta",
                              #"Maternal beta meta", "Maternal standard error meta", "Maternal p-value meta",
                              "Fetal p-value meta", "Maternal p-value meta", "P-value 2DF"),
      "Value" = c(beta_fetal_adj_1, sqrt(fetal_var_adj_1), results_fetal, 
                  beta_maternal_adj_1, sqrt(maternal_var_adj_1), results_maternal,
                  #beta_meta_fetal, sqrt(var_beta_meta_fetal), results_meta_fetal, 
                  #beta_meta_maternal, sqrt(var_beta_meta_maternal), results_meta_maternal, 
                  results_meta_fetal,results_meta_maternal,results_2df)
    )
  }
  
  )

}

shinyApp(ui = ui, server = server)