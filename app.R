#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
set.seed(200)
library(shiny)
library(ggpubr)
library(ggplot2)
library(ggrepel)
library(tidyverse) # for data manipulation
library(kableExtra) # for table printing
library(tmle)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("A shiny app for TMLE \n"),
    ## more buttons in: https://fontawesome.com/search?p=3&m=free
    actionButton("do", "Run tmle, click once and wait ...", icon("rocket"), 
                 style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
    # verbatimTextOutput("verb"),
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
      
        sidebarPanel(
          selectInput("A_vs_B", 
                      label = "Select the differences between treatments A and B: ", 
                      choices = c("A > B" = "A_greater_B", 
                                  "A == B" = "A_equal_B", 
                                  "A < B" = "A_less_B"), selected = "A_equal_B"
                      
          ),
            sliderInput("population",
                        label = "Number of patients:",
                        min = 50,
                        max = 200,
                        value = 100), 
          checkboxGroupInput("g.learners", 
                             label = "Select some learners to estimate E[A|W]: ", 
                             choices = c("Logit" = "SL.glmnet", 
                                         "Random Forest" = "SL.ranger", 
                                         "MARS" = "SL.earth"), selected = "SL.glmnet"
          ),
            checkboxGroupInput("Q.learners", 
                               label = "Select some learners to estimate E[Y|A,W]: ", 
                               choices = c("Logit" = "SL.glmnet_", 
                                           "Random Forest" = "SL.ranger_", 
                                           "MARS" = "SL.earth_", 
                                           "Logit.interaction" = "SL.glm.interaction_"), selected = "SL.glmnet_"
            )
        ), 
        # Show a plot of the generated distribution
        # mainPanel(
        #    plotOutput(outputId = "ps_plot"),
        #    fluidRow(column(10,dataTableOutput('coef_sl__model'))), 
        #    plotOutput(outputId = "sl_plot")
        # )
        
        mainPanel(
          tabsetPanel(
            tabPanel("Plots", plotOutput(outputId = "ps_plot"), 
                     plotOutput(outputId = "sl_plot")),
            tabPanel("SuperLearner", fluidRow(column(10,dataTableOutput('coef_sl__model'))))
          )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # output$verb <- eventReactive(input$do, {
  #   renderText({ "We are running tmle please wait ... " })
  # })
  
  # A_vs_B <- reactive({input$A_vs_B})
  valor <- reactive({if (input$A_vs_B == "A_greater_B") {
    5
  } else if (input$A_vs_B == "A_equal_B") {
    1
  } else {
    0.05
  }})
  
  n <- eventReactive(input$do, {input$population})
  W1 <- eventReactive(input$do, {rbinom(n(), size=1, prob=0.2)}) # binary confounder
  W2 <- eventReactive(input$do, {rbinom(n(), size=1, prob=0.5)}) # binary confounder
  W3 <- eventReactive(input$do, {round(runif(n(), min=2, max=7))}) # continuous confounder
  W4 <- eventReactive(input$do, {round(runif(n(), min=0, max=4))}) # continuous confounder
  A <- eventReactive(input$do, {rbinom(n(), size=1, prob= plogis(-0.2 + 0.2*W2() + log(0.1*W3()) + 0.3*W4() + 0.2*W1()*W4()))}) # binary treatment depends on confounders
  Y <- eventReactive(input$do, {rbinom(n(), size=1, prob= plogis(-1 + valor()*A() - 0.1*W1() + 0.2*W2() + 0.3*W3() - 0.1*W4() + sin(0.1*W2()*W4())))}) # binary outcome depends on confounders
  W <- eventReactive(input$do, {data.frame(W1(), W2(), W3(), W4())})
  W_A <- eventReactive(input$do, {data.frame(W(), A())})
  
  SL.glmnet_ <- function(Y, X, family, ...) {# load required packages
    require("glmnet")
    
    lrn_logit <- create.Learner("SL.glmnet")
    
    sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE),  #cvControl = list(stratifyCV = TRUE),
                                 # SL.library = c("SL.xgboost.mrna", "SL.xgboost.cnv"),
                                 SL.library = c(lrn_logit$names), ...)
    
    # pred is the predicted responses for newX (on the scale of the outcome)
    pred <-  sl___model$SL.predict
    # fit returns all objects needed for predict.SL.template
    fit <- list(object = sl___model)
    # declare class of fit for predict.SL.template
    class(fit) <- 'SL.template'
    # return a list with pred and fit
    out <- list(pred = pred, fit = fit)
    return(isolate(out))}
  SL.earth_ <- function(Y, X, family, ...) {
    require("earth")
    
    lrn_mars <- create.Learner("SL.earth")
    
    sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE),  #cvControl = list(stratifyCV = TRUE),
                                 SL.library = c(lrn_mars$names), ...)
    
    # pred is the predicted responses for newX (on the scale of the outcome)
    pred <-  sl___model$SL.predict
    # fit returns all objects needed for predict.SL.template
    fit <- list(object = sl___model)
    # declare class of fit for predict.SL.template
    class(fit) <- 'SL.template'
    # return a list with pred and fit
    out <- list(pred = pred, fit = fit)
    return(isolate(out))
  }
  SL.ranger_ <- function(Y, X, family, ...) {
    require("ranger")
    
    lrn_rf <- create.Learner("SL.ranger")
    
    sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE),  #cvControl = list(stratifyCV = TRUE),
                                 # SL.library = c("SL.rf.mrna", "SL.rf.cnv"),
                                 SL.library = c(lrn_rf$names), ...)
    
    # pred is the predicted responses for newX (on the scale of the outcome)
    pred <-  sl___model$SL.predict
    # fit returns all objects needed for predict.SL.template
    fit <- list(object = sl___model)
    # declare class of fit for predict.SL.template
    class(fit) <- 'SL.template'
    # return a list with pred and fit
    out <- list(pred = pred, fit = fit)
    return(isolate(out))
  }
  SL.glm.interaction_ <- function (Y, X, family, ...) {
    lrn_glm.inter <- create.Learner("SL.glm.interaction")
    
    sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE),  #cvControl = list(stratifyCV = TRUE),
                                 # SL.library = c("SL.rf.mrna", "SL.rf.cnv"),
                                 SL.library = c(lrn_glm.inter$names), ...)
    
    # pred is the predicted responses for newX (on the scale of the outcome)
    pred <-  sl___model$SL.predict
    # fit returns all objects needed for predict.SL.template
    fit <- list(object = sl___model)
    # declare class of fit for predict.SL.template
    class(fit) <- 'SL.template'
    # return a list with pred and fit
    out <- list(pred = pred, fit = fit)
    return(isolate(out))
  }
  
  # sl_ps <- reactive({
  #   CV.SuperLearner(Y = A(), X = W(), family = binomial(),
  #                   # For a real analysis we would use V = 10.
  #                   cvControl = list(V = 5, stratifyCV = TRUE), innerCvControl = list(list(V=5, stratifyCV = TRUE)),
  #                   parallel = "multicore", verbose = TRUE,
  #                   SL.library = gsub("_", "" ,input$g.learners) )
  # })
  
  sl_ps <- eventReactive(input$do, {
    CV.SuperLearner(Y = A(), X = W(), family = binomial(),
                    # For a real analysis we would use V = 10.
                    cvControl = list(V = 5, stratifyCV = TRUE), innerCvControl = list(list(V=5, stratifyCV = TRUE)),
                    parallel = "multicore", verbose = TRUE,
                    SL.library = input$g.learners )
  })
  
  # sl__model <- reactive({
  #   CV.SuperLearner(Y = Y(), X = W_A(), family = binomial(), control = list(saveFitLibrary = TRUE),
  #                   # For a real analysis we would use V = 10.
  #                   cvControl = list(V = 5, stratifyCV = TRUE), innerCvControl = list(list(V=3, stratifyCV = TRUE)),
  #                   parallel = "multicore", verbose = TRUE,
  #                   # SL.library = input$Q.learners )
  #                   SL.library = input$Q.learners )
  # })
  
  sl__model <- eventReactive(input$do, {
    CV.SuperLearner(Y = Y(), X = W_A(), family = binomial(), control = list(saveFitLibrary = TRUE),
                    # For a real analysis we would use V = 10.
                    cvControl = list(V = 5, stratifyCV = TRUE), innerCvControl = list(list(V=3, stratifyCV = TRUE)),
                    parallel = "multicore", verbose = TRUE,
                    # SL.library = input$Q.learners )
                    SL.library = input$Q.learners )
  })
  
  # summary(sl_ps)
  g1W <- eventReactive(input$do, {sl_ps()$SL.predict})
  
    output$ps_plot <- renderPlot({
      
      # PS plot
      ps <- data.frame(cbind(g1W=g1W(),treat=A()))
      ps_plot <- ggplot() +
        # Top
        geom_density(data=subset(ps, treat==0), aes(x = g1W, y = ..density..), fill="#69b3a2" ) +
        geom_label( aes(x=0.5, y=2.5, label="Treatment A"), color="#69b3a2") +
        # Bottom
        geom_density(data=subset(ps, treat==1), aes(x = g1W, y = -..density..), fill= "#404080") +
        geom_label( aes(x=0.5, y=-2.5, label="Treatment B"), color="#404080") +
        labs(x = "Propensity score estimated", y = "Density")
      
      ggarrange(ps_plot, plot(sl__model()) + theme_bw() + labs(title="E[Y|A,W]"), nrow = 1, ncol = 2, labels = c("A", "B"))
      
    })
    
    output$coef_sl__model <- renderDataTable({coef(sl__model())})
    
    output$sl_plot <- renderPlot({
      
      libs_final <- names(sl__model()$AllSL[[1]]$fitLibrary)
      tmle_comparison <- data.frame(matrix(data = 0, nrow = 0, ncol = length(sl__model()$AllSL) + 1))
      
      for (i in libs_final) {
        
        Q1W <- apply(data.frame(lapply(sl__model()$AllSL, function(x) SuperLearner::predict.SuperLearner(x$fitLibrary[[i]]$object, newdata = W_A() %>% mutate(A = 1), type = "response")$pred)), 1, mean)
        Q0W <- apply(data.frame(lapply(sl__model()$AllSL, function(x) SuperLearner::predict.SuperLearner(x$fitLibrary[[i]]$object, newdata = W_A() %>% mutate(A = 0), type = "response")$pred)), 1, mean)
        Q <- data.frame(Q0W,Q1W)
        
        ## ------------------------ ##
        ### ... g. tmle update ... ###
        ## ------------------------ ##
        result.Qcgc <- tmle(Y = Y(), A = A(), W = W(), Q = Q, g1W = g1W(), family = "binomial")
        
        ATE <- data.frame(result.Qcgc$estimates[2]); ATE[is.na(ATE)] <- 0; ATE <- ATE %>% tidyr::pivot_wider(names_from = ATE.CI, values_from = ATE.CI)
        RR <- data.frame(result.Qcgc$estimates[3]); RR[is.na(RR)] <- 0; RR <- RR %>% tidyr::pivot_wider(names_from = RR.CI, values_from = RR.CI)
        OR <- data.frame(result.Qcgc$estimates[4]); OR[is.na(OR)] <- 0; OR <- OR %>% tidyr::pivot_wider(names_from = OR.CI, values_from = OR.CI)
        ATT <- data.frame(result.Qcgc$estimates[6]); ATT[is.na(ATT)] <- 0; ATT <- ATT %>% tidyr::pivot_wider(names_from = ATT.CI, values_from = ATT.CI)
        ATC <- data.frame(result.Qcgc$estimates[7]); ATC[is.na(ATC)] <- 0; ATC <- ATC %>% tidyr::pivot_wider(names_from = ATC.CI, values_from = ATC.CI)
        
        tmle_comparison[i,1] <- i
        tmle_comparison[i,2:6] <- ATE
        tmle_comparison[i,7:12] <- RR
        tmle_comparison[i,13:18] <- OR
        tmle_comparison[i,19:24] <- ATT
        tmle_comparison[i,25:30] <- ATC
        tmle_comparison[i,31] <- apply(data.frame(lapply(sl__model()$AllSL, function(x) x$cvRisk)), 1, mean, na.rm =TRUE)[i]
        
        if (i == libs_final[length(libs_final)]) {
          colnames(tmle_comparison) <- c("learner_strategy", "ATE.psi", "ATE.var.psi", "ATE.pvalue", "ATE.CI.LO", "ATE.CI.UP", "RR.psi", "RR.pvalue", "RR.log.psi", "RR.var.log.psi", "RR.CI.LO", "RR.CI.UP",
                                         "OR.psi", "OR.pvalue", "OR.log.psi", "OR.var.log.psi", "OR.CI.LO", "OR.CI.UP","ATT.psi", "ATT.converged", "ATT.var.psi", "ATT.pvalue", "ATT.CI.LO", "ATT.CI.UP",
                                         "ATC.psi", "ATC.converged", "ATC.var.psi", "ATC.pvalue", "ATC.CI.LO", "ATC.CI.UP", "cvRisk_sl_model1")
        }
      }
      
      tmle_comparison$learner_strategy <- gsub("_All", "", tmle_comparison$learner_strategy)
      # tmle_comparison$learner_strategy <- gsub("Late", "Int", tmle_comparison$learner_strategy)
      
      ## superlearner
      df_super1 <- list()
      k <- 1
      for (i in sl__model()$AllSL) {
        for (j in libs_final) {
          df_super1[[j]] <- SuperLearner::predict.SuperLearner(i$fitLibrary[[j]]$object, newdata = W_A() %>% mutate(A = 1), type = "response")$pred * sl__model()$coef[k,j]
        }
        k <- k + 1
      }
      
      Q1W <- apply(data.frame(df_super1), 1, mean)
      
      df_super0 <- list()
      k <- 1
      for (i in sl__model()$AllSL) {
        for (j in libs_final) {
          df_super0[[j]] <- SuperLearner::predict.SuperLearner(i$fitLibrary[[j]]$object, newdata = W_A() %>% mutate(A = 0), type = "response")$pred * sl__model()$coef[k,j]
        }
        k <- k + 1
      }
      
      Q0W <- apply(data.frame(df_super0), 1, mean)
      
      Q <- data.frame(Q0W,Q1W)
      head(Q); tail(Q); dim(Q)
      
      result.Qcgc <- tmle(Y = Y(), A = A(), W = W(), Q = Q, g1W = g1W(), family = "binomial")
      
      # saveRDS(object = result.Qcgc, file = paste0(ruta, "os_lgg_data_var", nvar, "/tmle_result_var", nvar, ".rds"))
      
      ATE <- data.frame(result.Qcgc$estimates[2]); ATE[is.na(ATE)] <- 0; ATE <- ATE %>% tidyr::pivot_wider(names_from = ATE.CI, values_from = ATE.CI)
      RR <- data.frame(result.Qcgc$estimates[3]); RR[is.na(RR)] <- 0; RR <- RR %>% tidyr::pivot_wider(names_from = RR.CI, values_from = RR.CI)
      OR <- data.frame(result.Qcgc$estimates[4]); OR[is.na(OR)] <- 0; OR <- OR %>% tidyr::pivot_wider(names_from = OR.CI, values_from = OR.CI)
      ATT <- data.frame(result.Qcgc$estimates[6]); ATT[is.na(ATT)] <- 0; ATT <- ATT %>% tidyr::pivot_wider(names_from = ATT.CI, values_from = ATT.CI)
      ATC <- data.frame(result.Qcgc$estimates[7]); ATC[is.na(ATC)] <- 0; ATC <- ATC %>% tidyr::pivot_wider(names_from = ATC.CI, values_from = ATC.CI)
      
      tmle_comparison["SuperLearner",1] <- "SuperLearner"
      tmle_comparison["SuperLearner",2:6] <- ATE
      tmle_comparison["SuperLearner",7:12] <- RR
      tmle_comparison["SuperLearner",13:18] <- OR
      tmle_comparison["SuperLearner",19:24] <- ATT
      tmle_comparison["SuperLearner",25:30] <- ATC
      
      ## ------------------------------------------------------------ ##
      ### ... Graficos y m?s graficos de tmle y del superlearner ... ###
      ## ------------------------------------------------------------ ##
      
      ### ATE plot
      psi.ate <- ggplot(tmle_comparison, aes(x=reorder(learner_strategy, -ATE.psi), y=ATE.psi)) + 
        geom_line() + 
        geom_point() + 
        geom_hline(yintercept=0, linetype="dashed", color = "black") + 
        geom_text_repel(label = round(tmle_comparison$ATE.psi, digits = 4), size=3) + 
        geom_errorbar(aes(ymin=ATE.CI.LO,ymax=ATE.CI.UP), width=.2,
                      position=position_dodge(0.05)) + 
        theme(axis.text.x=element_blank(), axis.title.x = element_blank(), 
              text = element_text(size = 13))
      
      ### ATT plot
      psi.att <- ggplot(tmle_comparison, aes(x=reorder(learner_strategy, -ATE.psi), y=ATT.psi)) + 
        geom_line() + 
        geom_point() + 
        geom_hline(yintercept=0, linetype="dashed", color = "black") + 
        geom_text_repel(label = round(tmle_comparison$ATT.psi, digits = 4), size=3) + 
        geom_errorbar(aes(ymin=ATT.CI.LO,ymax=ATT.CI.UP), width=.2,
                      position=position_dodge(0.05)) + 
        theme(axis.text.x=element_blank(), axis.title.x = element_blank(), 
              text = element_text(size = 13))
      
      ### ATC plot
      psi.atc <- ggplot(tmle_comparison, aes(x=reorder(learner_strategy, -ATE.psi), y=ATC.psi)) + 
        geom_line() + 
        geom_point() + 
        geom_hline(yintercept=0, linetype="dashed", color = "black") + 
        geom_text_repel(label = round(tmle_comparison$ATC.psi, digits = 4), size=3) + 
        geom_errorbar(aes(ymin=ATC.CI.LO,ymax=ATC.CI.UP), width=.2,
                      position=position_dodge(0.05)) + 
        theme(axis.text.x = element_text(angle = 0, hjust=.75,vjust=0.95), 
              text = element_text(size = 13)) + xlab("Method")
      
      psi <- ggarrange(psi.ate, psi.att, psi.atc, # + rremove("x.text"),
                       labels = c("C", "D", "E"),
                       ncol = 1, nrow = 3)
      
      ### RR plot
      psi.rr <- ggplot(tmle_comparison, aes(x=reorder(learner_strategy, -ATE.psi), y=RR.psi)) + 
        geom_line() + 
        geom_point() + 
        geom_hline(yintercept=1, linetype="dashed", color = "black") + 
        geom_text_repel(label = round(tmle_comparison$RR.psi, digits = 4), size=3) + 
        geom_errorbar(aes(ymin=RR.CI.LO,ymax=RR.CI.UP), width=.2,
                      position=position_dodge(0.05)) + 
        theme(axis.text.x=element_blank(), axis.title.x = element_blank(), 
              text = element_text(size = 13))
      
      ### OR plot
      psi.or <- ggplot(tmle_comparison, aes(x=reorder(learner_strategy, -ATE.psi), y=OR.psi)) + 
        geom_line() + 
        geom_point() + 
        geom_hline(yintercept=1, linetype="dashed", color = "black") + 
        geom_text_repel(label = round(tmle_comparison$OR.psi, digits = 4), size=3) + 
        geom_errorbar(aes(ymin=OR.CI.LO,ymax=OR.CI.UP), width=.2,
                      position=position_dodge(0.05)) + 
        theme(axis.text.x = element_text(angle = 0, hjust=.75,vjust=0.95), 
              text = element_text(size = 13)) + xlab("Method")
      
      or.rr <- ggarrange(psi.rr, psi.or, # + rremove("x.text"),
                         labels = c("F", "G"),
                         ncol = 1, nrow = 2)
      
      ggarrange(psi, or.rr, ncol = 2)
      
      # ggarrange(p1, p2, ncol = 1, nrow = 2)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
