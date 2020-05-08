library(shiny)
library(devtools)
library(BiocManager)
options(repos = BiocManager::repositories())
library(GEOquery) 
library(R.utils)
library(reshape2)
library(ggplot2)
library(limma)
library(dplyr)


shinyServer(function(input, output) {
    clinical_outcome1 <-getGEO("GSE120396")
    clinical_outcome<- clinical_outcome1$GSE120396_series_matrix.txt.gz
    rejection_status1  <- clinical_outcome$characteristics_ch1.1
    rejection_status <- unlist(lapply( strsplit(as.character(rejection_status1), ": " ) , `[[` , 2)  )
    
    datadir = "/Users/apple/Documents/2020S1/DATA3888/LAB/GSE120396_RAW"
    fileNames <- list.files(datadir)
    gse = c()
    for(i in 1:length(fileNames)){
        temptable <- read.delim(file.path(datadir, fileNames[i]), header=TRUE)
        gse <- cbind(gse, temptable[,2])
        colnames(gse)[i] <- colnames(temptable)[2]
    }
    rownames(gse) = read.delim(file.path(datadir, fileNames[1]), header=TRUE)[,1]
    largevar = apply(gse, 1, var)
    ind = which(largevar > quantile(largevar, 0.9))
    X = as.matrix(t(gse[ind,]))
    y = rejection_status
    
    
    knn_acc <- reactive({
        cvK = 5  # number of CV folds
        cv_50acc5_knn = c()
        cv_acc_knn = c()
        
        n_sim = 25 ## number of repeats
        for (i in 1:n_sim) {
            
            cvSets = cvTools::cvFolds(nrow(X), cvK)  # permute all the data, into 5 folds
            cv_acc_knn = cv_acc_svm = cv_acc_rf = c()
            
            for (j in 1:cvK) {
                test_id = cvSets$subsets[cvSets$which == j]
                X_test = X[test_id, ]
                X_train = X[-test_id, ]
                y_test = y[test_id]
                y_train = y[-test_id]
                
                fit5 = class::knn(train = X_train, test = X_test, cl = y_train, k = 5)
                cv_acc_knn[j] = table(fit5, y_test) %>% diag %>% sum %>% `/`(length(y_test))
            }
            cv_50acc5_knn <- append(cv_50acc5_knn, mean(cv_acc_knn))
        }
        return(cv_50acc5_knn)
    })
    
    svm_acc <- reactive({
        cvK = 5  # number of CV folds
        cv_50acc5_svm = c()
        cv_acc_svm = c()
        
        n_sim = 25 ## number of repeats
        for (i in 1:n_sim) {
            
            cvSets = cvTools::cvFolds(nrow(X), cvK)  # permute all the data, into 5 folds
            cv_acc_svm = c()
            
            for (j in 1:cvK) {
                test_id = cvSets$subsets[cvSets$which == j]
                X_test = X[test_id, ]
                X_train = X[-test_id, ]
                y_test = y[test_id]
                y_train = y[-test_id]
                svm_res <- e1071::svm(x = X_train, y = as.factor(y_train))
                fit <- predict(svm_res, X_test)
                cv_acc_svm[j] = table(fit, y_test) %>% diag %>% sum %>% `/`(length(y_test))
            }
            cv_50acc5_svm <- append(cv_50acc5_svm, mean(cv_acc_svm))
        }
        return(cv_50acc5_svm)
    })
    
    rf_acc <- reactive({
        cvK = 5  # number of CV folds
        cv_50acc5_rf = c()
        cv_acc_rf = c()
        
        n_sim = 25 ## number of repeats
        for (i in 1:n_sim) {
            
            cvSets = cvTools::cvFolds(nrow(X), cvK)  # permute all the data, into 5 folds
            cv_acc_rf = c()
            
            for (j in 1:cvK) {
                test_id = cvSets$subsets[cvSets$which == j]
                X_test = X[test_id, ]
                X_train = X[-test_id, ]
                y_test = y[test_id]
                y_train = y[-test_id]
                rf_res <- randomForest::randomForest(x = X_train, y = as.factor(y_train))
                fit <- predict(rf_res, X_test)
                cv_acc_rf[j] = table(fit, y_test) %>% diag %>% sum %>% `/`(length(y_test))
            }
            cv_50acc5_rf <- append(cv_50acc5_rf, mean(cv_acc_rf))
        }
        return(cv_50acc5_rf)
    })
    cp_acc <- reactive({
      cvK = 5  # number of CV folds
      cv_50acc5_knn = cv_50acc5_svm = cv_50acc5_rf = c()
      cv_acc_knn = cv_acc_svm = cv_acc_rf = c()
      
      n_sim = 25 ## number of repeats
      for (i in 1:n_sim) {
        
        cvSets = cvTools::cvFolds(nrow(X), cvK)  # permute all the data, into 5 folds
        cv_acc_knn = cv_acc_svm = cv_acc_rf = c()
        
        for (j in 1:cvK) {
          test_id = cvSets$subsets[cvSets$which == j]
          X_test = X[test_id, ]
          X_train = X[-test_id, ]
          y_test = y[test_id]
          y_train = y[-test_id]
          
          ## KNN
          fit5 = class::knn(train = X_train, test = X_test, cl = y_train, k = 5)
          cv_acc_knn[j] = table(fit5, y_test) %>% diag %>% sum %>% `/`(length(y_test))
          
          ## SVM
          svm_res <- e1071::svm(x = X_train, y = as.factor(y_train))
          fit <- predict(svm_res, X_test)
          cv_acc_svm[j] = table(fit, y_test) %>% diag %>% sum %>% `/`(length(y_test))
          
          ## RandomForest
          rf_res <- randomForest::randomForest(x = X_train, y = as.factor(y_train))
          fit <- predict(rf_res, X_test)
          cv_acc_rf[j] = table(fit, y_test) %>% diag %>% sum %>% `/`(length(y_test))
        }
        cv_50acc5_knn <- append(cv_50acc5_knn, mean(cv_acc_knn))
        cv_50acc5_svm <- append(cv_50acc5_svm, mean(cv_acc_svm))
        cv_50acc5_rf <- append(cv_50acc5_rf, mean(cv_acc_rf))
      } 
      return(list(KNN = cv_50acc5_knn, RandomForest=cv_50acc5_rf, SVM=cv_50acc5_svm))
    })
    output$`Box-Plots` = renderPlot({
        comb <- switch(input$classifier,
               "KNN" = knn_acc(),
               "SVM" = svm_acc(),
               "RandomForest" = rf_acc(),
               "Comparison" = cp_acc())
        boxplot(comb, main="The boxplots of Prediction Accuracy",ylab="Accuracy")
    })
  
})
    


