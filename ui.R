
library(shiny)

shinyUI(
    fluidPage(
        headerPanel("Classifier accuracy for predicting kisney rejection"),
        sidebarLayout(
            sidebarPanel(
                selectInput("classifier","Select the classifier:", choices = c("KNN", "RandomForest", "SVM", "Comparison")),
            ),
        mainPanel(
            plotOutput(outputId = "Box-Plots")
            )
        )
    ))

