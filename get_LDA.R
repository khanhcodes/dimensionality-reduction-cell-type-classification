library(MASS)

## This script is to perform linear discriminant analysis (LDA) for dimensionality reudction
## and obtain the loadings of LDA to conduct projection 
## Return the list containg two objects: the output of LDA (%predictions) and the loadings (%loadings)

get_LDA <- function(scaled_training) {
  data <- read.csv(scaled_training)
  data$X <- NULL

  # Get the output of LDA
  LDA <- lda(cell_type ~ ., data=data)
  predictions_full <- predict(LDA)
  predictions <- predictions_full$x
  write.csv(predictions, "LDA_output_R.csv")

  # Get the loading matrix of LDA
  loadings <- round(LDA$scaling, 2)
  write.csv(loadings, "LDA_loadings_R.csv")

  newlist <- list(predictions = predictions,loadings = loadings)
  
  return(newlist)

}

#combo = get_LDA('/content/scaled_training_label_sample_official.csv')
#print(combo$predictions)
#print(combo$loadings)