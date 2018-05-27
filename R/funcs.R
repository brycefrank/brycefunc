#' Removes a column by name from a dataframe.
#'
#' @param data_frame A target dataframe
#' @param column_vector A string vector of column names to remove.
remove_col <- function(data_frame, column_vector) {
  data_frame[, !(names(data_frame) %in% column_vector)]
}

#' Retrieves the lm object from a leaps object by index.
#'
#' @param leaps_object The leaps object to retrieve from.
#' @param model_index The model index to retrieve.
#' @param x The predictors used to construct the leaps_object
#' @param y A vector of the predicted variable
select_leaps_model <- function(leaps_object, model_index, x, y) {
  leaps_frame <- x[, names(x)[leaps_object$which[model_index, ]]]
  leaps_frame$y <- y
  lm(y ~ . , data = leaps_frame)
}
