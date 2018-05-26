remove_col <- function(data_frame, column_vector) {
  data_frame[, !(names(data_frame) %in% column_vector)]
}


select_leaps_model <- function(leaps_object, model_index, x, y) {
  leaps_frame <- x[, names(x)[leaps_object$which[model_index, ]]]
  leaps_frame$y <- y
  lm(y ~ . , data = leaps_frame)
}
