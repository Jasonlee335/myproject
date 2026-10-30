required_packages <- c(
  "haven", "dplyr", "tidyr", "tibble", "survey", "ggplot2", "glmnet",
  "survival", "cmprsk", "sandwich", "lmtest", "mice", "splines",
  "patchwork", "svglite"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (getRversion() < "4.5.0") {
  warning("Archived analyses used R 4.5.0; this installation is older.")
}
if (length(missing_packages)) {
  stop("Install required packages before running: ",
       paste(missing_packages, collapse = ", "), call. = FALSE)
}
cat("All declared packages are available. R ", as.character(getRversion()), "\n", sep = "")
