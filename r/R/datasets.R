load_dataset <- function(name = NULL) {
  known <- c("athletics20", "worldathletics-mens-100m-2379")
  files <- c(
    athletics20 = "athletics20.csv",
    "worldathletics-mens-100m-2379" = "worldathletics-mens-100m-2379.csv"
  )
  if (is.null(name)) return(known)
  if (!is.character(name) || length(name) != 1L)
    stop("name must be a single string")
  if (!name %in% known)
    stop("unknown dataset '", name, "'; known names: ", paste(known, collapse = ", "))
  path <- system.file("extdata", files[[name]], package = "lex1dx", mustWork = TRUE)
  x <- utils::read.csv(path, header = TRUE)$time
  if (length(x) == 0L || !is.numeric(x) || any(!is.finite(x)))
    stop("dataset '", name, "' is not a finite numeric column")
  as.numeric(x)
}
