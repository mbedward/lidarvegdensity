.single_numeric <- function(x) {
  if (is.numeric(x) && length(x) == 1) {
    x
  } else {
    nm <- deparse(substitute(x))
    paste(nm, "should be a single numeric value")
  }
}

.single_integer <- function(x) {
  if (is.integer(x) && length(x) == 1) {
    x
  } else {
    nm <- deparse(substitute(x))
    paste(nm, "should be a single integer value")
  }
}

.gtzero <- function(x) {
  if (is.numeric(x) && all(x > 0)) {
    x
  } else {
    nm <- deparse(substitute(x))
    paste(nm, "should be greater than zero")
  }
}
