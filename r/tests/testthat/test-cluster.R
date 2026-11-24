builtin <- function() {
  # Same 800 points as cpp/tests/test_agree.cpp, stored as C++ IEEE bits
  # so the comparison is not sensitive to libm differences.
  hs <- readLines(testthat::test_path("builtin_x_bits.txt"), warn = FALSE)
  as.numeric(lex1dx:::.bits_to_double(hs))
}

load_golden <- function(path) {
  lines <- readLines(path, warn = FALSE)
  cases <- list()
  cur <- NULL
  flush <- function() {
    if (!is.null(cur)) cases[[length(cases) + 1L]] <<- cur
  }
  for (line in lines) {
    line <- trimws(line)
    if (!nzchar(line) || startsWith(line, "#")) next
    if (startsWith(line, "case=")) {
      flush()
      toks <- strsplit(line, " +")[[1]]
      kv <- strsplit(toks, "=", fixed = TRUE)
      meta <- setNames(vapply(kv, `[`, "", 2L), vapply(kv, `[`, "", 1L))
      cur <- list(mode = meta[["mode"]], k = as.integer(meta[["k"]]),
                  D = as.integer(meta[["D"]]), Dprime = as.integer(meta[["Dprime"]]))
    } else if (startsWith(line, "bounds=")) {
      cur$bounds <- as.integer(strsplit(sub("^bounds=", "", line), ",", fixed = TRUE)[[1]])
    } else if (startsWith(line, "sse_bits=")) {
      cur$sse <- lex1dx:::.bits_to_double(sub("^sse_bits=", "", line))
    } else if (startsWith(line, "centers_bits=")) {
      hs <- strsplit(sub("^centers_bits=", "", line), ",", fixed = TRUE)[[1]]
      cur$centers <- lex1dx:::.bits_to_double(hs)
    }
  }
  flush()
  cases
}

test_that("cluster is bit-identical to C++ on the builtin set", {
  x <- builtin()
  cases <- load_golden(testthat::test_path("builtin_cpp.txt"))
  expect_length(cases, 8L)
  for (c in cases) {
    r <- cluster(x, c$k, mode = c$mode)
    expect_identical(as.integer(r$bounds), c$bounds)
    expect_identical(r$sse, c$sse[[1]])
    expect_identical(as.numeric(r$centers), as.numeric(c$centers))
    expect_identical(as.integer(r$D), c$D)
    expect_identical(as.integer(r$Dprime), c$Dprime)
  }
})

test_that("empty input is rejected", {
  expect_error(cluster(numeric(), 2L), "empty")
})

test_that("NaN is rejected", {
  expect_error(cluster(c(1, NaN, 2), 2L), "NA or NaN")
})

test_that("Inf is rejected", {
  expect_error(cluster(c(1, Inf, 2), 2L), "Inf")
})
