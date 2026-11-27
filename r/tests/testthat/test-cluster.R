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

test_that("cluster matches C++ golden output on the builtin set", {
  x <- builtin()
  cases <- load_golden(testthat::test_path("builtin_cpp.txt"))
  expect_length(cases, 8L)
  strict_fp <- identical(Sys.getenv("LEX1DX_STRICT_FP"), "1")
  for (c in cases) {
    r <- cluster(x, c$k, mode = c$mode)
    expect_identical(as.integer(r$bounds), c$bounds)
    expect_identical(as.integer(r$D), c$D)
    expect_identical(as.integer(r$Dprime), c$Dprime)
    if (strict_fp) {
      expect_identical(r$sse, c$sse[[1]])
      expect_identical(as.numeric(r$centers), as.numeric(c$centers))
    } else {
      expect_equal(r$sse, c$sse[[1]], tolerance = 1e-12)
      expect_equal(as.numeric(r$centers), as.numeric(c$centers), tolerance = 1e-12)
    }
  }
})

test_that("default_M returns max(4096, 32 * k) and validates k", {
  expect_identical(default_M(1L), 4096L)
  expect_identical(default_M(64L), 4096L)
  expect_identical(default_M(256L), 8192L)
  expect_error(default_M(0L), "k must be")
})

test_that("k clamping and M override work", {
  x <- c(1, 1, 2, 2, 3)
  r <- cluster(x, 10L)
  expect_identical(length(r$sizes), 3L)
  expect_equal(r$sse, 0)
  r2 <- cluster(seq_len(20), 4L, M = 16L)
  expect_identical(length(r2$sizes), 4L)
})

test_that("invalid inputs are rejected", {
  expect_error(cluster(numeric(), 2L), "empty")
  expect_error(cluster(c(1, NaN, 2), 2L), "NA or NaN")
  expect_error(cluster(c(1, Inf, 2), 2L), "Inf")
  expect_error(cluster(c(1, 2, 3), 0L), "k must be")
  expect_error(cluster(c(1, 2, 3), 2L, mode = "bad"), "mode must be")
  expect_error(cluster(c(1, 2, 3), 2L, M = -1L), "M must be")
})

test_that("sig=5 and mode='sig5' keep k via the precision fallback", {
  x <- c(rep(0, 50), 1.000012 + 1e-10 * (0:49 - 25), 1.000017 + 1e-10 * (0:49 - 25), rep(2, 50))
  r5 <- cluster(x, 4L, sig = 5L, M = 16L)
  expect_identical(r5$sig, 6L)
  expect_identical(length(r5$sizes), 4L)
  expect_identical(as.integer(r5$Dprime), 4L)
  r_alias <- cluster(x, 4L, mode = "sig5")
  expect_identical(r_alias$sig, 6L)
  expect_identical(length(r_alias$sizes), 4L)
  expect_identical(as.integer(r_alias$bounds), as.integer(r5$bounds))
})

test_that("weights match expanding integer counts and bad weights are rejected", {
  vals <- c(1, 2, 3)
  wts <- c(1, 3, 1)
  expanded <- c(1, 2, 2, 2, 3)
  rw <- cluster(vals, 2L, weights = wts)
  re <- cluster(expanded, 2L)
  expect_equal(rw$sse, re$sse)
  expect_equal(as.numeric(rw$centers), as.numeric(re$centers))
  mass <- mapply(function(lo, hi) sum(wts[(lo + 1L):hi]), rw$bounds[-length(rw$bounds)], rw$bounds[-1])
  expect_equal(as.numeric(mass), as.numeric(re$sizes))
  expect_error(cluster(vals, 2L, weights = c(1, -1, 1)), "weights must be finite")
  expect_error(cluster(vals, 2L, weights = c(1, NaN, 1)), "NA or NaN")
  expect_error(cluster(vals, 2L, weights = c(1, 1)), "weights size must match")
  expect_error(cluster(vals, 2L, weights = c(0, 0, 0)), "positive finite total")
})

test_that("exact solvers agree on SSE and sig is validated", {
  set.seed(1)
  x <- rlnorm(2000)
  rw <- cluster(x, 8L, mode = "exact", solver = "wilber")
  rc <- cluster(x, 8L, mode = "exact", solver = "ckmeans")
  denom <- max(rw$sse, rc$sse)
  rel <- if (denom == 0) 0 else abs(rw$sse - rc$sse) / denom
  expect_lt(rel, 1e-12)
  expect_error(cluster(x, 4L, sig = 10L), "0..9")
  expect_error(cluster(x, 4L, sig = 5.5), "0..9")
  expect_error(cluster(x, 4L, solver = "other"), "solver must be")
  expect_error(cluster(x, 4L, mode = "bad"), "mode must be")
})

test_that("Dprime is nondecreasing in sig on a synthetic sample", {
  set.seed(1)
  x <- rlnorm(2000)
  prev <- 0L
  for (d in 1:9) {
    r <- cluster(x, 8L, sig = d)
    expect_identical(r$sig, as.integer(d))
    expect_true(as.integer(r$Dprime) >= prev)
    prev <- as.integer(r$Dprime)
  }
  r0 <- cluster(x, 8L, sig = 0L)
  rt <- cluster(x, 8L)
  expect_identical(r0$sig, 0L)
  expect_identical(as.integer(r0$bounds), as.integer(rt$bounds))
  expect_equal(r0$sse, rt$sse)
})

test_that("athletics20, k = 5: wilber and ckmeans give the same partition", {
  x <- load_dataset("athletics20")
  for (mode in c("exact", "tuned")) {
    rw <- cluster(x, 5L, mode = mode, solver = "wilber")
    rc <- cluster(x, 5L, mode = mode, solver = "ckmeans")
    expect_identical(as.integer(rw$bounds), as.integer(rc$bounds), info = mode)
    denom <- max(rw$sse, rc$sse)
    rel <- if (denom == 0) 0 else abs(rw$sse - rc$sse) / denom
    expect_lt(rel, 1e-12, label = mode)
    expect_identical(rw$solver_used, "wilber")
    expect_identical(rc$solver_used, "ckmeans")
  }
})

test_that("tolerance must be in (0, 1)", {
  x <- c(1, 2, 3, 4)
  expect_error(cluster(x, 2L, tolerance = 0), "tolerance must be in")
  expect_error(cluster(x, 2L, tolerance = 1), "tolerance must be in")
  expect_error(cluster(x, 2L, tolerance = -1e-11), "tolerance must be in")
  expect_error(cluster(x, 2L, tolerance = 2), "tolerance must be in")
  expect_silent(cluster(x, 2L, tolerance = 1e-11))
})

valid_partition <- function(r, kk) {
  expect_length(r$bounds, kk + 1L)
  expect_equal(length(r$sizes), kk)
  for (i in seq_len(kk)) {
    expect_lt(r$bounds[i], r$bounds[i + 1L])
    expect_gt(r$sizes[i], 0)
  }
}

test_that("wilber falls back on {0,1,2,3} k=3 and 0:255 k=27", {
  rw <- cluster(c(0, 1, 2, 3), 3L)
  rc <- cluster(c(0, 1, 2, 3), 3L, solver = "ckmeans")
  expect_equal(rw$sse, 0.5)
  expect_equal(rc$sse, 0.5)
  expect_identical(as.integer(rw$bounds), as.integer(rc$bounds))
  expect_identical(rw$solver_used, "ckmeans")
  expect_identical(rc$solver_used, "ckmeans")
  valid_partition(rw, 3L)

  x <- 0:255
  rw2 <- cluster(x, 27L)
  rc2 <- cluster(x, 27L, solver = "ckmeans")
  expect_identical(rw2$solver_used, "ckmeans")
  expect_identical(as.integer(rw2$bounds), as.integer(rc2$bounds))
  expect_equal(rw2$sse, rc2$sse)
  valid_partition(rw2, 27L)
})

test_that("skipped DP: solver_used mirrors the requested solver", {
  ties <- rep(2, 20)
  rw <- cluster(ties, 8L)
  rc <- cluster(ties, 8L, solver = "ckmeans")
  expect_equal(length(rw$sizes), 1L)
  expect_identical(rw$solver_used, "wilber")
  expect_identical(rc$solver_used, "ckmeans")
  xs <- c(1, 2, 3, 1, 2)
  cw <- cluster(xs, 10L)
  cck <- cluster(xs, 10L, solver = "ckmeans")
  expect_equal(length(cw$sizes), 3L)
  expect_identical(cw$solver_used, "wilber")
  expect_identical(cck$solver_used, "ckmeans")
})

