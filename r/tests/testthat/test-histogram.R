test_that("bin_table matches sizes and centers", {
  x <- load_dataset("athletics20")
  r <- cluster(x, 5L)
  tab <- bin_table(r, x)
  expect_equal(nrow(tab), length(r$sizes))
  expect_equal(sum(tab$n), length(x))
  expect_equal(tab$n, as.numeric(r$sizes))
  expect_equal(tab$mean, as.numeric(r$centers), tolerance = 1e-12)
  expect_equal(sum(tab$pct), 100, tolerance = 1e-9)
  sx <- sort(x)
  b <- as.integer(r$bounds)
  for (i in seq_len(nrow(tab))) {
    expect_equal(tab$lo[i], sx[b[i] + 1L])
    expect_equal(tab$hi[i], sx[b[i + 1L]])
  }
})

test_that("histogram writes a PNG", {
  x <- load_dataset("athletics20")
  r <- cluster(x, 5L)
  f <- tempfile(fileext = ".png")
  histogram(r, x, file = f, title = "athletics20")
  expect_true(file.exists(f))
  expect_gt(file.info(f)$size, 0)
  con <- file(f, "rb")
  magic <- readBin(con, what = "raw", n = 8)
  close(con)
  expect_identical(magic, as.raw(c(0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A)))
})
