test_that("load_dataset lists names and rejects unknowns", {
  expect_identical(load_dataset(), c("athletics20", "worldathletics-mens-100m-2379"))
  expect_error(load_dataset("missing"), "unknown dataset")
})

test_that("athletics20 has 20 finite sorted times and 12 distinct values", {
  x <- load_dataset("athletics20")
  expect_length(x, 20L)
  expect_true(is.numeric(x))
  expect_true(all(is.finite(x)))
  expect_true(all(diff(x) >= 0))
  expect_equal(x[1], 9.58)
  expect_equal(x[20], 9.84)
  expect_length(unique(x), 12L)
})

test_that("the 2379-row list starts with athletics20 and has 58 distinct times", {
  x <- load_dataset("worldathletics-mens-100m-2379")
  y <- load_dataset("athletics20")
  expect_length(x, 2379L)
  expect_true(all(is.finite(x)))
  expect_true(all(diff(x) >= 0))
  expect_equal(x[1:20], y)
  expect_length(unique(x), 58L)
})
