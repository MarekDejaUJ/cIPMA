test_that("the requirement at pi = 1 is the smallest predictor value among attaining cases", {
  X <- c(10, 25, 30, 42, 55, 61, 70, 85)
  Y <- c(5, 20, 35, 30, 60, 75, 70, 90)
  r <- cIPMA:::qnca_requirement(X, Y, 1, 60)
  expect_equal(r$status, "required")
  expect_equal(r$required, min(X[Y >= 60]))
})

test_that("a requirement at the lowest predictor value is not necessary", {
  X <- c(10, 25, 30, 42)
  Y <- c(50, 20, 35, 30)
  r <- cIPMA:::qnca_requirement(X, Y, 1, 40)
  expect_equal(r$status, "not necessary")
  expect_true(is.na(r$required))
})

test_that("a level above every outcome value is not attained", {
  r <- cIPMA:::qnca_requirement(c(1, 2, 3), c(10, 20, 30), 1, 31)
  expect_equal(r$status, "not attained")
})

test_that("the requirement rises or stays level as the tolerance falls", {
  set.seed(3)
  d <- qnca::generate_reverse_L(300)
  req <- vapply(c(1, 0.98, 0.95, 0.90), function(p) cIPMA:::qnca_requirement(d$X, d$Y, p, 50)$required, 1)
  expect_true(all(diff(req) >= -1e-12))
})

test_that("the failing count excludes the case at the required level", {
  f <- cIPMA:::qnca_failing(c(1, 2, 3, 4), 3)
  expect_equal(f[["count"]], 2)
  expect_equal(f[["percent"]], 50)
})

test_that("a seed leaves the random number state of the session unchanged", {
  set.seed(99)
  before <- .Random.seed
  cIPMA:::with_seed(1, stats::runif(3))
  expect_identical(.Random.seed, before)
})
