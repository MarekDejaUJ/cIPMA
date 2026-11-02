fixture <- function(level) {
  read.csv(test_path("fixtures", sprintf("cipma016_mobi_target%d.csv", level)))
}

for (level in c(60, 80)) {
  test_that(sprintf("pi = 1 reproduces the NCA results of cIPMA 0.1.6 at level %d", level), {
    res <- mobi_cipma(target_level = level, B = 0L)
    ref <- fixture(level)
    r <- res$cIPMA_results[match(ref$Construct, res$cIPMA_results$Construct), ]
    b <- res$Bottleneck_summary[match(ref$Construct, res$Bottleneck_summary$Construct), ]
    expect_equal(r$Importance, ref$Importance, tolerance = 1e-10)
    expect_equal(r$Performance, ref$Performance, tolerance = 1e-10)
    expect_equal(r$d, ref$NCA_d, tolerance = 1e-10)
    expect_identical(is.na(b$Required_Level), is.na(ref$Required_Level))
    ok <- !is.na(ref$Required_Level)
    expect_equal(round(b$Required_Level[ok], 3), ref$Required_Level[ok])
    # 0.1.6 compared the scores with the required level rounded to three decimals;
    # the difference in failing cases is the cases between the exact and the rounded level.
    scores <- res$Data_Rescaled
    eps <- sqrt(.Machine$double.eps)
    between <- vapply(seq_along(ref$Construct), function(i) {
      if (!ok[i]) return(0)
      x <- scores[[ref$Construct[i]]]
      sum(x >= b$Required_Level[i] - eps & x < ref$Required_Level[i] - eps)
    }, numeric(1))
    expect_equal(ref$Fail_Count - b$Fail_Count, between)
  })
}

test_that("pi = 1 agrees with a live run of the NCA package", {
  skip_if_not_installed("NCA")
  res <- mobi_cipma(target_level = 80, B = 0L)
  s <- res$Data_Rescaled
  for (cn in res$Predictors) {
    m <- NCA::nca_analysis(s, cn, "Loyalty", ceilings = "ce_fdh")
    d_nca <- m$summaries[[cn]]$params["Effect size", "ce_fdh"]
    expect_equal(res$cIPMA_results$d[res$cIPMA_results$Construct == cn], as.numeric(d_nca), tolerance = 1e-10)
  }
})

test_that("the failing shares at pi = 1 equal the percentile bottleneck table of NCA", {
  # The percentile bottleneck table of NCA is the source of the failing shares of
  # the original cIPMA (Hauff et al. 2024, Table 5).
  skip_if_not_installed("NCA")
  for (level in c(60, 80)) {
    res <- mobi_cipma(target_level = level, B = 0L)
    s <- res$Data_Rescaled
    for (cn in res$Predictors) {
      m <- NCA::nca_analysis(s, cn, "Loyalty", ceilings = "ce_fdh", bottleneck.x = "percentile",
                             bottleneck.y = "actual", steps = c(level, NA))
      bt <- m$bottlenecks$ce_fdh
      row <- which(abs(as.numeric(as.character(bt[[1]])) - level) < 1e-6)
      nca_pct <- suppressWarnings(as.numeric(as.character(bt[row, cn])))
      b <- res$Bottleneck_summary[res$Bottleneck_summary$Construct == cn, ]
      if (is.na(nca_pct)) {
        expect_equal(b$Fail_Count, 0)
      } else {
        expect_equal(round(b$Fail_Percentage, 1), nca_pct)
      }
    }
  }
})

test_that("the tolerance path covers 1.00 to 0.90 and moves in one direction", {
  res <- mobi_cipma(target_level = 80, pi = 0.95, B = 0L)
  p <- res$Tolerance_Path
  expect_equal(sort(unique(round(p$pi, 2)), decreasing = TRUE), round(seq(1, 0.90, by = -0.01), 2))
  for (cn in res$Predictors) {
    s <- p[p$Construct == cn, ]
    s <- s[order(-s$pi), ]
    expect_true(all(diff(s$d) >= -1e-9))
    req <- ifelse(s$Status == "not necessary", -1, s$Required_Level)
    expect_true(all(diff(req[s$Status != "not attained"]) >= -1e-9))
  }
  expect_equal(res$Bottleneck_summary$Required_Level,
               p$Required_Level[abs(p$pi - 0.95) < 1e-9][match(res$Bottleneck_summary$Construct,
                                                                p$Construct[abs(p$pi - 0.95) < 1e-9])])
})

test_that("necessity is decided at pi = 1 with the permutation test", {
  res <- mobi_cipma(target_level = 80, pi = 0.95, B = 19L)
  r <- res$cIPMA_results
  expect_true(all(r$p_value >= 1 / 20 & r$p_value <= 1))
  expect_identical(r$Is_Necessary, r$d >= 0.10 & r$p_value < 0.05)
  expect_true(all(r$d_pi >= r$d - 1e-9))
})

test_that("summaries, report and plots are produced", {
  res <- mobi_cipma(target_level = 80, pi = 0.95, B = 19L)
  for (lv in c("brief", "standard", "full")) {
    s <- summary(res, detail = lv)
    expect_s3_class(s, "summary.cipma")
    expect_output(print(s), "Combined map")
  }
  expect_output(print(res), "cIPMA: target Loyalty")
  expect_named(summary(res, detail = "full"), c("detail", "header", "main", "detail_table", "readings",
                                                "path_d", "path_required", "bottleneck"))
  f <- tempfile(fileext = ".md")
  lines <- cipma_report(res, file = f)
  expect_true(file.exists(f))
  expect_true(any(grepl("^## Tolerance path", lines)))
  for (type in c("map", "path", "bottleneck")) expect_s3_class(plot(res, type = type), "ggplot")
  expect_s3_class(plot(res, scale = "full"), "ggplot")
})

test_that("the diagnostics run on every predictor", {
  res <- mobi_cipma(target_level = 80, pi = 0.95, B = 0L)
  dg <- cipma_diagnostics(res, constructs = c("Image", "Value"), M = 9L, reps = 3L, seed = 2)
  expect_s3_class(dg, "cipma_diagnostics")
  expect_named(dg$results, c("Image", "Value"))
  expect_named(dg$results$Image, c("band", "outliers", "resolution", "probe"))
  expect_output(print(dg), "Reference band")
  lines <- cipma_report(res, diagnostics = dg)
  expect_true(any(grepl("^## QNCA diagnostics", lines)))
})

test_that("nca_rep is accepted with a deprecation warning", {
  m <- mobi_fit()
  expect_warning(cipma(m$boot, "Loyalty", m$data, m$scales, target_level = 80, nca_rep = 0,
                       pi_path = 1), "deprecated")
})

test_that("invalid tolerances are refused", {
  m <- mobi_fit()
  expect_error(cipma(m$boot, "Loyalty", m$data, m$scales, pi = 1.2, B = 0L), "pi")
  expect_error(cipma(m$boot, "Loyalty", m$data, m$scales, pi_path = c(1, 0), B = 0L), "pi_path")
})
