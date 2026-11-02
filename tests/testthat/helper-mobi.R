# The ECSI mobile phone model on the mobi data of seminr, shared by the tests.
mobi_fit <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      mobi <- seminr::mobi
      mm <- seminr::constructs(
        seminr::composite("Image", seminr::multi_items("IMAG", 1:5)),
        seminr::composite("Expectation", seminr::multi_items("CUEX", 1:3)),
        seminr::composite("Quality", seminr::multi_items("PERQ", 1:7)),
        seminr::composite("Value", seminr::multi_items("PERV", 1:2)),
        seminr::composite("Satisfaction", seminr::multi_items("CUSA", 1:3)),
        seminr::composite("Complaints", seminr::single_item("CUSCO")),
        seminr::composite("Loyalty", seminr::multi_items("CUSL", 1:3))
      )
      sm <- seminr::relationships(
        seminr::paths(from = "Image", to = c("Expectation", "Satisfaction", "Loyalty")),
        seminr::paths(from = "Expectation", to = c("Quality", "Value", "Satisfaction")),
        seminr::paths(from = "Quality", to = c("Value", "Satisfaction")),
        seminr::paths(from = "Value", to = "Satisfaction"),
        seminr::paths(from = "Satisfaction", to = c("Complaints", "Loyalty")),
        seminr::paths(from = "Complaints", to = "Loyalty")
      )
      pls <- suppressMessages(seminr::estimate_pls(mobi, mm, sm))
      boot <- suppressMessages(seminr::bootstrap_model(pls, nboot = 20, cores = 1, seed = 7))
      scales <- stats::setNames(rep(list(c(1, 10)), ncol(mobi)), names(mobi))
      cache <<- list(data = mobi, pls = pls, boot = boot, scales = scales)
    }
    cache
  }
})

mobi_cipma <- function(target_level = 80, pi = 1, B = 19L, ...) {
  m <- mobi_fit()
  cipma(m$boot, "Loyalty", m$data, m$scales, target_level = target_level, pi = pi,
        B = B, seed = 11, pls_model = m$pls, ...)
}
