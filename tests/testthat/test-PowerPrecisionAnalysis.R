# Estimation of minimum required sample size (MRSS)
test_that("Vectorized MRSS estimation for binomial proportion test is valid.", {
  stat0 <- c(0.00254452926208651, 0.0274217272536247, 0.514559169437911)
  n0 <- c(393, 9518, 12329)
  mrss_unconditional_round3 <- c(646033.442, 58373.382, 1468.654)
  mrss_conditional_round3 <- c(646036.188, 58376.127, 1471.399)
  mrss_asin_round3 <- c(645671.576, 58344.819, 1471.243)
  mrss_unconditional <- PowerAnalysis(
    stat0 = stat0, effect = 0.10, estimand = "mrss", resp_type = "binary",
    stat_type = "mean", effect_type = "relative", method = "unconditional"
  )
  mrss_conditional <- PowerAnalysis(
    stat0 = stat0, effect = 0.10, estimand = "mrss", resp_type = "binary",
    stat_type = "mean", effect_type = "relative", method = "conditional"
  )
  mrss_asin <- PowerAnalysis(
    stat0 = stat0, effect = 0.10, estimand = "mrss", resp_type = "binary",
    stat_type = "mean", effect_type = "relative", method = "asin"
  )
  expect_identical(round(mrss_unconditional, 3), mrss_unconditional_round3)
  expect_identical(round(mrss_conditional, 3), mrss_conditional_round3)
  expect_identical(round(mrss_asin, 3), mrss_asin_round3)
})

test_that("Vectorized MRSS estimation for z test is valid.", {
  # This example is similar to that for vectorized MRSS estimation for a
  # binomial proportion test, except that we introduce artificial standard
  # deviations (SDs), not the binomial proportion SDs.
  stat0 <- c(0.00254452926208651, 0.0274217272536247, 0.514559169437911)
  sd0 <- c(0.001, 0.01, 0.10)
  n0 <- c(393, 9518, 12329)
  mrss_unconditional_round6 <- c(242.450325, 208.760369, 59.288031)
  mrss_unconditional <- PowerAnalysis(
    stat0 = stat0, effect = 0.10, sd0 = sd0, n = n0, estimand = "mrss",
    resp_type = "continuous", stat_type = "mean", effect_type = "relative",
    method = "unconditional"
  )
  expect_identical(round(mrss_unconditional, 6), mrss_unconditional_round6)
})

# Estimation of minimum detectable effect (MDE)
test_that("Vectorized MDE estimation for binomial proportion test is valid.", {
  stat0 <- c(0.00254452926208651, 0.0274217272536247, 0.514559169437911)
  n0 <- c(393, 9518, 12329)
  mde_unconditional_round6 <- c(63.810855, 0.883672, 0.105841)
  mde_asin_round6 <- c(47.004082, 0.873652, 0.105945)
  sample_size <- DurationToSampleSize(
    duration = 1, n0 = n0, t0 = 4, arms = 2, traffic_prop = 0.85
  )
  mde_unconditional <- PowerAnalysis(
    stat0 = stat0, n = sample_size, estimand = "mde",
    resp_type = "binary", stat_type = "mean", effect_type = "relative",
    method = "unconditional"
  )
  mde_asin <- PowerAnalysis(
    stat0 = stat0, n = sample_size, estimand = "mde",
    resp_type = "binary", stat_type = "mean", effect_type = "relative",
    method = "asin"
  )
  expect_identical(round(mde_unconditional, 6), mde_unconditional_round6)
  expect_identical(round(mde_asin, 6), mde_asin_round6)
})

test_that("Vectorized MDE estimation for z test is valid.", {
  # This example is similar to that for vectorized MDE estimation for a
  # binomial proportion test, except that we introduce artificial standard
  # deviations (SDs), not the binomial proportion SDs.
  stat0 <- c(0.00254452926208651, 0.0274217272536247, 0.514559169437911)
  sd0 <- c(0.001, 0.01, 0.10)
  n0 <- c(393, 9518, 12329)
  mde_unconditional_round6 <- c(0.240963, 0.045435, 0.021274)
  sample_size <- DurationToSampleSize(
    duration = 1, n0 = n0, t0 = 4, arms = 2, traffic_prop = 0.85
  )
  mde_unconditional <- PowerAnalysis(
    stat0 = stat0, sd0 = sd0, n = sample_size, estimand = "mde",
    resp_type = "continuous", stat_type = "mean", effect_type = "relative",
    method = "unconditional"
  )
  expect_identical(round(mde_unconditional, 6), mde_unconditional_round6)
})
