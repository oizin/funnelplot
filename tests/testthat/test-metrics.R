context("test-metrics")

phat <- rep(0.5,20)
y <- c(rep(0,10),rep(1,10))

test_that("metrics are correct for noise models", {
  expect_equal(auc_roc(phat,y),0.5)
  expect_equal(accuracy(classify(phat),y),0.5)
})

phat <- c(rep(0,10),rep(1,10))
y <- c(rep(0,10),rep(1,10))

test_that("metrics are correct for oracle models", {
  expect_equal(auc_roc(phat,y),1.0)
  expect_equal(accuracy(classify(phat),y),1.0)
})
