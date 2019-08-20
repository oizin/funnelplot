context("test-eval")

df1 <- data.frame(y = c(0,0,0,1,1,0,0,0,1,1,1,0),unit = rep(c(0,1),each=6))
f1 <- funnel(y~1|unit,data=df1)

test_that("error if too few folds", {
  expect_s3_class(evalCasemixAdj(f1,folds = 2),"data.frame")
  expect_s3_class(evalCasemixAdj(f1,folds = 3),"data.frame")
  expect_s3_class(evalCasemixAdj(f1,folds = 4),"data.frame")
  expect_s3_class(evalCasemixAdj(f1,folds = 5),"data.frame")
  expect_error(evalCasemixAdj(f1,folds = 6))
})

df2 <- data.frame(y = rbinom(n=1e5,size=1,prob=0.5),unit = rep(c(0,1),each=1e5/2))
f2 <- funnel(y~1|unit,data=df2)
nfolds <- 2
res2 <- evalCasemixAdj(f2,folds = 2)

test_that("Check evalCasemixAdj returns values indicating noise", {
  expect_equal(res2$brier[nfolds+1],0.25,tolerance = 0.05)
  expect_equal(res2$accuracy[nfolds+1],0.5,tolerance = 0.05)
  expect_equal(res2$pseudoR2[nfolds+1],0.0,tolerance = 0.05)
  expect_equal(res2$auc_roc[nfolds+1],0.5,tolerance = 0.05)
  expect_equal(res2$no_info_rate[nfolds+1],0.5,tolerance = 0.05)
})
