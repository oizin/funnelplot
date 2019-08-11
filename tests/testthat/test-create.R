context("test-create")

df1 <- data.frame(y = c(0,0,0,1,1,0,0,0,1,1,1,0),
  unit = rep(0,1,each=6))
df2 <- data.frame(y = c(9,0,0,1,1,0,0,0,1,1,1,0),
  unit = rep(0,1,each=6))
df3 <- data.frame(y = c(NA,0,0,1,1,0,0,0,1,1,1,0),
  unit = rep(0,1,each=6))

test_that("outcome variable shold be binary only", {
  expect_type(funnel(y~1|unit,data=df1),"list")
  expect_error(funnel(y~1|unit,data=df2))
  expect_error(funnel(y~1|unit,data=df3))
})

test_that("formula is y~1|unit or y~covariates|unit", {
  expect_error(funnel(y~0|unit,data=df1))
  expect_error(funnel(y~1+2|unit,data=df1))
  expect_error(funnel(y~1+x|unit,data=df1))
})
