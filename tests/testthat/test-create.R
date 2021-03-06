context("test-create")

df1 <- data.frame(y = c(0,0,0,1,1,0,0,0,1,1,1,0),
  unit = rep(0,1,each=6))
df2 <- data.frame(y = c(9,0,0,1,1,0,0,0,1,1,1,0),
  unit = rep(0,1,each=6))
df3 <- data.frame(y = c(NA,0,0,1,1,0,0,0,1,1,1,0),
  unit = rep(0,1,each=6))
df4 <- data.frame(y = rbinom(n=1e2,size=1,prob=0.5),
                  x1 = rnorm(n=1e2),
                  x2 = rnorm(n=1e2),
                  unit = sample(c(1:10),size=1e2,replace=TRUE))

test_that("outcome variable shold be binary only", {
  expect_type(funnelModel(y~1|unit,data=df1),"list")
  expect_error(funnelModel(y~1|unit,data=df2))
  expect_error(funnelModel(y~1|unit,data=df3))
  expect_error(funnelModel(y~.|unit,data=df1))
})

test_that("formula is y~1|unit or y~covariates|unit", {
  expect_error(funnelModel(y~0|unit,data=df1))
  expect_error(funnelModel(y~1+2|unit,data=df1))
  expect_error(funnelModel(y~1+x|unit,data=df1))
})

test_that("formula works for types y ~ x1^2", {
  expect_error(funnelModel(y~0|unit,data=df1))
  expect_error(funnelModel(y~1+2|unit,data=df1))
  expect_type(funnelModel(y~ I(x1^2) |unit,data=df4),"list")
  expect_type(funnelModel(y~ x1^2 + x2*x1 |unit,data=df4),"list")
  expect_type(funnelModel(y~ x1^2 + splines::bs(x2) |unit,data=df4),"list")
})

test_that("check distTarget and pointTarget functions", {
  expect_type(funnelModel(y~ x1 + x2|unit,control = distTarget(),data=df4),"list")
  expect_type(funnelModel(y~ x1 + x2|unit,control = pointTarget(crtlOverDisp = TRUE),data=df4),"list")
  expect_type(funnelModel(y~ x1 + x2|unit,control = pointTarget(crtlOverDisp = TRUE, standardised = TRUE),data=df4),"list")
  expect_error(funnelModel(y~ x1 + x2|unit,control = pointTarget(multAdj = "bf"),data=df4))
  expect_error(funnelModel(y~ x1 + x2|unit,control = pointTarget(trim = 1.0),data=df4))
  expect_error(funnelModel(y~ x1 + x2|unit,control = distTarget(trim = 1.0),data=df4))
  expect_error(funnelModel(y~ x1 + x2|unit,control = distTarget(limits = 0.0),data=df4))
})


