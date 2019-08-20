context("test-plot")

df1 <- data.frame(y = c(0,0,0,1,1,0,0,0,1,1,1,0),unit = rep(c(0,1),each=6))
f1 <- funnel(y~1|unit,data=df1)
df2 <- data.frame(y = rbinom(n=1e2,size=1,prob=0.5), unit = rep(c(1:10),each=1e2/10))
f2 <- funnel(y~1|unit,data=df2)
df3 <- data.frame(y = rbinom(n=1e2,size=1,prob=0.5), unit = sample(c(1:10),size=1e2,replace=TRUE))
f3 <- funnel(y~1|unit,data=df3)

test_that("plotting as expected", {
  expect_error(plot(f1))
  expect_error(plot(f2))
  expect_s3_class(plot(f3),c("gg","ggplot"))
  expect_error(plot(f3,label="100"))
  expect_error(plot(f3,identify="100"))
})
