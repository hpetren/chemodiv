testSampData <- data.frame(compA = c(0.3,0.4,0.1,0.2),
                           compB = c(0.4,0.3,0.4,0.4),
                           compC = c(0.3,0.3,0.5,0.4))

testCompDis <- matrix(data = c(0,0.6,0.7,0.6,0,0.3,0.7,0.3,0), nrow = 3)
colnames(testCompDis) <- c("compA", "compB", "compC")
rownames(testCompDis) <- c("compA", "compB", "compC")


test_that("beta diversity is calculated", {
  expect_equal(ncol(calcBetaDiv(sampleData = testSampData)), 3)
  expect_equal(ncol(calcBetaDiv(sampleData = testSampData,
                                compDisMat = testCompDis, q = 3)), 3)
  expect_false(any(is.na(calcBetaDiv(sampleData = testSampData))))
  expect_false(any(is.na(calcBetaDiv(sampleData = testSampData,
                                     compDisMat = testCompDis, q = 3))))
})

test_that("wrong input is detected and gives error", {
  expect_error(calcBetaDiv(sampleData = testSampData,
                           q = -1))
})
