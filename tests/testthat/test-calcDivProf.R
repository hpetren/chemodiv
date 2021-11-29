testSampData <- data.frame(compA = c(0.3,0.4,0.1,0.2),
                           compB = c(0.4,0.3,0.4,0.4),
                           compC = c(0.3,0.3,0.5,0.4))
testCompDis <- matrix(data = c(0,0.6,0.7,0.6,0,0.3,0.7,0.3,0), nrow = 3)
colnames(testCompDis) <- c("compA", "compB", "compC")
rownames(testCompDis) <- c("compA", "compB", "compC")


test_that("a diversity profile with correct dim and no NA is generated", {
  expect_output(str(calcDivProf(sampleData = testSampData)), "List of 4")
  expect_output(str(calcDivProf(sampleData = testSampData,
                                compDisMat = testCompDis)), "List of 4")
  expect_equal(nrow(calcDivProf(sampleData = testSampData,
                                compDisMat = testCompDis)$divProf),
               nrow(testSampData))
  expect_equal(ncol(calcDivProf(sampleData = testSampData,
                                compDisMat = testCompDis,
                                qMin = 0, qMax = 3, step = 0.1)$divProf),
               length(seq(0, 3, by = 0.1)))
  expect_false(any(is.na(calcDivProf(sampleData = testSampData,
                                     compDisMat = testCompDis)$divProf)))

})

test_that("faulty input is detected and gives error", {
  expect_error(calcDivProf(sampleData = testSampData, step = 100))
  expect_error(calcDivProf(sampleData = testSampData, qMin = -1))
  expect_error(calcDivProf(sampleData = testSampData, qMin = 3, qMax = 1))
})
