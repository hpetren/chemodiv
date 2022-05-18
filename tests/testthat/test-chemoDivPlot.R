testSampData <- data.frame(compA = c(0.3,0.4,0.1,0.2),
                           compB = c(0.4,0.3,0.4,0.4),
                           compC = c(0.3,0.3,0.5,0.4))
testCompDis <- matrix(data = c(0,0.6,0.7,
                               0.6,0,0.3,
                               0.7,0.3,0), nrow = 3)
colnames(testCompDis) <- c("compA", "compB", "compC")
rownames(testCompDis) <- c("compA", "compB", "compC")
testSampDis <- sampDis(sampleData = testSampData,
                       compDisMat = testCompDis,
                       type = c("BrayCurtis", "GenUniFrac"))
groups <- c("I","I","II","II")

testDiv <- calcDiv(testSampData)
testDivProf <- calcDivProf(testSampData)

testChemoDivPlot <- suppressWarnings(chemoDivPlot(compDisMat = testCompDis,
                                                  divData = testDiv,
                                                  divProfData = testDivProf,
                                                  sampDisMat = testSampDis$BrayCurtis,
                                                  groupData = groups))

testChemoDivPlot2 <- suppressWarnings(chemoDivPlot(sampDisMat = testSampDis,
                                                   groupData = groups))

test_that("chemodiversity plots are outputted", {
  expect_match(typeof(testChemoDivPlot), "list")
  expect_equal(nrow(testChemoDivPlot), 2)
  expect_equal(ncol(testChemoDivPlot), 2)
})

test_that("chemoDivPlot can handle list as sampDisMat input", {
  expect_match(typeof(testChemoDivPlot2), "list")
  expect_equal(nrow(testChemoDivPlot2), 1)
  expect_equal(ncol(testChemoDivPlot2), 2)
})
