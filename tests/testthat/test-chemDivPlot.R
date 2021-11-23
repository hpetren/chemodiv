testSampData <- data.frame(compA = c(0.3,0.4,0.1,0.2),
                           compB = c(0.4,0.3,0.4,0.4),
                           compC = c(0.3,0.3,0.5,0.4))

testCompDis <- matrix(data = c(0,0.6,0.7,0.6,0,0.3,0.7,0.3,0), nrow = 3)
colnames(testCompDis) <- c("compA", "compB", "compC")
rownames(testCompDis) <- c("compA", "compB", "compC")

testSampDis <- sampleDis(sampleData = testSampData)

testDiv <- calcDiv(testSampData)

testDivProf <- calcDivProf(testSampData)

groups <- c("I","I","II","II")

testChemDivPlot <- chemDivPlot(compDisMat = testCompDis,
                               divData = testDiv,
                               divProfData = testDivProf,
                               sampleDisMat = testSampDis,
                               groupData = groups)


test_that("chemodiversity plots are outputted", {
  expect_match(typeof(testChemDivPlot), "list")
  expect_equal(nrow(testChemDivPlot), 2)
  expect_equal(ncol(testChemDivPlot), 2)
})


