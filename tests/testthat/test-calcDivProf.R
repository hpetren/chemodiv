testSampData <- data.frame(compA = c(0.3,0.4,0.1,0.2),
                           compB = c(0.4,0.3,0.4,0.4),
                           compC = c(0.3,0.3,0.5,0.4))
testCompDis <- matrix(data = c(0,0.6,0.7,
                               0.6,0,0.3,
                               0.7,0.3,0), nrow = 3)
colnames(testCompDis) <- c("compA", "compB", "compC")
rownames(testCompDis) <- c("compA", "compB", "compC")
testCompDis2 <- testCompDis
colnames(testCompDis2) <- c("compA", "compB", "compX")

test_that("a diversity profile with correct dim and no NA is generated", {
  expect_output(str(calcDivProf(sampleData = testSampData)), "List of 5")
  expect_output(str(calcDivProf(sampleData = testSampData,
                                compDisMat = testCompDis,
                                type = "FuncHillDiv")), "List of 5")
  expect_equal(nrow(calcDivProf(sampleData = testSampData)$divProf),
               nrow(testSampData))
  expect_equal(ncol(calcDivProf(sampleData = testSampData,
                                compDisMat = testCompDis,
                                type = "FuncHillDiv",
                                qMin = 0, qMax = 3, step = 0.1)$divProf),
               length(seq(0, 3, by = 0.1)))
  expect_false(any(is.na(calcDivProf(sampleData = testSampData,
                                     compDisMat = testCompDis,
                                     type = "FuncHillDiv")$divProf)))
})

test_that("wrong/non-logical input is detected and gives error/message", {
  expect_error(calcDivProf(sampleData = testSampData, step = 100))
  expect_error(calcDivProf(sampleData = testSampData, qMin = -1))
  expect_error(calcDivProf(sampleData = testSampData, qMin = 3, qMax = 1))
  expect_error(calcDivProf(sampleData = testSampData,
                           compDisMat = testCompDis,
                           type = "NotAnIndex"))
  expect_error(calcDivProf(sampleData = testSampData,
                           compDisMat = testCompDis,
                           type = c("HillDiv", "FuncHillDiv")))
  expect_error(calcDivProf(sampleData = testSampData,
                           type = "FuncHillDiv"))
  expect_error(calcDivProf(sampleData = testSampData,
                           compDisMat = testCompDis[3:1,3:1],
                           type = "FuncHillDiv"))
  expect_error(calcDivProf(sampleData = testSampData,
                           compDisMat = testCompDis2,
                           type = "FuncHillDiv"))
  expect_message(calcDivProf(sampleData = testSampData,
                             compDisMat = testCompDis,
                             type = "HillDiv"))
})
