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

test_that("beta diversity is calculated", {
  expect_equal(nrow(calcBetaDiv(sampleData = testSampData)), 1)
  expect_equal(nrow(calcBetaDiv(sampleData = testSampData,
                                compDisMat = testCompDis,
                                type = "FuncHillDiv",
                                q = 3)), 1)
  expect_equal(nrow(calcBetaDiv(sampleData = testSampData,
                                compDisMat = testCompDis,
                                type = c("HillDiv", "FuncHillDiv"),
                                q = 3)), 2)
  expect_false(any(is.na(calcBetaDiv(sampleData = testSampData))))
  expect_false(any(is.na(calcBetaDiv(sampleData = testSampData,
                                     compDisMat = testCompDis,
                                     type = "FuncHillDiv",
                                     q = 3))))
})

test_that("wrong/non-logical input is detected and gives error/message", {
  expect_error(calcBetaDiv(sampleData = testSampData,
                           q = -1))
  expect_error(calcBetaDiv(sampleData = testSampData,
                           type = "NotAnIndex"))
  expect_error(calcBetaDiv(sampleData = testSampData,
                           type = "FuncHillDiv"))
  expect_error(calcBetaDiv(sampleData = testSampData,
                           compDisMat = testCompDis[3:1,3:1],
                           type = c("FuncHillDiv")))
  expect_error(calcBetaDiv(sampleData = testSampData,
                           compDisMat = testCompDis2,
                           type = c("FuncHillDiv")))
  expect_message(calcBetaDiv(sampleData = testSampData,
                             compDisMat = testCompDis,
                             type = c("HillDiv")))
})
