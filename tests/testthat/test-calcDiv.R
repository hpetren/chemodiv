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

test_that("all diversity/evenness gives non-NA output", {
  expect_false(any(is.na(calcDiv(sampleData = testSampData,
                                 type = "HillDiv"))))
  expect_false(any(is.na(calcDiv(sampleData = testSampData,
                                 compDisMat = testCompDis,
                                 type = "FuncHillDiv"))))
  expect_false(any(is.na(calcDiv(sampleData = testSampData,
                                 type = "Shannon"))))
  expect_false(any(is.na(calcDiv(sampleData = testSampData,
                                 type = "Simpson"))))
  expect_false(any(is.na(calcDiv(sampleData = testSampData,
                                 type = "PielouEven"))))
  expect_false(any(is.na(calcDiv(sampleData = testSampData,
                                 type = "HillEven"))))
  expect_false(any(is.na(calcDiv(sampleData = testSampData,
                                 type = "RaoQ",
                                 compDisMat = testCompDis))))
})

test_that("all indices can be calculated simultaneously", {
  expect_equal(ncol(calcDiv(sampleData = testSampData,
                            compDisMat = testCompDis,
                            type = c("HillDiv", "FuncHillDiv", "Shannon",
                                     "Simpson", "PielouEven", "HillEven",
                                     "RaoQ"),
                            q = 1.5)), 7)
})

test_that("wrong/non-logical input is detected and gives error/message", {
  expect_error(calcDiv(sampleData = testSampData,
                       type = "HillDiv", q = -1))
  expect_error(calcDiv(sampleData = testSampData,
                       type = "NotAnIndex"))
  expect_error(calcDiv(sampleData = testSampData,
                       type = "FuncHillDiv"))
  expect_error(calcDiv(sampleData = testSampData,
                       compDisMat = testCompDis[3:1,3:1],
                       type = "FuncHillDiv"))
  expect_error(calcDiv(sampleData = testSampData,
                       compDisMat = testCompDis2,
                       type = "FuncHillDiv"))
  expect_message(calcDiv(sampleData = testSampData,
                         compDisMat = testCompDis,
                         type = c("HillDiv", "Shannon",
                                  "Simpson", "PielouEven", "HillEven")))
})

test_that("FuncHillDiv and Raos Q in utils.R produce correct output", {
  expect_equal(round(as.numeric(calcDiv(sampleData = testSampData[1,],
                                        compDisMat = testCompDis,
                                        type = "FuncHillDiv",
                                        q = 1)), 3), 3.169)
  expect_equal(as.numeric(calcDiv(sampleData = testSampData[1,],
                                  compDisMat = testCompDis,
                                  type = "FuncHillDiv",
                                  q = 0)), 3.2)
  expect_equal(as.numeric(calcDiv(sampleData = testSampData[1,],
                                  compDisMat = testCompDis,
                                  type = "RaoQ")), 0.342)
})
