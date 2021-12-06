testSampData <- data.frame(compA = c(0.3,0.4,0.1,0.2),
                           compB = c(0.4,0.3,0.4,0.4),
                           compC = c(0.3,0.3,0.5,0.4))
testCompDis <- matrix(data = c(0,0.6,0.7,0.6,0,0.3,0.7,0.3,0), nrow = 3)
colnames(testCompDis) <- c("compA", "compB", "compC")
rownames(testCompDis) <- c("compA", "compB", "compC")


test_that("all diversity/evenness gives non-NA output", {
  expect_false(any(is.na(calcDiv(sampleData = testSampData,
                                 type = "HillDiv"))))
  expect_false(any(is.na(calcDiv(sampleData = testSampData,
                                 type = "FuncHillDiv",
                                 compDisMat = testCompDis))))
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
  expect_false(any(is.na(calcDiv(sampleData = testSampData,
                                 type = "FAD",
                                 compDisMat = testCompDis))))
})

test_that("all indices can be calculated simultaneously", {
  expect_equal(ncol(calcDiv(sampleData = testSampData,
                       compDisMat = testCompDis,
                       type = c("HillDiv", "FuncHillDiv", "Shannon",
                                "Simpson", "PielouEven", "HillEven", "RaoQ"),
                       q = 1.5)), 7)
})

test_that("wrong input is detected and gives error", {
  expect_error(calcDiv(sampleData = testSampData,
                       type = "HillDiv", q = -1))
  expect_error(calcDiv(sampleData = testSampData,
                       type = "NotAnIndex"))
  expect_error(calcDiv(sampleData = testSampData,
                       type = "FuncHillDiv"))
})


