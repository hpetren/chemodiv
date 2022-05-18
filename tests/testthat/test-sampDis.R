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

test_that("Bray-Curtis and UniFracs work", {
  expect_equal(nrow(sampDis(testSampData)$BrayCurtis),
               nrow(testSampData))
  expect_equal(nrow(sampDis(testSampData,
                            testCompDis,
                            type = "GenUniFrac")$GenUniFrac),
               nrow(testSampData))
  expect_output(str(sampDis(testSampData,
                            testCompDis,
                            type = c("BrayCurtis","GenUniFrac"))),
               "List of 2")
})

test_that("error actions and messages work", {
  expect_message(sampDis(testSampData*2))
  expect_message(sampDis(testSampData, testCompDis, type = "BrayCurtis"))
  expect_error(sampDis(testSampData, testCompDis[3:1,3:1]))
  expect_error(sampDis(testSampData, testCompDis2))
  expect_error(sampDis(testSampData, testCompDis, type = "wrong"))
  expect_error(sampDis(testSampData, type = "GenUniFrac"))
})
