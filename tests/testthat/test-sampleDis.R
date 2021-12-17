testSampData <- data.frame(compA = c(0.3,0.4,0.1,0.2),
                           compB = c(0.4,0.3,0.4,0.4),
                           compC = c(0.3,0.3,0.5,0.4))
testCompDis <- matrix(data = c(0,0.6,0.7,0.6,0,0.3,0.7,0.3,0), nrow = 3)
colnames(testCompDis) <- c("compA", "compB", "compC")
rownames(testCompDis) <- c("compA", "compB", "compC")
testCompDis2 <- testCompDis
testCompDis3 <- testCompDis
colnames(testCompDis2) <- c("compX", "compB", "compC")
rownames(testCompDis3) <- c("compA", "compY", "compC")


test_that("Bray-Curtis and UniFracs work", {
  expect_equal(nrow(sampleDis(testSampData)$BrayCurtis),
               nrow(testSampData))
  expect_equal(nrow(sampleDis(testSampData,
                              testCompDis,
                              type = "GenUniFrac")$GenUniFrac),
               nrow(testSampData))
  expect_output(str(sampleDis(testSampData,
                              testCompDis,
                              type = c("BrayCurtis","GenUniFrac"))),
               "List of 2")
})

test_that("error actions and messages work", {
  expect_message(sampleDis(testSampData*2))
  expect_message(sampleDis(testSampData, testCompDis, type = "BrayCurtis"))
  expect_error(sampleDis(testSampData, testCompDis2))
  expect_error(sampleDis(testSampData, testCompDis3))
  expect_error(sampleDis(testSampData, testCompDis, type = "wrong"))
  expect_error(sampleDis(testSampData, type = "GenUniFrac"))
})
