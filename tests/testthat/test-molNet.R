testCompDis <- matrix(data = c(0,0.6,0.7,
                               0.6,0,0.3,
                               0.7,0.3,0), nrow = 3)
colnames(testCompDis) <- c("compA", "compB", "compC")
rownames(testCompDis) <- c("compA", "compB", "compC")
testCompDis2 <- testCompDis
colnames(testCompDis2)[1] <- "compX"
testNpcTable <- data.frame(compound = c("compA", "compB", "compC"),
                           pathway = c("Path1", "Path1", "Path2"))

test_that("network object and properties are calculated", {
  expect_equal(length(molNet(testCompDis)), 4)
  expect_true(is.list(molNet(testCompDis)$networkObject))
  expect_true(is.numeric(molNet(testCompDis, testNpcTable)$nNpcPathways))
})

test_that("faulty input is detected and gives error", {
  expect_error(molNet(testCompDis, cutOff = "minPathway"))
  expect_error(molNet(testCompDis2))
})
