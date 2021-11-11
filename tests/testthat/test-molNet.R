testCompDis <- matrix(data = c(0,0.6,0.7,0.6,0,0.3,0.7,0.3,0), nrow = 3)
colnames(testCompDis) <- c("compA", "compB", "compC")
rownames(testCompDis) <- c("compA", "compB", "compC")

test_that("network object and properties are calculated", {
  expect_equal(length(molNet(testCompDis)), 5)
  expect_true(is.list(molNet(testCompDis)$networkObject))
})
