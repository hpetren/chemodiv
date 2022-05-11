testCompDis <- matrix(data = c(0,0.6,0.7,
                               0.6,0,0.3,
                               0.7,0.3,0), nrow = 3)
colnames(testCompDis) <- c("compA", "compB", "compC")
rownames(testCompDis) <- c("compA", "compB", "compC")
group <- c("I","I","II","II")
testSampData <- data.frame(compA = c(0.3,0.4,0.1,0.2),
                           compB = c(0.4,0.3,0.4,0.4),
                           compC = c(0.3,0.3,0.5,0.4))
testNpcTable <- data.frame(compound = c("compA", "compB", "compC"),
                           pathway = c("Polyketides",
                                       "Polyketides",
                                       "Carbohydrates"))

testMolNet <- molNet(testCompDis)

test_that("network plots are generated", {
  expect_match(typeof(molNetPlot(testSampData,
                                 testMolNet$networkObject)), "list")
  expect_match(typeof(molNetPlot(testSampData,
                                 testMolNet$networkObject,
                                 groupData = group)), "list")
  expect_match(typeof(molNetPlot(testSampData,
                                 testMolNet$networkObject,
                                 groupData = group,
                                 npcTable = testNpcTable,)), "list")
  expect_match(typeof(molNetPlot(testSampData,
                                 testMolNet$networkObject,
                                 npcTable = testNpcTable)), "list")
  expect_match(typeof(molNetPlot(testSampData,
                                 testMolNet$networkObject,
                                 plotNames = TRUE)), "list")
  expect_match(typeof(molNetPlot(testSampData,
                                 testMolNet$networkObject,
                                 npcTable = testNpcTable,
                                 plotNames = TRUE)), "list")
})

test_that("faulty input is detected and gives error", {
  expect_error(molNetPlot(testSampData,
                          testMolNet$networkObject,
                          groupData = group,
                          plotNames = TRUE))
})
