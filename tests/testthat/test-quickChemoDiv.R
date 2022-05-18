testSampData <- data.frame(limonene = c(0.3,0.4,0.1,0.2),
                           benzaldehyde = c(0.4,0.3,0.4,0.4),
                           Unknown1 = c(0.3,0.3,0.5,0.4))
groups <- c("I","I","II","II")

# Only testing without compound data
testQuickChemoDiv1 <- quickChemoDiv(sampleData = testSampData,
                                    output = "data")

testQuickChemoDiv2 <- suppressWarnings(quickChemoDiv(sampleData = testSampData,
                                                     groupData = groups,
                                                     output = "plots"))

test_that("quickChemoDiv works with data output", {
  expect_match(typeof(testQuickChemoDiv1), "list")
  expect_equal(length(testQuickChemoDiv1), 3)
})

test_that("quickChemoDiv works with plot output", {
  expect_match(typeof(testQuickChemoDiv2), "list")
  expect_equal(nrow(testQuickChemoDiv2), 2)
  expect_equal(ncol(testQuickChemoDiv2), 2)
})

test_that("wrong output argument gives error", {
  expect_error(quickChemoDiv(compoundData = testCompData,
                             sampleData = testSampData,
                             groupData = groups,
                             output = "wrong"))
  expect_error(quickChemoDiv(compoundData = testCompData,
                             sampleData = testSampData,
                             groupData = groups,
                             output = c("plots", "data")))
})
