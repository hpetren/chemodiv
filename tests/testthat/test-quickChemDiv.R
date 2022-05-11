testSampData <- data.frame(limonene = c(0.3,0.4,0.1,0.2),
                           benzaldehyde = c(0.4,0.3,0.4,0.4),
                           Unknown1 = c(0.3,0.3,0.5,0.4))
groups <- c("I","I","II","II")

# Only testing without compound data
testQuickChemDiv1 <- quickChemDiv(sampleData = testSampData,
                                  output = "data")

testQuickChemDiv2 <- suppressWarnings(quickChemDiv(sampleData = testSampData,
                                                   groupData = groups,
                                                   output = "plots"))

test_that("quickChemDiv works with data output", {
  expect_match(typeof(testQuickChemDiv1), "list")
  expect_equal(length(testQuickChemDiv1), 3)
})

test_that("quickChemDiv works with plot output", {
  expect_match(typeof(testQuickChemDiv2), "list")
  expect_equal(nrow(testQuickChemDiv2), 2)
  expect_equal(ncol(testQuickChemDiv2), 2)
})

test_that("wrong output argument gives error", {
  expect_error(quickChemDiv(compoundData = testCompData,
                            sampleData = testSampData,
                            groupData = groups,
                            output = "wrong"))
  expect_error(quickChemDiv(compoundData = testCompData,
                            sampleData = testSampData,
                            groupData = groups,
                            output = c("plots", "data")))
})
