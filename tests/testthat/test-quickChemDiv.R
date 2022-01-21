testCompData <- data.frame(compounds = c("limonene",
                                         "benzaldehyde",
                                         "Unknown1"),
                           smiles = c("CC1=CCC(CC1)C(=C)C",
                                      "C1=CC=C(C=C1)C=O",
                                      NA),
                           inchikey = c("XMGQYMWWDOXHJM-UHFFFAOYSA-N",
                                        "HUMNYLRZRPPJDN-UHFFFAOYSA-N",
                                        NA))
testSampData <- data.frame(limonene = c(0.3,0.4,0.1,0.2),
                           benzaldehyde = c(0.4,0.3,0.4,0.4),
                           Unknown1 = c(0.3,0.3,0.5,0.4))
groups <- c("I","I","II","II")

testQuickChemDiv1 <- quickChemDiv(compoundData = testCompData,
                                  sampleData = testSampData,
                                  groupData = groups,
                                  output = "data")

testQuickChemDiv2 <- quickChemDiv(compoundData = testCompData,
                                  sampleData = testSampData,
                                  groupData = groups,
                                  output = "plots")

test_that("quichChemDiv works with data output", {
  expect_match(typeof(testQuickChemDiv1), "list")
  expect_equal(length(testQuickChemDiv1), 4)
})

test_that("quichChemDiv works with plot output", {
  expect_match(typeof(testQuickChemDiv2), "list")
  expect_equal(nrow(testQuickChemDiv2), 2)
  expect_equal(ncol(testQuickChemDiv2), 2)
})

test_that("wrong output argument gives error", {
  expect_error(quickChemDiv(compoundData = testCompData,
                            sampleData = testSampData,
                            groupData = groups,
                            output = "wrong"))
})
