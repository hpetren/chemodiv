testCompData <- data.frame(compound = c("limonene",
                                        "benzaldehyde"),
                           smiles = c("CC1=CCC(CC1)C(=C)C",
                                      "C1=CC=C(C=C1)C=O"),
                           inchikey = c("XMGQYMWWDOXHJM-UHFFFAOYSA-N",
                                        "HUMNYLRZRPPJDN-UHFFFAOYSA-N"))
testSampData <- data.frame(limonene = c(0.3,0.4,0.1),
                           benzaldehyde = c(0.7,0.6,0.9))

test_that("Formatting ok message is printed", {
  expect_message(chemoDivCheck(testSampData, testCompData),
                 "The two datasets")
})

test_that("Formatting problem messages are printed", {
  expect_error(chemoDivCheck(as.matrix(testSampData), testCompData),
               "sampleData should")
  expect_error(chemoDivCheck(testSampData, as.matrix(testCompData)),
               "compoundData should")
  expect_message(chemoDivCheck(testSampData, testCompData[,1:2]),
                 "compoundData should include")
  expect_message(chemoDivCheck(testSampData * 2, testCompData),
                 "Not all row sums")
  expect_message(chemoDivCheck(testSampData, testCompData[c(2,1),]),
                 "The name and order")
  expect_message(chemoDivCheck(testSampData[,c(2,1)], testCompData),
                 "The name and order")
})
