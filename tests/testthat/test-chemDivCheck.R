testCompData <- data.frame(compound = c("limonene",
                                        "benzaldehyde"),
                           smiles = c("CC1=CCC(CC1)C(=C)C",
                                      "C1=CC=C(C=C1)C=O"),
                           inchikey = c("XMGQYMWWDOXHJM-UHFFFAOYSA-N",
                                        "HUMNYLRZRPPJDN-UHFFFAOYSA-N"))

testSampData <- data.frame(limonene = c(0.3,0.4,0.1),
                           benzaldehyde = c(0.7,0.6,0.9))

test_that("Formatting ok message is printed", {
  expect_message(chemDivCheck(testCompData, testSampData),
                 "The two datasets")
})

test_that("Formatting problem messages are printed", {
  expect_error(chemDivCheck(testCompData, as.matrix(testSampData)),
               "sampleData should")
  expect_error(chemDivCheck(as.matrix(testCompData), testSampData),
               "compoundData should")
  expect_message(chemDivCheck(testCompData[,1:2], testSampData),
                 "compoundData should include")
  expect_message(chemDivCheck(testCompData, testSampData * 2),
                 "Not all row sums")
  expect_message(chemDivCheck(testCompData[c(2,1),], testSampData),
                 "The name and order")
  expect_message(chemDivCheck(testCompData, testSampData[,c(2,1)]),
                 "The name and order")
})


