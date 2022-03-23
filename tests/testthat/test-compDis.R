testCompData <- data.frame(compounds = c("limonene",
                                         "benzaldehyde",
                                         "Unknown1"),
                           smiles = c("CC1=CCC(CC1)C(=C)C",
                                      "C1=CC=C(C=C1)C=O",
                                      NA),
                           inchikey = c("XMGQYMWWDOXHJM-UHFFFAOYSA-N",
                                        "HUMNYLRZRPPJDN-UHFFFAOYSA-N",
                                        NA))

compDisRes <- suppressWarnings(compDis(testCompData,
                                       type = c("NPClassifier",
                                                "PubChemFingerprint",
                                                "fMCS")))

test_that("all three types of compDis and their mean works", {
  expect_output(str(compDisRes), "List of 4")
  expect_equal(ncol(compDisRes$meanDisMat), 3)
  expect_equal(nrow(compDisRes$meanDisMat), 3)
})

test_that("function stops if no type is correct", {
  expect_error(compDis(testCompData, type = "wrong"))
})
