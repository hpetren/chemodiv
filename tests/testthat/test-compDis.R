testCompData <- data.frame(compound = c("limonene",
                                        "benzaldehyde",
                                        "Unknown1",
                                        "faulty"),
                           smiles = c("CC1=CCC(CC1)C(=C)C",
                                      "C1=CC=C(C=C1)C=O",
                                      NA,
                                      "NOTSMILES"),
                           inchikey = c("XMGQYMWWDOXHJM-UHFFFAOYSA-N",
                                        "HUMNYLRZRPPJDN-UHFFFAOYSA-N",
                                        NA,
                                        "NOTINCHIKEY"))

testNpcTable <- data.frame(compound = c("limonene",
                                        "benzaldehyde"),
                           smiles = c("CC1=CCC(CC1)C(=C)C",
                                      "C1=CC=C(C=C1)C=O"),
                           inchikey = c("XMGQYMWWDOXHJM-UHFFFAOYSA-N",
                                        "HUMNYLRZRPPJDN-UHFFFAOYSA-N"),
                           pathway = c("Terpenoids",
                                       "Shikimates and Phenylpropanoids"),
                           superclass = c("Monoterpenoids",NA),
                           class = c("Menthane monoterpenoids",NA),
                           class2 = c("Monocyclic monoterpenoids",NA))

# Skip test that uses internet resources
test_that("all three types of compDis and their mean works", {
  skip_on_cran()
  skip_if_offline()
  expect_output(str(suppressWarnings(compDis(testCompData,
                                             type = c("NPClassifier",
                                                      "PubChemFingerprint",
                                                      "fMCS")))),
                "List of 4")
})

test_that("function stops if no type is correct", {
  expect_error(compDis(testCompData, type = "wrong"))
})

test_that("function works with NPClassifier and table", {
  expect_output(str(suppressWarnings(compDis(testCompData[1:2,],
                                             type = "NPClassifier",
                                             npcTable = testNpcTable))),
                "List of 1")
})
