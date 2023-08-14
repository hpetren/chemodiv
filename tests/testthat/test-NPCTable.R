testCompData <- data.frame(compound = c("limonene",
                                        "unknown",
                                        "faulty",
                                        "empty"),
                           smiles = c("CC1=CCC(CC1)C(=C)C",
                                      NA,
                                      "NOTSMILES",
                                      ""),
                           inchikey = c("XMGQYMWWDOXHJM-UHFFFAOYSA-N",
                                        NA,
                                        "NOTINCHIKEY",
                                        ""))

# Skip test that uses internet resources
test_that("NPC-classification is generated", {
  skip_on_cran()
  skip_if_offline()
  expect_equal(nrow(NPCTable(testCompData[1,])), nrow(testCompData[1,]))
  expect_message(NPCTable(testCompData), "Is the SMILES correct?")
})

# Skip test that uses internet resources
test_that("warnings and messages work", {
  skip_on_cran()
  skip_if_offline()
  expect_error(NPCTable(testCompData[2,]))
  expect_error(NPCTable(testCompData[3,]))
  expect_error(NPCTable(testCompData[4,]))
})
