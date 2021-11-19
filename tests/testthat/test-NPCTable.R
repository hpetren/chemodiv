testCompData <- data.frame(compounds = c("limonene", "unknown"),
                           smiles = c("CC1=CCC(CC1)C(=C)C", "NOTSMILES"),
                           inchikey = c("XMGQYMWWDOXHJM-UHFFFAOYSA-N",
                                        "NOTINCHIKEY"))


test_that("NPCTable works", {
  expect_equal(nrow(NPCTable(testCompData[1,])), nrow(testCompData[1,]))
  expect_error(NPCTable(testCompData[2,]))
  expect_message(NPCTable(testCompData), "Is the SMILES correct?")
})
