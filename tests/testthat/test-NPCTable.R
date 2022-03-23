testCompData <- data.frame(compounds = c("limonene", "unknown"),
                           smiles = c("CC1=CCC(CC1)C(=C)C", "NOTSMILES"),
                           inchikey = c("XMGQYMWWDOXHJM-UHFFFAOYSA-N",
                                        "NOTINCHIKEY"))

test_that("NPC-classification is generated", {
  expect_equal(nrow(NPCTable(testCompData[1,])), nrow(testCompData[1,]))
  expect_message(NPCTable(testCompData), "Is the SMILES correct?")
})

test_that("warnings and messages work", {
  expect_error(NPCTable(testCompData[2,]))
})
