testCompData1 <- data.frame(compounds = c("limonene"),
                            smiles = c("CC1=CCC(CC1)C(=C)C"),
                            inchikey = c("XMGQYMWWDOXHJM-UHFFFAOYSA-N"))

testCompData2 <- data.frame(compounds = c("limonene"),
                            smiles = c("NOTASMILES"),
                            inchikey = c("XMGQYMWWDOXHJM-UHFFFAOYSA-N"))

test_that("NPCTable works", {
  expect_equal(nrow(NPCTable(testCompData1)), nrow(testCompData1))
  expect_warning(NPCTable(testCompData2), "Is the SMILES correct?")
})
