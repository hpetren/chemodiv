minCompData1 <- data.frame(compounds = c("limonene"),
                           smiles = c("CC1=CCC(CC1)C(=C)C"),
                           inchikey = c("XMGQYMWWDOXHJM-UHFFFAOYSA-N"))

minCompData2 <- data.frame(compounds = c("limonene"),
                           smiles = c("NOTASMILES"),
                           inchikey = c("XMGQYMWWDOXHJM-UHFFFAOYSA-N"))

test_that("NPCTable works", {
  expect_equal(nrow(NPCTable(minCompData1)), nrow(minCompData1))
  expect_warning(NPCTable(minCompData2), "Is the SMILES correct?")
})
