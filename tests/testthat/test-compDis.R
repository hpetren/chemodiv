testCompData1 <- data.frame(compounds = c("limonene",
                                          "benzaldehyde",
                                          "Unknown1"),
                            smiles = c("CC1=CCC(CC1)C(=C)C",
                                       "C1=CC=C(C=C1)C=O",
                                       NA),
                            inchikey = c("XMGQYMWWDOXHJM-UHFFFAOYSA-N",
                                         "HUMNYLRZRPPJDN-UHFFFAOYSA-N",
                                         NA))

compDisRes <- compDis(testCompData1[1:2,])

# Only doing one test here, as this is the function that
# takes most time to run (so not wanting to make it too complicated)
test_that("all three types of compDis and their mean works", {
  expect_output(str(compDisRes), "List of 4")
})

# So stop() works correctly
test_that("function stopped if no type is correct", {
  expect_error(compDis(testCompData1[1:2,], type = "wrong"))
})




# So vegan throws a warning with empty compounds. This is no problem
# since I add 1 or mean values. But it makes it so that the test_that
# below (currently not commented) throws it's own warning. I still think
# the test passes (see manual), and also check() passes. So the test is
# ok and I don't have to suppress the warning I think (this can be done
# with suppressWarnings()). But maybe to be clear one could suppress the
# vegdist warning and write own one? Seems there are pros and cons with that:
# https://www.r-bloggers.com/2012/05/a-warning-about-warning/
# Not run as of 2021-11-15, due to changing how compDis() shows
# the get_cid() closed connection warning

# test_that("message is printed with unknown compounds", {
#   expect_message(compDis(testCompData1))
# })

