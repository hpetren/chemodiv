library(chemdiv)

test_that("testing the test_function", {
  expect_equal(nchar(test_function("Some input")), 26)
})
