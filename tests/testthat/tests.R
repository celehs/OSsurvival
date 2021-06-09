test_that("Kern function works", {
  expect_equal(as.numeric(round(Kern.FUN(1,1,1),1)), 0.4)
})
