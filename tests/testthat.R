Sys.setenv("R_TESTS" = "")
if (require(testthat)) {
  test_local()
}