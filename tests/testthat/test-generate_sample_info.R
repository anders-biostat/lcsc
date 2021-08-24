test_that("generate_sample_info(cells) returns a tibble", {
  
  cells_test <- readRDS("cells_test.rds")
  
  sample_info <- generate_sample_info(cells_test)
  
  expect_type(sample_info, "list")
})


test_that("generate_sample_info(cells) returns a tibble with #row = #samples", {
  
  cells_test <- readRDS("cells_test.rds")
  
  sample_info <- generate_sample_info(cells_test)
  
  samples <- unique(cells_test$sample)
  
  expect_equal(nrow(sample_info), length(samples))
})



