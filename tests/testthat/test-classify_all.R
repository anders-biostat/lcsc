test_that("classify_all returns a list", {
  
  annotation_test <- readRDS("annotation.rds")
  
  class_all_samples <- classify_all(annotation_test)
  
  expect_type(class_all_samples, "list")
})

test_that("classify_all returns a list containing only annotated samples", {
  
  annotation_test <- readRDS("annotation.rds")
  
  count = 0
  for (entry in annotation_test) {
    if (length(entry)!=0) {
      count = count + 1
    }
  }
  
  class_all_samples <- classify_all(annotation_test)
  
  expect_equal(length(class_all_samples), count)
})



