test_that("qs2 preserves genomic objects, tables, NULL values and attributes", {
  source(file.path(dcs_project_root(), "scripts", "dominant_serialization.R"))
  path <- tempfile(fileext=".qs"); on.exit(unlink(path))
  object <- list(ranges=GenomicRanges::GRangesList(TX1=GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(10,20), strand="+")),
    table=data.table::data.table(sample=c("A","B"), count=c(0,NA_real_)),
    matrix=matrix(1:6,2,dimnames=list(c("a","b"),c("x","y","z"))), empty=NULL)
  dominant_save_object(object,path)
  expect_equal(dominant_read_object(path),object)
  expect_true(file.exists(path))
})

test_that("serialization thread override is validated", {
  source(file.path(dcs_project_root(), "scripts", "dominant_serialization.R"))
  withr::local_envvar(DOMINANT_QS_THREADS="1")
  expect_identical(dominant_serialization_threads(),1L)
  Sys.setenv(DOMINANT_QS_THREADS="4")
  expect_identical(dominant_serialization_threads(),4L)
  for (value in c("0","-1","1.5","no","Inf","2147483648")) {
    Sys.setenv(DOMINANT_QS_THREADS=value)
    expect_error(dominant_serialization_threads(),"positive integer")
  }
})

test_that("external QS sibling preference respects freshness", {
  source(file.path(dcs_project_root(), "scripts", "dominant_serialization.R"))
  folder <- tempfile();dir.create(folder);on.exit(unlink(folder,recursive=TRUE))
  rds <- file.path(folder,"input.rds");qs <- file.path(folder,"input.qs")
  saveRDS(1,rds);dominant_save_object(2,qs)
  Sys.setFileTime(qs,Sys.time()-60)
  expect_identical(dominant_preferred_object_file(rds),rds)
  Sys.setFileTime(qs,Sys.time()+60)
  expect_identical(dominant_preferred_object_file(rds),qs)
  unlink(rds)
  expect_identical(dominant_preferred_object_file(rds),qs)
})
