source(file.path(dcs_project_root(),'scripts','dominant_observatory_urls.R'))

observatory_fixture <- function() {
  api <- dominant_observatory_api()
  api$observatory_state_from_inputs('all_samples-Homo_sapiens',c('tissue','cell_line'),'browser',
    list(gene='ATF4-ENSG00000128272',tx='ENST00000674920',frames_type='columns',go=TRUE),
    list(index=c('1','2','3'),active_selection_id='3',
      labels=list('1'='Control: WT','2'='Case: stress','3'='All merged'),
      plot_selections=list('1'=c('SRR1','SRR2'),'2'='SRR3','3'=character()),
      data_table_selections=list('1'=c('SRR1','SRR2'),'2'='SRR3','3'=character())))
}

test_that('old empty All merged guard blocks readiness and normalization restores it', {
  api <- dominant_observatory_api()
  original <- observatory_fixture()
  url <- paste0('https://example.org/apps/ribo/#Observatory?obs_state=',api$make_observatory_url_state_param(original))
  before <- api$parse_observatory_url_hash(sub('^[^#]*','',url))
  resolved <- function(s) RiboCrypt:::observatory_resolve_library_selections(
    s$selections$data_table_selections,c('SRR1','SRR2','SRR3','UNREQUESTED'),
    s$selections$active_selection_id,s$selections$labels,s$selections$index)
  expect_false(RiboCrypt:::observatory_browser_ready_to_kickoff(before,before$browser,resolved(before)))
  updated <- dominant_observatory_upgrade_url(url)
  expect_match(updated,'https://example.org/apps/ribo/#Observatory?obs_state=',fixed=TRUE)
  after <- api$parse_observatory_url_hash(sub('^[^#]*','',updated))
  expect_identical(after$selections$index,c('1','2'))
  expect_identical(after$selections$active_selection_id,'1')
  expect_identical(after$selections$data_table_selections[['1']],c('SRR1','SRR2'))
  expect_identical(after$selections$data_table_selections[['2']],'SRR3')
  expect_true(RiboCrypt:::observatory_browser_ready_to_kickoff(after,after$browser,resolved(after)))
  expect_identical(dominant_observatory_upgrade_url(updated),updated)
})

test_that('missing matched cohorts do not autostart the entire collection', {
  state <- observatory_fixture()
  state$selections$runs <- character()
  for (field in c('p','d')) state$selections[[field]] <- list('1'=integer(),'2'=integer(),'3'=integer())
  new <- dominant_observatory_normalize_state(state)
  expect_false(new$browser$go)
  expect_length(new$selections$order,1L)
  expect_identical(new$selections$labels[[1]],'No matched case/control subset')
  expect_length(new$selections$runs,0L)
})

test_that('populated groups including All merged and partial-pooling cohorts survive', {
  api <- dominant_observatory_api()
  state <- observatory_fixture()
  state$selections$p[['3']] <- state$selections$d[['3']] <- c(1L,2L,3L)
  new <- dominant_observatory_normalize_state(state)
  expect_identical(new$selections$order,c('1','2','3'))
  expect_identical(new$selections$active,'3')
  expanded <- api$observatory_expand_selections(new$selections)
  expect_identical(expanded$data_table_selections[['3']],c('SRR1','SRR2','SRR3'))
})

test_that('stored table URLs are normalized without changing review labels', {
  api <- dominant_observatory_api()
  old <- paste0('https://ribocrypt.org/#Observatory?obs_state=',api$make_observatory_url_state_param(observatory_fixture()))
  dt <- data.table::data.table(observatory_url=c(old,'',NA_character_),
    partial_pooling_observatory_url=c(old,'bad-link',''),manual_review_status=c('positive','negative',''))
  new <- dominant_observatory_upgrade_table(dt)
  expect_identical(new$manual_review_status,c('positive','negative',''))
  expect_identical(new$observatory_url[2:3],c('',NA_character_))
  expect_identical(new$partial_pooling_observatory_url[2],'bad-link')
  expect_identical(new$observatory_url[1],dominant_observatory_upgrade_url(old))
})
