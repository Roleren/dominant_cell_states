# Use the local development backend's URL codec and keep saved cohorts explicit.
dominant_observatory_api <- local({
  api <- NULL
  function() {
    if (!is.null(api)) return(api)
    repo <- Sys.getenv('RIBOCRYPT_REPO',unset=path.expand('~/livemount/forks/RiboCrypt'))
    if (!file.exists(file.path(repo,'DESCRIPTION'))) stop('Local RiboCrypt checkout required for Observatory URLs.')
    devtools::load_all(repo,quiet=TRUE)
    names <- c('observatory_state_from_inputs','observatory_expand_selections',
      'make_observatory_url_state_param','parse_observatory_url_hash','observatory_validate_url_state')
    api <<- stats::setNames(lapply(names,getFromNamespace,ns='RiboCrypt'),names)
    api
  }
})

dominant_observatory_normalize_state <- function(state) {
  api <- dominant_observatory_api()
  selections <- api$observatory_expand_selections(state$selections)
  ids <- selections$index
  populated <- vapply(ids,function(id) {
    length(selections$plot_selections[[id]])>0L || length(selections$data_table_selections[[id]])>0L
  },logical(1))
  guard <- vapply(ids,function(id) identical(selections$labels[[id]],'All merged'),logical(1)) & !populated
  if (any(populated)) {
    keep <- ids[!guard]
    for (field in c('plot_selections','data_table_selections','labels')) selections[[field]] <- selections[[field]][keep]
    selections$index <- keep
    if (!selections$active_selection_id %in% ids[populated]) selections$active_selection_id <- ids[populated][1L]
  } else {
    # No materialized cohort: open the selected transcript for manual selection.
    # Autostart must not substitute the entire collection for missing evidence.
    state$browser$go <- FALSE
    id <- ids[1L]
    selections <- list(index=id,active_selection_id=id,
      plot_selections=stats::setNames(list(character()),id),
      data_table_selections=stats::setNames(list(character()),id),
      labels=stats::setNames(list('No matched case/control subset'),id))
  }
  result <- api$observatory_state_from_inputs(state$exp,state$color_by,state$view,state$browser,selections)
  api$observatory_validate_url_state(result)
}

dominant_observatory_state_from_inputs <- function(selected_experiment,color_by,
                                                   view='browser',browser=list(),selections=list()) {
  api <- dominant_observatory_api()
  dominant_observatory_normalize_state(api$observatory_state_from_inputs(
    selected_experiment,color_by,view,browser,selections))
}

dominant_observatory_state_param <- function(state) {
  api <- dominant_observatory_api()
  api$make_observatory_url_state_param(dominant_observatory_normalize_state(state))
}

dominant_observatory_upgrade_url <- function(url) {
  if (length(url)!=1L || is.na(url) || !grepl('#Observatory?obs_state=',url,fixed=TRUE)) return(url)
  api <- dominant_observatory_api()
  parsed <- api$parse_observatory_url_hash(sub('^[^#]*','',url))
  if (is.null(parsed)) return(url)
  state <- dominant_observatory_state_from_inputs(parsed$exp,parsed$color_by,parsed$view,parsed$browser,parsed$selections)
  paste0(sub('#.*$','',url),'#Observatory?obs_state=',dominant_observatory_state_param(state))
}

dominant_observatory_upgrade_table <- function(dt) {
  for (field in grep('observatory_url$',names(dt),value=TRUE)) {
    old <- as.character(dt[[field]])
    unique_urls <- unique(old[!is.na(old) & nzchar(old)])
    upgraded <- vapply(unique_urls,dominant_observatory_upgrade_url,character(1))
    valid <- !is.na(old) & nzchar(old)
    old[valid] <- upgraded[match(old[valid],unique_urls)]
    data.table::set(dt,j=field,value=old)
  }
  dt
}
