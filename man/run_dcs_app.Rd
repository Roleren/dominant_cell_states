\name{run_dcs_app}
\alias{run_dcs_app}
\title{Run the RDG atlas browser Shiny app}
\usage{
run_dcs_app(
  root = dcs_project_root(),
  host = Sys.getenv("RDG_ATLAS_BROWSER_HOST", "127.0.0.1"),
  port = {
     value <- Sys.getenv("RDG_ATLAS_BROWSER_PORT", "")
     if (nzchar(value)) as.integer(value) else NULL
  },
  launch.browser = interactive(),
  ...
)
}
\arguments{
\item{root}{Dominant-cell-states project root. Defaults to
\code{dcs_project_root()}.}

\item{host}{Host interface passed to \code{shiny::runApp()}.}

\item{port}{Port passed to \code{shiny::runApp()}. \code{NULL} lets Shiny
choose.}

\item{launch.browser}{Whether Shiny should open a browser.}

\item{...}{Additional arguments passed to \code{shiny::runApp()}.}
}
\value{
The value returned by \code{shiny::runApp()}.
}
\description{
Launches the dominant-cell-state ribosome decision graph atlas browser using
the local analysis checkout as the data root. Generated app inputs are read
from \code{results/}.
}
