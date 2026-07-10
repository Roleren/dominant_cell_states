\name{dcs_project_root}
\alias{dcs_project_root}
\alias{dcs_source_path}
\alias{dcs_results_dir}
\alias{dcs_result_path}
\title{Dominant-cell-states project paths}
\usage{
dcs_project_root(start = getwd())

dcs_source_path(..., root = dcs_project_root())

dcs_results_dir(root = dcs_project_root(), create = FALSE)

dcs_result_path(..., root = dcs_project_root(), create_parent = FALSE)
}
\arguments{
\item{start}{Directory where root discovery starts.}

\item{...}{Path components.}

\item{root}{Dominant-cell-states project root.}

\item{create}{Whether to create the results directory.}

\item{create_parent}{Whether to create the parent directory of the result path.}
}
\value{
A normalized project root or constructed path.
}
\description{
Helpers for finding the standalone dominant-cell-states checkout and for
constructing source and generated-result paths. Generated outputs should be
placed under \code{results/}.
}
