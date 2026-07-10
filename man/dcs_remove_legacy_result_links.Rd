\name{dcs_remove_legacy_result_links}
\alias{dcs_remove_legacy_result_links}
\title{Remove stale root-level result links}
\usage{
dcs_remove_legacy_result_links(root = dcs_project_root())
}
\arguments{
\item{root}{Dominant-cell-states project root.}
}
\value{
The removed link paths, invisibly.
}
\description{
Removes stale compatibility symlinks that point from the repository root back
into \code{results/}. Current pipeline runs should write directly under
\code{results/} and should not recreate these links.
}
