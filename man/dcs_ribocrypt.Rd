\name{dcs_ribocrypt_repo}
\alias{dcs_ribocrypt_repo}
\alias{dcs_load_ribocrypt}
\title{Find and load local RiboCrypt}
\usage{
dcs_ribocrypt_repo(start = getwd(), required = TRUE)

dcs_load_ribocrypt(repo = dcs_ribocrypt_repo(required = TRUE), quiet = TRUE)
}
\arguments{
\item{start}{Directory where repository discovery starts.}

\item{required}{Whether failure to find RiboCrypt should stop.}

\item{repo}{Path to the local RiboCrypt checkout.}

\item{quiet}{Passed to \code{devtools::load_all()}.}
}
\value{
The RiboCrypt repository path, invisibly for \code{dcs_load_ribocrypt()}.
}
\description{
These helpers enforce use of a local RiboCrypt checkout, normally supplied by
the \code{RIBOCRYPT_REPO} environment variable.
}
