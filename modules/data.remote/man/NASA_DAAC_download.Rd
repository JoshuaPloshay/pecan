% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NASA_DAAC_download.R
\name{NASA_DAAC_download}
\alias{NASA_DAAC_download}
\title{Parallel download data from the NASA ORNL DAAC server given period, spatial bounding box, and data DOI.}
\usage{
NASA_DAAC_download(
  ul_lat,
  ul_lon,
  lr_lat,
  lr_lon,
  ncore = 1,
  from,
  to,
  outdir,
  doi,
  netrc_file = NULL,
  just_path = FALSE
)
}
\arguments{
\item{ul_lat}{Numeric: upper left latitude.}

\item{ul_lon}{Numeric: upper left longitude.}

\item{lr_lat}{Numeric: lower right latitude.}

\item{lr_lon}{Numeric: lower right longitude.}

\item{ncore}{Numeric: numbers of core to be used if the maximum core}

\item{from}{Character: date from which the data search starts. In the form
"yyyy-mm-dd".}

\item{to}{Character: date on which the data search end. In the form
"yyyy-mm-dd".}

\item{outdir}{Character: path of the directory in which to save the
downloaded files.Default to the working directory. If it doesn't exist it
will be created. Ignored if \code{just_path=TRUE}}

\item{doi}{Character: data DOI on the NASA DAAC server, it can be obtained directly from the NASA ORNL DAAC data portal (e.g., GEDI L4A through https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=2056).}

\item{netrc_file}{Character: path to the credential file, default is NULL.}

\item{just_path}{Boolean: if we just want the metadata and URL or proceed the actual download.}
}
\value{
A list containing meta data for each data downloaded.
}
\description{
Parallel download data from the NASA ORNL DAAC server given period, spatial bounding box, and data DOI.
}
\examples{
\dontrun{
ul_lat <- 35
ul_lon <- -121
lr_lat <- 33
lr_lon <- -117
from <- "2022-02-23"
to <- "2022-05-30"
doi <- "10.3334/ORNLDAAC/2183"
metadata <- NASA_DAAC_download(ul_lat = ul_lat, 
                               ul_lon = ul_lon, 
                               lr_lat = lr_lat, 
                               lr_lon = lr_lon, 
                               from = from, 
                               to = to, 
                               doi = doi,
                               just_path = T)
}
}
\author{
Dongchen Zhang
}
