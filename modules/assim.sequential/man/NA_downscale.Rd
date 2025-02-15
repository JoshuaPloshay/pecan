% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downscale_function.R
\name{NA_downscale}
\alias{NA_downscale}
\title{North America Downscale Function}
\usage{
NA_downscale(data, coords, date, C_pool, covariates)
}
\arguments{
\item{data}{In quotes, file path for .rds containing ensemble data.}

\item{coords}{In quotes, file path for .csv file containing the site coordinates, columns named "lon" and "lat".}

\item{date}{In quotes, if SDA site run, format is yyyy/mm/dd, if NEON, yyyy-mm-dd. Restricted to years within file supplied to 'data'.}

\item{C_pool}{In quotes, carbon pool of interest. Name must match carbon pool name found within file supplied to 'data'.}

\item{covariates}{SpatRaster stack, used as predictors in randomForest. Layers within stack should be named. Recommended that this stack be generated using 'covariates' instructions in assim.sequential/inst folder}
}
\value{
It returns the `downscale_output` list containing lists for the training and testing data sets, models, and predicted maps for each ensemble member.
}
\description{
This function uses the randomForest model.
}
\details{
This function will downscale forecast data to unmodeled locations using covariates and site locations
}
\author{
Joshua Ploshay
}
