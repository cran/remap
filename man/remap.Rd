% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remap.R
\name{remap}
\alias{remap}
\title{Build separate models for mapping multiple regions.}
\usage{
remap(
  data,
  regions,
  region_id,
  model_function,
  buffer,
  min_n = 1,
  distances,
  cores = 1,
  progress = FALSE,
  ...
)
}
\arguments{
\item{data}{An sf data frame with point geometry.}

\item{regions}{An sf dataframe with polygon or multipolygon geometry.}

\item{region_id}{Optional name of column in 'regions' that contains the id
that each region belongs to (no quotes). If null, it will be assumed that
each row of 'regions' is its own region.}

\item{model_function}{A function that can take a subset of 'data' and
output a model that can be used to predict new values when passed to generic
function predict().}

\item{buffer}{The length of the buffer zone around each region in km where
observations are included in the data used to build models for each region.
(Can be a named vector with different values for each unique 'region_id' in
'region'.)}

\item{min_n}{The minimum number of observations to use when building a model.
If there are not enough observations in the region and buffer, then the
closest min_n observations are used. \code{min_n} must be at least 1.}

\item{distances}{An optional matrix of distances between 'data' and 'regions'
generated by \code{redist()} function (calculated internally if not
provided). Note that unless you know that you have min_n within a certain
distance, no max_dist parameter should be used in \code{redist()}.}

\item{cores}{Number of cores for parallel computing. 'cores' above
default of 1 will require more memory.}

\item{progress}{If true, a text progress bar is printed to the console.
(Progress bar only appears if 'cores' = 1.)}

\item{...}{Extra arguments to pass to 'model_function' function.}
}
\value{
A \emph{remap} S3 object containing:
\describe{
  \item{\emph{models}}{A list of models containing a model output by
  'model_function' for each region.}
  \item{\emph{regions}}{'regions' object passed to the function (used for
  prediction). The first column is 'region_id' or the row number of 'regions'
  if 'region_id is missing. The second column is the region geometry.}
  \item{\emph{call}}{Shows the parameters that were passed to the function.}
}
}
\description{
Separate models are built for each given region and combined into one S3
object that can be used to predict on new data using generic function
predict().
}
\details{
If a model fails for a region, a warning is given but the modeling
process will continue.

A description of the methodology can be found in Wagstaff and Bean (2023)
"remap: Regionalized Models with Spatially Smooth Predictions"
<doi:10.32614/RJ-2023-004>.
}
\examples{
library(remap)
data(utsnow)
data(utws)

# We will keep these examples simple by only modeling non-zero values of
# snow water equivalent (WESD)

utsnz <- utsnow[utsnow$WESD > 0, ]

# Build a remap model with lm that has formula WESD ~ ELEVATION
# The buffer to collect data around each region is 30km
# The minimum number of observations per region is 10
remap_model <- remap(data = utsnz,
                     regions = utws,
                     region_id = HUC2,
                     model_function = lm,
                     formula = log(WESD) ~ ELEVATION,
                     buffer = 20,
                     min_n = 10,
                     progress = TRUE)

# Resubstitution predictions
remap_preds <- exp(predict(remap_model, utsnz, smooth = 10))
head(remap_preds)

}
\seealso{
\code{\link{predict.remap}} - used for predicting on new data.
  \code{\link{redist}} - used for pre-computing distances.
}
