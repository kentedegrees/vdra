\name{pmn}
\alias{pmn}
\title{
	PopMedNet Simulator
}
\description{
  This function is intended to act as a proxy for PopMedNet when developing code to run on PopMedNet, or when testing out the distributed regression programs provided with this package.  When used, it is expected that this function, the analysis center, and the data partner(s) will each be run in their own instances of R.  The analysis center and the data partner(s) will share a common directory where the subdirectories dp0, dp1, ... will be stored.  The directory dp0 is the monitor folder for the analysis center.  With the exception of 2-party regression, it is assumed that the data partner using dp1 is also the data partner which holds the response.  In the case of 2-party regression, the analysis center holds the response.
}
\usage{
  pmn(numDataPartners, directory = getwd())
}
\arguments{
  \item{numDataPartners}{the number of data partners involved in the multiple regression.  If a data partner is also acting as the analysis center, then that data partner is not counted.}
	\item{directory}{the directory where the directories dp0, dp1, ... are located, which are used by the analysis center and data partner(s) to save data and receive data from each other.}
}
\seealso{
  \code{\link{AnalysisCenter.2Party}}, \code{\link{AnalysisCenter.3Party}}, \code{\link{AnalysisCenter.KParty}}
}
