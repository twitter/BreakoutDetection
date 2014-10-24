\name{breakout}
\alias{breakout}
\title{
	Breakout Detector: Energy Divisive with Medians
}
\description{
	A technique for robustly, i.e., in the presence of anomalies, detecting single or multiple change points in univariate time series.
}
\usage{
	breakout(Z, min.size = 30, method = 'amoc', ...)
}

\arguments{
	\item{Z}{
		The input time series. This is either a numeric vector or a data.frame which has 'timestamp' and 'count' components.
	}
	\item{min.size}{
		The minimum number of observations between change points.
	}
	\item{method}{
		Method must be one of either 'amoc' (At Most One Change) or 'multi' (Multiple Changes). For 'amoc' at most one change point 
		location will be returned.
	}
	\item{...}{
		See the details section for information about additional arguments.
	}
}

\details{
	The additional arguments that can be supplied depend upon whether single ('amoc') or multiple ('multi') change point 
	analysis is being performed.
	
	In both cases the following arguments are accepted:
	\itemize{
		\item \code{plot:} A flag indicating if a plot of the time series along with the estimated change point, indicated by vertical lines, should also be returned.
		\item \code{xlab:} X-axis label if plotting is enabled.
		\item \code{ylab:} Y-axis label if plotting is enabled.
		\item \code{title:} Plot title if plotting is enabled.
	}

	For singe change analysis, the following arguments are accepted:
	\itemize{
		\item \code{alpha:} 
		The alpha parameter used to weight the distance between observations. This is a real number in the interval (0,2]. The default value is \code{alpha=2}.
		\item \code{exact:}
		This flag is for selecting the use of true medians (\code{TRUE}) or approximate medians (\code{FALSE}) when determining change points. The default value is \code{exact=TRUE}.	
		\item \code{sig.lvl:}
		Once a change point is found its statistical significance is determined through a hypothesis test. \code{sig.lvl} specifies the significance for the hypothesis test. The default value is \code{sig.lvl=0.05}.	
		\item \code{nperm:}
		The number of permutations to perform in order to obtain an approximate p-value. If \code{nperm=0} then then permutation test is not
		performed. The default value is \code{nperm=0}.
		}
	
	For multiple change analysis, the following arguments are accepted:
	\itemize{
	\item \code{degree:}
		The degree of the penalization polynomial. \code{degree} can take the values 0, 1, and 2. The default value is \code{degree=1}.
	
	\item \code{beta:}
		A real numbered constant used to further control the amount of penalization. This is the default form of penalization, if neither (or both) \code{beta} or (and) \code{percent} are supplied this argument will be used. The default value is \code{beta=0.008}.
	
	\item \code{percent:}
		A real numbered constant used to further control the amount of penalization. This value specifies the minimum percent change in the goodness of fit statistic to consider adding an additional change point. A value of 0.25 corresponds to a 25\% increase. \code{percent} doesn't have a default value.
		}
	
}

\value{
	The returned value is a list with the following components.
	\item{loc}{
		The estimated change point location(s).
	}
	\item{time}{
		The amount of required processing time, in seconds.
	}
	\item{pval}{
		The approximate p-value obtained from the permutation test. When \code{nperm=0} or \code{method='multi'}, \code{pval=NA}.
	}
	\item{plot}{
		A ggplot graphical object if plotting was requested by the user. The supplied image is of the input time series along with 
		the estimated change point location(s).
	}
}

\references{
Nicholas A. James, Arun Kejariwal, David S. Matteson, "Leveraging Cloud Data to Mitigate User Experience from 'Breaking Bad': The Twitter Approach, 2014
}

\examples{
	library(BreakoutDetection)
	data(Scribe)
	res = breakout(Scribe, min.size=24, method='multi', beta=.001, degree=1, plot=TRUE)
	res$plot
}
