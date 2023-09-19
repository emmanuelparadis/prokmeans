\name{kmeansConstr}
\alias{kmeansConstr}
\title{K-Means With Constrained Observations}
\description{
  Performs k-means classification with contrained memberships of
  some observations which are not reclassified but included in the
  calculations of the parameters (multivariate center of each class).
}
\usage{
kmeansConstr(x, cluster, K, iter.lim = 100, restart.max = 10,
             diag4restart = 10, quiet = FALSE)
}
\arguments{
  \item{x}{a numeric matrix (or vector which is converted as a
    one-column matrix).}
  \item{cluster}{a vector of class membership for each observation in
    \code{x} (see details).}
  \item{K}{an integer giving the number of groups (or clusters); should
    not be less than \code{max(cluster)}.}
  \item{iter.lim}{an integer giving the maximum number of iterations.}
  \item{restart.max}{the maximum number of restarts (see details).}
  \item{diag4restart}{the number of consecutive times when the number of
    reclassified observations is not decreasing before the whole
    procedure is restarted (see details).}
  \item{quiet}{a logical value specifying whether the progress of the
    computations should be printed or not.}
}
\details{
  The vector \code{cluster} should be with integer values giving the
  memberships of constrained observations (using the values 1, 2, and so
  on). The observations with unknown membership can be \code{NA}, 0, or
  any negative value.

  The value given to \code{K} can be greater than the number of classes
  observed in \code{cluster}.

  The algorithm is very close to the classical k-means procedure: the
  observations with unknown membership are assigned to one of the
  \code{K} classes randomly, the (multivariate) centers of each class
  are computed, and these observations are reassigned to their
  respective closest class. These steps are repeated \code{iter.lim}
  times or until no observation is reclassified.

  The number of reclassified observations (i.e., those that changed
  membership at the end of the above steps) is monitored at each
  iteration: if this number is not decreasing during \code{diag4restart}
  consecutive iterations, the whole procedure is restarted with new
  random class assignments. This restart is done at most
  \code{restart.max} times.

  If some observations are still reclassified after \code{iter.lim}
  iterations but no restart has been done, it certainly means that the
  algorithm has not yet converged fully and it may be useful to increase
  the value of this argument.
}
\value{a vector of integers giving the final group assignments.}
\author{Emmanel Paradis}
\seealso{\code{\link{pkm}}}
\keyword{models}
