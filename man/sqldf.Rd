\name{sqldf}
\alias{sqldf}
\title{Helper function for read.csv.sql}
\usage{
  sqldf(x, stringsAsFactors = FALSE, row.names = FALSE,
    envir = parent.frame(),
    method = getOption("sqldf.method"),
    file.format = list(), dbname,
    drv = getOption("sqldf.driver"), user, password = "",
    host = "localhost", port, dll = getOption("sqldf.dll"),
    connection = getOption("sqldf.connection"),
    verbose = isTRUE(getOption("sqldf.verbose")))
}
\arguments{
  \item{x}{As in \code{sqldf}.}

  \item{stringsAsFactors}{As in \code{sqldf}.}

  \item{row.names}{As in \code{sqldf}.}

  \item{envir}{As in \code{sqldf}.}

  \item{method}{As in \code{sqldf}.}

  \item{file.format}{As in \code{sqldf}.}

  \item{dbname}{As in \code{sqldf}.}

  \item{drv}{As in \code{sqldf}.}

  \item{user}{As in \code{sqldf}.}

  \item{password}{As in \code{sqldf}.}

  \item{host}{As in \code{sqldf}.}

  \item{port}{As in \code{sqldf}.}

  \item{dll}{As in \code{sqldf}.}

  \item{connection}{As in \code{sqldf}.}

  \item{verbose}{As in \code{sqldf}.

  For arguments, value, and other information, see
  \code{sqldf} - this function is a direct copy of that
  function.}
}
\description{
  Used internally by read.csv.sql, which drives the
  \code{makeDb} function. Not necessary to call this
  function directly when using the tornado package.
}
\references{
  http://cran.r-project.org/web/packages/sqldf/sqldf.pdf
}
\seealso{
  \code{makeDb}
}

