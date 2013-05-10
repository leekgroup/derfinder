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
\description{
  Used internally by read.csv.sql, which drives the
  \code{makeDb} function. Not necessary to call this
  function directly when using the tornado package.
}
\details{
  For arguments, value, and other information, see
  \code{sqldf} - this function is a direct copy of that
  function.
}
\references{
  http://cran.r-project.org/web/packages/sqldf/sqldf.pdf
}
\seealso{
  \code{makeDb}
}

