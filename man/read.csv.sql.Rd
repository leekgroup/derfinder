\name{read.csv.sql}
\alias{read.csv.sql}
\title{Main workhorse of makeDb}
\usage{
  read.csv.sql(file, sql = "select * from file",
    header = TRUE, sep = ",", row.names, eol, skip, filter,
    nrows, field.types, comment.char, dbname = tempfile(),
    drv = "SQLite", ...)
}
\arguments{
  \item{file}{A file path or a URL (beginning with
  \code{http://} or \code{ftp://}). If the \code{filter}
  argument is used and no file is to be input to the filter
  then \code{file} can be omitted, \code{NULL}, \code{NA}
  or \code{""}. The \code{textfile} argument from
  \code{makeDb} is used here.}

  \item{sql}{character string holding an SQL statement. The
  table representing the file should be referred to as
  \code{file}.}

  \item{header}{as in \code{read.csv}}

  \item{sep}{as in \code{read.csv}}

  \item{row.names}{as in \code{read.csv}}

  \item{eol}{Character that ends lines}

  \item{skip}{Skip indicated number of lines in input
  file.}

  \item{filter}{see \code{read.csv.sql} in \code{sqldf}}

  \item{nrows}{Number of rows used to determine column
  types. It defaults to 50. Using -1 causes it to use all
  rows for determining column types. This argument is
  rarely needed.}

  \item{field.types}{A list whose names are the column
  names and whose contents are the SQLite types (not the R
  class names) of the columns. Specifying these types
  improves how fast it takes. Unless speed is very
  important this argument is not normally used.}

  \item{comment.char}{If specified this character and
  anything following it on any line of the input will be
  ignored.}

  \item{dbname}{As in \code{sqldf} except that the default
  is \code{tempfile()}. Specifying NULL will put the
  database in memory which may improve speed but will limit
  the size of the database by the available memory. When
  using tornado, \code{tempfile()} is not used:
  \code{dbname} must be provided as an argument to
  \code{makeDb}.}

  \item{drv}{ignored: the only drive used can be SQLite.}

  \item{\dots}{arguments passed to \code{sqldf}.}
}
\value{
  If the sql statement is a select statement then a data
  frame is returned. In tornado, this is never the case.
}
\description{
  Edited version of \code{read.csv.sql} in the \code{sqldf}
  package (\code{tcltk} functionality removed), used mainly
  to create a SQLite database from a text file. The
  function \code{makeDb} wraps this function in the tornado
  package, so no need to execute this directly.
}
\details{
  This function is entirely wrapped by \code{makeDb} --
  experienced users may use this for debugging, but
  otherwise not usually necessary to call directly.
}
\references{
  http://cran.r-project.org/web/packages/sqldf/sqldf.pdf
}
\seealso{
  \code{\link{makeDb}}
}

