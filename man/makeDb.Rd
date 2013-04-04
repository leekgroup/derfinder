\name{makeDb}
\alias{makeDb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ create SQLite database from text file
}
\description{
Dumps the contents of a table (saved as a text file) into a SQLite database, performing some filtering along the way.
}
\usage{
makeDb(dbfile, textfile, tablename, sep = "\t", cutoff = 5)
}
\arguments{
  \item{dbfile}{
Character string giving the file name/location of the database to be created.  Generally ends in \code{.db}. 
}
  \item{textfile}{
The text file containing the table to be dumped into \code{dbfile}.
}
  \item{tablename}{
Character string containing name to give the table inside \code{dbfile}.
}
  \item{sep}{
The separator used in \code{textfile}.  The tornado pipeline creates tab-separated text files, so \code{"\t"} is the default.
}
  \item{cutoff}{
Rows in \code{textfile} must have at least one entry (not counting the first column, which is assumed to hold genomic position) greater than \code{cutoff} to be included in \code{dbfile}.
}
}
\value{
No return, but writes the file \code{dbfile} containing table \code{tablename} by filtering \code{textfile} according to \code{cutoff}.
}
\author{
Alyssa Frazee
}
\note{
The workhorse of this function is a modified version of \code{\link{read.csv.sql}}, found in the \code{sqldf} package.
}
