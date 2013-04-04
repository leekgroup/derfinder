\name{getSvs}
\alias{getSvs}
\title{
do SVA on matrix stored in database
}
\description{
performs surrogate variable analysis on a matrix, usually base-pair by sample coverage, so the surrogate variables can be adjusted for in later analysis steps.  SVA is meant to provide a method to adjust for unknown confounders.
}
\usage{
getSvs(dbfile, tablename, group, chunksize = 1e+05, colsubset = c(-1))
}
\arguments{
  \item{dbfile}{
Name/location (as character string) of database (usually ".db") file containing nucleotide by sample coverage.  
}
  \item{tablename}{
Name of the table the database contains
}
  \item{group}{
a 0/1 vector grouping the samples (columns) in the database. 
}
  \item{chunksize}{
How many rows of the database should be processed at a time? 
}
  \item{colsubset}{
Optional vector of column indices of the input file that denote samples you wish to include in analysis. Should NEVER include 1 (genomic position).
}
}
\details{
}
\value{
An object of class sva - see help files for sva.
}
\references{
}
\author{
Jeff Leek
}
\note{
}

