#'Main workhorse of makeDb
#'
#'Edited version of \code{read.csv.sql} in the \code{sqldf} package
#'(\code{tcltk} functionality removed), used mainly to create a SQLite database
#'from a text file. The function \code{makeDb} wraps this function in the
#'tornado package, so no need to execute this directly.
#'
#'This function is entirely wrapped by \code{makeDb} -- experienced users may
#'use this for debugging, but otherwise not usually necessary to call directly.
#'
#'@param file A file path or a URL (beginning with \code{http://} or
#'\code{ftp://}). If the \code{filter} argument is used and no file is to be
#'input to the filter then \code{file} can be omitted, \code{NULL}, \code{NA}
#'or \code{""}. The \code{textfile} argument from \code{makeDb} is used here.
#'@param sql character string holding an SQL statement. The table representing
#'the file should be referred to as \code{file}.
#'@param header as in \code{read.csv}
#'@param sep as in \code{read.csv}
#'@param row.names as in \code{read.csv}
#'@param eol Character that ends lines
#'@param skip Skip indicated number of lines in input file.
#'@param filter see \code{read.csv.sql} in \code{sqldf}
#'@param nrows Number of rows used to determine column types. It defaults to
#'50. Using -1 causes it to use all rows for determining column types. This
#'argument is rarely needed.
#'@param field.types A list whose names are the column names and whose contents
#'are the SQLite types (not the R class names) of the columns. Specifying these
#'types improves how fast it takes. Unless speed is very important this
#'argument is not normally used.
#'@param comment.char If specified this character and anything following it on
#'any line of the input will be ignored.
#'@param dbname As in \code{sqldf} except that the default is
#'\code{tempfile()}. Specifying NULL will put the database in memory which may
#'improve speed but will limit the size of the database by the available
#'memory. When using tornado, \code{tempfile()} is not used: \code{dbname} must
#'be provided as an argument to \code{makeDb}.
#'@param drv ignored: the only drive used can be SQLite.
#'@param \dots arguments passed to \code{sqldf}.
#'@return If the sql statement is a select statement then a data frame is
#'returned. In tornado, this is never the case.
#'@seealso \code{\link{makeDb}}
#'@references http://cran.r-project.org/web/packages/sqldf/sqldf.pdf
#'@export

read.csv.sql <- function (file, sql = "select * from file", header = TRUE, sep = ",", row.names, eol, skip, filter, nrows, field.types, comment.char, dbname = tempfile(), drv = "SQLite", ...) {
	require(proto)
    file.format <- list(header = header, sep = sep)
    if (!missing(eol)) 
        file.format <- append(file.format, list(eol = eol))
    if (!missing(row.names)) 
        file.format <- append(file.format, list(row.names = row.names))
    if (!missing(skip)) 
        file.format <- append(file.format, list(skip = skip))
    if (!missing(filter)) 
        file.format <- append(file.format, list(filter = filter))
    if (!missing(nrows)) 
        file.format <- append(file.format, list(nrows = nrows))
    if (!missing(field.types)) 
        file.format <- append(file.format, list(field.types = field.types))
    if (!missing(comment.char)) 
        file.format <- append(file.format, list(comment.char = comment.char))
    pf <- parent.frame()
    if (missing(file) || is.null(file) || is.na(file)) 
        file <- ""
    tf <- NULL
    if (substring(file, 1, 7) == "http://" || substring(file, 
        1, 6) == "ftp://") {
        tf <- tempfile()
        on.exit(unlink(tf), add = TRUE)
        download.file(file, tf, mode = "wb")
        file <- tf
    }
    p <- proto(pf, file = file(file))
    p <- do.call(proto, list(pf, file = file(file)))
    sqldf(sql, envir = p, file.format = file.format, dbname = dbname, 
        drv = drv, ...)
}
