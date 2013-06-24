#'Helper function for read.csv.sql
#'
#'The same as \code{strapply} in the \code{gsubfn} package, but with
#'\code{tcltk} capabilities removed.
#'
#'See details in \code{gsubfn} package.
#'
#'@param X list or (atomic) vector of character strings to be used
#'@param pattern character string containing a regular (or character string for
#'\code{fixed=TRUE} to be matched in the given character vector.
#'@param FUN a function, formula, character string, list or proto object to be
#'applied to each element of \code{X}. See discussion in \code{\link{gsubfn}.}
#'@param backref see \code{gsubfn}
#'@param \dots optional arguments to \code{gsubfn}
#'@param empty If there is no match to a string return this value.
#'@param ignore.case If TRUE then case is ignored in the \code{pattern}
#'argument.
#'@param perl If TRUE then \code{engine="R"} is used with perl regular
#'expressions.  It is required to keep this argument at TRUE, since \code{tcl}
#'engine capabilities have been removed from this function.
#'@param engine Should always be set to \code{"R"}, since the \code{tcl} engine
#'is not available in the tornado package.
#'@param simplify logical or function. If logical, should the result be
#'simplified to a vector or matrix, as in \code{sapply} if possible? If
#'function, that function is applied to the result with each component of the
#'result passed as a separate argument. Typically if the form is used it will
#'typically be specified as rbind.
#'@param USE.NAMES logical; if \code{TRUE} and if \code{X} is character, use
#'\code{X} as 'names' for the result unless it had names already.
#'@param combine combine is a function applied to the components of the result
#'of \code{FUN}. The default is \code{"c"}. \code{"list"} is another common
#'choice. The default may change to be \code{"list"} in the future.
#'@return A list of character strings
#'@note Does not need to be used directly in tornado; \code{makeDb} wraps this
#'entirely.
#'@author G. Grothendieck
#'@export
#'@seealso \code{\link{makeDb}}
#'@references http://cran.r-project.org/web/packages/gsubfn/gsubfn.pdf

strapply <- function (X, pattern, FUN = function(x, ...) x, backref = NULL, ..., empty = NULL, ignore.case = FALSE, perl = TRUE, engine = "R", simplify = FALSE, USE.NAMES = FALSE, combine = c) {
    combine <- match.funfn(combine)
    stopifnot(!missing(pattern))
    pattern <- as.character(pattern)
    if (is.proto(FUN) || perl) 
        engine <- "R"
    if (identical(engine, "R")) 
        return(ostrapply(X = X, ignore.case = ignore.case, pattern = pattern, 
            FUN = FUN, backref = backref, ..., empty = empty, 
            perl = perl, simplify = simplify, USE.NAMES = USE.NAMES, 
            combine = combine))
    if (is.proto(FUN)) {
    }
    else if (is.character(FUN)) {
        FUN.orig <- FUN
        FUN <- function(...) FUN.orig
    }
    else if (is.list(FUN)) {
        values.replacement <- FUN
        names.replacement <- names(FUN)
        FUN <- function(...) {
            idx <- match(..1, names.replacement, nomatch = match("", 
                names.replacement, nomatch = 0))
            if (idx > 0) 
                values.replacement[[idx]]
            else ..1
        }
    }
    else {
        FUN <- match.funfn(FUN)
    }
    ff <- function(x) {
        s <- strapply1(x, pattern, backref, ignore.case)
        if (length(s) == 0 && !is.null(empty)) 
            s <- matrix(empty, 1)
        L <- lapply(seq_len(ncol(s)), function(j) {
            combine(do.call(FUN, as.list(s[, j])))
        })
        do.call("c", L)
    }
    result <- sapply(X, ff, simplify = is.logical(simplify) && 
        simplify, USE.NAMES = USE.NAMES)
    if (is.logical(simplify)) 
        result
    else do.call(match.funfn(simplify), result)
}
