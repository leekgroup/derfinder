#'Helper function for read.csv.sql
#'
#'Same as the \code{strapply} or \code{ostrapply} functions in \code{gsubfn}
#'package, but with \code{tcltk} capabilities removed.
#'
#'See \code{strapply} in \code{gsubfn}.  In the tornado package this function
#'should rarely be called on its own.  It is an internal process of
#'\code{read.csv.sql}, which is an internal process of \code{makeDb}.
#'
#'@param X list or (atomic) vector of character strings to be used
#'@param pattern character string containing a regular (or character string for
#'\code{fixed=TRUE} to be matched in the given character vector.
#'@param FUN a function, formula, character string, list or proto object to be
#'applied to each element of \code{X}. See discussion in \code{\link{gsubfn}.}
#'@param ignore.case If \code{TRUE} then case is ignored in the \code{pattern}
#'argument.
#'@param \dots optional arguments to \code{gsubfn}.
#'@param empty If there is no match to a string return this value.
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
#'@return A list of character strings.
#'@author G. Grothendieck
#'@export
#'@seealso \code{\link{makeDb}},\code{link{read.csv.sql}}
#'@references http://cran.r-project.org/web/packages/gsubfn/gsubfn.pdf

ostrapply <- function (X, pattern, FUN = function(x, ...) x, ignore.case = FALSE, ..., empty = NULL, simplify = FALSE, USE.NAMES = FALSE, combine = c) {
    here <- environment()
    combine <- match.funfn(combine)
    if (is.character(FUN)) {
        FUN.orig <- FUN
        FUN <- function(...) FUN.orig
    }
    else if (is.list(FUN)) {
        values.replacement <- FUN
        names.replacement <- names(FUN)
        here$FUN <- function(...) {
            idx <- match(..1, names.replacement, nomatch = match("", 
                names.replacement, nomatch = 0))
            if (idx > 0) 
                values.replacement[[idx]]
            else ..1
        }
    }
    p <- if (is.proto(FUN)) {
        FUN$X <- X
        FUN$pattern <- pattern
        FUN$simplify <- simplify
        FUN$USE.NAMES <- USE.NAMES
        FUN$combine <- combine
        proto(pre = function(this) {
            this$first <- TRUE
            this$v <- NULL
            if (!is.null(FUN[["pre"]])) 
                FUN$pre()
        }, fun = function(this, ...) {
            FUN$count <- this$count
            this$v <- if (this$first) 
                combine(FUN$fun(...))
            else c(this$v, combine(FUN$fun(...)))
            this$first <- FALSE
        }, post = function(this) {
            if (this$first) 
                this$v <- NULL
            if (!is.null(FUN[["post"]])) 
                FUN$post()
        }, )
    }
    else {
        FUN <- match.funfn(FUN)
        proto(pre = function(this) {
            this$first <- TRUE
            this$v <- NULL
        }, fun = function(this, ...) {
            this$v <- if (this$first) 
                combine(FUN(...))
            else c(this$v, combine(FUN(...)))
            this$first <- FALSE
        }, post = function(this) {
            if (this$first) 
                this$v <- NULL
        })
    }
    ff <- function(x, ...) {
        gsubfn(pattern, p, x, ignore.case = ignore.case, 
            ...)
        p$v
    }
    result <- sapply(X, ff, ..., simplify = isTRUE(simplify), 
        USE.NAMES = USE.NAMES)
    if (is.logical(simplify)) 
        result
    else {
        do.call(match.funfn(simplify), result)
    }
}
