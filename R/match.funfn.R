#'Generic extended version of R match.fun
#'
#'A generic \code{match.fun}.
#'
#'The default method is the same as \code{match.fun} and the \code{formula}
#'method is the same as \code{as.function.formula}.  This function can be used
#'within the body of a function to convert a function specification whether its
#'a function, character string or formula into an actual function.
#'
#'@param FUN Function, character name of function or formula describing
#'function.
#'@param descend logical; control whether to search past non-function objects.
#'@return Returns a function.
#'@note This is exactly the same as \code{match.funfn} in the \code{gsubfn}
#'package.  Do not load the \code{gsubfn} package to use this function, as
#'\code{gsubfn} loads \code{tcltk}, which is not advised.
#'@author G. Grothendieck
#'@export

match.funfn <- function (FUN, descend = TRUE) {
    if (is.function(FUN)) 
        return(FUN)
    if (inherits(FUN, "formula")) 
        return(as.function.formula(FUN))
    if (!(is.character(FUN) && length(FUN) == 1 || is.symbol(FUN))) {
        FUN <- eval.parent(substitute(substitute(FUN)))
        if (!is.symbol(FUN)) 
            stop(gettextf("'%s' is not a function, character or symbol", 
                deparse(FUN)), domain = NA)
    }
    envir <- parent.frame(2)
    if (descend) 
        FUN <- get(as.character(FUN), mode = "function", envir = envir)
    else {
        FUN <- get(as.character(FUN), mode = "any", envir = envir)
        if (!is.function(FUN)) 
            stop(gettextf("found non-function '%s'", FUN), domain = NA)
    }
    return(FUN)
}
