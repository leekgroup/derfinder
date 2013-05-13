#'Helper function for read.csv.sql
#'
#'see \code{gsubfn} in \code{sqldf} package - this function is equivalent, but
#'functionality requiring \code{tcltk} has been removed.
#'
#'If \code{replacement} is a string then it acts like \code{gsub}. If
#'\code{replacement} is a function then each matched string is passed to the
#'replacement function and the output of that function replaces the matched
#'string in the result. The ﬁrst argument to the replacement function is the
#'matched string and subsequent arguments are the backreferences, if any. If
#'\code{replacement} is a list then the result of the regular expression match
#'is, in turn, matched against the names of that list and the value
#'corresponding to the ﬁrst name in the list that is match is returned. If
#'there are no names matching then the ﬁrst unnamed component is returned and
#'if there are no matches then the string to be matched is returned. If
#'\code{backref} is not speciﬁed or is specified and is positive then the
#'entire match is used to lookup the value in the list whereas if
#'\code{backref} is negative then the identiﬁed backreference is used. If
#'\code{replacement} is a formula instead of a function then a one line
#'function is created whose body is the right hand side of the formula and
#'whose arguments are the left hand side separated by + signs (or any other
#'valid operator). The environment of the function is the environment of the
#'formula. If the arguments are omitted then the free variables found on the
#'right hand side are used in the order encountered. 0 can be used to indicate
#'no arguments. \code{letters}, \code{LETTERS} and \code{pi} are never
#'automatically used as arguments. If \code{replacement} is a proto object then
#'it should have a \code{fun} method which is like the replacement function
#'except its ﬁrst argument is the object and the remaining arguments are as in
#'the replacement function and are affected by \code{backref} in the same way.
#'\code{gsubfn} automatically inserts the named arguments in the call to
#'\code{gsubfn} into the proto object and also maintains a \code{count}
#'variable which counts matches within strings. The user may optionally specify
#'\code{pre} and \code{post} methods in the proto object which are ﬁred at the
#'beginning and end of each string (not each match). They each take one
#'argument, the object. Note that if \code{backref} is non-negative then
#'internally the pattern will be parenthesized. A utility function \code{cat0}
#'is available. They are like \code{\link{cat}} and \code{\link{paste}} except
#'that their default sep value is "".
#'
#'@param pattern Same as \code{pattern} in \code{gsub}
#'@param replacement A character string, function, list, formula or proto
#'object. See Details.
#'@param x Same as \code{x} in \code{gsub}
#'@param backref Number of backreferences to be passed to function. If zero or
#'positive the match is passed as the ﬁrst argument to the replacement function
#'followed by the indicated number of backreferences as subsequent arguments.
#'If negative then only the that number of backreferences are passed but the
#'match itself is not. If omitted it will be determined automatically, i.e. it
#'will be 0 if there are no backreferences and otherwise it will equal negative
#'the number of back references. It determines this by counting the number of
#'non-escaped left parentheses in the pattern. Also if the function contains an
#'ampersand as an argument then \code{backref} will be taken as non-negative
#'and the ampersand argument will get the full match.
#'@param USE.NAMES See \code{USE.NAMES} in \code{sapply}
#'@param ignore.case If \code{TRUE} then case is ignored in the \code{pattern}
#'argument.
#'@param env Environment in which to evaluate the replacement function.
#'Normally this is left at its default value.
#'@param \dots Other \code{gsub} arguments
#'@return as in \code{gsub}
#'@author G. Grothendieck
#'@export
#'@seealso \code{\link{strapply}}

gsubfn <- function (pattern, replacement, x, backref, USE.NAMES = FALSE, ignore.case = FALSE, env = parent.frame(), ...) {
	engine <- "R"
    R.engine <- identical(engine, "R")
    here <- environment()
    if (missing(replacement)) 
        here$replacement <- function(...) eval(parse(text = paste(..., 
            sep = "")), env)
    if (is.character(replacement)) {
        if (R.engine) 
            return(base::gsub(pattern, replacement, x, ...))
        else { stop("trying to use tcltk.  don't do that.")        }
    }
    if (is.list(replacement)) {
        values.replacement <- replacement
        names.replacement <- names(replacement)
        here$replacement <- function(...) {
            idx <- match(..1, names.replacement, nomatch = match("", 
                names.replacement, nomatch = 0))
            if (idx > 0) 
                values.replacement[[idx]]
            else ..1
        }
    }
    if (missing(pattern)) 
        pattern <- "[$]([[:alpha:]][[:alnum:].]*)|`([^`]+)`"
    pattern <- as.character(pattern)
    e <- NULL
    if (!inherits(replacement, "formula") && !is.function(replacement)) {
        e <- replacement
        e$pattern <- pattern
        e$x <- x
        e$backref <- if (!missing(backref)) 
            backref
        e$USE.NAMES <- USE.NAMES
        e$env <- env
        dots <- list(...)
        if (!is.null(names(dots))) {
            nam <- names(dots)
            for (n in nam[nam != ""]) assign(n, dots[[n]], e)
        }
        e$replacement <- function(this, ...) {
            this$count <- this$count + 1
            this$match <- c(...)
            this$fun(...)
        }
        here$replacement <- e$replacement
    }
    here$replacement <- match.funfn(replacement)
    if (missing(backref) || is.null(backref)) {
        noparen <- base::gsub("\\\\.", "", pattern)
        noparen <- base::gsub("\\[[^\\]]*\\]", "", noparen)
        backref <- -nchar(base::gsub("[^(]", "", noparen))
    }
    if (names(formals(here$replacement))[[1]] == "&") {
        backref <- abs(backref)
        if (!is.null(e)) 
            e$backref <- backref
    }
    j <- (identical(engine, "R") && !is.null(backref) && backref >= 
        0) + abs(backref)
    i <- if (!R.engine && backref >= 0) 
        0
    else 1
    j <- max(i, j)
    stopifnot(is.character(pattern), is.character(x), is.function(replacement))
    force(env)
    gsub.function <- function(x) {
        if (R.engine && !is.null(backref) && backref >= 0) {
            pattern <- paste("(", pattern, ")", sep = "")
        }
        if (!is.null(e)) {
            e$count <- 0
            if ("pre" %in% ls(e)) 
                e$pre()
        }
        repl <- function(i, j) {
            rs <- paste("\\", seq(i, j), collapse = "\002", sep = "")
            rs <- paste("\001\002", rs, "\001", sep = "")
            if (R.engine) 
                tryCatch(base::gsub(pattern, rs, x, ignore.case = ignore.case, 
                  ...), error = function(x) if (j > i) 
                  repl(i, j - 1)
                else stop(x))
            else { stop("trying to use tcltk. don't do that (2)")}
        }
        z <- repl(i, j)
        z <- strsplit(z, "\001")[[1]]
        f <- function(s) {
            if (nchar(s) > 0 && substring(s, 1, 1) == "\002") {
                s <- sub("\002$", "\002\002", s)
                L <- as.list(strsplit(s, "\002")[[1]][-1])
                do.call(replacement, L)
            }
            else s
        }
        z <- paste(sapply(z, f), collapse = "")
        if (!is.null(e) && "post" %in% ls(e)) 
            e$post()
        z
    }
    sapply(x, gsub.function, USE.NAMES = USE.NAMES)
}
