#'Get the last element
#'
#'The function returns the last element in a given object (any object with a \code{tail} method
#'will work).
#'
#'
#'@param x any object with a \code{tail} method
#'@return same as \code{x[length(x)]}
#'@note Helper function for \code{getParams}
#'@author Alyssa Frazee
#'@export
#'@examples
#'
#'x = c(1:20)
#'last(x) == 20
#'
last <- function(x) return(tail(x, n=1))
