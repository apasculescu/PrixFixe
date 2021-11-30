NA2 <- function(x, replacement)
  {
    x[is.na(x)] <- replacement
    x
  }
Inf2 <- function(x, replacement)
  {
    x[is.infinite(x)] <- replacement
    x
  }
#posInf2 <- function(x, replacement)
#  {
#    x[x == +Inf] <- replacement
#  }
#negInf2 <- function(x, replacement)
#  {
#    x[x == -Inf] <- replacement
#  }
NaN2 <- function(x, replacement)
  {
    x[is.nan(x)] <- replacement
    x
  }


TODAY <- function() format(Sys.Date(), "%F")


inf2NA <- function(x)
  {
    warning("'inf2NA(x)' is a legacy function. consider using Inf2(x, NA) as a replacement.",
            call. = TRUE,
            immediate. = FALSE)
    Inf2(x, NA)
  }


"%ni%" <- Negate("%in%")


## treats x and y as character vectors (through coersion if needed), thus the sorting is "lexicographic".
## NOTE THAT THIS MEANS (if decreasing = FALSE) 50 < 8.
## delegates to paste().
## ... arguments are passed as such to paste().
## the only really useful ones are sep, which are hard-coded with defaults here.
## recycles the shorter of the two vectors with a warning if x and y aren't the same length.
pairAndSort <- function(x, y = NULL, sep = " ", decreasing = FALSE, ...)
  {
    if(is.null(y)) {
      stopifnot(is.matrix(x) || is.data.frame(x))
      y <- x[,2]
      x <- x[,1] ## overwrites the original 'x', which is why the 'y' assignment must come first.
    }
    
    if(length(x) != length(y)) {
      warning("x and y are not equal in length, recycling the shorter.")
    }

    x <- as.character(x)
    y <- as.character(y)

    if(decreasing == FALSE) {
      paste(pmin(x, y), pmax(x, y), sep = sep, ...)
    } else {
      paste(pmax(x, y), pmin(x, y), sep = sep, ...)
    }
  }


linearRescale <- function(x, original_range = range(x, na.rm = TRUE), scaled_range = c(0, 1))
  {
    if(any(is.na(x))) warning("NAs detected.")

    (x - min(original_range)) * abs(diff(scaled_range) / diff(original_range)) + min(scaled_range)
  }



rotateMatrix <- function(m, nturns = 1, clockwise = TRUE)
  {
    stopifnot(is.matrix(m))
    stopifnot(nturns >= 0)
    stopifnot(is.logical(clockwise))

    nturns <- nturns %% 4 ## modulo operator, nturns is now 0, 1, 2, or 3.
    if(nturns == 0) return(m) ## nturns can now be only 1, 2, or 3.
    ## 3 turns in one direction is the same as 1 turn in the opposite direction.
    if(nturns == 3) {
      nturns == 1
      clockwise = !clockwise
    }

    clock90 <- function(m) t(m[nrow(m):1,])
    counterclock90 <- function(m) t(m)[ncol(m):1,]
    if(clockwise)
      for(i in 1:nturns) m <- clock90(m)
    else
      for(i in 1:nturns) m <- counterclock90(m)

    return(m)
  }
