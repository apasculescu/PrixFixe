## no GPF dependencies.
initUTIL <- function()
  {
    PUBLIC <- list()


    PUBLIC$makeUnique <- function(x, sep = "_")
      {
        x <- as.character(x)
        counts <- tapply(x, x, length)
        counts_gt_1 <- counts[counts > 1]
        for(duplicated_string in names(counts_gt_1)) {
          ix <- which(x == duplicated_string)
          x[ix] <- sprintf("%s%s%s", x[ix], sep, 1:counts_gt_1[duplicated_string])
        }
        return(x)
      }


    return(PUBLIC)
  }
