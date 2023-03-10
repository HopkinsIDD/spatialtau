##' Generalized version of \code{get.pi}
##'
##' Generalized version of the \code{get.pi} function that takes in an arbitrary function and
##' returns the probability that a point within a particular range of a point of interest shares the relationship
##' specified by the passed in function with that point.
##'
##' @param posmat a matrix with columns x, y and any other named columns
##'    columns needed by fun
##' @param fun a function that takes in two rows of posmat and returns:
##' \enumerate{
##'      \item  for pairs included in the numerator and denominator
##'      \item for pairs that should only be included in the denominator
##'      \item for pairs that should be ignored all together}
##' Note that names from \code{posmat} are not perserved in calls to \code{fun}, so the columns of the matrix should be
##' referenced numerically
##' so this is not available to the fun
##' @param r the series of spatial distances (or there maximums) we are
##'          interested in
##' @param r.low the low end of each range, 0 by default
##'
##' @return  pi value for each distance range that we look at. Where:
##'
##' \deqn{ \pi(d_1,d_2) = \frac{\sum \boldsymbol{1} (d_{ij} \in (d_1,d_2)) \boldsymbol{1} (f(i,j)=1) }{\sum \sum \boldsymbol{1} (d_{ij} \in (d_1,d_2)) \boldsymbol{1} (f(i,j) \in \{1,2\}) }}
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.pi
##' @family spatialtau
##'
##' @example example/get_pi.R
##'
get.pi <- function(posmat,
                   fun,
                   r = 1,
                   r.low=rep(0,length(r))) {

    xcol <-  which(colnames(posmat)=="x")
    ycol <- which(colnames(posmat)=="y")

    #check that both columns exist
    if (length(xcol)!=1 & length(ycol)!=1) {
        stop("unique x and y columns must be defined")
    }

    rc <- .Call("get_pi",
                posmat,
                fun,
                r,
                r.low,
                1:nrow(posmat),
                xcol,
                ycol)
    return(rc)
}

##'
##' Optimized version of \code{get.pi} for typed data.
##'
##' Version of the \code{get.pi} function that is optimized for statically typed data. That is
##' data where we are interested in the probability of points within some distance of points of
##' typeA are of typeB.
##'
##' @param posmat a matrix with columns type, x, y
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##'
##' @return pi values for all the distances we looked at
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.pi
##'
##' @example example/get_pi_typed.R
##'
get.pi.typed <- function(posmat,
                         typeA = -1,
                         typeB = -1,
                         r=1,
                         r.low=rep(0,length(r))) {
  
  return(.C("get_pi_typed",
            as.integer(posmat[,"type"]),
            as.double(posmat[,"x"]),
            as.double(posmat[,"y"]),
            as.integer(nrow(posmat)),
            as.integer(typeA),
            as.integer(typeB),
            as.double(r.low),
            as.double(r),
            as.integer(length(r)),
            as.integer(1:nrow(posmat)),
            rc=double(length(r))
  )$rc)
}


##'
##' Optimized version of \code{get.pi} for typed data.
##'
##' Version of the \code{get.pi} function that is optimized for statically typed data. That is
##' data where we are interested in the probability of points within some distance of points of
##' typeA are of typeB.
##'
##' @param posmat a matrix with columns type, x, y, and weight
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##'
##' @return pi values for all the distances we looked at
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.pi
##'
##' @example example/get_pi_typed.R
##'
get.pi.typed.wts <- function(posmat,
                         typeA = -1,
                         typeB = -1,
                         r=1,
                         r.low=rep(0,length(r))) {

    return(.C("get_pi_typed_wts",
              as.integer(posmat[,"type"]),
              as.double(posmat[,"x"]),
              as.double(posmat[,"y"]),
              as.double(posmat[,"weight"]),
              as.integer(nrow(posmat)),
              as.integer(typeA),
              as.integer(typeB),
              as.double(r.low),
              as.double(r),
              as.integer(length(r)),
              as.integer(1:nrow(posmat)),
              rc=double(length(r))
              )$rc)
}


##'
##' Optimized version of \code{get.pi} for typed data.
##'
##' Version of the \code{get.pi} function that is optimized for statically typed data. That is
##' data where we are interested in the probability of points within some distance of points of
##' typeA are of typeB.
##'
##' @param posmat a matrix with columns type, x, y, and weight
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##'
##' @return pi values for all the distances we looked at
##'
##' @author Shaun Truelove
##'
##' @family get.pi
##'
##' @example example/get_pi_typed.R
##'
get.pi.clustsurvey <- function(posmat,
                             r=1,
                             r.low=rep(0,length(r))) {
  
  return(.C("get_pi_clustsurvey",
            as.double(posmat[,"p"]),
            as.double(posmat[,"x"]),
            as.double(posmat[,"y"]),
            as.integer(posmat[,"s"]),
            as.integer(nrow(posmat)),
            as.double(r.low),
            as.double(r),
            as.integer(length(r)),
            as.integer(1:nrow(posmat)),
            rc=double(length(r))
  )$rc)
}

##'
##' Optimized version of \code{get.pi} for typed data.
##'
##' Version of the \code{get.pi} function that is optimized for statically typed data. That is
##' data where we are interested in the probability of points within some distance of points of
##' typeA are of typeB.
##'
##' @param posmat a matrix with columns type, x, y, and weight
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##'
##' @return pi values for all the distances we looked at
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.pi
##'
##' @example example/get_pi_typed.R
##'
get.pi.clustsurvey.wts <- function(posmat,
                             r=1,
                             r.low=rep(0,length(r))) {
  
  return(.C("get_pi_clustsurvey_wts",
            as.double(posmat[,"p"]),
            as.double(posmat[,"x"]),
            as.double(posmat[,"y"]),
            as.integer(posmat[,"s"]),
            as.double(posmat[,"weight"]),
            as.integer(nrow(posmat)),
            as.double(r.low),
            as.double(r),
            as.integer(length(r)),
            as.integer(1:nrow(posmat)),
            rc=double(length(r))
  )$rc)
}


##' Calculate bootstrapped confidence intervals for \code{get.pi} values.
##'
##' Wrapper to \code{get.pi.bootstrap} that takes care of calculateing the
##' confience interavals based on the bootstrapped values..
##'
##'
##' @param posmat a matrix with columns type, x and y
##' @param fun the function to decide relationships
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range. 0 by default
##' @param boot.iter the number of bootstrap iterations
##' @param ci.low the low end of the ci...0.025 by default
##' @param ci.high the high end of the ci...0.975 by default
##'
##' @return a matrix with a row for the high and low values and
##'     a column per distance
##'
##' @author Justin Lessler
##'
##' @family get.pi
##'
##' @example example/get_pi_ci.R
##'
get.pi.ci <- function(posmat,
                      fun,
                      r=1,
                      r.low=rep(0,length(r)),
                      boot.iter = 1000,
                      ci.low=0.025,
                      ci.high=0.975) {
    boots <- get.pi.bootstrap(posmat, fun, r, r.low, boot.iter)

    rc <- matrix(nrow=2,ncol=ncol(boots))

    rownames(rc) <- c(ci.low,ci.high)

    for (i in 1:ncol(rc)) {
        rc[,i] <- quantile(boots[,i], probs=c(ci.low, ci.high))
    }

    return(rc)
}


##' Calculate bootstrapped confidence intervals for \code{get.pi.typed} values.
##'
##' Wrapper to \code{get.pi.typed.bootstrap} that takes care of calculateing the
##' confience interavals based on the bootstrapped values..
##'
##'
##' @param posmat a matrix with columns type, x and y
##' @param fun the function to decide relationships
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range. 0 by default
##' @param boot.iter the number of bootstrap iterations
##' @param ci.low the low end of the ci...0.025 by default
##' @param ci.high the high end of the ci...0.975 by default
##'
##' @return a matrix with a row for the high and low values and
##'     a column per distance
##'
##' @author Justin Lessler
##'
##' @family get.pi
##'
##' @example example/get_pi_ci.R
##'
get.pi.typed.ci <- function(posmat,
                      typeA = -1,
                      typeB = -1,
                      r=1,
                      r.low=rep(0,length(r)),
                      boot.iter = 1000,
                      ci.low=0.025,
                      ci.high=0.975) {
  boots <- get.pi.typed.bootstrap(posmat, typeA = -1, typeB = -1, r, r.low, boot.iter)
  
  rc <- matrix(nrow=2,ncol=ncol(boots))
  
  rownames(rc) <- c(ci.low,ci.high)
  
  for (i in 1:ncol(rc)) {
    rc[,i] <- quantile(boots[,i], probs=c(ci.low, ci.high))
  }
  
  return(rc)
}


##' Calculate bootstrapped confidence intervals for \code{get.pi.typed.wts} values.
##'
##' Wrapper to \code{get.pi.typed.wts.bootstrap} that takes care of calculateing the
##' confience interavals based on the bootstrapped values..
##'
##'
##' @param posmat a matrix with columns type, x and y
##' @param fun the function to decide relationships
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range. 0 by default
##' @param boot.iter the number of bootstrap iterations
##' @param ci.low the low end of the ci...0.025 by default
##' @param ci.high the high end of the ci...0.975 by default
##'
##' @return a matrix with a row for the high and low values and
##'     a column per distance
##'
##' @author Justin Lessler
##'
##' @family get.pi
##'
##' @example example/get_pi_ci.R
##'
get.pi.typed.wts.ci <- function(posmat,
                            typeA = -1,
                            typeB = -1,
                            r=1,
                            r.low=rep(0,length(r)),
                            boot.iter = 1000,
                            ci.low=0.025,
                            ci.high=0.975) {
  boots <- get.pi.typed.wts.bootstrap(posmat, typeA = -1, typeB = -1, r, r.low, boot.iter)
  
  rc <- matrix(nrow=2,ncol=ncol(boots))
  
  rownames(rc) <- c(ci.low,ci.high)
  
  for (i in 1:ncol(rc)) {
    rc[,i] <- quantile(boots[,i], probs=c(ci.low, ci.high))
  }
  
  return(rc)
}



##' Bootstrap \code{get.pi} values.
##'
##' Runs \code{get.pi} on multiple bootstraps of the data. Is formulated
##' such that the relationships between
##' points and themselves will not be calculated.
##'
##' @param posmat a matrix with columns type, x and y
##' @param fun the function to decide relationships
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range. 0 by default
##' @param boot.iter the number of bootstrap iterationsxs
##'
##' @return pi values for all the distances we looked at
##'
##' @note In each bootstrap iteration N observations are drawn from the existing data with replacement. To avoid errors in
##' inference resulting from the same observatin being compared with itself in the bootstrapped data set, original indices
##' are perserved, and pairs of points in the bootstrapped dataset with the same original index are ignored.
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.pi
##'
##' @example example/get_pi_bootstrap.R
##'
get.pi.bootstrap <- function(posmat,
                             fun,
                             r=1,
                             r.low=rep(0,length(r)),
                             boot.iter = 500) {


    xcol <-  which(colnames(posmat)=="x")
    ycol <- which(colnames(posmat)=="y")

    #check that both columns exist
    if (length(xcol)!=1 & length(ycol)!=1) {
        stop("unique x and y columns must be defined")
    }

    rc <- matrix(nrow=boot.iter, ncol=length(r))
    for (i in 1:boot.iter) {
        inds <- sample(nrow(posmat), replace=T)
        rc[i,] <- .Call("get_pi",
                        posmat[inds,],
                        fun,
                        r,
                        r.low,
                        inds,
                        xcol,
                        ycol)
    }
    return(rc)
}


##' runs bootstrapping on \code{get.pi.typed}
##'
##' Bootstraps typed pi values. Makes sure distances between a sample and
##' another draw of itself are left out
##'
##' @param boot.iter the number of bootstrap iterations
##' @param posmat a matrix with columns type, x and y
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##'
##' @return pi values for all the distances we looked at
##'
##' @family get.pi
##'
##' @example example/get_pi_typed_bootstrap.R
##'
get.pi.typed.bootstrap <- function(posmat,
                   typeA = -1,
                   typeB = -1,
                   r=1,
                   r.low=rep(0,length(r)),
                   boot.iter) {


  rc <- matrix(nrow=boot.iter, ncol=length(r))
  for (i in 1:boot.iter) {
    inds <- sample(nrow(posmat), replace=T)
    rc[i,] <- .C("get_pi_typed",
                 as.integer(posmat[inds,"type"]),
                 as.double(posmat[inds,"x"]),
                 as.double(posmat[inds,"y"]),
                 as.integer(nrow(posmat)),
                 as.integer(typeA),
                 as.integer(typeB),
                 as.double(r.low),
                 as.double(r),
                 as.integer(length(r)),
                 as.integer(inds),
                 rc=double(length(r))
                 )$rc
  }
  return(rc)
}


##' runs bootstrapping on \code{get.pi.typed.wts}
##'
##' Bootstraps typed pi values. Makes sure distances between a sample and
##' another draw of itself are left out
##'
##' @param boot.iter the number of bootstrap iterations
##' @param posmat a matrix with columns type, x and y
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##'
##' @return pi values for all the distances we looked at
##'
##' @family get.pi
##'
##' @example example/get_pi_typed_bootstrap.R
##'
get.pi.typed.wts.bootstrap <- function(posmat,
                                   typeA = -1,
                                   typeB = -1,
                                   r=1,
                                   r.low=rep(0,length(r)),
                                   boot.iter) {
  
  rc <- matrix(nrow=boot.iter, ncol=length(r))
  for (i in 1:boot.iter) {
    inds <- sample(nrow(posmat), replace=T)
    rc[i,] <- .C("get_pi_typed_wts",
                 as.integer(posmat[inds,"type"]),
                 as.double(posmat[inds,"x"]),
                 as.double(posmat[inds,"y"]),
                 as.double(posmat[inds,"weight"]),
                 as.integer(nrow(posmat)),
                 as.integer(typeA),
                 as.integer(typeB),
                 as.double(r.low),
                 as.double(r),
                 as.integer(length(r)),
                 as.integer(inds),
                 rc=double(length(r))
    )$rc
  }
  return(rc)
}




##' get the null distribution of the \code{get.pi} function
##'
##' Does permutations to calculate the null distribution of get pi
##' if there were no spatial dependence. Randomly reassigns coordinates
##' to each observation permutations times
##'
##' @param posmat a matrix with columns type, x and y
##' @param fun the function to evaluate
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##' @param permutations the number of permute iterations
##'
##' @return pi values for all the distances we looked at
##'
##' @family get.pi
##'
##' @example example/get_pi_permute.R
##'
get.pi.permute <- function(posmat,
                           fun,
                           r=1,
                           r.low=rep(0,length(r)),
                           permutations) {


    xcol <-  which(colnames(posmat)=="x")
    ycol <- which(colnames(posmat)=="y")

    #check that both columns exist
    if (length(xcol)!=1 & length(ycol)!=1) {
        stop("unique x and y columns must be defined")
    }

  rc <- matrix(nrow=permutations, ncol=length(r))
  for (i in 1:permutations) {
    inds <- sample(nrow(posmat))#, replace=T)
    tmp.posmat <- posmat
    tmp.posmat[,"x"] <- posmat[inds,"x"]
    tmp.posmat[,"y"] <- posmat[inds,"y"]
    rc[i,] <- .Call("get_pi",
                    tmp.posmat,
                    fun,
                    r,
                    r.low,
                    1:nrow(posmat),
                    xcol,
                    ycol)
  }
  return(rc)
}


##' get the null distribution of the get.pi.typed function
##'
##' Does permutations to calculate the null distribution of get pi
##' if there were no spatial dependence. Randomly reassigns coordinates
##' to each observation permutations times
##'
##' @param posmat a matrix with columns type, x and y
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##' @param permutations the number of permute iterations
##'
##' @return pi values for all the distances we looked at
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.pi
##'
##' @example example/get_pi_typed_permute.R
##'
get.pi.typed.permute <- function(posmat,
                                 typeA = -1,
                                 typeB = -1,
                                 r=1,
                                 r.low=rep(0,length(r)),
                                 permutations) {

    xcol <-  which(colnames(posmat)=="x")
    ycol <- which(colnames(posmat)=="y")

    #check that both columns exist
    if (length(xcol)!=1 & length(ycol)!=1) {
        stop("unique x and y columns must be defined")
    }


  rc <- matrix(nrow=permutations, ncol=length(r))
  for (i in 1:permutations) {
    inds <- sample(nrow(posmat))#, replace=T)
    rc[i,] <- .C("get_pi_typed",
                 as.integer(posmat[,"type"]),
                 as.double(posmat[inds,"x"]),
                 as.double(posmat[inds,"y"]),
                 as.integer(nrow(posmat)),
                 as.integer(typeA),
                 as.integer(typeB),
                 as.double(r.low),
                 as.double(r),
                 as.integer(length(r)),
                 as.integer(1:nrow(posmat)),
                 rc=double(length(r))
                 )$rc
  }
  return(rc)
}



##' get the null distribution of the \code{get.pi.typed.wts} function
##'
##' Does permutations to calculate the null distribution of get pi
##' if there were no spatial dependence. Randomly reassigns coordinates
##' to each observation permutations times
##'
##' @param posmat a matrix with columns type, x and y
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##' @param permutations the number of permute iterations
##'
##' @return pi values for all the distances we looked at
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.pi
##'
##' @example example/get_pi_typed_permute.R
##'
get.pi.typed.wts.permute <- function(posmat,
                                 typeA = -1,
                                 typeB = -1,
                                 r=1,
                                 r.low=rep(0,length(r)),
                                 permutations) {
  
  xcol <-  which(colnames(posmat)=="x")
  ycol <- which(colnames(posmat)=="y")
  
  #check that both columns exist
  if (length(xcol)!=1 & length(ycol)!=1) {
    stop("unique x and y columns must be defined")
  }
  
  rc <- matrix(nrow=permutations, ncol=length(r))
  for (i in 1:permutations) {
    inds <- sample(nrow(posmat))#, replace=T)
    rc[i,] <- .C("get_pi_typed_wts",
                 as.integer(posmat[,"type"]),
                 as.double(posmat[inds,"x"]),
                 as.double(posmat[inds,"y"]),
                 as.double(posmat[inds,"weight"]),
                 as.integer(nrow(posmat)),
                 as.integer(typeA),
                 as.integer(typeB),
                 as.double(r.low),
                 as.double(r),
                 as.integer(length(r)),
                 as.integer(1:nrow(posmat)),
                 rc=double(length(r))
    )$rc
  }
  return(rc)
}


##' generalized version of \code{get.tau}
##'
##' returns the relative probability that a points at some distance
##' from an index point share some relationship with that point versus
##' the probability any point shares that relationship with that point.
##'
##' @param posmat a matrix with columns x, y and any other named columns
##'    columns needed by fun
##' @param fun a function that takes in two rows of posmat and returns:
##' \enumerate{
##'      \item  for pairs included in the numerator and denominator
##'      \item for pairs that should only be included in the denominator
##'      \item for pairs that should be ignored all together}
##' Note that names from \code{posmat} are not perserved in calls to
##' \code{fun}, so the columns of the matrix should be referenced numerically
##' so this is not available to fun
##' @param r the series of spatial distances (or there maximums) we are
##'          interested in
##' @param r.low the low end of each range, 0 by default
##'
##' @return The tau value for each distance we look at. Where:
##'
##' \code{tau = get.pi(posmat, fun, r, r.low)/get.pi(posmat,fun,0,infinity)}
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.tau
##' @family spatialtau
##'
##' @example example/get_tau.R
get.tau <- function(posmat,
                   fun,
                   r = 1,
                   r.low=rep(0,length(r))) {

    xcol <-  which(colnames(posmat)=="x")
    ycol <- which(colnames(posmat)=="y")

    #check that both columns exist
    if (length(xcol)!=1 & length(ycol)!=1) {
        stop("unique x and y columns must be defined")
    }

    rc <- .Call("get_tau",
                posmat,
                fun,
                r,
                r.low,
                1:nrow(posmat),
                xcol,
                ycol)
    return(rc)
}


##' Optimizewd version of \code{get.tau} for typed data
##'
##' Version of th e \code{get.tau} function that is optimized for
##' statically typed data. That is dat where we want the relationship between
##' points of type A and points of type B
##'
##' @param posmat a matrix with columns type, x, y
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##'
##' @return tau values for all the distances we looked at
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.tau
##'
##' @example example/get_tau_typed.R
##'
get.tau.typed <- function(posmat,
                          typeA = -1,
                          typeB = -1,
                          r=1,
                          r.low=rep(0,length(r))) {
  return(.C("get_tau_typed",
            as.integer(posmat[,"type"]),
            as.double(posmat[,"x"]),
            as.double(posmat[,"y"]),
            as.integer(nrow(posmat)),
            as.integer(typeA),
            as.integer(typeB),
            as.double(r.low),
            as.double(r),
            as.integer(length(r)),
            as.integer(1:nrow(posmat)),
            rc=double(length(r))
  )$rc)
}



##' Optimizewd version of \code{get.tau} for typed data
##'
##' Version of th e \code{get.tau} function that is optimized for
##' statically typed data. That is dat where we want the relationship between
##' points of type A and points of type B
##'
##' @param posmat a matrix with columns type, x, y, and weight
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##'
##' @return tau values for all the distances we looked at
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.tau
##'
##' @example example/get_tau_typed.R
##'
get.tau.typed.wts <- function(posmat,
                          typeA = -1,
                          typeB = -1,
                          r=1,
                          r.low=rep(0,length(r))) {
  return(.C("get_tau_typed_wts",
            as.integer(posmat[,"type"]),
            as.double(posmat[,"x"]),
            as.double(posmat[,"y"]),
            as.double(posmat[,"weight"]),
            as.integer(nrow(posmat)),
            as.integer(typeA),
            as.integer(typeB),
            as.double(r.low),
            as.double(r),
            as.integer(length(r)),
            as.integer(1:nrow(posmat)),
            rc=double(length(r))
            )$rc)
}


##' Optimizewd version of \code{get.tau} for typed data
##'
##' Version of th e \code{get.tau} function that is optimized for
##' statically typed data. That is dat where we want the relationship between
##' points of type A and points of type B
##'
##' @param posmat a matrix with columns type, x, y, and weight
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##'
##' @return tau values for all the distances we looked at
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.tau
##'
##' @example example/get_tau_typed.R
##'
get.tau.clustsurvey <- function(posmat,
                                    r=1,
                                    r.low=rep(0,length(r))) {
  return(.C("get_tau_clustsurvey",
            as.integer(posmat[,"p"]),
            as.double(posmat[,"x"]),
            as.double(posmat[,"y"]),
            as.double(posmat[,"s"]),
            as.integer(nrow(posmat)),
            as.double(r.low),
            as.double(r),
            as.integer(length(r)),
            as.integer(1:nrow(posmat)),
            rc=double(length(r))
  )$rc)
}

##' Optimizewd version of \code{get.tau} for typed data
##'
##' Version of th e \code{get.tau} function that is optimized for
##' statically typed data. That is dat where we want the relationship between
##' points of type A and points of type B
##'
##' @param posmat a matrix with columns type, x, y, and weight
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##'
##' @return tau values for all the distances we looked at
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.tau
##'
##' @example example/get_tau_typed.R
##'
get.tau.clustsurvey.wts <- function(posmat,
                              r=1,
                              r.low=rep(0,length(r))) {
  return(.C("get_tau_clustsurvey_wts",
            as.integer(posmat[,"p"]),
            as.double(posmat[,"x"]),
            as.double(posmat[,"y"]),
            as.double(posmat[,"s"]),
            as.double(posmat[,"weight"]),
            as.integer(nrow(posmat)),
            as.double(r.low),
            as.double(r),
            as.integer(length(r)),
            as.integer(1:nrow(posmat)),
            rc=double(length(r))
  )$rc)
}

##' Bootstrap confidence interval for the \code{get.tau} values
##'
##' Wrapper to \code{get.tau.bootstrap} that takes care of calulating
##' the confidence intervals based on the bootstrapped values
##'
##' @param posmat a matrix appropriate for input to \code{get.tau}
##' @param fun a function appropriate as input to \code{get.pi}
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##' @param boot.iter the number of bootstrap iterations
##' @param ci.low the low end of the ci...0.025 by default
##' @param ci.high the high end of the ci...0.975 by default
##'
##' @return tau values for all the distances examined
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.tau
##'
##' @example example/get_tau_ci.R
##'
get.tau.ci <- function(posmat,
                       fun,
                       r=1,
                       r.low=rep(0,length(r)),
                       boot.iter = 1000,
                       ci.low=0.025,
                       ci.high=0.975) {
    boots <- get.tau.bootstrap(posmat, fun, r, r.low, boot.iter)

    rc <- matrix(nrow=2, ncol=ncol(boots))

    rownames(rc) <- c(ci.low,ci.high)

    for (i in 1:ncol(rc)) {
        rc[,i] <- quantile(boots[,i], probs=c(ci.low, ci.high))
    }

    return(rc)
}



##' Bootstrap confidence interval for the \code{get.tau.typed} values
##'
##' Wrapper to \code{get.tau.bootstrap} that takes care of calulating
##' the confidence intervals based on the bootstrapped values
##'
##' @param posmat a matrix appropriate for input to \code{get.tau}
##' @param fun a function appropriate as input to \code{get.pi}
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##' @param boot.iter the number of bootstrap iterations
##' @param ci.low the low end of the ci...0.025 by default
##' @param ci.high the high end of the ci...0.975 by default
##'
##' @return tau values for all the distances examined
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.tau
##'
##' @example example/get_tau_ci.R
##'
get.tau.typed.ci <- function(posmat,
                       typeA = -1,
                       typeB = -1,
                       r=1,
                       r.low=rep(0,length(r)),
                       boot.iter = 1000,
                       ci.low=0.025,
                       ci.high=0.975) {
  boots <- get.tau.typed.bootstrap(posmat, typeA = -1, typeB = -1, r, r.low, boot.iter)
  
  rc <- matrix(nrow=2, ncol=ncol(boots))
  
  rownames(rc) <- c(ci.low,ci.high)
  
  for (i in 1:ncol(rc)) {
    rc[,i] <- quantile(boots[,i], probs=c(ci.low, ci.high))
  }
  
  return(rc)
}




##' Bootstrap confidence interval for the \code{get.tau.typed.wts} values
##'
##' Wrapper to \code{get.tau.bootstrap} that takes care of calulating
##' the confidence intervals based on the bootstrapped values
##'
##' @param posmat a matrix appropriate for input to \code{get.tau}
##' @param fun a function appropriate as input to \code{get.pi}
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##' @param boot.iter the number of bootstrap iterations
##' @param ci.low the low end of the ci...0.025 by default
##' @param ci.high the high end of the ci...0.975 by default
##'
##' @return tau values for all the distances examined
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.tau
##'
##' @example example/get_tau_ci.R
##'
get.tau.typed.wts.ci <- function(posmat,
                             typeA = -1,
                             typeB = -1,
                             r=1,
                             r.low=rep(0,length(r)),
                             boot.iter = 1000,
                             ci.low=0.025,
                             ci.high=0.975) {
  boots <- get.tau.typed.wts.bootstrap(posmat, typeA = -1, typeB = -1, r, r.low, boot.iter)
  
  rc <- matrix(nrow=2, ncol=ncol(boots))
  
  rownames(rc) <- c(ci.low,ci.high)
  
  for (i in 1:ncol(rc)) {
    rc[,i] <- quantile(boots[,i], probs=c(ci.low, ci.high))
  }
  
  return(rc)
}




##' Bootstrap \code{get.tau} values.
##'
##' Runs \code{get.tau} on multiple bootstraps of the data. Is formulated
##' such that the relationship between points and themselves will not be
##' calculated
##'
##'
##' @param posmat a matrix appropriate for input to \code{get.tau}
##' @param fun a function appropriate as input to \code{get.pi}
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##' @param boot.iter the number of bootstrap iterations
##'
##' @return tau values for all the distances we looked at
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.tau
##'
##' @example example/get_tau_bootstrap.R
##'
get.tau.bootstrap <- function(posmat,
                              fun,
                              r=1,
                              r.low=rep(0,length(r)),
                              boot.iter) {


    xcol <-  which(colnames(posmat)=="x")
    ycol <- which(colnames(posmat)=="y")

    #check that both columns exist
    if (length(xcol)!=1 & length(ycol)!=1) {
        stop("unique x and y columns must be defined")
    }

  rc <- matrix(nrow=boot.iter, ncol=length(r))
  for (i in 1:boot.iter) {
      inds <- sample(nrow(posmat), replace=T)
      rc[i,] <- .Call("get_tau",
                      posmat[inds,],
                      fun,
                      r,
                      r.low,
                      inds,
                      xcol,
                      ycol)
  }
  return(rc)
}



##' runs bootstrapping for \code{get.tau.typed}
##'
##' @param posmat a matrix with columns type, x and y
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##' @param boot.iter the number of bootstrap iterations
##'
##' @return tau values for all the distances we looked at
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.tau
##'
##' @example example/get_tau_typed_bootstrap.R
##'
get.tau.typed.bootstrap <- function(posmat,
                                    typeA = -1,
                                    typeB = -1,
                                    r=1,
                                    r.low=rep(0,length(r)),
                                    boot.iter) {


  rc <- matrix(nrow=boot.iter, ncol=length(r))
  for (i in 1:boot.iter) {
    inds <- sample(nrow(posmat), replace=T)
    rc[i,] <- .C("get_tau_typed",
                 as.integer(posmat[inds,"type"]),
                 as.double(posmat[inds,"x"]),
                 as.double(posmat[inds,"y"]),
                 as.integer(nrow(posmat)),
                 as.integer(typeA),
                 as.integer(typeB),
                 as.double(r.low),
                 as.double(r),
                 as.integer(length(r)),
                 as.integer(inds),
                 rc=double(length(r))
                 )$rc
  }
  return(rc)
}



##' runs bootstrapping for \code{get.tau.typed.wts}
##'
##' @param posmat a matrix with columns type, x and y
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##' @param boot.iter the number of bootstrap iterations
##'
##' @return tau values for all the distances we looked at
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.tau
##'
##' @example example/get_tau_typed_bootstrap.R
##'
get.tau.typed.wts.bootstrap <- function(posmat,
                                    typeA = -1,
                                    typeB = -1,
                                    r=1,
                                    r.low=rep(0,length(r)),
                                    boot.iter) {
  
  rc <- matrix(nrow=boot.iter, ncol=length(r))
  for (i in 1:boot.iter) {
    inds <- sample(nrow(posmat), replace=T)
    rc[i,] <- .C("get_tau_typed_wts",
                 as.integer(posmat[inds,"type"]),
                 as.double(posmat[inds,"x"]),
                 as.double(posmat[inds,"y"]),
                 as.double(posmat[inds,"weight"]),
                 as.integer(nrow(posmat)),
                 as.integer(typeA),
                 as.integer(typeB),
                 as.double(r.low),
                 as.double(r),
                 as.integer(length(r)),
                 as.integer(inds),
                 rc=double(length(r))
    )$rc
  }
  return(rc)
}


##' get the null distribution of the \code{get.tau} function
##'
##' Does permutations to calculate the null distribution of get pi
##' if there were no spatial dependence. Randomly reassigns coordinates
##' to each observation permutations times
##'
##' @param posmat a matrix appropriate for input to \code{get.tau}
##' @param fun a function appropriate for input to \code{get.tau}
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##' @param permutations the number of permute iterations
##'
##' @return tau values for all the distances we looked at
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.tau
##'
##' @example example/get_tau_permute.R
##'
get.tau.permute <- function(posmat,
                            fun,
                            r=1,
                            r.low=rep(0,length(r)),
                            permutations) {


    xcol <-  which(colnames(posmat)=="x")
    ycol <- which(colnames(posmat)=="y")

    #check that both columns exist
    if (length(xcol)!=1 & length(ycol)!=1) {
        stop("unique x and y columns must be defined")
    }

    rc <- matrix(nrow=permutations, ncol=length(r))
    for (i in 1:permutations) {
        inds <- sample(nrow(posmat))#, replace=T)
        tmp.posmat <- posmat
        tmp.posmat[,"x"] <- posmat[inds,"x"]
        tmp.posmat[,"y"] <- posmat[inds,"y"]
        rc[i,] <- .Call("get_tau",
                        tmp.posmat,
                        fun,
                        r,
                        r.low,
                        1:nrow(posmat),
                        xcol,
                        ycol)
    }

    return(rc)
}



##' get the null distribution for the \code{get.tau.typed} function
##'
##'
##' @param posmat a matrix with columns type, x and y
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##' @param permutations the number of permute iterations
##'
##' @return a matrix with permutation tau values for each distance specified
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.tau
##'
##' @example example/get_tau_typed_permute.R
##'
get.tau.typed.permute <- function(posmat,
                                  typeA = -1,
                                  typeB = -1,
                                  r=1,
                                  r.low=rep(0,length(r)),
                                  permutations) {

    xcol <-  which(colnames(posmat)=="x")
    ycol <- which(colnames(posmat)=="y")

    #check that both columns exist
    if (length(xcol)!=1 & length(ycol)!=1) {
        stop("unique x and y columns must be defined")
    }

  rc <- matrix(nrow=permutations, ncol=length(r))
  for (i in 1:permutations) {
    inds <- sample(nrow(posmat))#, replace=T)
    rc[i,] <- .C("get_tau_typed",
                 as.integer(posmat[,"type"]),
                 as.double(posmat[inds,"x"]),
                 as.double(posmat[inds,"y"]),
                 as.integer(nrow(posmat)),
                 as.integer(typeA),
                 as.integer(typeB),
                 as.double(r.low),
                 as.double(r),
                 as.integer(length(r)),
                 #as.integer(inds),
                 as.integer(1:nrow(posmat)),
                 rc=double(length(r))
                 )$rc
  }
  return(rc)
}





##' get the null distribution for the \code{get.tau.typed} function
##'
##'
##' @param posmat a matrix with columns type, x and y
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##' @param permutations the number of permute iterations
##'
##' @return a matrix with permutation tau values for each distance specified
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.tau
##'
##' @example example/get_tau_typed_permute.R
##'
get.tau.typed.wts.permute <- function(posmat,
                                  typeA = -1,
                                  typeB = -1,
                                  r=1,
                                  r.low=rep(0,length(r)),
                                  permutations) {
  
  xcol <-  which(colnames(posmat)=="x")
  ycol <- which(colnames(posmat)=="y")
  
  #check that both columns exist
  if (length(xcol)!=1 & length(ycol)!=1) {
    stop("unique x and y columns must be defined")
  }
  
  rc <- matrix(nrow=permutations, ncol=length(r))
  for (i in 1:permutations) {
    inds <- sample(nrow(posmat))#, replace=T)
    rc[i,] <- .C("get_tau_typed_wts",
                 as.integer(posmat[,"type"]),
                 as.double(posmat[inds,"x"]),
                 as.double(posmat[inds,"y"]),
                 as.double(posmat[inds,"weight"]),
                 as.integer(nrow(posmat)),
                 as.integer(typeA),
                 as.integer(typeB),
                 as.double(r.low),
                 as.double(r),
                 as.integer(length(r)),
                 #as.integer(inds),
                 as.integer(1:nrow(posmat)),
                 rc=double(length(r))
    )$rc
  }
  return(rc)
}










NULL

##' @name DengueSimulationR01
##' @title Simulated dataset of dengue transmission with basic reproductive number of 1
##' @format Matrix with five columns representing the X and Y coordinates of infected individuals, the time of infection, the genotype of the infecting pathogen and the serotype of the infecting pathogen.
##' @description Dataset simulated using an agent based model with a spatially heterogeneous population structure. Infectious agents were introduced resulting in agent to agent transmission. The distance between successive cases in a transmission chain were randomly drawn from a uniform distribution U(0,100). Each infectious agent resulted in a single transmission to another agent after a delay of 15 days, reflecting the generation time of dengue. There are 11 transmission chains, each with a different genotype. The genotypes are subdivided into four serotypes. 
##' @docType data
##' @usage DengueSimulationR01
##' @author Justin Lessler and Henrik Salje

NULL

##' @name DengueSimulationR02
##' @title Simulated dataset of dengue cases with basic reproductive number of 2
##' @format Matrix with five columns representing the X and Y coordinates of infected individuals, the time of infection, the genotype of the infecting pathogen and the serotype of the infecting pathogen.
##' @description Dataset simulated using an agent based model with a spatially heterogeneous population structure. Infectious agents were introduced resulting in agent to agent transmission. The distance between successive cases in a transmission chain were randomly drawn from a uniform distribution U(0,100). Each infectious agent resulted in transmissions to two other agents after a delay of 15 days, reflecting the generation time of dengue. There are 11 transmission chains, each with a different genotype. The genotypes are subdivided into four serotypes. 
##' @docType data
##' @usage DengueSimulationR01
##' @author Justin Lessler and Henrik Salje

NULL

