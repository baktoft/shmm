# Spatial hidden Markov model (SHMM)
#    Copyright (C) 2015-2016  Martin Waever Pedersen, mawp@dtu.dk or wpsgodd@gmail.com
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


#' @name fit.shmm
#' @title Fit a continuous-time surplus production model to data.
#' @details Fit a continuous-time surplus production model to data.
#' @param inp List of input variables as output by check.inp.
#' @param dbg Debugging option. Will print out runtime information useful for debugging if set to 1. Will print even more if set to 2.
#' @return List containing results.
#' @export
#' @examples
#' rep <- fit.shmm(inp)
#' @import TMB
fit.shmm_new <- function(inp, dbg=0){
    ## FIXME: Reorder arrays to 'good' order before passing to TMB
    inp$datlik <- reorder_time(inp$datlik, last=TRUE)
    # Check input list
    inp <- check.inp_new(inp)
    rep <- list(inp=inp)
    
    tic <- Sys.time()
    cat('\nCreate objective function...')
    obj <- make.obj(inp) # Smoothing not performed here

    if (inp$do.estimation){
        cat('\nEstimating parameters...')
        rep$opt <- try(nlminb(obj$par, obj$fn, obj$gr))
    } else {
        cat('\nEvaluating for fixed parameters...')
        nouse <- obj$fn() # Evaluate fn for initial parameters
    }
    #rep$do.smoo <- inp$do.smoo
    TMB::FreeADFun(obj) ## Free AD tapes
    obj$env$type <- setdiff(obj$env$type, "ADFun") ## And prevent retaping (!)
    cat('\nSave report...')
    rep$reportraw <- obj$report()
    ## FIXME: Reorder arrays back to 'bad' order
    rep$reportraw <- reorder_time(rep$reportraw, last=FALSE)
    rep$report <- list()
    for (nm in names(rep$reportraw)){
        vecdistr <- rep$reportraw[[nm]]
        if (length(grep('out', nm)) > 0){
            nmpl <- sub('out', '', nm)
            if (length(dim(vecdistr)) == 2){
                rep$report[[nmpl]] <- format.distr(vecdistr, rep)
                #rep$report[[nm]] <- NULL
            }
        } else {
            rep$report[[nm]] <- vecdistr
        }
    }
    
    if (inp$do.sd.report){
        cat('\nCalculating uncertainties...')
        rep <- try(TMB::sdreport(obj))
    }

    rep$obj <- obj
    cat('\nCalculating tracks...')
    rep <- get.tracks(rep)

    rep$computing.time <- difftime(Sys.time(), tic, unit='secs')
    
    if(!is.null(rep)){
        class(rep) <- "shmmcls"
    }
    cat('\n')
    return(rep)
}

## Kasper: All 3D arrays in original code are badly ordered as (time, x, y). On TMB side, all arrays should ideally have the time dimension slowest running i.e., come *last*.
## > lapply(inp$datlik,dim)
## $sst
## [1] 217 440 200
## $bathy
## [1] 217 440 200
## $sun
## [1] 217 440 200
## $release
## [1]   1 440 200
## $recapture
## [1]   1 440 200
## $all
## [1]   217 88000
## The following function solves this problem. Example:
## > tmp <- inp$datlik
## > tmp <- reorder_time(tmp,TRUE)  ## Change to good order
## > tmp <- reorder_time(tmp,FALSE) ## Change back to bad order
## > identical(tmp, inp$datlik)
## [1] TRUE
reorder_time <- function(x, last=TRUE) {
    if (is.list(x))
        return (lapply(x, reorder_time, last=last))
    if (is.null(dim(x)))
        return (x)
    p <- 1:length(dim(x))
    if (last)
        p <- c(p[-1], p[1])
    else
        p <- c(p[length(p)], p[-length(p)])
    aperm(x, p)
}

#' @name make.obj
#' @title Create TMB obj using TMB::MakeADFun and squelch screen printing.
#' @param inp List of input variables as output by check.inp.
#' @param phase Estimation phase, integer.
#' @return List to be used as data input to TMB.
#' @export
#' @import TMB
make.obj <- function(inp, phase=1){
    datlist <- make.datlist(inp)
    parlist <- make.parlist(inp)
    # This one requires memory if not using atomic
    obj <- TMB::MakeADFun(data=datlist,
                          parameters=parlist,
                          random=NULL,
                          DLL=inp$scriptname,
                          checkParameterOrder=FALSE,
                          map=inp$map)
    # Make TMB quiet
    TMB:::config(trace.optimize=0, DLL=inp$scriptname)
    verbose <- inp$verbose
    if (is.null(verbose)) verbose <- TRUE
    obj$env$tracemgc <- verbose
    obj$env$inner.control$trace <- verbose
    obj$env$silent <- ! verbose
    obj$fn(obj$par)
    return(obj)
}


#' @name make.datlist
#' @title Create data list used as input to TMB::MakeADFun.
#' @param inp List of input variables as output by check.inp.
#' @return List to be used as data input to TMB::MakeADFun.
#' @export
make.datlist <- function(inp){
    datlist <- list(doviterbi=inp$do.viterbi,
                    datlik=inp$datlik$all,
                    solvetype=inp$solvetype,
                    ns=inp$ns,
                    iobs=inp$iobs,
                    #isave=inp$isave,
                    #dosmoo=inp$do.smoo,
                    #doviterbi=inp$do.viterbi,
                    I=inp$gen$I,
                    dt=inp$dt,
                    m=inp$gen$m,
                    Sns=inp$gen$Sns,
                    Sew=inp$gen$Sew,
                    Sx=inp$gen$Sx,
                    Sy=inp$gen$Sy,
                    lgam=inp$gen$lgam)
    ## Add 'rawdat' flag + rawdat if present
    datlist$rawdat <- !is.null(inp$rawdat)
    datlist <- c(datlist, inp$rawdat)
    return(datlist)
}


#' @name make.parlist
#' @title Create parameter list used as input to TMB::MakeADFun.
#' @param inp List of input variables as output by check.inp.
#' @return List to be used as parameter input to TMB::MakeADFun.
#' @export
make.parlist <- function(inp){
    if (is.list(inp$parlist))
        return (inp$parlist)
    parlist <- list(logDx=log(inp$gen$Dx),
                    logDy=log(inp$gen$Dy))
    if (!is.null(inp$rawdat)) {
        k <- length(inp$rawdat$obs[[1]])
        parlist$logsd <- rep(0, k)
        parlist$bias <- rep(0, k)
    }
    return(parlist)
}

