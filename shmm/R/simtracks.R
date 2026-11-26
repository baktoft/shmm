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


#' @name sim.tracks
#' @title Simulate tracks from the posterior distribution of the track.
#' @details Draws complete tracks from p(x_1:ns | y_1:ns) by forward-filtering
#'     backward-sampling. The forward filter has already been run by
#'     fit.shmm_new() and is available as rep$report$phi; the backward pass
#'     draws the last state from it and steps backwards drawing
#'
#'         x_t | x_t+1, y_1:t  ~  phi_t(i) * P(i -> x_t+1),
#'
#'     where P = expm(G*dt). The draws are exact, and unlike the mean track in
#'     rep$tracks$mean each one is a coherent trajectory, so quantities such as
#'     occupation time or the extremes of the track can be computed per draw and
#'     summarised over draws.
#'
#'     method="clean" returns those daily positions and nothing else. Land
#'     constrains the fitted model through the emission, where the con mask
#'     makes the likelihood zero on land, but not through the generator, which
#'     diffuses over the full rectangle. Two consecutive daily positions can
#'     therefore both be in water with the straight line between them crossing a
#'     peninsula.
#'
#'     method="bridge" adds m sub-daily positions per day for the steps where
#'     that happens, drawn from a diffusion bridge of the same fitted generator
#'     restricted to the water cells, by the same backward sampling. It takes
#'     the emission mask at face value, as saying the fish is in water rather
#'     than merely observed there. The sub-daily positions are model draws and
#'     cannot be on land, though the straight line drawn between two of them can
#'     still clip a coast; land.crossed() measures how much.
#'
#'     Either way this repairs the drawn path, not the inference: the daily
#'     positions are drawn under the unmasked generator. Nothing derived from
#'     the sub-daily points should be quoted, and path length is not a usable
#'     quantity at any resolution.
#'
#'     The draws are conditional on the estimated parameters.
#' @param rep Result of running fit.shmm_new.
#' @param nsim Number of tracks to simulate.
#' @param method "clean" for the daily positions, "bridge" to also fill in
#'     sub-daily positions where a step crosses land.
#' @param m Sub-daily steps per day, method="bridge" only.
#' @param addJitter Spread positions uniformly within their grid cell,
#'     method="clean" only.
#' @param land Logical matrix (nx, ny) marking land. Defaults to the mask the
#'     model itself uses, con != 0 from rep$inp$rawdat, which is TRUE on water.
#'     check.inp_new() overwrites rep$inp$land with an empty transposed matrix,
#'     so that is not usable; if there is no rawdat, pass land = inp$land from
#'     the unchecked inp, where TRUE is water, the opposite of find.land().
#' @param wet Is the land argument TRUE on water rather than on land?
#' @param snap Move daily positions that fall on land to the nearest water cell.
#' @param tol Tolerance of the land test, in grid cells.
#' @param seed Optional seed.
#' @return A list with X and Y (one column per track), the cell indices ix and
#'     iy of the daily positions, daily (the rows of X and Y holding them), and
#'     for method="bridge" the log probability logp of each repaired step under
#'     the water-masked diffusion, so steps the model can barely make while
#'     staying in water are visible.
#' @export
#' @examples
#' sim <- sim.tracks(rep, nsim = 100, seed = 1)
#' quantile(colSums(sim$X < -20), c(0.025, 0.5, 0.975))
#' bri <- sim.tracks(rep, nsim = 100, method = "bridge", seed = 1)
sim.tracks <- function(rep, nsim=100, method=c('clean', 'bridge'), m=24,
                       addJitter=FALSE, land=NULL, wet=TRUE, snap=TRUE,
                       tol=1/32, seed=NULL){
    method <- match.arg(method)
    if (addJitter && method == 'bridge'){
        warning('addJitter is ignored for method="bridge"')
    }
    if (!is.null(seed)) set.seed(seed)
    grid <- rep$inp$grid
    phi <- rep$report$phi
    if (is.null(phi)){
        stop('rep$report$phi not found, was the model fitted with fit.shmm_new()?')
    }
    ns <- dim(phi)[1]
    trans <- get.transition(rep)
    Ax <- trans$Ax
    Ay <- trans$Ay

    ix <- matrix(0L, ns, nsim)
    iy <- matrix(0L, ns, nsim)
    last <- rcat2(phi[ns, , ], nsim)
    ix[ns, ] <- last$ix
    iy[ns, ] <- last$iy
    for (t in (ns-1):1){
        phit <- phi[t, , ]
        wy <- Ay[, iy[t+1, ], drop=FALSE]
        wx <- Ax[, ix[t+1, ], drop=FALSE] * (phit %*% wy)
        ix[t, ] <- rcat(wx)
        iy[t, ] <- rcat(t(phit[ix[t, ], , drop=FALSE]) * wy)
    }

    if (method == 'bridge'){
        return(bridge.tracks(rep, ix, iy, m, land, wet, snap, tol))
    }
    X <- matrix(grid$xx[ix], ns, nsim)
    Y <- matrix(grid$yy[iy], ns, nsim)
    if (addJitter){
        X <- X + stats::runif(length(X), -0.5, 0.5) * grid$dx
        Y <- Y + stats::runif(length(Y), -0.5, 0.5) * grid$dy
    }
    out <- list(X=X, Y=Y, ix=ix, iy=iy, daily=seq_len(ns), m=1L,
                date=rep$inp$date, nsim=nsim, method=method)
    class(out) <- 'shmmtracks'
    return(out)
}


bridge.tracks <- function(rep, ix, iy, m, land, wet, snap, tol){
    grid <- rep$inp$grid
    isWet <- water.mask(rep, land, wet)
    ns <- nrow(ix)
    nsim <- ncol(ix)

    nsnapped <- 0
    onland <- which(!isWet[cbind(as.vector(ix), as.vector(iy))])
    if (snap && length(onland) > 0){
        wcell <- which(isWet, arr.ind=TRUE)
        for (k in onland){
            j <- which.min((wcell[, 1] - ix[k])^2 + (wcell[, 2] - iy[k])^2)
            ix[k] <- wcell[j, 1]
            iy[k] <- wcell[j, 2]
        }
        nsnapped <- length(onland)
    }
    X <- matrix(grid$xx[ix], ns, nsim)
    Y <- matrix(grid$yy[iy], ns, nsim)

    bad <- which(!clear(isWet, grid, as.vector(X[-ns, ]), as.vector(Y[-ns, ]),
                        as.vector(X[-1, ]), as.vector(Y[-1, ]), tol))
    gen <- water.generator(isWet, grid, get.D(rep))
    fix <- vector('list', length(bad))
    logp <- rep(NA_real_, length(bad))
    nfailed <- 0
    for (n in seq_along(bad)){
        k <- bad[n]
        t <- (k - 1) %% (ns - 1) + 1
        s <- (k - 1) %/% (ns - 1) + 1
        br <- bridge.seg(gen, rep$inp$dt, m, gen$id[ix[t, s], iy[t, s]],
                         gen$id[ix[t+1, s], iy[t+1, s]])
        if (is.null(br)){
            nfailed <- nfailed + 1
        } else {
            fix[[n]] <- br$path
            logp[n] <- br$logp
        }
    }

    nfine <- (ns - 1) * m + 1
    fx <- matrix(0, nfine, nsim)
    fy <- matrix(0, nfine, nsim)
    for (j in seq_len(m)){
        rows <- seq(j, by=m, length.out=ns-1)
        u <- (j - 1) / m
        fx[rows, ] <- X[-ns, ] + u * (X[-1, ] - X[-ns, ])
        fy[rows, ] <- Y[-ns, ] + u * (Y[-1, ] - Y[-ns, ])
    }
    fx[nfine, ] <- X[ns, ]
    fy[nfine, ] <- Y[ns, ]
    for (n in seq_along(bad)){
        if (is.null(fix[[n]])) next
        k <- bad[n]
        t <- (k - 1) %% (ns - 1) + 1
        s <- (k - 1) %/% (ns - 1) + 1
        rows <- (t - 1) * m + seq_len(m)
        fx[rows, s] <- grid$xx[gen$cx[fix[[n]]]]
        fy[rows, s] <- grid$yy[gen$cy[fix[[n]]]]
    }
    out <- list(X=fx, Y=fy, ix=ix, iy=iy, daily=seq(1, by=m, length.out=ns),
                m=m, date=rep$inp$date, nsim=nsim, method='bridge', logp=logp,
                nbridged=length(bad) - nfailed, nfailed=nfailed,
                nsnapped=nsnapped, nseg=(ns - 1) * nsim, tol=tol)
    class(out) <- 'shmmtracks'
    return(out)
}


#' @name land.crossed
#' @title How much land a set of drawn tracks passes over.
#' @param rep Result of running fit.shmm_new.
#' @param X Matrix of longitudes, one column per track.
#' @param Y Matrix of latitudes.
#' @param land Logical matrix (nx, ny) marking land. Defaults as in sim.tracks.
#' @param wet Is the land argument TRUE on water rather than on land?
#' @param dens Number of samples per drawn segment.
#' @return List with the number of drawn segments, how many cross land, the
#'     total land distance crossed in grid cells, and the worst single crossing.
#' @export
land.crossed <- function(rep, X, Y, land=NULL, wet=TRUE, dens=200){
    grid <- rep$inp$grid
    isWet <- water.mask(rep, land, wet)
    n <- nrow(X)
    x0 <- as.vector(X[-n, ]); x1 <- as.vector(X[-1, ])
    y0 <- as.vector(Y[-n, ]); y1 <- as.vector(Y[-1, ])
    len <- sqrt(((x1 - x0) / grid$dx)^2 + ((y1 - y0) / grid$dy)^2)
    hits <- numeric(length(x0))
    for (u in seq(0, 1, length.out=dens)){
        i <- clamp(round((x0 + u * (x1 - x0) - grid$xx[1]) / grid$dx) + 1, grid$nx)
        j <- clamp(round((y0 + u * (y1 - y0) - grid$yy[1]) / grid$dy) + 1, grid$ny)
        hits <- hits + !isWet[cbind(i, j)]
    }
    over <- hits / dens * len
    return(list(nseg=length(over), ncross=sum(over > 0), cells=sum(over),
                worst=max(over)))
}


#' @name get.transition
#' @title One-step transition matrices of the fitted model.
#' @details Returns Ax and Ay such that the full transition matrix is their
#'     Kronecker product, P[(i,j),(k,l)] = Ax[i,k] * Ay[j,l], mirroring the
#'     separable matrix exponential used on the TMB side.
#' @param rep Result of running fit.shmm_new.
#' @return List with elements Ax and Ay.
#' @export
get.transition <- function(rep){
    inp <- rep$inp
    if (!is.null(inp$solvetype) && inp$solvetype != 'separable'){
        warning('get.transition() assumes solvetype "separable".')
    }
    D <- get.D(rep)
    return(list(Ax=expm.sym(D$Dx * inp$dt * inp$gen$Sx),
                Ay=expm.sym(D$Dy * inp$dt * inp$gen$Sy)))
}


#' @name get.D
#' @title Extract estimated diffusivities from a fit.
#' @details Reads them from obj$env$last.par, the parameters obj$report() used,
#'     so the transition matrix is consistent with the reported filter.
#' @param rep Result of running fit.shmm_new.
#' @return List with elements Dx and Dy.
#' @export
get.D <- function(rep){
    par <- rep$obj$env$last.par
    opt <- if (inherits(rep$opt, 'try-error')) NULL else rep$opt$par
    if (!is.null(par) && !is.null(opt)){
        if (identical(names(par), names(opt)) && max(abs(par - opt)) > 1e-6){
            warning('reported distributions are not at the optimum, sim.tracks() uses the reported ones')
        }
    }
    if (is.null(par)) par <- opt
    return(list(Dx=inp.or.par(par, 'logDx', rep$inp$gen$Dx),
                Dy=inp.or.par(par, 'logDy', rep$inp$gen$Dy)))
}


#' @name get.sd
#' @title Movement standard deviation of the fitted model, in grid units.
#' @param rep Result of running fit.shmm_new.
#' @return List with elements sdx and sdy.
#' @export
get.sd <- function(rep){
    D <- get.D(rep)
    return(list(sdx=D2sd(D$Dx * rep$inp$dt, rep$inp$grid$dx),
                sdy=D2sd(D$Dy * rep$inp$dt, rep$inp$grid$dy)))
}


#' @name plotshmm.simtracks
#' @title Plot simulated posterior tracks.
#' @param rep Result of running fit.shmm_new.
#' @param sim Output of sim.tracks.
#' @param nplot Number of tracks to draw. The default keeps the plot readable
#'     rather than a solid block; the title says how many of how many are
#'     shown. Use nplot = Inf to draw all of them.
#' @param add.map Add map?
#' @param add Add to existing plot?
#' @return Nothing.
#' @export
plotshmm.simtracks <- function(rep, sim, nplot=25, add.map=TRUE, add=FALSE){
    plot.shmmtracks(sim, rep=rep, nplot=nplot, add.map=add.map, add=add)
}


#' @name plot.shmmtracks
#' @title Plot simulated posterior tracks.
#' @param x Output of sim.tracks.
#' @param rep Result of running fit.shmm_new, to add the mean track. Optional.
#' @param nplot Number of tracks to draw. The default keeps the plot readable
#'     rather than a solid block; the title says how many of how many are
#'     shown. Use nplot = Inf to draw all of them.
#' @param add.map Add map?
#' @param add Add to existing plot?
#' @param ... Ignored.
#' @return Nothing.
#' @export
plot.shmmtracks <- function(x, rep=NULL, nplot=25, add.map=TRUE, add=FALSE, ...){
    n <- ncol(x$X)
    ii <- seq_len(min(nplot, n))
    if (!add){
        main <- if (length(ii) < n) paste(length(ii), 'of', n, 'posterior tracks')
                else paste(n, 'posterior tracks')
        plot(range(x$X), range(x$Y), type='n', xlab='Longitude', ylab='Latitude',
             main=main)
        if (add.map) maps::map(add=TRUE, fill=TRUE, col='grey80')
    }
    for (i in ii) lines(x$X[, i], x$Y[, i], col=grDevices::rgb(0, 0, 0, 0.2))
    if (!is.null(rep$tracks$mean)){
        lines(rep$tracks$mean$X, rep$tracks$mean$Y, col='blue', lwd=2)
    }
    invisible(NULL)
}


inp.or.par <- function(par, nm, default){
    if (!is.null(par) && nm %in% names(par)) return(exp(unname(par[nm])))
    return(default)
}



# The water mask. check.inp_new() overwrites inp$land with an empty transposed
# matrix, so the default is the mask the model itself uses: con is zero exactly
# on land.
water.mask <- function(rep, land=NULL, wet=TRUE){
    grid <- rep$inp$grid
    if (is.null(land)){
        con <- rep$inp$rawdat$con
        if (!is.null(con) && identical(dim(con), c(grid$nx, grid$ny))){
            return(con != 0)
        }
        land <- rep$inp$land
    }
    if (!is.matrix(land) || !identical(dim(land), c(grid$nx, grid$ny))){
        stop('land must be an (nx, ny) = ', grid$nx, ' x ', grid$ny,
             ' matrix, got ', paste(dim(land), collapse=' x '),
             '. Pass land = inp$land from the unchecked inp.')
    }
    isWet <- if (wet) land == TRUE else land == FALSE
    if (!any(isWet)) stop('no water cells, check the wet argument')
    return(isWet)
}


clamp <- function(i, n) pmin(pmax(i, 1), n)


expm.sym <- function(A){
    ei <- eigen(A, symmetric=TRUE)
    B <- ei$vectors %*% (exp(ei$values) * t(ei$vectors))
    B[B < 0] <- 0
    return(B / rowSums(B))
}


# Sample one row index per column of a matrix of non-negative weights
rcat <- function(w){
    cw <- apply(w, 2, cumsum)
    tot <- cw[nrow(cw), ]
    if (any(!is.finite(tot)) || any(tot <= 0)){
        stop('degenerate weights in backward sampling')
    }
    u <- stats::runif(ncol(cw)) * tot
    return(colSums(cw < rep(u, each=nrow(cw))) + 1L)
}


rcat2 <- function(p, nsim){
    j <- sample.int(length(p), nsim, replace=TRUE, prob=as.vector(p))
    return(list(ix=((j-1) %% nrow(p)) + 1L, iy=((j-1) %/% nrow(p)) + 1L))
}


# Does the straight line stay in water? Vectorised over segments.
clear <- function(isWet, grid, x0, y0, x1, y1, tol=1/32){
    len <- sqrt(((x1 - x0) / grid$dx)^2 + ((y1 - y0) / grid$dy)^2)
    nchk <- max(2, ceiling(max(len) / tol) + 1)
    good <- rep(TRUE, length(x0))
    for (u in seq(0, 1, length.out=nchk)){
        i <- clamp(round((x0 + u * (x1 - x0) - grid$xx[1]) / grid$dx) + 1, grid$nx)
        j <- clamp(round((y0 + u * (y1 - y0) - grid$yy[1]) / grid$dy) + 1, grid$ny)
        good <- good & isWet[cbind(i, j)]
    }
    return(good)
}


# Fitted generator restricted to the water cells. The diagonal is minus the row
# sum, so the coast is reflecting and probability is conserved.
water.generator <- function(isWet, grid, D){
    nx <- grid$nx
    ny <- grid$ny
    id <- matrix(0L, nx, ny)
    w <- which(isWet)
    id[w] <- seq_along(w)
    nw <- length(w)
    cx <- ((w - 1) %% nx) + 1
    cy <- ((w - 1) %/% nx) + 1
    edge <- function(dx, dy, rate){
        to <- cbind(cx + dx, cy + dy)
        ok <- to[, 1] >= 1 & to[, 1] <= nx & to[, 2] >= 1 & to[, 2] <= ny
        j <- rep(0L, nw)
        j[ok] <- id[to[ok, , drop=FALSE]]
        k <- which(j > 0)
        cbind(k, j[k], rate)
    }
    E <- rbind(edge(1, 0, D$Dx), edge(-1, 0, D$Dx),
               edge(0, 1, D$Dy), edge(0, -1, D$Dy))
    G <- Matrix::sparseMatrix(i=E[, 1], j=E[, 2], x=E[, 3], dims=c(nw, nw))
    Matrix::diag(G) <- -Matrix::rowSums(G)
    return(list(G=G, id=id, cx=cx, cy=cy, nw=nw, lam=max(-Matrix::diag(G))))
}


# exp(G*h) %*% v by uniformization, without forming the matrix
expm.act <- function(gen, h, v){
    lam <- gen$lam
    nterm <- ceiling(lam * h + 4 * sqrt(lam * h) + 5)
    term <- v
    out <- v
    for (k in 1:nterm){
        term <- (as.vector(gen$G %*% term) / lam + term) * (lam * h / k)
        out <- out + term
    }
    return(out * exp(-lam * h))
}


# Diffusion bridge from cell ia to cell ib in m steps through water, by the same
# forward-backward sampling used for the daily positions.
bridge.seg <- function(gen, dt, m, ia, ib){
    if (ia < 1 || ib < 1) return(NULL)
    h <- dt / m
    r <- vector('list', m - 1)
    v <- numeric(gen$nw)
    v[ib] <- 1
    logscale <- 0
    for (k in (m-1):1){
        v <- expm.act(gen, h, v)
        sc <- max(v)
        if (!is.finite(sc) || sc <= 0) return(NULL)
        v <- v / sc                     # rescale, the bridge underflows otherwise
        logscale <- logscale + log(sc)
        r[[k]] <- v
    }
    path <- integer(m)
    path[1] <- ia
    logp <- NA_real_
    for (k in 1:(m-1)){
        q <- expm.act(gen, h, unit.vec(gen$nw, path[k]))
        wgt <- q * r[[k]]
        tot <- sum(wgt)
        if (!is.finite(tot) || tot <= 0) return(NULL)
        if (k == 1) logp <- log(tot) + logscale
        path[k+1] <- sample.int(gen$nw, 1, prob=wgt)
    }
    return(list(path=path, logp=logp))
}


unit.vec <- function(n, i){
    v <- numeric(n)
    v[i] <- 1
    return(v)
}
