library(shmm)

## Hack 1: Get a dummy 'inp' object that 'fit.shmm' will accept.
## Release and re-capture are obtained from the first and last datalikelihood image (FIXME)
dat <- readRDS("hba2/inp_34861.rds")
inp <- calc.data.likelihood(dat)
rmna <- function(x){x[!is.finite(x)] <- 0; x}
inp$datlik$sst <- rmna(inp$datlik$sst) ## NOTE
inp$datlik$all <- rmna(inp$datlik$all) ## NOTE

## Hack 2: Add 'rawdat' to 'inp'. (cpp code will then bypass the old 'datalikelihood' format - except 1st and last frame...)

## Load raw data (I think I created this myself via a script from hba...)
load("dat.RData")
names(dat)
##  [1] "Dawn"         "Dusk"         "sst"          "sun_dat_dawn" "sun_dat_dusk"
##  [6] "sst_env"      "land"         "psat"         "tagDate"      "popDate"     
## [11] "tagLon"       "tagLat"       "popLon"       "popLat"       "lon_range"   
## [16] "lat_range"

## Reshape 1 (messy - cleanup)
get_maps <- function(dat) {
    maps <- dat[sapply(dat, function(x) length(dim(x)) == 3)]
    maps <- shmm:::reorder_time(maps)
    ans <- array(0, c(dim(maps[[1]]), length(maps)))
    for (i in 1:length(maps)) ans[,,,i] <- maps[[i]]
    dimnames(ans) <- list(x=NULL,y=NULL,time=NULL,datatype=names(maps))
    ans
}
map <- get_maps(dat)
## > dim(map)
## [1] 440 200 217   3
dimnames(map)
## Land (constraint)
con <- !is.na(map[,,100,"sst_env"])
## Get corresponding observations
obs <- cbind(dat$Dawn,dat$Dusk,dat$sst)
## select datatypes
select <- 1:3
obs <- obs[,select,drop=FALSE]
map <- map[,,,select,drop=FALSE]
## Get constraint
con <- !dat$land
rawdat <- list(map=map,obs=obs,con=con)

## Reshape to lists (now used by the cpp)
rawdat$map <- lapply(1:217, function(i)rawdat$map[,,i,])
rawdat$obs <- lapply(1:217, function(i)rawdat$obs[i,])
## Non-missing indicator
rawdat$map_keep <- lapply(rawdat$map, function(x)is.finite(x))
rawdat$obs_keep <- lapply(rawdat$obs, function(x)is.finite(x))
## Remove NAs
rmna <- function(x){x[!is.finite(x)] <- 0; x}
rawdat$map <- lapply(rawdat$map, rmna)
rawdat$obs <- lapply(rawdat$obs, rmna)

## FINAL RAW DATA LAYOUT:
## * 'map' and 'map_keep' are lists of conformable 3D arrays (x,y,datatype)
## * 'obs' and 'obs_keep' are lists of conformable 1D vectors (datatype)
## * 'con' (constraint) is a 2D array (x,y)
typeof(rawdat)
## [1] "list"
names(rawdat)
## [1] "map"      "obs"      "con"      "map_keep" "obs_keep"

## Set rawdat
inp$rawdat <- rawdat

## Remove redundant copies from workspace
rm(list=c("dat","map","rawdat")); gc()
##TMB::config(optimize.instantly=0)

system.time(rep <- fit.shmm(inp, dbg = TRUE))
## Calculating tracks...
##    user  system elapsed 
## 833.363   5.169 561.514 

## - 2 diffusion parameters
## - 3 sd parameters (one for each datatype)
## - 3 bias parameters (one for each datatype)
rep$opt$par
##        logDx        logDy        logsd        logsd        logsd         bias 
##  2.998408855  1.716505790 -4.204136762 -4.237383030 -1.883358934  0.011129956 
##         bias         bias 
## -0.007268652 -0.430622567
rep$opt$objective
## [1] -743.5178

## In this case track is crossing land. There are two possible explanations:
## 1. The animal moves too fast compared to the solver timestep => Solution: Use the other solver (uniformization) which has the constraint built into the generator.
## 2. The track is estimated as the 'snapshot mean' which can in theory be on land. => Solution: Use viterbi to get most probable track
image(inp$grid$xx,inp$grid$yy,con)
points(rep$tracks$mean$X,rep$tracks$mean$Y,type="b")
