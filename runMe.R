

# devtools::install_github('kaskr/adcomp/TMB', ref='v1.7.3') # virker :-)
# devtools::install_github('kaskr/adcomp/TMB', ref='v1.7.10') # virker :-)
# devtools::install_github('kaskr/adcomp/TMB', ref='v1.7.11') # virker :-)
# devtools::install_github('kaskr/adcomp/TMB', ref='v1.7.12') # FAIL!!!!
# devtools::install_github('kaskr/adcomp/TMB', ref='v1.7.15') # FAIL!!! virker ikke - kan ikke installere pakken, men hvis pakken allerede er inst virker det...
devtools::install('C:/Users/hbak/Desktop/shmm/shmm')

rm(list=ls())
library(shmm)

# # test doesn't work because error in land-stuff - uncomment line 81 and comment line 80 in checkinp.R to make test work.
# rep <- test.shmm(nobs=50, do.viterbi=0)
# plot(rep$tracks$mean$X, rep$tracks$mean$Y, type="l")
# lines(rep$inp$true$X, rep$inp$true$Y, col="red")

