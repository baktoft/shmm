remotes::install_github('baktoft/shmm', ref='dev_2021')

rm(list=ls())
library(shmm)

# # test doesn't work because error in land-stuff - uncomment line 81 and comment line 80 in checkinp.R to make test work.
rep <- test.shmm(nobs=100, do.viterbi=0)
plot(rep$tracks$mean$X, rep$tracks$mean$Y, type="l")
lines(rep$inp$true$X, rep$inp$true$Y, col="red")




