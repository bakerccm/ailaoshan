# R code to run final occupancy models

library("here")
library("R2jags")

# get file names from command line arguments
args <- commandArgs(trailingOnly = TRUE)

data.filename <- args[1]    # file containing jags.data
jags.filename <- args[2]    # file containing jags model
output.filename <- args[3]  # file to save model output to as data stream

# read in prepared data
load(file = data.filename)

########################################################################################
# get start time
starttime <- Sys.time()
print(paste("Starting:", starttime))

########################################################################################
# run model
# note use of do.call with jags.parallel to avoid issue with mcmc.settings not being found
# note also that z.start and num.species need to be included in jags.data$model.data for this to work
model.output <- do.call(
    jags.parallel,
    args = list(
        model.file = jags.filename,
        data = jags.data$model.data, inits = jags.data$inits, parameters.to.save = jags.data$params,
        n.iter = jags.data$mcmc.settings$ni, n.thin = jags.data$mcmc.settings$nt,
        n.burnin = jags.data$mcmc.settings$nb, n.chains = jags.data$mcmc.settings$nc,
        working.directory = getwd()
    )
)

# save model results as data stream
saveRDS(model.output, file = output.filename)

########################################################################################
# get end time
endtime <- Sys.time()
print(paste("Finished:", endtime))
print(paste("Elapsed:", endtime - starttime))

########################################################################################

# session info
writeLines(capture.output(sessionInfo()), here("sessioninfo", "Ailaoshan_model_parallel_final.sessioninfo.txt"))
