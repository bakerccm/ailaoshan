# R code for calculating posterior mean jaccard matrices from occupancy model output

library("here")
library("tidyverse")
library("R2jags")

########################################################################################
# get file names from command line arguments
args <- commandArgs(trailingOnly = TRUE)

modeldata.rdata.filename <- args[1]    # file containing model data
# modeldata.rdata.filename <- here("rdata", "Ailaoshan_model_data_final_LSU.rdata")
# modeldata.rdata.filename <- here("rdata", "Ailaoshan_model_data_final_SSU.rdata")

modeloutput.rds.filename <- args[2]    # file containing modelling results
# modeloutput.rds.filename <- here("rds", "Ailaoshan_model_output_final_LSU.rds")
# modeloutput.rds.filename <- here("rds", "Ailaoshan_model_output_final_SSU.rds")

site_distances.rds.filename <- args[3]    # file name to save site jaccard distances to
# site_distances.rds.filename <- here("rds", "Ailaoshan_final_LSU_jaccard_site.rds")
# site_distances.rds.filename <- here("rds", "Ailaoshan_final_SSU_jaccard_site.rds")

species_distances.rds.filename <- args[4]    # file name to save species jaccard distances to
# species_distances.rds.filename <- here("rds", "Ailaoshan_final_LSU_jaccard_species.rds")
# species_distances.rds.filename <- here("rds", "Ailaoshan_final_SSU_jaccard_species.rds")

########################################################################################
# get model data
# (take polygon and OTU labels from jags.data$model.data$y or jags.data$model.data$z.start as stored in this rdata file)
load(modeldata.rdata.filename)

# get modelling output
model.output <- readRDS(file = modeloutput.rds.filename)

########################################################################################
# get PolygonIDs (i.e. polygons with location information only) and species to calculate Jaccard distances for

all.PolygonIDs <- dimnames(jags.data$model.data$y)[[1]]

real.PolygonID.indexes <- which(!is.na(as.numeric(all.PolygonIDs)))
names(real.PolygonID.indexes) <- all.PolygonIDs[real.PolygonID.indexes]

OTU.names <- dimnames(jags.data$model.data$y)[[3]]

num.interations <- dim(model.output$BUGSoutput$sims.list$z)[1]

########################################################################################
# get iterations to use for averaging

    # use to do only a subset of iterations
    # set.seed(42)
    # num.interations <- 100
    # random.iterations <- sample(x = 1:350000, size = num.interations, replace = FALSE)

    num.interations <- dim(model.output$BUGSoutput$sims.list$z)[1]

########################################################################################
# functions to calculate pairwise jaccard distances

jaccard.site <- function (iter, site.i, site.j) {
    shared.species <- sum(model.output$BUGSoutput$sims.list$z[iter,site.i,] * model.output$BUGSoutput$sims.list$z[iter,site.j,])
    total.species <- sum(model.output$BUGSoutput$sims.list$z[iter,site.i,]) + sum(model.output$BUGSoutput$sims.list$z[iter,site.j,]) - shared.species
    return(shared.species/total.species)
}

jaccard.species <- function (iter, species.i, species.j) {
    # note: performs calculation on 'real' polygon IDs, i.e. excludes those that just represent a ranger without location information
    shared.sites <- sum(model.output$BUGSoutput$sims.list$z[iter,real.PolygonID.indexes,species.i] * model.output$BUGSoutput$sims.list$z[iter,real.PolygonID.indexes,species.j])
    total.sites <- sum(model.output$BUGSoutput$sims.list$z[iter,real.PolygonID.indexes,species.i]) + sum(model.output$BUGSoutput$sims.list$z[iter,real.PolygonID.indexes,species.j]) - shared.sites
    return(shared.sites/total.sites)
}

########################################################################################
# jaccard similarities for sites

    jaccard.site.avgs <- array(NA, dim = c(length(real.PolygonID.indexes), length(real.PolygonID.indexes)), dimnames = list(names(real.PolygonID.indexes), names(real.PolygonID.indexes)))

    # only performs calculations for 'real' polygon IDs, i.e. excludes those that just represent a ranger without location information
    for (i in seq_along(real.PolygonID.indexes)) {
        for (j in 1:(i-1)) {
            jaccard.site.dists <- sapply(1:num.interations, function (X) jaccard.site(X, real.PolygonID.indexes[i], real.PolygonID.indexes[j]) )
            # jaccard.site.dists <- sapply(1:num.interations, function (X) jaccard.site(random.iterations[X], real.PolygonID.indexes[i], real.PolygonID.indexes[j]) ) # use for a subset of iterations
            jaccard.site.avgs[i,j] <- mean(jaccard.site.dists)
        }
        jaccard.site.avgs[i,i] <- 1
        jaccard.site.avgs[upper.tri(jaccard.site.avgs)] <- t(jaccard.site.avgs)[upper.tri(jaccard.site.avgs)]
    }

# jaccard similarities for species

    jaccard.species.avgs <- array(NA, dim = c(length(OTU.names), length(OTU.names)), dimnames = list(OTU.names, OTU.names))

    for (i in seq_along(OTU.names)) {
        for (j in 1:(i-1)) {
            jaccard.species.dists <- sapply(1:num.interations, function (X) jaccard.species(X, i, j) )
            # jaccard.species.dists <- sapply(1:num.interations, function (X) jaccard.species(random.iterations[X], i, j) ) # use for a subset of iterations
            jaccard.species.avgs[i,j] <- mean(jaccard.species.dists)
        }
        jaccard.species.avgs[i,i] <- 1
        jaccard.species.avgs[upper.tri(jaccard.species.avgs)] <- t(jaccard.species.avgs)[upper.tri(jaccard.species.avgs)]
    }

########################################################################################
# output results
saveRDS(jaccard.site.avgs, file = site_distances.rds.filename)
saveRDS(jaccard.species.avgs, file = species_distances.rds.filename)

# session info
writeLines(capture.output(sessionInfo()), here("sessioninfo", "Ailaoshan_posterior_jaccard.sessioninfo.txt"))
