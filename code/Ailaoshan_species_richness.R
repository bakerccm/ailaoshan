# R code to examine difference between species richness and community occupancy

library("here")
library("tidyverse")
library("R2jags")

########################################################################################
# input and output file names

    filenames <- list(

        modeldata.LSU.rdata = here("rdata", "Ailaoshan_model_data_final_LSU.rdata"), # rdata file containing data for SSU model
        modeldata.both.rdata = here("rdata", "Ailaoshan_model_data_final.rdata"), # rdata file containing data for both LSU and SSU models

        modeloutput.rds = here("rds", "Ailaoshan_model_output_final_LSU.rds") # file containing modelling results

    )

########################################################################################
# get data

    # file containing LSU model data
    # (take polygon and OTU labels from jags.data$model.data$y or jags.data$model.data$z.start as stored in this rdata file)
    load(filenames$modeldata.LSU.rdata)

    # file containing model data prior to pulling out only the data required for the models
    # just get leech.augmented from this .rdata file for the unscaled covariate values
    both = new.env()
    load(filenames$modeldata.both.rdata, envir = both)
    leech.augmented = both$leech.augmented
    rm(both)

    # modelling output
    model.output <- readRDS(file = filenames$modeloutput.rds)

########################################################################################
# number of MCMC samples

    # sims.list has model output arranged by variable
    # and appears to be a rearranged version of sims.matrix
    nsamp <- nrow(model.output$BUGSoutput$sims.list[[1]]) # nrow should be the same for any item in the list

    nsites <- dim(model.output$BUGSoutput$sims.list$z)[2]
    nspec <- dim(model.output$BUGSoutput$sims.list$z)[3]

    nmammals <- sum(jags.data$model.data$g == 1)
    nfrogs <- sum(jags.data$model.data$g == 2)
########################################################################################
# occupancy estimates

# elevation

    # predictor values
    # jags.data$model.data$occ[1,] is elev
    # jags.data$model.data$occ[2,] is reserve

    # jags.data$model.data$occ["elev",]
    # jags.data$model.data$occ["reserve",]
    # note that these do have colnames

    community.pred <- rep(NA, nsites)
    names(community.pred) <- colnames(jags.data$model.data$occ)

    # make sites with NA = zero
    # omit this if you just want to exclude those points
        jags.data$model.data$occ[is.na(jags.data$model.data$occ)] <- 0

    # posterior mean community occupancy per site
    for(site in 1:nsites){
            mammals <- plogis(model.output$BUGSoutput$sims.list$mu.eta[,1,1] + model.output$BUGSoutput$sims.list$mu.beta[,1] * jags.data$model.data$occ["elev",site] + model.output$BUGSoutput$sims.list$mu.beta[,2] * jags.data$model.data$occ["reserve",site])
            frogs <- plogis(model.output$BUGSoutput$sims.list$mu.eta[,1,2] + model.output$BUGSoutput$sims.list$mu.beta[,1] * jags.data$model.data$occ["elev",site] + model.output$BUGSoutput$sims.list$mu.beta[,2] * jags.data$model.data$occ["reserve",site])
            community <- (mammals * nmammals/nspec) + (frogs * nfrogs/nspec)
            community.pred[site] <- mean(community)
    }

    Nsite <- model.output$BUGSoutput$summary[grep("^Nsite", rownames(model.output$BUGSoutput$summary)),]

    plot(Nsite[,"mean"], community.pred)

    cor(Nsite[,"mean"], community.pred, use="complete.obs")
    # [1] 0.9345153


    sites.occupied.LSU002 <- apply(jags.data$model.data$y[,,"LSU002"],MAR =1 ,FUN = function (X) sum(X,na.rm=TRUE)) > 0


###
    z.mean.LSU002 <- apply(model.output$BUGSoutput$sims.list$z[,,1] ,MAR =2 ,mean)

    data.frame(sites.occupied = sites.occupied.LSU002, z.mean = z.mean.LSU002) %>% arrange(sites.occupied)
###

# session info
writeLines(capture.output(sessionInfo()), here("sessioninfo", "Ailaoshan_species_richness.sessioninfo.txt"))
