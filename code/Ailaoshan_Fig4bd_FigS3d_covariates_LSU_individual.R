# R code for generating occupancy and detection predictions for individual species
# LSU dataset

library("here")
library("tidyverse")
library("R2jags")

########################################################################################
# input and output file names

    filenames <- list(

        modeldata.LSU.rdata = here("rdata", "Ailaoshan_model_data_final_LSU.rdata"), # rdata file containing data for LSU model
        modeldata.both.rdata = here("rdata", "Ailaoshan_model_data_final.rdata"), # rdata file containing data for both LSU and SSU models

        modeloutput.rds = here("rds", "Ailaoshan_model_output_final_LSU.rds"), # file containing LSU modelling results

        output.rdata = here("rdata", "Ailaoshan_occupancy_covariates_LSU_individual_predictions.rdata") # file name to save predictions to

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

########################################################################################
# extract info from model data

    num.species <- jags.data$model.data$num.species

    species.groups <- tibble(OTU = names(jags.data$model.data$g), group.no = jags.data$model.data$g) %>%
            mutate(group.name = recode(group.no, "1" = "mammals/birds", "2" = "amphibians/reptiles"))

    species.names <- leech.augmented$LSU %>% select(OTU,consensus.short, consensus.class) %>% distinct() %>% mutate(consensus.class = tolower(consensus.class))

########################################################################################
# range of unscaled and equivalent scaled values to generate predictions for

    pred.vals <- list(
        elev = data.frame(
            scaled = seq(min(jags.data$model.data$occ["elev",], na.rm = TRUE), max(jags.data$model.data$occ["elev",], na.rm = TRUE), length.out = 500),
            unscaled = seq(min(leech.augmented$LSU$elevation_median, na.rm = TRUE), max(leech.augmented$LSU$elevation_median, na.rm = TRUE), length.out = 500)),
        reserve = data.frame(
            scaled = seq(min(jags.data$model.data$occ["reserve",], na.rm = TRUE), max(jags.data$model.data$occ["reserve",], na.rm = TRUE), length.out = 500),
            unscaled = seq(min(leech.augmented$LSU$distance_to_nature_reserve_boundary, na.rm = TRUE), max(leech.augmented$LSU$distance_to_nature_reserve_boundary, na.rm = TRUE), length.out = 500)),
        numleeches = seq(1, 100, length.out = 100)
    )

# lists of empty arrays to fill with predictions

    predictions <- list(
        elev = array(NA, dim = c(500, num.species), dimnames = list(1:500, colnames(jags.data$model.data$z.start))),
        reserve = array(NA, dim = c(500, num.species), dimnames = list(1:500, colnames(jags.data$model.data$z.start))),
        numleeches = array(NA, dim = c(100, num.species), dimnames = list(1:100, colnames(jags.data$model.data$z.start)))
    )

# generate predictions
    for(k in 1:num.species){
        # empty arrays to fill with predictions for species k
            elev.pred <- array(NA, dim = c(nrow(pred.vals$elev), nsamp)) # empty array to be filled
            reserve.pred <- array(NA, dim = c(nrow(pred.vals$reserve), nsamp)) # empty array to be filled
            numleeches.pred <- array(NA, dim = c(length(pred.vals$numleeches), nsamp)) # empty array to be filled
        # generate predictions based on the mcmc iterations
            for(s in 1:nsamp){
                elev.pred[,s] <- plogis(model.output$BUGSoutput$sims.list$beta0[s,k] + model.output$BUGSoutput$sims.list$beta[s,1,k] * pred.vals$elev$scaled) # psi ~ elev for species k
                reserve.pred[,s] <- plogis(model.output$BUGSoutput$sims.list$beta0[s,k] + model.output$BUGSoutput$sims.list$beta[s,2,k] * pred.vals$reserve$scaled) # psi ~ reserve for species k
                numleeches.pred[,s] <- 1 - ((1 - plogis(model.output$BUGSoutput$sims.list$gamma0[s,k])) ^ (pred.vals$numleeches/100))     # p ~ numleeches for species k
            }
        #calculate means and store in list of prediction arrays
            predictions$elev[,k] <- apply(elev.pred, 1, mean)
            predictions$reserve[,k] <- apply(reserve.pred, 1, mean)
            predictions$numleeches[,k] <- apply(numleeches.pred, 1, mean)
        rm(elev.pred, reserve.pred, numleeches.pred)
    }

# convert arrays to dataframes

    LSU.individual.predictions <- list(
        elev = predictions$elev %>% as.data.frame() %>%
            mutate(elev.scaled = pred.vals$elev$scaled, elev.unscaled = pred.vals$elev$unscaled) %>%
            pivot_longer(-c(elev.scaled, elev.unscaled), names_to = "OTU", values_to="estimated occupancy"),
        reserve = predictions$reserve %>% as.data.frame() %>%
            mutate(reserve.scaled = pred.vals$reserve$scaled, reserve.unscaled = pred.vals$reserve$unscaled) %>%
            pivot_longer(-c(reserve.scaled, reserve.unscaled), names_to = "OTU", values_to="estimated occupancy"),
        numleeches = predictions$numleeches %>% as.data.frame() %>%
            mutate(leeches = pred.vals$numleeches) %>%
            pivot_longer(-leeches, names_to = "OTU", values_to="estimated detection")
    )

    LSU.individual.predictions <- lapply(LSU.individual.predictions, function (X) {
        X %>% left_join(species.groups %>% select(OTU, group.name), by = "OTU") %>%
            left_join(species.names, by = "OTU")
    })

########################################################################################
# save predictions to file

    save(LSU.individual.predictions, file = filenames$output.rdata)

########################################################################################

# session info
writeLines(capture.output(sessionInfo()), here("sessioninfo", "Ailaoshan_Fig4bd_FigS3d_covariates_LSU_individual.sessioninfo.txt"))
