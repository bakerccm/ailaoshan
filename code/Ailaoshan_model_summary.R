# R code for postprocessing occupancy modelling results

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

modelsummary.rds.filename <- args[3]    # file name to save results summary to
# modelsummary.rds.filename <- here("rds", "Ailaoshan_model_summary_final_LSU.rds")
# modelsummary.rds.filename <- here("rds", "Ailaoshan_model_summary_final_SSU.rds")

########################################################################################
# get model data
# (take polygon and OTU labels from jags.data$model.data$y or jags.data$model.data$z.start as stored in this rdata file)
load(modeldata.rdata.filename)

########################################################################################
# get modelling output
model.output <- readRDS(file = modeloutput.rds.filename)

########################################################################################
# parse model estimates from JAGS output

# note that
#   rownames(jags.data$model.data$z.start) == dimnames(jags.data$model.data$y)[[1]] # Polygon_ID labels
#   colnames(jags.data$model.data$z.start) == dimnames(jags.data$model.data$y)[[3]] # OTU labels

# Nsite varies by site
    Nsite.rows <- grep("Nsite\\[", rownames(model.output$BUGSoutput$summary))
    Nsite.output <- model.output$BUGSoutput$summary[Nsite.rows,]
    Nsite.rownames <- rownames(Nsite.output)
    Nsite.output <- as_tibble(Nsite.output) %>%
        # store results more efficiently since these have to be integer
            mutate(n.eff = as.integer(n.eff)) %>%
        # parse rownames into indexes representing Polygon_ID
            mutate(Polygon_ID_ = Nsite.rownames) %>%
            mutate(Polygon_ID_ = gsub("Nsite\\[","",Polygon_ID_)) %>% mutate(Polygon_ID_ = gsub("\\]","",Polygon_ID_)) %>%
            mutate(Polygon_ID_ = as.numeric(Polygon_ID_)) %>%
        # convert indexes back to our original Polygon_ID labels
            mutate(Polygon_ID = rownames(jags.data$model.data$z.start)[Polygon_ID_]) %>%
            select(Polygon_ID, everything(), -Polygon_ID_)

# occ varies by species
    estocc.rows <- grep("estimated.occupancy\\[", rownames(model.output$BUGSoutput$summary))
    estocc.output <- model.output$BUGSoutput$summary[estocc.rows,]
    estocc.rownames <- rownames(estocc.output)
    estocc.output <- as_tibble(estocc.output) %>%
        # store results more efficiently since these have to be integer
            mutate(n.eff = as.integer(n.eff)) %>%
        # parse rownames into indexes representing OTU
            mutate(OTU_ = estocc.rownames) %>%
            mutate(OTU_ = gsub("estimated.occupancy\\[","",OTU_)) %>% mutate(OTU_ = gsub("\\]","",OTU_)) %>%
            mutate(OTU_ = as.numeric(OTU_)) %>%
        # convert indexes back to our original OTU labels
            mutate(OTU = colnames(jags.data$model.data$z.start)[OTU_]) %>%
            select(OTU, everything(), -OTU_)

# z varies by species and site
    z.rows <- grep("z\\[", rownames(model.output$BUGSoutput$summary))
    z.output <- model.output$BUGSoutput$summary[z.rows,]
    z.rownames <- rownames(z.output)
    z.output <- as_tibble(z.output) %>%
        # store results more efficiently since these have to be integer
            mutate(`2.5%`=as.integer(`2.5%`), `25%`=as.integer(`25%`), `50%`=as.integer(`50%`),
               `75%`=as.integer(`75%`), `97.5%`=as.integer(`97.5%`), n.eff = as.integer(n.eff)) %>%
        # parse rownames into indexes representing Polygon_ID and OTU
            mutate(var = z.rownames) %>%
            mutate(var = gsub("z\\[","",var)) %>% mutate(var = gsub("\\]","",var)) %>%
            separate(var, sep=',', into=c("Polygon_ID_","OTU_")) %>%
            mutate(Polygon_ID_ = as.numeric(Polygon_ID_), OTU_ = as.numeric(OTU_)) %>%
        # convert indexes back to our original Polygon_ID and OTU labels
            mutate(Polygon_ID = rownames(jags.data$model.data$z.start)[Polygon_ID_]) %>%
            mutate(OTU = colnames(jags.data$model.data$z.start)[OTU_]) %>%
            select(Polygon_ID, OTU, everything(), -Polygon_ID_, -OTU_)

# beta0 varies by species
# beta0 is the species-specific constant for occupancy
# converting to probability scale should give occupancy estimate (=psi) for constant (i.e. mean) values
# for environmental covariates included in model, since predictors are scaled and centered
    beta0.rows <- grep("beta0\\[", rownames(model.output$BUGSoutput$summary))
    beta0.output <- model.output$BUGSoutput$summary[beta0.rows,]
    beta0.rownames <- rownames(beta0.output)
    beta0.output <- as_tibble(beta0.output) %>%
        # store results more efficiently since these have to be integer
            mutate(n.eff = as.integer(n.eff)) %>%
        # parse rownames into indexes representing OTU
            mutate(OTU_ = beta0.rownames) %>%
            mutate(OTU_ = gsub("beta0\\[","",OTU_)) %>% mutate(OTU_ = gsub("\\]","",OTU_)) %>%
            mutate(OTU_ = as.numeric(OTU_)) %>%
        # convert indexes back to our original OTU labels
            mutate(OTU = colnames(jags.data$model.data$z.start)[OTU_]) %>%
            select(-OTU_) %>%
        # add columns for beta0 values converted to probability scale (plogis = inverse logit)
            mutate(prob_mean = plogis(mean), `prob_2.5%` = plogis(`2.5%`), `prob_50%` = plogis(`50%`), `prob_97.5%` = plogis(`97.5%`)) %>%
            select(OTU, prob_mean, `prob_2.5%`, `prob_50%`, `prob_97.5%`, everything())

# beta varies by species
# beta is indexed as beta[x,k] where x={1,2,3} denotes which slope and k denotes the species
    beta.rows <- grep("beta\\[", rownames(model.output$BUGSoutput$summary))
    beta.output <- model.output$BUGSoutput$summary[beta.rows,]
    beta.rownames <- rownames(beta.output)
    beta.output <- as_tibble(beta.output) %>%
        # store results more efficiently since these have to be integer
            mutate(n.eff = as.integer(n.eff)) %>%
        # parse rownames into indexes representing coefficient and OTU
            mutate(var = beta.rownames) %>%
            mutate(var = gsub("beta\\[","",var)) %>% mutate(var = gsub("\\]","",var)) %>%
            separate(var, sep=',', into=c("occupancy.covariate_","OTU_")) %>%
            mutate(occupancy.covariate_ = as.numeric(occupancy.covariate_), OTU_ = as.numeric(OTU_)) %>%
        # convert OTU index back to our original OTU labels and occupancy.covariate index back to original covariate labels
            mutate(OTU = colnames(jags.data$model.data$z.start)[OTU_]) %>%
            mutate(occupancy.covariate = rownames(jags.data$model.data$occ)[occupancy.covariate_]) %>%
            select(OTU, occupancy.covariate, everything(), -OTU_, -occupancy.covariate_)
    beta.output <- split(beta.output, beta.output$occupancy.covariate) # split into list of dataframes, one for each occupancy covariate

# gamma0 varies by species
# gamma0 is the species-specific constant for detection probability per 100 leeches on the logit scale
# converting to probability scale should give estimate of detection probability per 100 leeches (=r)
    gamma0.rows <- grep("gamma0\\[", rownames(model.output$BUGSoutput$summary))
    gamma0.output <- model.output$BUGSoutput$summary[gamma0.rows,]
    gamma0.rownames <- rownames(gamma0.output)
    gamma0.output <- as_tibble(gamma0.output) %>%
        # store results more efficiently since these have to be integer
            mutate(n.eff = as.integer(n.eff)) %>%
        # parse rownames into indexes representing OTU
            mutate(OTU_ = gamma0.rownames) %>%
            mutate(OTU_ = gsub("gamma0\\[","",OTU_)) %>% mutate(OTU_ = gsub("\\]","",OTU_)) %>%
            mutate(OTU_ = as.numeric(OTU_)) %>%
        # convert indexes back to our original OTU labels
            mutate(OTU = colnames(jags.data$model.data$z.start)[OTU_]) %>%
            select(-OTU_) %>%
        # add columns for gamma0 values converted to probability scale (plogis = inverse logit)
            mutate(prob_mean = plogis(mean), `prob_2.5%` = plogis(`2.5%`), `prob_50%` = plogis(`50%`), `prob_97.5%` = plogis(`97.5%`)) %>%
            select(OTU, prob_mean, `prob_2.5%`, `prob_50%`, `prob_97.5%`, everything())

    # mu.beta varies by taxonomic group
    mu.beta.rows <- grep("mu.beta\\[", rownames(model.output$BUGSoutput$summary))
    mu.beta.output <- model.output$BUGSoutput$summary[mu.beta.rows,]
    mu.beta.rownames <- rownames(mu.beta.output)
    mu.beta.output <- as_tibble(mu.beta.output) %>%
        # store results more efficiently since these have to be integer
            mutate(n.eff = as.integer(n.eff)) %>%
        mutate(covariate = mu.beta.rownames) %>%
        select(covariate, everything())

########################################################################################
# saves processed model output as a list
model.summary <- list(
    Nsite.output = Nsite.output,
    estocc.output = estocc.output,
    z.output = z.output,
    beta0.output = beta0.output,
    beta.output = beta.output,
    gamma0.output = gamma0.output,
    mu.beta.output = mu.beta.output
)
saveRDS(model.summary, file = modelsummary.rds.filename)

########################################################################################

# session info
writeLines(capture.output(sessionInfo()), here("sessioninfo", "Ailaoshan_model_summary.sessioninfo.txt"))
