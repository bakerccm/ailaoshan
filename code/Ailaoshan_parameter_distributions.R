# R code to get parameter posterior distributions from initial occupancy model run to help constrain priors in Bayesian variable selection

library("here")
library("tidyverse")
library("R2jags")

# get file names from command line arguments

    args <- commandArgs(trailingOnly = TRUE)
    input.filename <- args[1]    # file containing model results
    output.filename <- args[2]    # file name to save results summary to

# read in model data

    model.output <- readRDS(file = input.filename)

# convert to mcmc so coda can read it

    output.mcmc <- as.mcmc(model.output)

# get sims and put in data.frame

    # get variable names

        mu.varnames <- varnames(output.mcmc) %>% grep("^mu.", ., value = TRUE)
        sigma.varnames <- varnames(output.mcmc) %>% grep("^sigma.", ., value = TRUE)

    # get sims

        mu.sims <- lapply(mu.varnames, function (X) output.mcmc[,X] %>% unlist())
        names(mu.sims) <- mu.varnames
        mu.sims <- do.call(cbind, mu.sims) %>% as.data.frame() %>% pivot_longer(everything(), names_to = "variable", values_to = "sims")
        mu.sims.summary <- mu.sims %>% group_by(variable) %>% summarize(mean = mean(sims), sd = sd(sims), .groups = "drop")

        sigma.sims <- lapply(sigma.varnames, function (X) output.mcmc[,X] %>% unlist())
        names(sigma.sims) <- sigma.varnames
        sigma.sims <- do.call(cbind, sigma.sims) %>% as.data.frame() %>% pivot_longer(everything(), names_to = "variable", values_to = "sims")
        sigma.sims.summary <- sigma.sims %>%
            mutate(sims = sims * (2 * rbinom(n(), size = 1, prob = 0.5) - 1)) %>%
            group_by(variable) %>% summarize(mean = mean(sims), sd = sd(sims), .groups = "drop")

# send output to file
    sink(output.filename)
    print(mu.sims.summary, n = Inf) # set n = Inf to print all rows of tibble
    print(sigma.sims.summary, n = Inf) # set n = Inf to print all rows of tibble
    sink()

# session info
writeLines(capture.output(sessionInfo()), here("sessioninfo", "Ailaoshan_parameter_distributions.sessioninfo.txt"))
