# R code to get posterior probabilities of different models after Bayesian variable selection

library("here")
library("tidyverse")
library("R2jags")

########################################################################################
# get filenames from command line arguments

args <- commandArgs(trailingOnly = TRUE)

model.output.filename <- args[1]    # .rds file containing bernoulli model output
model.probabilities.filename <- args[2]    # .txt file to write model probabilities table to

########################################################################################

model.output <- readRDS(file = here("rds", model.output.filename)) # filename from command line args

# convert to mcmc so coda can read it
output.mcmc <- as.mcmc(model.output)

sims = model.output$BUGSoutput$sims.list$delta.beta
colnames(sims) <- c("elev", "tpi", "road", "stream", "reserve")
delta.beta.tbl = as_tibble(sims)

# inspect covariances of delta.beta's
# (delta.beta.cov <- cov(delta.beta.tbl))

# calculate table with probabilities of different models
delta.beta.model.results <- delta.beta.tbl %>%
    mutate(model_code = paste0(elev, tpi, road, stream, reserve)) %>%
    mutate(model_string = paste0(
        ifelse(elev == 1, "elev + ", ""),
        ifelse(tpi == 1, "tpi + ", ""),
        ifelse(road == 1, "road + ", ""),
        ifelse(stream == 1, "stream + ", ""),
        ifelse(reserve == 1, "reserve + ", ""))) %>%
    mutate(model_string = str_replace(model_string, " \\+ $", "")) %>%
    select(model_code, model_string) %>% group_by(model_code, model_string) %>% tally() %>% ungroup() %>%
    arrange(desc(n)) %>% mutate(prob = n/sum(n)) %>% select(-n)

sink(file = here("modelprobs", model.probabilities.filename)) # filename from command line args
print(delta.beta.model.results, n = Inf)
sink()

# session info
writeLines(capture.output(sessionInfo()), here("sessioninfo", "Ailaoshan_model_probabilities.sessioninfo.txt"))
