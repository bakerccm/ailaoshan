# R code for generating Rhat plots to assess convergence

library("here")
library("tidyverse")

########################################################################################
# get model summary filenames

    LSU.modelsummary.rds.filename <- here("rds", "Ailaoshan_model_summary_final_LSU.rds")

    SSU.modelsummary.rds.filename <- here("rds", "Ailaoshan_model_summary_final_SSU.rds")

########################################################################################
# get model summaries

    LSU.model.summary <- readRDS(file = LSU.modelsummary.rds.filename)

    SSU.model.summary <- readRDS(file = SSU.modelsummary.rds.filename)

    # join summaries from two models together
    Nsite.output <- list(LSU  = LSU.model.summary$Nsite.output, SSU  = SSU.model.summary$Nsite.output)
    estocc.output <- list(LSU  = LSU.model.summary$estocc.output, SSU  = SSU.model.summary$estocc.output)
    z.output <- list(LSU  = LSU.model.summary$z.output, SSU  = SSU.model.summary$z.output)
    beta0.output <- list(LSU  = LSU.model.summary$beta0.output, SSU  = SSU.model.summary$beta0.output)
    beta.output <- list(LSU  = LSU.model.summary$beta.output, SSU  = SSU.model.summary$beta.output)
    gamma0.output <- list(LSU  = LSU.model.summary$gamma0.output, SSU  = SSU.model.summary$gamma0.output)

    # remove individual model summaries
    rm(LSU.model.summary, SSU.model.summary)

########################################################################################
## check model convergence

    # Rhat values should be close to 1 for convergence, and they are

    # Rhat due to Gelman and Rubin 1992 (aka Brooks-Gelman-Rubin statistic; see also Gelman et al 2014)

    # LSU

        LSU_Rhat <- bind_rows(
            Nsite = Nsite.output$LSU %>% select(Rhat),
            estocc = estocc.output$LSU %>% select(Rhat),
            z = z.output$LSU %>% select(Rhat),
            beta0 = beta0.output$LSU %>% select(Rhat),
            beta1 = beta.output$LSU$elev %>% select(Rhat),
            beta2 = beta.output$LSU$reserve %>% select(Rhat),
            gamma0.output = gamma0.output$LSU %>% select(Rhat),
            .id = "model.parameter"
        )

        LSU_Rhat %>%
            ggplot() + geom_histogram(aes(x=Rhat), bins = 12) + facet_wrap(~model.parameter, scales = "free")

        ggsave(filename = here("pdfs", "Ailaoshan_model_Rhat_LSU.pdf"), useDingbats = FALSE)

    # SSU

        SSU_Rhat <- bind_rows(
            Nsite = Nsite.output$SSU %>% select(Rhat),
            estocc = estocc.output$SSU %>% select(Rhat),
            z = z.output$SSU %>% select(Rhat),
            beta0 = beta0.output$SSU %>% select(Rhat),
            beta1 = beta.output$SSU$elev %>% select(Rhat),
            gamma0.output = gamma0.output$SSU %>% select(Rhat),
            .id = "model.parameter"
        )

        SSU_Rhat %>%
            ggplot() + geom_histogram(aes(x=Rhat), bins = 12) + facet_wrap(~model.parameter, scales = "free")

        ggsave(filename = here("pdfs", "Ailaoshan_model_Rhat_SSU.pdf"), useDingbats = FALSE)

########################################################################################

# session info
writeLines(capture.output(sessionInfo()), here("sessioninfo", "Ailaoshan_model_Rhat.sessioninfo.txt"))
