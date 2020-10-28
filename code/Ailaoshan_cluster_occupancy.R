# R code for getting species-level occupancies for each cluster of sites

# This file take inputs <modeldata.rdata>, <modeloutput.rds> and <siteclusters.rds>
# and saves estimated occupancy for each species in each cluster as <output.rds>

library("here")
library("tidyverse")
library("R2jags")

# some useful tools for tracking memory usage
# see http://adv-r.had.co.nz/memory.html
# library("pryr")

########################################################################################
# get file names from command line arguments
args <- commandArgs(trailingOnly = TRUE)

modeldata.rdata.filename <- args[1]  # file containing model data
modeloutput.rds.filename <- args[2]  # file containing modelling results
siteclusters.rds.filename <- args[3]  # file containing site clusters
output.rds.filename <- args[4]  # file name to save site jaccard distances to

########################################################################################
# model input data
# (take polygon and OTU labels from jags.data$model.data$y or jags.data$model.data$z.start as stored in this rdata file)

    load(modeldata.rdata.filename)

# modelling output

    model.output <- readRDS(file = modeloutput.rds.filename)

# site clusters

    site.clusters <- readRDS(file = siteclusters.rds.filename)
    site.clusters <- site.clusters %>% mutate(jaccard.cut = as.numeric(jaccard.cut))

    n.clusters <- site.clusters$jaccard.cut %>% unique() %>% length()

########################################################################################
# extract model sims and discard other model output to save memory

    sims.list <- model.output$BUGSoutput$sims.list$z
    # this take almost no additional memory:
    # mem_change(sims.list <- model.output$BUGSoutput$sims.list$z)
    # >792 B

    rm(model.output)
    # this frees up a lot of space:
    # mem_change(rm(model.output))
    # >-73.4 GB

########################################################################################
# attach PolygonIDs and OTUs from model input data to sims.list as dimnames

    # jags.data$model.data$y is indexed as (PolygonIDs, replicates, species)
    # dim(jags.data$model.data$y)
    # [1] 209  40  59

    # sims.list is indexed as (iterations, PolygonIDs, species)
    # dim(sims.list)
    # [1] 350000    209     59

    n.iterations <- dim(sims.list)[1] # 350,000 iterations

    n.species <- dim(sims.list)[3]

    dimnames(sims.list) <- list(
        1:n.iterations,
        dimnames(jags.data$model.data$y)[[1]],
        dimnames(jags.data$model.data$y)[[3]]
    )

########################################################################################
# which sites are in each cluster?

    cluster.PolygonIDs <- lapply(1:n.clusters, function (X) site.clusters$Polygon_ID[site.clusters$jaccard.cut == X])

# which columns in sims.list correspond to these sites?

    cluster.indexes <- lapply(1:n.clusters, function (X) which(dimnames(sims.list)[[2]] %in% cluster.PolygonIDs[[X]]))

# how many sites are in each cluster?

    cluster.sitecounts <- sapply(1:n.clusters, function (X) sum(site.clusters$jaccard.cut == X))

    cluster.sitecounts

########################################################################################
# for each MCMC iteration, get fraction of sites in each cluster occupied by each species

    cluster.occ <- array(NA, dim = c(n.iterations, n.clusters, n.species))

    dimnames(cluster.occ) <- list(1:n.iterations, 1:n.clusters, dimnames(sims.list)[[3]])

    for (cluster in 1:n.clusters) {
        for (species in 1:n.species) {
            cluster.occ[,cluster,species] <- (sims.list[,cluster.indexes[[cluster]],species] / cluster.sitecounts[cluster]) %>% apply(MARGIN = 1, sum)
        }
    }

    # inspect results
    # cluster.occ[1:10,,1:2]

########################################################################################
# summarize results and save to file

    # get means and 95% Bayesian confidence intervals
    cluster.occ.summary = list(
        mean = apply(cluster.occ, MARGIN = c(2,3), mean),
        BCI = apply(cluster.occ, MARGIN = c(2,3), function (X) quantile(X, probs = c(0.025,0.975)))
    )

    # dim(cluster.occ.summary$mean)
    # [1]  3 59

    # dim(cluster.occ.summary$BCI)
    # [1]  2  3 59

    # convert summary results to list of data.frames
    cluster.occ.summary.dflist <- list(
        mean = cluster.occ.summary$mean[,] %>%
            as.data.frame() %>%
            rownames_to_column(var = "cluster") %>%
            pivot_longer(!cluster, names_to = "OTU", values_to = "mean"),
        "2.5%" = cluster.occ.summary$BCI["2.5%",,] %>%
            as.data.frame() %>%
            rownames_to_column(var = "cluster") %>%
            pivot_longer(!cluster, names_to = "OTU", values_to = "2.5%"),
        "97.5%" = cluster.occ.summary$BCI["97.5%",,] %>%
            as.data.frame() %>%
            rownames_to_column(var = "cluster") %>%
            pivot_longer(!cluster, names_to = "OTU", values_to = "97.5%")
    )

    # join data.frames ready for export
    cluster.occ.summary.df <- cluster.occ.summary.dflist[["mean"]] %>%
        full_join(cluster.occ.summary.dflist[["2.5%"]], by = c("cluster", "OTU")) %>%
        full_join(cluster.occ.summary.dflist[["97.5%"]], by = c("cluster", "OTU")) %>%
        select(OTU, cluster, everything())

    # save output
    saveRDS(cluster.occ.summary.df, file = output.rds.filename)

########################################################################################

# session info
writeLines(capture.output(sessionInfo()), here("sessioninfo", "Ailaoshan_cluster_occupancy.sessioninfo.txt"))
