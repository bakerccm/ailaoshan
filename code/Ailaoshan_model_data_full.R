# R code to prepare to run initial occupancy models to gauge values for parameters

library("here")
library("tidyverse")

# get raw data

    load(file=here("rdata", "Ailaoshan_OTU_table.rdata"))

# augment read data and convert to presence-absence

    leech.augmented <- sapply(c("LSU", "SSU"), simplify = FALSE, FUN = function (X) {
        leech %>%
            # add extra rows for Polygon_IDs with missing data
                bind_rows(leech.supplement) %>%
            # split by dataset
                filter(dataset == X) %>%
            # convert reads to presence-absence
                mutate(present = as.integer(ifelse(reads > 0, 1, 0))) %>%
                select(-reads) %>%
            # add counts of other taxa
                group_by(Polygon_ID, Lab_ID) %>%
                mutate(other = sum(present) - present) %>%
                ungroup() %>%
            # scale variables
                mutate(leech_qty_scaled = scale(leech_qty))
    })

# make empty lists of arrays to be filled with data

    # get max dimensions for sizing data arrays

        # number of sites i.e. Polygon_IDs for each dataset
        num.sites <- sapply(names(leech.augmented), simplify = FALSE,
            FUN = function (X) leech.augmented[[X]] %>% select(Polygon_ID) %>% distinct() %>% nrow() )

        # (maximum) number of replicates per Polygon_ID in each dataset
        num.reps <- sapply(names(leech.augmented), simplify = FALSE,
            FUN = function (X) leech.augmented[[X]] %>% group_by(Polygon_ID, OTU) %>% tally() %>% group_by(Polygon_ID) %>% summarize(max.n = max(n), .groups = "drop") %>% pull(max.n) %>% max() )

        # total number of observed species for each dataset
        num.species <- sapply(names(leech.augmented), simplify = FALSE,
            FUN = function (X) leech.augmented[[X]] %>% select(OTU) %>% distinct() %>% nrow() )

    # make empty data arrays to be filled later

        # species presences for each site, to be indexed as [num.sites, num.reps, num.species]
        y <- sapply(names(leech.augmented), simplify = FALSE, FUN = function (X) array(NA, dim=c(num.sites[[X]], num.reps[[X]], num.species[[X]])) )

        # number of leeches per replicate, to be indexed as [num.sites, num.reps]
        numleeches <- sapply(names(leech.augmented), simplify = FALSE, FUN = function (X) array(NA, dim=c(num.sites[[X]], num.reps[[X]])) )

        # group for each species (mammals/birds, amphibians/reptiles)
        g <- sapply(names(leech.augmented), simplify = FALSE, FUN = function (X) rep(NA, num.species[[X]]) )

    # check data dimensions
    # (note that num.sites now includes Polygon_IDs that never had any reads but also extra Polygon_IDs
    # derived from Ranger_IDs, so exceeds number of actual polygons
        unlist(num.sites) #LSU: 209, SSU: 209
        unlist(num.reps) #LSU: 40, SSU: 38
        unlist(num.species) #LSU: 59, SSU: 72

# prepare data

    # names for ensuring that arrays are aligned

        # species
        OTU.labels <- sapply(names(leech.augmented), simplify = FALSE,
            FUN = function (X) leech.augmented[[X]] %>% select(OTU) %>% distinct() %>% arrange(OTU) %>% pull(OTU) )

        # sites
        polygon.labels <- sapply(names(leech.augmented), simplify = FALSE,
            FUN = function (X) leech.augmented[[X]] %>% select(Polygon_ID) %>% distinct() %>% arrange(Polygon_ID) %>% pull(Polygon_ID) )

        # replicates
        replicate.labels <- sapply(names(leech.augmented), simplify = FALSE,
            FUN = function (X) leech.augmented[[X]] %>% select(replicate_no) %>% distinct() %>% arrange(replicate_no) %>% pull(replicate_no) )

    num.species.groups <- list()    # 2 -- i.e. Mammals/Birds, Amphibians/Reptiles
    species.groups <- list()        # i.e. Mammals/Birds, Amphibians/Reptiles

    for (i in names(leech.augmented)) {

        # fill y[[i]] with data
            temp <- leech.augmented[[i]] %>%
                select(Polygon_ID, OTU, replicate_no, present) %>%
                pivot_wider(names_from = replicate_no, values_from = present)
            for (k in seq_along(OTU.labels[[i]])) {
                # this is just the wide data from before
                    temp.k <- temp %>% filter(OTU == OTU.labels[[i]][k]) # filter to the kth OTU
                # arrange rows and then remove Polygon_ID column
                    temp.k <- temp.k[match(polygon.labels[[i]], temp.k$Polygon_ID),] %>% select(starts_with("replicate_"))
                # arrange columns
                    temp.k <- temp.k[,match(replicate.labels[[i]], names(temp.k))]
                # convert to matrix
                    y[[i]][,,k] <- temp.k %>% data.matrix()
            }
            rm(temp, k, temp.k)
            dimnames(y[[i]]) <- list(polygon.labels[[i]], replicate.labels[[i]], OTU.labels[[i]])

        # fill numleeches[[i]] with data
            temp <- leech.augmented[[i]] %>%
                select(dataset, Polygon_ID, replicate_no, leech_qty) %>%
                distinct() %>% pivot_wider(names_from = replicate_no, values_from = leech_qty) ### N.B. using unscaled leech_qty ###
            temp <- temp[match(polygon.labels[[i]], temp$Polygon_ID),]
            temp <- temp %>% select(starts_with("replicate_"))
            temp <- temp[,match(replicate.labels[[i]], names(temp))]
            numleeches[[i]] <- temp %>% data.matrix()
            rm(temp)
            dimnames(numleeches[[i]]) <- list(polygon.labels[[i]], replicate.labels[[i]])

        # fill g[[i]] with species group data in numeric form (i.e. group that each species belongs to: mammals/birds, amphibians/reptiles)
            num.species.groups[[i]] <- 2
            species.groups[[i]] <- list(Mammals = 1, Birds = 1, Amphibians = 2, Reptiles = 2)
            temp <- leech.augmented[[i]] %>% select(dataset, OTU, consensus.class) %>% distinct()
            temp <- temp[match(OTU.labels[[i]], temp$OTU),]
            g[[i]] <- temp %>%
                mutate(g = recode(consensus.class, !!!species.groups[[i]])) %>%
                mutate(g = as.integer(g)) %>% pull(g)
            rm(temp)
            names(g[[i]]) <- OTU.labels[[i]]
    }

    rm(i, OTU.labels, polygon.labels, replicate.labels)

# observed occurrence values as starting values for z

    # note that these vary by dataset

    # note also that we have some site/species combinations that don't even have one replicate
    #    any(is.na(y$LSU[,1,]))
    #    any(is.na(y$SSU[,1,]))

    z.start <- sapply(names(leech.augmented), simplify = FALSE, function (X) {
        suppressWarnings( # for the 'all NA' sites, max returns -Inf with a warning
            zst <- apply(y[[X]], c(1,3), max, na.rm = TRUE) # observed occurrences as starting values for z
        )
        zst[!is.finite(zst)] <- 1 # performs replacement for values of -Inf (arising from sites with all NAs, i.e. extra Polygon_IDs without any replicates) as well as NA
        return(zst)
    })

# extract and scale per-site occupancy covariates

    site.covariates.scaled <- sapply(names(leech.augmented), simplify = FALSE, function (X) {
        env.data = leech.augmented[[X]] %>%
            select(Polygon_ID, elevation_median, tpi_median, distance_to_road_median, distance_to_stream_median, distance_to_nature_reserve_boundary) %>%
            distinct()
        site.covs = tibble(Polygon_ID = y[[X]] %>% rownames()) %>% left_join(env.data, by = "Polygon_ID") # ensure row order is same as y
        site.covs.scaled <- scale(site.covs[,-1])
        rownames(site.covs.scaled) <- site.covs$Polygon_ID
        return(site.covs.scaled)
    })

# MCMC settings

    mcmc.settings <- list(
        ni = 50000,
        nt = 1,     # no thinning
        nb = 10000,
        nc = 3
    )

# package data for JAGS

    model.data <- sapply(names(leech.augmented), simplify = FALSE, function (X) {
        list(y = y[[X]], num.sites = num.sites[[X]], num.reps = num.reps[[X]], num.species = num.species[[X]],
            # occupancy covariates
                occ = rbind(elev = site.covariates.scaled[[X]][,"elevation_median"],
                    tpi = site.covariates.scaled[[X]][,"tpi_median"],
                    road = site.covariates.scaled[[X]][,"distance_to_road_median"],
                    stream = site.covariates.scaled[[X]][,"distance_to_stream_median"],
                    reserve = site.covariates.scaled[[X]][,"distance_to_nature_reserve_boundary"]),
                    occupancy.slopes = 5,  # number of occupancy slope coefficients
            # sampling covariates
                numleeches = numleeches[[X]],
            # other
                z.start = z.start[[X]],
                g = g[[X]], # group for each species (1 = Mammals/Birds, 2 = Amphibians/Reptiles)
                num.species.groups = num.species.groups[[X]] # number of different species groups, i.e. 2
        )
    })

# initial values

    # note that z.start is evaluated at runtime, meaning that it will be the object stored in the model data
    inits <- function() list(z = z.start)

# parameters to monitor

    params <- c("delta.beta", "delta.gamma", "mu.beta", "sigma.beta", "mu.eta", "sigma.eta")

# save data to files for modelling

    datasets <- c("LSU", "SSU")

    for (d in datasets) {
        jags.data <- list(
            model.data = model.data[[d]],
            inits = inits,
            params = params,
            mcmc.settings = mcmc.settings
        )
        jags.data.filename <- here("rdata", paste0("Ailaoshan_model_data_full_", d, ".rdata"))
        save(jags.data, file=jags.data.filename)
    }

# save prepared MCMC data

    save(inits, leech, leech.supplement, leech.augmented, mcmc.settings, model.data, numleeches,
         params, site.covariates.scaled, y, z.start, file = here("rdata", "Ailaoshan_model_data_full.rdata"))

# session info
writeLines(capture.output(sessionInfo()), here("sessioninfo", "Ailaoshan_model_data_full.sessioninfo.txt"))
