# R code to draw Fig 2ab

library("here")
library("tidyverse")
library("vegan")
library("iNEXT")

load(file=here("rdata","Ailaoshan_OTU_table.rdata"))

# plot breakdown of OTUs by taxonomic class

    leech %>%
        distinct(dataset, consensus.class, consensus.short) %>%
        group_by(dataset, consensus.class) %>%
        mutate(consensus.class = tolower(consensus.class)) %>%
        ggplot(aes(x = consensus.class, fill=consensus.class)) + stat_count(show.legend = FALSE) + facet_wrap("dataset") +
        labs(y = "number of species") +
        theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x = element_blank())

    ggsave(filename = here("figures","Fig2a_taxonomic_breakdown.pdf"), width = 4, height=3)

    # altenative figure for SAS talk

        # leech %>%
        #     filter(dataset == "SSU") %>% select(OTU) %>% distinct() %>% nrow() # 72 OTUs in SSU dataset

        # leech %>%
        #     filter(dataset == "SSU") %>%
        #     distinct(dataset, consensus.class, consensus.short) %>%
        #     group_by(dataset, consensus.class) %>%
        #     ggplot(aes(x = consensus.class, fill=consensus.class)) + stat_count(show.legend = FALSE) +
        #     labs(x = "Class", y = "Number of species") +
        #     theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x = element_blank(),
        #         plot.background = element_rect(fill = "#EBEBEB",color = "white"), axis.ticks = element_blank())

        # ggsave(filename = here("figures","SSU_taxonomic_breakdown.pdf"), width = 2.7, height=2.7, useDingbats=FALSE)

    # alternative figures for vISEC2020 talk

        # leech %>%
        #     select(dataset, OTU) %>% distinct() %>%
        #     group_by(dataset) %>% tally() # 59 (LSU) and 72 (SSU) OTUs

        # leech %>%
        #     distinct(dataset, consensus.class, consensus.short) %>%
        #     group_by(dataset, consensus.class) %>%
        #     ggplot(aes(x = consensus.class, fill=consensus.class)) + stat_count(show.legend = FALSE) + facet_wrap("dataset") +
        #     labs(x = "Class", y = "Number of species") +
        #     theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x = element_blank(),
        #         plot.background = element_rect(fill = "#EBEBEB",color = "white"), axis.ticks = element_blank())

        # ggsave(filename = here("figures","taxonomic_breakdown.pdf"), width = 5, height=2.7, useDingbats=FALSE)

# rarefaction and extrapolation using iNext

    # collate LabID-wise incidence data in wide format

        incidence.Lab_ID <- list(
            LSU = leech %>% filter(dataset == "LSU") %>%
                    mutate(LabID.incidence = ifelse(reads > 0, 1, 0)) %>%
                    select(Lab_ID, OTU, LabID.incidence) %>%
                    spread(key = OTU, value = LabID.incidence) %>%
                    column_to_rownames("Lab_ID") %>%
                    as.matrix() %>% t(),
            SSU = leech %>% filter(dataset == "SSU") %>%
                    mutate(LabID.incidence = ifelse(reads > 0, 1, 0)) %>%
                    select(Lab_ID, OTU, LabID.incidence) %>%
                    spread(key = OTU, value = LabID.incidence) %>%
                    column_to_rownames("Lab_ID") %>%
                    as.matrix() %>% t(),
            combined = leech %>%
                    # summarize by consensus.short to avoid double counting OTUs when both datasets included
                        group_by(Lab_ID, consensus.short) %>%
                        summarize(LabID.incidence = ifelse(sum(reads) > 0, 1, 0), .groups = "drop") %>%
                    # convert to wide format data
                        select(Lab_ID, consensus.short, LabID.incidence) %>%
                        spread(key = consensus.short, value = LabID.incidence, fill = 0) %>%
                        column_to_rownames("Lab_ID") %>%
                        as.matrix() %>% t()
        )

    # perform iNext calculations

        iNext.LabID <- iNEXT(incidence.Lab_ID, q=0, datatype="incidence_raw", endpoint = 1500, knots = 150)

    # asymptotic diversity estimates

        iNext.LabID$AsyEst

    # default ggplot for iNext objects
    #   ggiNEXT(iNext.LabID, type=1)

    # convert iNext to dataframe and separate out elements for plotting
    # (easier to do this than to use default and override settings to match other plots in manuscript)

        iNext.LabID.df <- fortify(iNext.LabID, type=1) %>%
            rename(dataset = site) %>%
            mutate(dataset = factor(dataset, c("combined", "SSU", "LSU")))

        iNext.LabID.df.point <- iNext.LabID.df %>% filter(method == "observed")

        iNext.LabID.df.line <- iNext.LabID.df %>% filter(method != "observed") %>%
            mutate(method = factor(method, c("interpolated", "extrapolated"), c("interpolation", "extrapolation")))

    # draw iNext results using ggplot

        ggplot(iNext.LabID.df, aes(x = x, y = y, col = dataset)) +
            geom_point(data = iNext.LabID.df.point) +
            geom_line(data = iNext.LabID.df.line, aes(linetype = method)) +
            geom_ribbon(aes(ymin = y.lwr, ymax = y.upr, fill = dataset), alpha = 0.2, col = NA) +
            theme(legend.justification = c(1,0), legend.position = c(0.96,0.05), legend.key = element_blank()) +
            labs(x = "number of replicates", y = "species richness") + guides(linetype = FALSE)

        ggsave(here("figures","Fig2b_LabID_rarefaction.pdf"), width = 3.5, height = 3.5)

# save session info
writeLines(capture.output(sessionInfo()), here("sessioninfo", "Ailaoshan_Fig2ab_taxonomy_rarefaction.sessioninfo.txt"))
