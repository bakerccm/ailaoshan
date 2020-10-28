# snakefile for workflow from Ailaoshan study
# Ji et al. (2020) "Measuring protected-area outcomes with leech iDNA: large-scale quantification of vertebrate biodiversity in Ailaoshan reserve"

configfile: 'config/config.yaml'

# all output files
rule all:
    input:
        'fasta/FileS3_LSU_representative_sequences.fasta', 'fasta/FileS4_SSU_representative_sequences.fasta', 'figures/Fig1a_Ailaoshan_map_continent.pdf', 'figures/Fig1b_Ailaoshan_site_colour.png', 'figures/Fig1b_Ailaoshan_site_grey.png', 'figures/Fig2a_taxonomic_breakdown.pdf', 'figures/Fig2b_LabID_rarefaction.pdf', 'figures/Fig2c_occupancy_detection_labels.pdf', 'figures/Fig2c_occupancy_detection.pdf', 'figures/Fig3a_LSU_observed_richness.pdf', 'figures/Fig3b_SSU_observed_richness.pdf', 'figures/Fig3c_LSU_estimated_richness.pdf', 'figures/Fig3d_SSU_estimated_richness.pdf', 'figures/Fig3e_LSU_histogram.pdf', 'figures/Fig3e_LSU_scatter.pdf', 'figures/Fig3f_SSU_histogram.pdf', 'figures/Fig3f_SSU_scatter.pdf', 'figures/Fig4_occupancy_covariates.pdf', 'figures/Fig5a_nmds_sites_LSU_labels.pdf', 'figures/Fig5a_nmds_sites_LSU.pdf', 'figures/Fig5b_nmds_sites_SSU_labels.pdf', 'figures/Fig5b_nmds_sites_SSU.pdf', 'figures/Fig5c_map_sites_LSU.pdf', 'figures/Fig5d_map_sites_SSU.pdf', 'figures/Fig6ab_cluster_occupancy.pdf', 'figures/FigS10_SSU_cluster_occupancy_10kg.pdf', 'figures/FigS10_SSU_cluster_occupancy_domestic.pdf', 'figures/FigS10_SSU_cluster_occupancy.pdf', 'figures/FigS1a_elevation.pdf', 'figures/FigS1b_elevation.pdf', 'figures/FigS1c_TPI.pdf', 'figures/FigS1d_TPI.pdf', 'figures/FigS1e_road.pdf', 'figures/FigS1f_road.pdf', 'figures/FigS1g_stream.pdf', 'figures/FigS1h_stream.pdf', 'figures/FigS1i_reserve.pdf', 'figures/FigS1j_reserve.pdf', 'figures/FigS3a_LSU_richhist.pdf', 'figures/FigS3b_SSU_richhist.pdf', 'figures/FigS3cdef_detection_covariates.pdf', 'figures/FigS4a_obsvrich_LabID.pdf', 'figures/FigS4b_richness_numleeches.pdf', 'figures/FigS4c_obsvrich_PolygonID.pdf', 'figures/FigS4d_observed_richness_sampling.pdf', 'figures/FigS4e_estimated_richness_PolygonID.pdf', 'figures/FigS4f_estimated_richness_sampling.pdf', 'figures/FigS5_LSU_elev_10kg.pdf', 'figures/FigS5_LSU_elev_domestic.pdf', 'figures/FigS5_LSU_elev.pdf', 'figures/FigS6_LSU_reserve_10kg.pdf', 'figures/FigS6_LSU_reserve_domestic.pdf', 'figures/FigS6_LSU_reserve.pdf', 'figures/FigS7_SSU_elev_10kg.pdf', 'figures/FigS7_SSU_elev_domestic.pdf', 'figures/FigS7_SSU_elev.pdf', 'figures/FigS8a_tree_sites_LSU.pdf', 'figures/FigS8b_tree_sites_SSU.pdf', 'figures/FigS9_LSU_cluster_occupancy_10kg.pdf', 'figures/FigS9_LSU_cluster_occupancy_domestic.pdf', 'figures/FigS9_LSU_cluster_occupancy.pdf', 'rdata/Ailaoshan_gis.rdata', 'modelprobs/Ailaoshan_model_output_bernoulli_LSU_probabilities.txt', 'modelprobs/Ailaoshan_model_output_bernoulli_SSU_probabilities.txt', 'modelslopes/Ailaoshan_model_output_full_LSU.txt', 'modelslopes/Ailaoshan_model_output_full_SSU.txt', 'pdfs/Ailaoshan_model_output_bernoulli_LSU.pdf', 'pdfs/Ailaoshan_model_output_bernoulli_SSU.pdf', 'pdfs/Ailaoshan_model_output_full_LSU.pdf', 'pdfs/Ailaoshan_model_output_full_SSU.pdf', 'pdfs/Ailaoshan_model_Rhat_LSU.pdf', 'pdfs/Ailaoshan_model_Rhat_SSU.pdf', 'preOTU_networks/preOTU_network_bipartite_weighted.pdf', 'preOTU_networks/preOTU_network_LSU_weighted.pdf', 'preOTU_networks/preOTU_network_SSU_weighted.pdf', 'rdata/Ailaoshan_dataexport.rdata', 'rdata/Ailaoshan_environmental.rdata', 'rdata/Ailaoshan_IUCNdata.rdata', 'rdata/Ailaoshan_model_data_bernoulli_LSU.rdata', 'rdata/Ailaoshan_model_data_bernoulli_SSU.rdata', 'rdata/Ailaoshan_model_data_bernoulli.rdata', 'rdata/Ailaoshan_model_data_final_LSU.rdata', 'rdata/Ailaoshan_model_data_final_SSU.rdata', 'rdata/Ailaoshan_model_data_final.rdata', 'rdata/Ailaoshan_model_data_full_LSU.rdata', 'rdata/Ailaoshan_model_data_full_SSU.rdata', 'rdata/Ailaoshan_model_data_full.rdata', 'rdata/Ailaoshan_occupancy_covariates_LSU_community_predictions.rdata', 'rdata/Ailaoshan_occupancy_covariates_LSU_individual_predictions.rdata', 'rdata/Ailaoshan_occupancy_covariates_SSU_community_predictions.rdata', 'rdata/Ailaoshan_occupancy_covariates_SSU_individual_predictions.rdata', 'rdata/Ailaoshan_OTU_table_preOTUs.rdata', 'rdata/Ailaoshan_OTU_table.rdata', 'rdata/Ailaoshan_preOTUs.rdata', 'rds/Ailaoshan_final_LSU_cluster_occupancy.rds', 'rds/Ailaoshan_final_LSU_jaccard_site.rds', 'rds/Ailaoshan_final_LSU_jaccard_species.rds', 'rds/Ailaoshan_final_SSU_cluster_occupancy.rds', 'rds/Ailaoshan_final_SSU_jaccard_site.rds', 'rds/Ailaoshan_final_SSU_jaccard_species.rds', 'rds/Ailaoshan_jaccard_sites_cut_LSU.rds', 'rds/Ailaoshan_jaccard_sites_cut_SSU.rds', 'rds/Ailaoshan_model_output_bernoulli_LSU.rds', 'rds/Ailaoshan_model_output_bernoulli_SSU.rds', 'rds/Ailaoshan_model_output_final_LSU.rds', 'rds/Ailaoshan_model_output_final_SSU.rds', 'rds/Ailaoshan_model_output_full_LSU.rds', 'rds/Ailaoshan_model_output_full_SSU.rds', 'rds/Ailaoshan_model_summary_final_LSU.rds', 'rds/Ailaoshan_model_summary_final_SSU.rds', 'tables/leeches_supplement_file_S2.xlsx', 'tables/leeches_supplement_file_S5.xlsx', 'tables/Table1_environmental_covariates_summary.csv', 'tables/Table3_IUCN_OTUs_latex.txt'

# process preOTU data
rule preOTUs:
    input:
        'data/OTUs_protax_outputs_12S_16S_weighted_unweighted_20190624.xlsx',
        'data/Ailaoshan_leech_2016_collections_20180406.xlsx',
        'data/2016_AilaoshanLeeches_info.xlsx'
    output:
        'rdata/Ailaoshan_preOTUs.rdata'
    log:
        'sessioninfo/Ailaoshan_preOTUs.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_preOTUs.R'

# draw networks to visualize correlations among preOTUs
rule preOTU_networks:
    input:
        'rdata/Ailaoshan_preOTUs.rdata'
    output:
        'preOTU_networks/preOTU_network_LSU_weighted.pdf',
        'preOTU_networks/preOTU_network_SSU_weighted.pdf',
        'preOTU_networks/preOTU_network_bipartite_weighted.pdf'
    log:
        'sessioninfo/Ailaoshan_preOTU_networks.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_preOTU_networks.R'

# process site environmental metadata
rule environmental:
    input:
        'data/environmental_variables_20180801.xlsx'
    output:
        'rdata/Ailaoshan_environmental.rdata'
    log:
        'sessioninfo/Ailaoshan_environmental.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_environmental.R'

# process final OTU data
rule OTU_table:
    input:
        'data/final taxonomy assignments/OTUs_protax_outputs_12S_16S_weighted_unweighted_20190624_20190830_20190924.xlsx',
        'data/Ailaoshan_leech_2016_collections_20180406.xlsx',
        'data/2016_AilaoshanLeeches_info.xlsx',
        'rdata/Ailaoshan_environmental.rdata',
        'data/pantheria/PanTHERIA_1-0_WR05_Aug2008.txt'
    output:
        'rdata/Ailaoshan_OTU_table.rdata',
        'rdata/Ailaoshan_OTU_table_preOTUs.rdata'
    log:
        'sessioninfo/Ailaoshan_OTU_table.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_OTU_table.R'

# get IUCN Redlist information
rule IUCNdata:
    input:
        'rdata/Ailaoshan_OTU_table.rdata'
    output:
        'rdata/Ailaoshan_IUCNdata.rdata'
    log:
        'sessioninfo/Ailaoshan_IUCNdata.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_IUCNdata.R {config[iucn_key]}'

# prepare data for drawing maps of Ailaoshan
# N.B. all shapefile files should be located in folder data/gis/ailaoshan_polygons
rule gis:
    input:
        'data/gis/ailaoshan_polygons/study_area.shp',
        'data/gis/srtm_data/n23_e100_1arc_v3.tif',
        'data/gis/srtm_data/n23_e101_1arc_v3.tif',
        'data/gis/srtm_data/n24_e100_1arc_v3.tif',
        'data/gis/srtm_data/n24_e101_1arc_v3.tif',
        'data/gis/srtm_data/n25_e100_1arc_v3.tif',
        'data/gis/srtm_data/n25_e101_1arc_v3.tif'
    output:
        'rdata/Ailaoshan_gis.rdata'
    log:
        'sessioninfo/Ailaoshan_gis.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_gis.R'

# draw Fig 1a
rule Fig1a_map_continent:
    input:
        'rdata/Ailaoshan_gis.rdata'
    output:
        'figures/Fig1a_Ailaoshan_map_continent.pdf'
    log:
        'sessioninfo/Ailaoshan_Fig1a_map_continent.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_Fig1a_map_continent.R'

# draw Fig 1b
rule Fig1b_map_3D:
    input:
        'rdata/Ailaoshan_gis.rdata'
    output:
        'figures/Fig1b_Ailaoshan_site_colour.png',
        'figures/Fig1b_Ailaoshan_site_grey.png'
    log:
        'sessioninfo/Ailaoshan_Fig1b_map_3D.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_Fig1b_map_3D.R'

# prepare Table 1 and draw Fig S1
rule Table1_FigS1_covariates:
    input:
        'rdata/Ailaoshan_environmental.rdata',
        'rdata/Ailaoshan_gis.rdata'
    output:
        'figures/FigS1a_elevation.pdf',
        'figures/FigS1c_TPI.pdf',
        'figures/FigS1e_road.pdf',
        'figures/FigS1g_stream.pdf',
        'figures/FigS1i_reserve.pdf',
        'figures/FigS1b_elevation.pdf',
        'figures/FigS1d_TPI.pdf',
        'figures/FigS1f_road.pdf',
        'figures/FigS1h_stream.pdf',
        'figures/FigS1j_reserve.pdf',
        'tables/Table1_environmental_covariates_summary.csv'
    log:
        'sessioninfo/Ailaoshan_Table1_FigS1_covariates.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_Table1_FigS1_covariates.R'

# examine descriptive statistics on dataset
# no file output: use interactively
rule descriptive:
    input:
        'rdata/Ailaoshan_OTU_table.rdata'
    log:
        'sessioninfo/Ailaoshan_descriptive.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_descriptive.R'

# examine missing data
# no file output: use interactively
rule missingdata:
    input:
        'rdata/Ailaoshan_OTU_table.rdata'
    log:
        'sessioninfo/Ailaoshan_missingdata.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_missingdata.R'

# draw Fig 2ab
rule Fig2ab_taxonomy_rarefaction:
    input:
        'rdata/Ailaoshan_OTU_table.rdata'
    output:
        'figures/Fig2a_taxonomic_breakdown.pdf',
        'figures/Fig2b_LabID_rarefaction.pdf'
    log:
        'sessioninfo/Ailaoshan_Fig2ab_taxonomy_rarefaction.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_Fig2ab_taxonomy_rarefaction.R'

# prepare data for occupancy modelling with all occupancy covariates
rule model_data_full:
    input:
        'rdata/Ailaoshan_OTU_table.rdata'
    output:
        'rdata/Ailaoshan_model_data_full_LSU.rdata',
        'rdata/Ailaoshan_model_data_full_SSU.rdata',
        'rdata/Ailaoshan_model_data_full.rdata'
    log:
        'sessioninfo/Ailaoshan_model_data_full.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_model_data_full.R'

# prepare data for occupancy modelling with Bayesian variable selection
rule model_data_bernoulli:
    input:
        'rdata/Ailaoshan_OTU_table.rdata'
    output:
        'rdata/Ailaoshan_model_data_bernoulli_LSU.rdata',
        'rdata/Ailaoshan_model_data_bernoulli_SSU.rdata',
        'rdata/Ailaoshan_model_data_bernoulli.rdata'
    log:
        'sessioninfo/Ailaoshan_model_data_bernoulli.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_model_data_bernoulli.R'

# prepare data for occupancy modelling with final variable selections
rule model_data_final:
    input:
        'rdata/Ailaoshan_OTU_table.rdata'
    output:
        'rdata/Ailaoshan_model_data_final_LSU.rdata',
        'rdata/Ailaoshan_model_data_final_SSU.rdata',
        'rdata/Ailaoshan_model_data_final.rdata'
    log:
        'sessioninfo/Ailaoshan_model_data_final.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_model_data_final.R'

# run occupancy model with all occupancy covariates
# alter model specification as desired by modifying dummy variables in input.jags
rule run_model_full:
    input:
        data = 'rdata/Ailaoshan_model_data_full_{dataset}.rdata',
        jags = 'jags/Ailaoshan_model_full.jags'
    output:
        rds = 'rds/Ailaoshan_model_output_full_{dataset}.rds',
        pdf = 'pdfs/Ailaoshan_model_output_full_{dataset}.pdf'
    log:
        'sessioninfo/Ailaoshan_model_parallel_full_{dataset}.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_model_parallel.R {input.data} {input.jags} {output.rds} {output.pdf}'

# summarize parameter posterior distributions from full model to help constrain priors in Bayesian variable selection
rule parameter_distributions:
    input:
        'rds/Ailaoshan_model_output_full_{dataset}.rds'
    output:
        'modelslopes/Ailaoshan_model_output_full_{dataset}.txt'
    log:
        'sessioninfo/Ailaoshan_parameter_distributions_{dataset}.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_parameter_distributions.R {input} {output}'

# run occupancy model with Bayesian variable selection
rule run_model_bernoulli:
    input:
        data = 'rdata/Ailaoshan_model_data_bernoulli_{dataset}.rdata',
        jags = 'jags/Ailaoshan_model_bernoulli_{dataset}.jags'
    output:
        rds = 'rds/Ailaoshan_model_output_bernoulli_{dataset}.rds',
        pdf = 'pdfs/Ailaoshan_model_output_bernoulli_{dataset}.pdf'
    log:
        'sessioninfo/Ailaoshan_model_parallel_bernoulli_{dataset}.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_model_parallel.R {input.data} {input.jags} {output.rds} {output.pdf}'

# get posterior probabilities for different model specifications
rule model_probabilities:
    input:
        'rds/Ailaoshan_model_output_bernoulli_{dataset}.rds'
    output:
        'modelprobs/Ailaoshan_model_output_bernoulli_{dataset}_probabilities.txt'
    log:
        'sessioninfo/Ailaoshan_model_probabilities_{dataset}.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_model_probabilities.R {input} {output}'

# run occupancy model with final variable selections
rule run_model_final:
    input:
        data = 'rdata/Ailaoshan_model_data_final_{dataset}.rdata',
        jags = 'jags/Ailaoshan_model_final_{dataset}.jags'
    output:
        'rds/Ailaoshan_model_output_final_{dataset}.rds'
    log:
        'sessioninfo/Ailaoshan_model_parallel_final_{dataset}.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_model_parallel_final.R {input.data} {input.jags} {output}'

# extract posterior summary information from final model
rule model_summary:
    input:
        rdata = 'rdata/Ailaoshan_model_data_final_{dataset}.rdata',
        rds = 'rds/Ailaoshan_model_output_final_{dataset}.rds'
    output:
        'rds/Ailaoshan_model_summary_final_{dataset}.rds'
    log:
        'sessioninfo/Ailaoshan_model_summary_{dataset}.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_model_summary.R {input.rdata} {input.rds} {output}'

# examine Rhat to help assess convergence
rule model_Rhat:
    input:
        'rds/Ailaoshan_model_summary_final_LSU.rds',
        'rds/Ailaoshan_model_summary_final_SSU.rds'
    output:
        'pdfs/Ailaoshan_model_Rhat_LSU.pdf',
        'pdfs/Ailaoshan_model_Rhat_SSU.pdf'
    log:
        'sessioninfo/Ailaoshan_model_Rhat.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_model_Rhat.R'

# compare species richness to community mean occupancy
# no file output: use interactively
rule species_richness:
    input:
        'rdata/Ailaoshan_model_data_final_LSU.rdata',
        'rdata/Ailaoshan_model_data_final.rdata',
        'rds/Ailaoshan_model_output_final_LSU.rds'
    log:
        'sessioninfo/Ailaoshan_species_richness.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_species_richness.R'

# draw Fig 2c
rule Fig2c_occupancy_detection:
    input:
        'rdata/Ailaoshan_OTU_table.rdata',
        'rds/Ailaoshan_model_summary_final_LSU.rds',
        'rds/Ailaoshan_model_summary_final_SSU.rds'
    output:
        'figures/Fig2c_occupancy_detection.pdf',
        'figures/Fig2c_occupancy_detection_labels.pdf'
    log:
        'sessioninfo/Ailaoshan_Fig2c_occupancy_detection.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_Fig2c_occupancy_detection.R'

# draw Fig 3 and Fig S3ab
rule Fig3_FigS3ab_occupancy:
    input:
        'rdata/Ailaoshan_OTU_table.rdata',
        'rdata/Ailaoshan_environmental.rdata',
        'rdata/Ailaoshan_gis.rdata',
        'rds/Ailaoshan_model_summary_final_LSU.rds',
        'rds/Ailaoshan_model_summary_final_SSU.rds'
    output:
        'figures/FigS3a_LSU_richhist.pdf',
        'figures/FigS3b_SSU_richhist.pdf',
        'figures/Fig3a_LSU_observed_richness.pdf',
        'figures/Fig3b_SSU_observed_richness.pdf',
        'figures/Fig3c_LSU_estimated_richness.pdf',
        'figures/Fig3d_SSU_estimated_richness.pdf',
        'figures/Fig3e_LSU_scatter.pdf',
        'figures/Fig3e_LSU_histogram.pdf',
        'figures/Fig3f_SSU_scatter.pdf',
        'figures/Fig3f_SSU_histogram.pdf'
    log:
        'sessioninfo/Ailaoshan_Fig3_FigS3ab_occupancy.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_Fig3_FigS3ab_occupancy.R'

# draw Fig S4
rule FigS4_richness:
    input:
        'rdata/Ailaoshan_OTU_table.rdata',
        'rds/Ailaoshan_model_summary_final_LSU.rds',
        'rds/Ailaoshan_model_summary_final_SSU.rds'
    output:
        'figures/FigS4a_obsvrich_LabID.pdf',
        'figures/FigS4b_richness_numleeches.pdf',
        'figures/FigS4c_obsvrich_PolygonID.pdf',
        'figures/FigS4d_observed_richness_sampling.pdf',
        'figures/FigS4e_estimated_richness_PolygonID.pdf',
        'figures/FigS4f_estimated_richness_sampling.pdf'
    log:
        'sessioninfo/Ailaoshan_FigS4_richness.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_FigS4_richness.R'

# compare results between datasets using coinertia etc
# no file output: use interactively
rule coinertia_etc:
    input:
        'rdata/Ailaoshan_OTU_table.rdata',
        'rds/Ailaoshan_model_summary_final_LSU.rds',
        'rds/Ailaoshan_model_summary_final_SSU.rds'
    log:
        'sessioninfo/Ailaoshan_coinertia_etc.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_coinertia_etc.R'

# prepare predictions for Fig 4ac and Fig S3c
rule Fig4ac_FigS3c_covariates_LSU_community:
    input:
        'rdata/Ailaoshan_model_data_final_LSU.rdata',
        'rdata/Ailaoshan_model_data_final.rdata',
        'rds/Ailaoshan_model_output_final_LSU.rds'
    output:
        'rdata/Ailaoshan_occupancy_covariates_LSU_community_predictions.rdata'
    log:
        'sessioninfo/Ailaoshan_Fig4ac_FigS3c_covariates_LSU_community.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_Fig4ac_FigS3c_covariates_LSU_community.R'

# prepare predictions for Fig 4bd and Fig S3d
rule Fig4bd_FigS3d_covariates_LSU_individual:
    input:
        'rdata/Ailaoshan_model_data_final_LSU.rdata',
        'rdata/Ailaoshan_model_data_final.rdata',
        'rds/Ailaoshan_model_output_final_LSU.rds'
    output:
        'rdata/Ailaoshan_occupancy_covariates_LSU_individual_predictions.rdata'
    log:
        'sessioninfo/Ailaoshan_Fig4bd_FigS3d_covariates_LSU_individual.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_Fig4bd_FigS3d_covariates_LSU_individual.R'

# prepare predictions for Fig 4e and Fig S3e
rule Fig4e_FigS3e_covariates_SSU_community:
    input:
        'rdata/Ailaoshan_model_data_final_SSU.rdata',
        'rdata/Ailaoshan_model_data_final.rdata',
        'rds/Ailaoshan_model_output_final_SSU.rds'
    output:
        'rdata/Ailaoshan_occupancy_covariates_SSU_community_predictions.rdata'
    log:
        'sessioninfo/Ailaoshan_Fig4e_FigS3e_covariates_SSU_community.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_Fig4e_FigS3e_covariates_SSU_community.R'

# prepare predictions for Fig 4f and Fig S3f
rule Fig4f_FigS3f_covariates_SSU_individual:
    input:
        'rdata/Ailaoshan_model_data_final_SSU.rdata',
        'rdata/Ailaoshan_model_data_final.rdata',
        'rds/Ailaoshan_model_output_final_SSU.rds'
    output:
        'rdata/Ailaoshan_occupancy_covariates_SSU_individual_predictions.rdata'
    log:
        'sessioninfo/Ailaoshan_Fig4f_FigS3f_covariates_SSU_individual.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_Fig4f_FigS3f_covariates_SSU_individual.R'

# draw Fig 4 and Fig S3cdef
rule Fig4_FigS3cdef_covariates:
    input:
        'rdata/Ailaoshan_occupancy_covariates_LSU_community_predictions.rdata',
        'rdata/Ailaoshan_occupancy_covariates_LSU_individual_predictions.rdata',
        'rdata/Ailaoshan_occupancy_covariates_SSU_community_predictions.rdata',
        'rdata/Ailaoshan_occupancy_covariates_SSU_individual_predictions.rdata'
    output:
        'figures/Fig4_occupancy_covariates.pdf',
        'figures/FigS3cdef_detection_covariates.pdf'
    log:
        'sessioninfo/Ailaoshan_Fig4_FigS3cdef_covariates.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_Fig4_FigS3cdef_covariates.R'

# draw Fig S5, Fig S6 and Fig S7
rule FigS5_FigS6_FigS7_occupancy_slopes:
    input:
        'rdata/Ailaoshan_OTU_table.rdata',
        'rdata/Ailaoshan_IUCNdata.rdata',
        'rds/Ailaoshan_model_summary_final_LSU.rds',
        'rds/Ailaoshan_model_summary_final_SSU.rds'
    output:
        'figures/FigS5_LSU_elev.pdf',
        'figures/FigS5_LSU_elev_10kg.pdf',
        'figures/FigS5_LSU_elev_domestic.pdf',
        'figures/FigS6_LSU_reserve.pdf',
        'figures/FigS6_LSU_reserve_10kg.pdf',
        'figures/FigS6_LSU_reserve_domestic.pdf',
        'figures/FigS7_SSU_elev.pdf',
        'figures/FigS7_SSU_elev_10kg.pdf',
        'figures/FigS7_SSU_elev_domestic.pdf'
    log:
        'sessioninfo/Ailaoshan_FigS5_FigS6_FigS7_occupancy_slopes.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_FigS5_FigS6_FigS7_occupancy_slopes.R'

# calculate posterior mean jaccard matrix
rule posterior_jaccard:
    input:
        rdata = 'rdata/Ailaoshan_model_data_final_{dataset}.rdata',
        rds = 'rds/Ailaoshan_model_output_final_{dataset}.rds'
    output:
        site = 'rds/Ailaoshan_final_{dataset}_jaccard_site.rds',
        species = 'rds/Ailaoshan_final_{dataset}_jaccard_species.rds'
    log:
        'sessioninfo/Ailaoshan_posterior_jaccard_{dataset}.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_posterior_jaccard.R {input.rdata} {input.rds} {output.site} {output.species}'

# draw Fig 5 and Fig S8
rule Fig5_FigS8_betadiv:
    input:
        'rdata/Ailaoshan_OTU_table.rdata',
        'rdata/Ailaoshan_environmental.rdata',
        'rdata/Ailaoshan_gis.rdata',
        'rds/Ailaoshan_final_LSU_jaccard_site.rds',
        'rds/Ailaoshan_final_SSU_jaccard_site.rds'
    output:
        'figures/FigS8a_tree_sites_LSU.pdf',
        'figures/FigS8b_tree_sites_SSU.pdf',
        'figures/Fig5a_nmds_sites_LSU.pdf',
        'figures/Fig5a_nmds_sites_LSU_labels.pdf',
        'figures/Fig5b_nmds_sites_SSU.pdf',
        'figures/Fig5b_nmds_sites_SSU_labels.pdf',
        'figures/Fig5c_map_sites_LSU.pdf',
        'figures/Fig5d_map_sites_SSU.pdf',
        'rds/Ailaoshan_jaccard_sites_cut_LSU.rds',
        'rds/Ailaoshan_jaccard_sites_cut_SSU.rds'
    log:
        'sessioninfo/Ailaoshan_Fig5_FigS8_betadiv.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_Fig5_FigS8_betadiv.R'

# calculate occupancy for each species in each site cluster
rule cluster_occupancy:
    input:
        rdata = 'rdata/Ailaoshan_model_data_final_{dataset}.rdata',
        rds = 'rds/Ailaoshan_model_output_final_{dataset}.rds',
        cluster = 'rds/Ailaoshan_jaccard_sites_cut_{dataset}.rds'
    output:
        'rds/Ailaoshan_final_{dataset}_cluster_occupancy.rds'
    log:
        'sessioninfo/Ailaoshan_cluster_occupancy_{dataset}.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_cluster_occupancy.R {input.rdata} {input.rds} {input.cluster} {output}'

# draw Fig 6, Fig S9 and Fig S10
rule Fig6_FigS9_FigS10_cluster_occupancy:
    input:
        'rdata/Ailaoshan_OTU_table.rdata',
        'rds/Ailaoshan_final_LSU_cluster_occupancy.rds',
        'rds/Ailaoshan_final_SSU_cluster_occupancy.rds',
        'rds/Ailaoshan_model_summary_final_LSU.rds',
        'rds/Ailaoshan_model_summary_final_SSU.rds'
    output:
        'figures/Fig6ab_cluster_occupancy.pdf',
        'figures/FigS9_LSU_cluster_occupancy.pdf',
        'figures/FigS9_LSU_cluster_occupancy_domestic.pdf',
        'figures/FigS9_LSU_cluster_occupancy_10kg.pdf',
        'figures/FigS10_SSU_cluster_occupancy.pdf',
        'figures/FigS10_SSU_cluster_occupancy_domestic.pdf',
        'figures/FigS10_SSU_cluster_occupancy_10kg.pdf'
    log:
        'sessioninfo/Ailaoshan_Fig6_FigS9_FigS10_cluster_occupancy.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_Fig6_FigS9_FigS10_cluster_occupancy.R'

# export tables for publication
rule exportdata:
    input:
        'rdata/Ailaoshan_OTU_table.rdata',
        'rdata/Ailaoshan_IUCNdata.rdata',
        'rds/Ailaoshan_model_summary_final_LSU.rds',
        'rds/Ailaoshan_model_summary_final_SSU.rds'
    output:
        'tables/Table3_IUCN_OTUs_latex.txt',
        'tables/leeches_supplement_file_S2.xlsx',
        'tables/leeches_supplement_file_S5.xlsx',
        'rdata/Ailaoshan_dataexport.rdata'
    log:
        'sessioninfo/Ailaoshan_exportdata.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_exportdata.R'

# export fasta files for publication
rule exportfasta:
    input:
        'rdata/Ailaoshan_OTU_table_preOTUs.rdata',
        'data/fasta/LSU_16S_otu_table_swarm_lulu_20190624.fas',
        'data/fasta/SSU_12S_otu_table_swarm_lulu_20190624.fas'
    output:
        'fasta/FileS3_LSU_representative_sequences.fasta',
        'fasta/FileS4_SSU_representative_sequences.fasta'
    log:
        'sessioninfo/Ailaoshan_exportfasta.sessioninfo.txt'
    shell:
        'Rscript code/Ailaoshan_exportfasta.R'
