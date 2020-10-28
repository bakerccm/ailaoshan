# Analysis code for Ailaoshan study

## Citation

Yinqiu JI\*, Christopher CM BAKER\*,\*\*, Viorel D POPESCU, Jiaxin WANG, Chunying WU, Zhengyang WANG, Yuanheng LI, Lin WANG, Chaolang HUA, Zhongxing YANG, Chunyan YANG, Charles CY XU, Qingzhong WEN, Naomi E PIERCE\*\* and Douglas W YU\*\*. Measuring protected-area outcomes with leech iDNA: large-scale quantification of vertebrate biodiversity in Ailaoshan reserve.

\* - equal contributions<br/>
\*\* - corresponding authors

A preprint of this study is available from bioRxiv at https://www.biorxiv.org/content/10.1101/2020.02.10.941336v1

## Workflow and repository structure

The workflow for this study is documented in the `snakefile` located in the root of the repository, and illustrated in the `snakemake_*.pdf` files.

Note that some steps in this workflow are computationally intensive and/or have high memory requirements, especially model estimation and calculations on the posterior sample.

 - `/data` contains the input data for this study. Note that `/data/gis` contains digital elevation model files downloaded from USGS Earth Explorer, and `/data/pantheria` contains data on mammals downloaded from the PanTHERIA database (see readme in that folder for details).

 - `/code` contains the R code for our analysis. `/jags` contains the JAGS code for the multi-species occupancy models.

 - `/config` contains `config.yaml` which may be used to store user-specific parameters such as the API key used by `/code/Ailaoshan_IUCNdata.R` to access IUCN data (the API key is not included in this repository, and you should supply your own if you want to use this code).

 - Intermediate and final output files get saved to `/preOTU_networks`, `/figures`, `/tables`, `/fasta`, `/modelprobs`, `/modelslopes`, `/pdfs`, `/rdata` and `/rds`. These files are not included in this repository but can all be re-generated using the input data and code provided here.

 - `/sessioninfo` contains information about the R packages loaded when running the code in this workflow. This may be particularly helpful if you need to identify the specific version of any package that we used.

## Further information

- Please contact (Chris Baker)[https://github.com/bakerccm] for any enquiries regarding this repository.
