# Overall and serotype-specific trends in dengue virus transmission in Mexico, 2016-2023

This repository includes .rds files which contain the data used to run the catalytic models used in the paper "Overall and serotype-specific trends in dengue virus transmission in Mexico, 2016-2023", an R script ("Code to run models.r") showing how to run the models for states in Mexico, and .STAN files of the four models used/discussed in the paper. 
These models are:

* Model 1, the non-serotype-specific time-constant model from Imai et al. 2016 [1]
* Model 2, an adapted version of the non-serotype-specific time-varying model from O'Driscoll et al. 2019 [2], modified to allow differential reporting rates for children under 15 years old and a time-varying historic FOI
* Model 3, as above but with a constant historic FOI
* Model 4, a serotype-specific time-varying model developed for this paper

The models are run using the "rstan" package, and functions from the "tidyverse" package are needed to perform some of the data manipulation required to format the case and population data to then run the models.

Provided within the "Figure Code" folder are the .R files that were used to generate the figures presented in the Results section of the paper, Fig1-7. Whilst the aggregated case data required to run the models has been provided, the raw case data used to generate this file, which is needed for the data analysis parts of the figures (Fig 1, Fig 2, and parts of Fig 3 and Fig 6) needs to be accessed online for the years 2020-2023 (available at: https://www.gob.mx/salud/documentos/datos-abiertos-bases-historicas-de-enfermedades-transmitidas-por-vector) and by request from INAI (Instituto Nacional de Acceso a la Información, folio 330026924000362) for years 2016-2019.

[1] Imai N, Dorigatti I, Cauchemez S, Ferguson NM. Estimating Dengue Transmission Intensity from Case-Notification Data from Multiple Countries. PLoS Negl Trop Dis. 2016;10(7):e0004833. Epub 20160711. doi: 10.1371/journal.pntd.0004833. PubMed PMID: 27399793; PubMed Central PMCID: PMCPMC4939939.
[2] O'Driscoll M, Imai N, Ferguson NM, Hadinegoro SR, Satari HI, Tam CC, et al. Spatiotemporal variability in dengue transmission intensity in Jakarta, Indonesia. PLoS Negl Trop Dis. 2020;14(3):e0008102. Epub 20200306. doi: 10.1371/journal.pntd.0008102. PubMed PMID: 32142516; PubMed Central PMCID: PMCPMC7080271.
