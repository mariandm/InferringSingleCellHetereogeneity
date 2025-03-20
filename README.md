# *Code for:* Inferring single-cell heterogeneity of bacteriophage life-history traits from population-scale dynamics
*By Marian Dominguez-Mirazo, 2025*

## General Description

In this repository you can find all code needed to replicate the figures and analyses shown in the paper. All code is written in Matlab R2024_a.

## Folder content description

### ./Data

Contains all original experimental data.

- **TimeseriesVirus.csv :** Population-level timeseries of free infectious particles (PFU/ml). The first column represents sampling time, the next four columns are individual replicates. The data is used for population-level trait estimation in Fig 1. 

- **SingleCellData.csv :** Raw data from the single-cell lysis detection protocol. The first column is the experimental replicate. The second column presents the identification number for the plaque assay. The rest of the columns are the number of plaques detected at a specific plaque assay in the multiple sampling timepoints ranging from 4 to 12 hr, and annotated in the column name. The data is used for Figures 2 and 3. 

- **PlasticAdhesion.csv :** Data from the plastic adhesion experiment. The data are number of plaques detected in a plaque assay after allowing virions in a 96-well plate for x number of hours. The number of hours is described in the header. The first column represents immediate plaque assays without adding to well plate. Each individual row represents a replicate (6 replicates). 


### ./Scripts

Contains all code divided into three folders. Code should be ran from inside the folder containing the script.

#### **Functions :** 

Recurring functions used throughout the code. 

- **ODE_SEnIV_sink.m :** ODE function for numerical integration of the dynamical model of lytic infection. See methods for mathematical description. 

- **gamma_fromstats.m :** Calculate gamma shape and scale given a mean and coefficient of variation

- **lognormal_fromstats.m :** Calculate lognormal mu and sigma given a mean and coefficient of variation

- **lhsu.m :** Lating Hypercube Sampling from uniform distribution. Budiman (2003)

- **get_CDF.m :** Get CDF probabilities at specific timepoints for the specified distribution model and parameter set

- **Likelihood_CDF.m :** Calculate likelihood for single-cell experimental observations given a latent period distribution model

- **ML_CI.m :** Calculate CDF parameter confidence intervals based on likelihood

- **burstmean_perrep :** Calculate average burst size per experimental replicate across timepoints

- **burst_bestfit :** Randomly sample a parameter space to find the best patameter combination that describes the effective burst size data given a model of latent period to burst size relationship

- **burst_kannoly.m :** Model of latent period to burst size relationship. Kannoly et al. (2023), Microbiology Spectrum

- **burst_linear.m :** Linear model of latent period to burst size relationship

- **burst_mm.m :** Hill function model of latent period to burst size relationship

- **burst_RSME.m :** Calculate the root squared mean error (RSME) for single-cell experimental observations given a latent period to burst size model and a set of parameters

- **calculate_betaeff.m :** Calculate expected beta effective given a latent period to burst size model and corresponding parameters, and a latent period distribution. See methods for mathematical description

- **simulateSingleCellExp.m :** Simulate the single-cell lysis detection protocol

- **random_ads.m :** Simulate adsorption times in single-cell lysis detection protocol

- **pdf_ads.m :** Compute an PDF for adsorption times given initial cell and phage densities, adsorption rate, and co-incubation times

#### **IntermediateScripts :**

Scripts that generate intermediate files. Code may either take to long to create, so an intermediate file has been created, or are recurring computations. The output of the scripts in this folder are stored in ./IntermediateFiles and bare the same name as the script that produces them. 


- **SingleCellInfections.m :** Transform single-cell raw data into infected cell and lysed cell counts

- **visualCI_popLevel.m :** Calculate 95% confidence intervals for population-level timeseries fit using MCMC chains. Required for Fig1. 

- **MLE_CDF.m :** Find Maximum Likelihood Estimate parameters for 4 latent period CDF models using single-cell data. Required for Fig2. 

- **burstMean_perReplicate.m :** Calculate the cumulative effective burst size from raw single cell data. Required for Fig3. 

- **LS_Burst.m :** Find the parameters that minimize the squared error for a latent period to burst size function that best describe the single-cell cumulative effective burst size data. Required for Fig3. 

- **bootstrap_burstmodel.m :** Bootstrap single-cell cumulative effective burst size data and fit to calculate 95% CI. Required for Fig3. 

- **GMM_burst.m :** GMM clustering to predict individual bursts that ocurred at half hour intervals. Required for Fig3. 

- **burstDistribution_predicted.m :** Calculate the predicted burst size distribution based on the latent period distribution and latent period to burst size relationship. Required for Fig3. 

- **visualCI_burstfit.m :** Calculate 95% confidence intervals for cumulative effective burst size data. Required for Fig3. 

- **CDFAccuracy.m :** Simulates single cell experiments for a range of latent period distributions and use it to calculate MLE latent period parameters

- **incubationError.m :** Calculate the observed latent period mean and CV when including adsorption times. 

#### **FigureScripts :**

- **Figure1CD_PopLevel_Fit_PDF.m :** Produces Figure 1 panels C and D of manuscript. Population-level fit of viral time series (C), and predicted latent period distribution (D)

- **Figure2CD_SingleCell_CDF_PDF.m :** Produces Figure 2 panels C and D. Latent period CDF fit to single-cell experimental data (C), and predicted latent period distribution from population-level and single-cell data (D). 

- **Figure3_BurstSizeFunction.m :** Produces all panels of Figure 3. Cumulative effective burst size fit (A), predicted latent period to burst size relationship (B), and expected burst size distribution (C). 

- **SupplFig1_EtoCV.m :** Supplementary Figure 1: Relationship between number of E compartments and CV of latent period distribution in population-level dynamical model. 

- **SupplFig2_ChainDiagnostic.m :** Supplementary Figure 2: MCMC chain convergence diagnostics. 

- **SupplFig3_CDF_multimodel.m :** Supplementary Figure 3: Latent period CDF fit to single-cell experimental data to four CDF models. 

- **SupplFig4_ParameterEstimateComparison.m :** Supplementary Figure 4: Comparison of parameter estimated (Latent period mean, CV, burst size) using the population-level and single-cell data.

- **SupplFig5_ProportionInfected.m :** Supplementary Figure 5: Proportion of infected cells through sampling time points in single-cell lysis detection protocol. 

- **SupplFig6_CDFAccuracy.m :** Supplementary Figure 6: Latent period prediction accuracy when simulating the single-cell lysis detection protocol across a range of latent period distributions. 

- **SupplFig7_IncubationError.m :** Supplementary Figure 7: Effect of incubation times in observed latent period mean and CV. 

- **SupplFig8_adhesion.m :** Supplementary Figure 8: Time-dependent adhesion of viral particles to plastic surface of well.

- **SupplFig9_BurstSizeFitModels.m :**  Supplementary Figure 9: Cumulative effective burst size fit for three different latent period to burst size relationship models. 
 
### ./IntermediateFiles
Contains all intermediate files. The code that produces the intermediate files is found in ./Scripts/IntermediateScripts and has the same name as the file it produces. 

### ./Figures
Contains all figures for the manuscript and supplementary information. The code for all Figures is found in ./Scripts/FigureScripts and the script has the same name as the figure it produces. 