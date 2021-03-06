# Bayesian-Approach-to-EVT

## Abstract

In this work is developed an extreme value analysis for maximum daily wind gusts in a Bayesian framework.
The first part of the work focuses on simplistic inferential procedures, investigating both the two main approaches for extreme data selection and analysis, i.e. the BlockMaximaApproach and the PeaksOverThresholdApproach, only for specific sites. On the other hand, in the second part, a hierarchical random effects model is built in order to better incorporate in the analysis the structural complexity of the data. In particular, it attempts to identify site and seasonal effects for the marginal densities of daily maximum wind gusts, as well as the serial dependence at each location through a first order Markov chain model. 


## Models

* ***GEV*** **model**: 

<p align="center"><img src="https://github.com/MatteoPierdomenico/BayesianApproachToEVT/blob/main/Report/Immagini/GEV_model.png" width="65%" height=65%">

* ***GPD*** **model**: 
<p align="center"><img src="https://github.com/MatteoPierdomenico/BayesianApproachToEVT/blob/main/Report/Immagini/GPD_model.png" width="60%" height="60%">

* ***Hierarchical random effects*** **model**: 
<p align="center"><img src="https://github.com/MatteoPierdomenico/BayesianApproachToEVT/blob/main/Report/Immagini/Likelihood.png" width="50%" height="50%">

<p align="center"><img src="https://github.com/MatteoPierdomenico/BayesianApproachToEVT/blob/main/Report/Immagini/Hmodel.png" width="50%" height="50%">

To better understand the mathematical analysis inherent to those models, see the *FabbrucciPierdomenicoRandazzo_report.pdf* file in the *Report* folder.


## File system

The data analyzed are collected in the **Ireland_daily_from1990.csv** file, in the folder *Dataset*, together with a little description of the source in the **Dataset.txt** file.

Codes are organized in the main folder in the following way: 

* The *BlockMaximaApproach* folder containing: 
  * **mcmc_gev_stan.R** file: it is a simple one site analysis with non-informative priors adopting the GEV model written in the **Gev.stan** file,
  * **gev_with_prior_elicitation.R** file: it considers the extrapolation, directly from a portion of data, and use of "informative" priors for the GEV model, it calls the **GEV2.stan** file;

* The *PeaksOverThresholdApproach* folder: 
  * **mcmc_gpd_jags.R** file: it corresponds to the Extreme Value analysis of our data through a Generalized Pareto distribution model, in which also the threshold is an unknown parameter of the model. The file calls the **GPD.stan** file present in the folder;

* The *HierarchicalRandomEffectsModel* folder, divided in three subfolders: 
  * the main one is the *Random effects model for sites and months* subfolder, where there is the complete analysis of the datasat written both in STAN and NIMBLE, using a GPD     to model extremes and considering a first order Markov chain to model the serial dependence in the daily wind gusts observations within every month, namely:
    * **FullHierarchicalModel_Nimble.R**: the NIMBLE version of our Hierarchical model with site and month random effects,
    * **Hierarchical_Model_Full_stan.R**: the STAN version of our Hierarchical model with site and month random effects which calls the  
      **Hierarchical_Model_Full.stan**  file;  
  
  * the other 2 subfolders correspond to preliminary steps in order to reach the more complex and complete hierarchical model for sites and months effects. In particular:
    * in the *Model the time dependence for one site* subfolder, it is implemented a preliminary simple GPD model for single sites observations, as the one in the *PeaksOverThresholdApproach* folder, with the serial dependence modelled by the first order Markov chain; 
    * in the *Random effects model for site* folder, there are both the STAN and NIMBLE implementation of a slightly more complex model, the last preliminary step, which analyzes the whole dataset, considering just the random effect for the site variation.
    
    
    
## Some results:
* Parameters' trace plots of the **GEV** model for the counties *Kerry* and *Westmeath* respectively:
<img src="https://github.com/MatteoPierdomenico/BayesianApproachToEVT/blob/main/Report/Immagini/traces.png" width="70%" height="70%">

* Return levels for *Kerry* and *Westmeath* from the **GPD** model analysis:
<img src="https://github.com/MatteoPierdomenico/BayesianApproachToEVT/blob/main/Report/Immagini/Ret_lev_GPD2.png" width="50%" height="50%">

* Parameters trace plots of the **hierarchical random effects** model for the counties *Kerry* and *Westmeath* in *July* and *December*:
<img src="https://github.com/MatteoPierdomenico/BayesianApproachToEVT/blob/main/Report/Immagini/Trace_Plots_HM.png" width="70%" height="70%">

* Return levels for *Kerry* and *Westmeath* in *December* and *July* from the **hierarchical random effects** model:
<img src="https://github.com/MatteoPierdomenico/BayesianApproachToEVT/blob/main/Report/Immagini/Return_levels_HM.png" width="50%" height="50%">

## Acknowledgments

Bayesian statistics project written by Davide Fabbrucci, Matteo Pierdomenico and Giacomo Randazzo for the Bayesian statistics 2020/21 course, professor Alessandra Guglielmi, Politecnico di Milano.
