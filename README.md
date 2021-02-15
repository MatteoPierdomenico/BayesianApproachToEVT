# Bayesian-Approach-to-EVT
Bayesian statistics project written by Davide Fabbrucci, Matteo Pierdomenico and Giacomo Randazzo; Bayesian statistics 2020/21 course, Politecnico di Milano.
Our work aims to study Extreme daily wind gust in 5 different counties of Ireland adopting a Bayesian approach. In this repository you can find the codes we implemented in order to accomplish our goal.

## Abstract

In this work is developed an extreme value analysis for maximum daily wind gusts in a Bayesian framework.
The first part of the work focuses on simplistic inferential procedures, investigating both the two main approaches for extreme data selection and analysis, i.e. the BlockMaximaApproach and the PeaksOverThresholdApproach, only for specific sites. On the other hand, in the second part, a hierarchical random effects model is built in order to better incorporate in the analysis the structural complexity of the data. In particular, it attempts to identify site and seasonal effects for the marginal densities of daily maximum wind gusts, as well as the serial dependence at each location through a first order Markov chain model. 

### Data
![alt text](https://github.com/MatteoPierdomenico/BayesianApproachToEVT/blob/main/Report/Immagini/Data.png)

## Models

* ***GEV*** **model**: 

![alt text](https://github.com/MatteoPierdomenico/BayesianApproachToEVT/tree/main/Report/Immagini/GEV_model.png)

* ***GPD*** **model**: 

![alt text](https://github.com/MatteoPierdomenico/BayesianApproachToEVT/tree/main/Report/Immagini/GPD_model.png)

* ***Hierarchical random effects*** **model**: 

![alt text](https://github.com/MatteoPierdomenico/BayesianApproachToEVT/tree/main/Report/Immagini/Likelihood.png)

![alt text](https://github.com/MatteoPierdomenico/BayesianApproachToEVT/tree/main/Report/Immagini/Hmodel.png)



## File system

The data analyzed are collected in the **Ireland_daily_from1990.csv** file, in the folder *Dataset*, together with a little description of the source in the **Dataset.txt** file.

In the *Report* folder is contained the *FabbrucciPierdomenicoRandazzo_report.pdf*, containing all the analysis of our work.

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
