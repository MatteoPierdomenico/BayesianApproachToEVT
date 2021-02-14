# Bayesian-Approach-to-EVT
Bayesian statistics project written by Davide Fabbrucci, Matteo Pierdomenico and Giacomo Randazzo; Bayesian statistics 2020/21 course, Politecnico di Milano.
Our work aims to study Extreme daily wind gust in 5 different counties of Ireland adopting a Bayesian approach. In this repository you can find the codes we implemented in order to accomplish our goal.

## Abstract

In this work is developed an extreme value analysis for maximum daily wind gusts in a Bayesian framework.
The first part of the work focuses on simplistic inferential procedures, investigating both the two main approaches for extreme data selection and analysis, i.e. the BlockMaximaApproach and the PeaksOverThresholdApproach, only for specific sites. On the other hand, in the second part, a hierarchical random effects model is built in order to better incorporate in the analysis the structural complexity of the data. In particular, it attempts to identify site and seasonal effects for the marginal densities of daily maximum wind gusts, as well as the serial dependence at each location through a first order Markov chain model. 


## File system

The data analyzed are collected in a .csv file, in the folder "Dataset", together with a little description of the source in the .txt file.
Codes are organized in the main folder in the following way: the BlockMaximaApproach folder, where there 
are two .R files that refers to the two .stan files, one of them is a simple analysis with non-informative priors adopting the GEV model, the other considers the extrapolation, directly from a portion of data, and use of "informative" priors for the GEV model;
in the PeaksOverThresholdApproach folder there is a .R file together with a .jags file that correspond to the Extreme Value Theory analysis of our data through a Generalized Pareto distribution model, in which also the threshold is an unknown parameter of the model;
in the HierarchicalRandomEffectsModel folder there are three subfolders: the main one is the "Random effects model for sites and months" subfolder, where there is the complete analysis of the datasat written both in STAN and NIMBLE, using a GPD to model extremes and considering a first order Markov chain to model the serial dependence in the daily wind gusts observations within every month; the other 2 subfolders correspond to preliminary steps in order to reach the more complex and complete hierarchical model for sites and months effects. In particular, in the "Model the time dependence for one site" subfolder, it is implemented a preliminary simple GPD model for single sites observations, as the one in the PeaksOverThresholdApproach folder, with the serial dependence modelled by the first order Markov chain; while in the "Random effects model for site" folder, there are both the STAN and NIMBLE implementation of a slightly more complex model, the last preliminary step, which analyzes the whole dataset, considering just the random effect for the site variation.
