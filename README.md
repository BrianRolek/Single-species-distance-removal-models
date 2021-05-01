# Single-species-distance-removal-models
This repository contains files and data for analyses within Rolek, B.W., D.J. Harrison, D.W. Linden, C.S. Loftin, P.B. Wood. 2021. Associations between breeding conifer-associated birds, forestry treatments, years-since-harvest, and vegetation in regenerating stands.

R scripts include abundance models with detection probabilities from distance and time-removal sampling for avian point counts. We use R as an interface to implement models in Just Another Gibb's Sampler (JAGS). 

Simplified models (both Poisson and zero-inflated Poisson) that lack covariates are provided in the file 01-analysis_basic-models-pois-ZIP.R

We provide two Rmarkdown documents for reproduceable workflows. One Rmarkdown document provides code for analyses of abundance in response to vegetaion characteristics. The other Rmarkdown document provides code for analyses of abundance in response to forestry treatments.


