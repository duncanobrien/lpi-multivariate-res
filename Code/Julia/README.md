## Recommended code workflow

[*ODE_preamble.jl*](ODE_preamble.jl) - preamble script listing the necessary files and loading supporting scripts.

[*simulate_communities.jl*](simulate_communities.jl) - simulate multispecies communities using the interaction matrices and starting abundances generated by [*save_communities.R*](https://github.com/duncanobrien/lpi-multivariate-res/tree/main/Code/R/save_communities.R).

## Supporting functions
[*generalised_lotka_volterra_invasive.jl*](generalised_lotka_volterra_invasive.jl) - collection of supporting functions to fit threshold generalised additive models.

[*generalised_lotka_volterra_stochastic.jl*](generalised_lotka_volterra_stochastic.jl) - a modfication of the `ewsnet_predict()` function from the [EWSmethods R package](https://www.authorea.com/doi/full/10.22541/au.166801190.00303336) to use the Impulse version of the EWSNet machine learning model. These weights are provided in this [repository](https://github.com/duncanobrien/ews-assessments/tree/main/python/weights/Pretrained).

[*extract_ews_pred_fn.R*](extract_ews_pred_fn.R) - function to extract the prediction made by each EWS indicator for each EWS computation method.