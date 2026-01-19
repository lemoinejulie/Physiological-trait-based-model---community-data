# Physiological-trait-based-model---community-data
# Description

This physiological model is based on a trait-based approach and is constructed using equations derived from experimental data as described in Lemoine et al., 2025. The model is designed to go beyond the classical descriptive analysis of planktonic communities by providing a functional perspective.

Using the outputs of the model, you can estimate the respiration rate, grazing pressure, and fecal pellet production of the total community, or for each zooplankton PFT (Plankton Functional Type). The main aim is to provide a tool to better understand trophic interactions, such as bottom-up or top-down controls, mediated by predator-prey dynamics. At a broader scale, the model helps quantify the contribution of the mesozooplankton community to the global carbon budget.

The model was developed for quantitative imaging data, such as the Tara Oceans dataset, but it can also be applied to other datasets, provided certain data requirements are met (see below).

# Data Requirements

To use the model correctly, datasets must include:

A complete size spectrum, including both prey and predators at the same station.

Note: It is also possible to simulate prey using Chl-a and POC data (from in situ measurements or even satellite data). For details, see the GitHub repository Physiological-trait-based-model — timeseries_mortality
.

Temperature data corresponding to each station.

# Available Scripts

Two types of scripts are provided depending on the scale of analysis:

1. Single-Station Analysis

Allows the study of dynamics at a specific station.

station_script_diagnostic

station_function_diagnostic

station_plot_diagnostic

2. Multi-Station Analysis

Applies the model to multiple stations simultaneously for comparative studies.

global_script_diagnostic

global_function_diagnostic

global_plot_diagnostic

The outputs and plots in these scripts are non-exhaustive and can be modified, extended, or adapted for specific research questions or datasets.

# Contact

For more information, help with running the scripts, or advice on adapting the scripts to your dataset, please contact:
julie.lemoine@imev-mer.fr

Please cite Lemoine et al., [in prep] when using this model.

# Reference
