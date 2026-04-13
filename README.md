# π p → η(′) π p in the double-Regge region

https://inspirehep.net/literature/1859521

## Overview

This repository contains Julia code for modeling and fitting `π p → η(′) π p` in the double-Regge region, together with analysis scripts for fits, projections, asymmetries, toy studies, and plotting.

## Scripts

### Fit setup and fitting

- [publish_settings.jl](scripts/publish_settings.jl): creates fit output directories and writes a `settings.toml` file for a chosen model tag and exchange set.
- [fit_modelS.jl](scripts/fit_modelS.jl): loads fit settings, builds the selected model, performs the extended negative-log-likelihood fit, and writes `fit-results.toml`.
- [constrained_projections.jl](scripts/constrained_projections.jl): reads a saved fit result, projects the fitted model onto constrained partial waves, and pulls the solution through the mass bins.
- [constrained_projections_random_start.jl](scripts/constrained_projections_random_start.jl): repeats constrained partial-wave projections from randomized initial conditions to study local minima and projection stability.

### Model projections and diagnostics

- [project_model.jl](scripts/project_model.jl): evaluates a fitted model in mass bins, computes partial-wave projections, forward/backward intensities, and diagram-interference contribution matrices, and writes summary plots.
- [project_symmetric_model.jl](scripts/project_symmetric_model.jl): studies a simplified symmetric two-exchange model and compares its projected odd/even-wave content with the data.
- [alpha_prime_deps.jl](scripts/alpha_prime_deps.jl): scans the Regge slope parameter `α′` and plots how the `cosθ` dependence of the model changes.
- [etapi_mc.jl](scripts/etapi_mc.jl): generates phase-space Monte Carlo in the `ηπ` channel and compares several model-weighted kinematic distributions to data-driven weights.

### Asymmetry studies

- [phi_asymmetry_data.jl](scripts/phi_asymmetry_data.jl): computes and plots `ϕ`-asymmetry observables directly from the reconstructed data amplitudes.
- [phi_asymmetry_model.jl](scripts/phi_asymmetry_model.jl): computes and plots `ϕ`-asymmetry observables for individual exchange contributions in the model.

### Data preparation and auxiliary studies

- [shift_etappi.jl](scripts/shift_etappi.jl): copies the raw `η′π` input tables into the shifted dataset and applies the phase convention change for the odd-`M` phase files.
- [saving_bootstrap.jl](scripts/saving_bootstrap.jl): saves central `cosθ` distributions and bootstrap uncertainty bands from the data into `data/exp_pro/`, then plots the resulting intervals.
- [against_dimas_data.jl](scripts/against_dimas_data.jl): reads four-vector toy events from a ROOT file, reconstructs observables, and compares them with saved partial-wave projections.

## Pipeline

- [pipeline/fit.jl](pipeline/fit.jl): pipeline-oriented fitting entrypoint that reads a TOML config, fits the selected model, and writes `fit-results.toml`.
- [pipeline/plot_metrics.jl](pipeline/plot_metrics.jl): post-fit plotting entrypoint that reads fit results and produces model diagnostics such as intensities, asymmetries, and contribution plots.
- [pipeline/fit_settings.toml](pipeline/fit_settings.toml): example pipeline configuration describing the system, fit range, exchanges, and initial parameters.
- [pipeline/fit-results_pipeline.toml](pipeline/fit-results_pipeline.toml): saved example fit result for the pipeline workflow.

## Notes

- Most scripts assume the project is activated with DrWatson and that the expected input data already exist under `data/`.
- Several scripts are exploratory analysis jobs and can run for a long time, especially on a fresh Julia environment because of precompilation.
