# π p → η(′) π p in the double-Regge region

https://inspirehep.net/literature/1859521

## Overview

This repository contains Julia code for modeling and fitting `π p → η(′) π p` in the double-Regge region, together with analysis scripts for fits, projections, asymmetries, toy studies, and plotting.

## Scripts

Cost guide:

- `Light`: quick utility or small diagnostic; suitable for routine sanity checks.
- `Moderate`: noticeable compute or plotting work, but usually still interactive.
- `Heavy`: fit, large integration/projection sweep, bootstrap, or large Monte Carlo job.
- `Very heavy`: long-running exploratory job; expect batch-style runtime.

### Fit setup and fitting

- [publish_settings.jl](scripts/publish_settings.jl): `Light`. Creates the fit output directories for the hard-coded tag in the script and writes the corresponding `settings.toml`.
- [fit_modelS.jl](scripts/fit_modelS.jl): `Heavy`. Loads a saved fit configuration, builds the selected model, runs the extended negative-log-likelihood fit, and writes `fit-results.toml`.
- [constrained_projections.jl](scripts/constrained_projections.jl): `Heavy`. Reads a saved fit result, projects the fitted model onto constrained partial waves, tracks ambiguities, and writes projection summaries.
- [constrained_projections_random_start.jl](scripts/constrained_projections_random_start.jl): `Very heavy`. Repeats constrained partial-wave projections from randomized initial conditions to study local minima; intended as a long exploratory batch job.

### Model projections and diagnostics

- [project_model.jl](scripts/project_model.jl): `Heavy`. Evaluates a saved fitted model in mass bins, computes partial-wave projections, forward/backward intensities, interference contributions, and several diagnostic plots; also expects `pdftk` for the final PDF merge.
- [project_symmetric_model.jl](scripts/project_symmetric_model.jl): `Moderate`. Studies a simplified two-exchange symmetric model and compares its projected odd/even-wave content with the `ηπ` data.
- [alpha_prime_deps.jl](scripts/alpha_prime_deps.jl): `Moderate`. Scans a shared Regge slope parameter `α′` and plots how the `cosθ` dependence of the model changes.
- [etapi_mc.jl](scripts/etapi_mc.jl): `Heavy`. Generates `10^6` phase-space Monte Carlo points in the `ηπ` channel and compares several model-weighted kinematic distributions to data-driven weights.

### Asymmetry studies

- [phi_asymmetry_data.jl](scripts/phi_asymmetry_data.jl): `Light`. Computes and plots `ϕ`-asymmetry observables directly from the reconstructed data amplitudes.
- [phi_asymmetry_model.jl](scripts/phi_asymmetry_model.jl): `Light`. Computes and plots `ϕ`-asymmetry observables for individual exchange contributions in a fixed `ηπ` model setup.

### Data preparation and auxiliary studies

- [shift_etappi.jl](scripts/shift_etappi.jl): `Light`. Copies the raw `η′π` input tables into `data/exp_raw/PLB_shifted` and applies the phase-convention shift for the odd-`M` phase files.
- [saving_bootstrap.jl](scripts/saving_bootstrap.jl): `Heavy`. Saves central `cosθ` distributions and bootstrap uncertainty bands from the data into `data/exp_pro/`, then plots the resulting intervals.
- [export_default_model_crosscheck.jl](scripts/export_default_model_crosscheck.jl): `Light`. Re-evaluates the preserved default-model fit results at a few fixed phase-space points and rewrites `test/fixtures/default_model_crosscheck.json`.
- [against_dimas_data.jl](scripts/against_dimas_data.jl): `Moderate`. Reads four-vector toy events from a ROOT file, reconstructs observables, and compares them with saved projection inputs such as `data/exp_pro/main_point.txt`.

## Pipeline

- [pipeline/fit.jl](pipeline/fit.jl): `Heavy`. Pipeline-oriented fitting entrypoint that reads a TOML config, fits the selected model, and writes `fit-results.toml`.
- [pipeline/plot_metrics.jl](pipeline/plot_metrics.jl): `Heavy`. Post-fit plotting entrypoint that reads fit results and produces model diagnostics such as intensities, asymmetries, contribution plots, and combined PDFs.
- [pipeline/fit_settings.toml](pipeline/fit_settings.toml): example pipeline configuration describing the system, fit range, exchanges, and initial parameters.
- [pipeline/fit-results_pipeline.toml](pipeline/fit-results_pipeline.toml): saved example fit result for the pipeline workflow.

## Notes

- Most scripts assume the project is activated with DrWatson and that the expected input data already exist under `data/`.
- The scripts marked `Heavy` or `Very heavy` are not intended as quick smoke tests; they are better treated as batch jobs, especially on a fresh Julia environment with precompilation overhead.
