# π p → η(′) π p in the double-Regge region
 
https://inspirehep.net/literature/1859521


## Content

### Plotting data, Exploration of the model

- [against_dimas_data.jl](scripts/against_dimas_data.jl)
- [against_henris_data.jl](scripts/against_henris_data.jl)
- [shift_etappi.jl](scripts/shift_etappi.jl)

- [alpha_prime_deps.jl](scripts/alpha_prime_deps.jl)
- [phi_asymmetry_data.jl](scripts/phi_asymmetry_data.jl)
- [phi_asymmetry_model.jl](scripts/phi_asymmetry_model.jl)

### Constrained Partial Waves

- [constrained_projections.jl](scripts/constrained_projections.jl)
- [constrained_projections_random_start.jl](scripts/constrained_projections_random_start.jl)
- [project_model.jl](scripts/project_model.jl)
- [project_symmetric_model.jl](scripts/project_symmetric_model.jl)

### Modeling

- [fit_modelD.jl](scripts/fit_modelD.jl): 
- [fit_modelS.jl](scripts/fit_modelS.jl):
- [etapi_mc.jl](scripts/etapi_mc.jl): generate phase space MC, plot with model weights 

- [bootstrap_reconstruction.jl](scripts/bootstrap_reconstruction.jl)
- [saving_bootstrap.jl](scripts/saving_bootstrap.jl): bootstrap data, save results

### Technical

- [ELLH.jl](scripts/ELLH.jl): prototype likelihood computation with `Rank1Matrix` object
- [code_warntype_investigation.jl](scripts/code_warntype_investigation.jl): MWE of a typo problem
- [publish_settings.jl](scripts/publish_settings.jl): create folders, save settings
