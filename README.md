# π p → η(′) π p in the double-Regge region
 
https://inspirehep.net/literature/1859521

## Content

- [fit_modelS.jl](scripts/fit_modelS.jl): fit the model with extended nll
- [constrained_projections.jl](scripts/constrained_projections.jl): extended-nll PW analysis of a model

## Phi (s2) dependence

- [phi_asymmetry_data.jl](scripts/phi_asymmetry_data.jl): investigate phi asymmetry in data
- [phi_asymmetry_model.jl](scripts/phi_asymmetry_model.jl): investigate phi asymmetry in model

### Technical

- [publish_settings.jl](scripts/publish_settings.jl): create folders, save settings
- [saving_bootstrap.jl](scripts/saving_bootstrap.jl): bootstrap data, save results in `data/exp_pro/`
- [against_dimas_data.jl](scripts/against_dimas_data.jl): compare partial waves with toys from four-vector ROOT file


(to be fixed next)
- [project_model.jl](scripts/project_model.jl)
- [project_symmetric_model.jl](scripts/project_symmetric_model.jl)


## Scripts to be fixed for new API

### Plotting data, Exploration of the model

- [shift_etappi.jl](scripts/shift_etappi.jl)
- [alpha_prime_deps.jl](scripts/alpha_prime_deps.jl)

### Constrained Partial Waves


- [constrained_projections_random_start.jl](scripts/constrained_projections_random_start.jl)

### Modeling

- [etapi_mc.jl](scripts/etapi_mc.jl): generate phase space MC, plot with model weights 
- [bootstrap_reconstruction.jl](scripts/bootstrap_reconstruction.jl)
