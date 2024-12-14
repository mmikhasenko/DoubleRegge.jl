# π p → η(′) π p in the double-Regge region
 
https://inspirehep.net/literature/1859521

## Content

- [saving_bootstrap.jl](scripts/saving_bootstrap.jl): bootstrap data, save results in `data/exp_pro/`
- [against_dimas_data.jl](scripts/against_dimas_data.jl): open Dima's ROOT file, compare with our data
- [fit_modelS.jl](scripts/fit_modelS.jl): fit the model with extended nll
- [constrained_projections.jl](scripts/constrained_projections.jl): extended-nll PW analysis of a model

(to be fixed next)
- [project_model.jl](scripts/project_model.jl)
- [project_symmetric_model.jl](scripts/project_symmetric_model.jl)

### Plotting data, Exploration of the model

- [against_henris_data.jl](scripts/against_henris_data.jl): open Dima's and Henri's ROOT file, compare with our data
- [shift_etappi.jl](scripts/shift_etappi.jl)

- [alpha_prime_deps.jl](scripts/alpha_prime_deps.jl)
- [phi_asymmetry_data.jl](scripts/phi_asymmetry_data.jl)
- [phi_asymmetry_model.jl](scripts/phi_asymmetry_model.jl)

### Constrained Partial Waves

- [constrained_projections_random_start.jl](scripts/constrained_projections_random_start.jl)

### Modeling

- [etapi_mc.jl](scripts/etapi_mc.jl): generate phase space MC, plot with model weights 
- [bootstrap_reconstruction.jl](scripts/bootstrap_reconstruction.jl)

### Technical

- [ELLH.jl](scripts/ELLH.jl): prototype likelihood computation with `Rank1Matrix` object
- [code_warntype_investigation.jl](scripts/code_warntype_investigation.jl): MWE of a typo problem
- [publish_settings.jl](scripts/publish_settings.jl): create folders, save settings
