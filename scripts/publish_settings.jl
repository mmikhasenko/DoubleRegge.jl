using DrWatson
@quickactivate "DoubleRegge"
# 
using TOML
using DoubleRegge

# # # # # # # # # # # # # # # # # # # # 
# 
tag = "etappi_a2Po-a2f2-f2f2-PoPo_opposite-sign"
# 
# # # # # # # # # # # # # # # # # # # # 

settings = Dict(
    "system" => "compass_η′π",
    "pathtodata" => joinpath("data", "exp_raw", "PLB_shifted"),
    "fitrange" => [2.4, 3.0],
    "t2" => -0.2,
    "exchanges" => [1, 2, 4, 5],
    "initial_pars" => [0.2, 1.2, -9.0, -0.5],
    "scale_α" => 0.8,
)

for folder in ([plotsfolder(), fitsfolder()])
    if !(tag ∈ readdir(folder))
        println("--> Create ", joinpath(folder, tag))
        mkdir(joinpath(folder, tag))
    end
    println("--> Finished")
end

output_name = fitsfolder(tag, "settings.toml")
open(output_name, "w") do io
    TOML.print(io, settings)
end


# bottom-Po model same sign
# settings = Dict(
#     "system" => "compass_ηπ",
#     "pathtodata" => joinpath("data","exp_raw","PLB_shifted"),
#     "fitrange" => [2.4, 3.0],
#     "t2" => -0.2,
#     "tag" => "a2Po-f2Po-PoPo_oppoite-sign",
#     "exchanges" => [1,3,5],
#     "initial_pars" => [0.7, 0.7, 0.0 ],
#     "s2_shift" => 30.0,
#     "scale_α" => 0.8,
# )

# # bottom-Po model opposite sign
# settings = Dict(
#     "system" => "compass_ηπ",
#     "pathtodata" => joinpath("data","exp_raw","PLB_shifted"),
#     "fitrange" => [2.4, 3.0],
#     "t2" => -0.2,
#     "tag" => "a2Po-f2Po-PoPo_opposite-sign",
#     "exchanges" => [1,3,5],
#     "initial_pars" => [0.7, -0.7, 0.0 ],
#     "s2_shift" => 30.0,
#     "scale_α" => 0.8,
# )

# cesar model with f2/f2: same sign
# settings = Dict(
#     "system" => "compass_ηπ",
#     "pathtodata" => joinpath("data","exp_raw","PLB_shifted"),
#     "fitrange" => [2.4, 3.0],
#     "t2" => -0.2,
#     "tag" => "a2Po-f2f2-PoPo_same-sign",
#     "exchanges" => [1,4,5],
#     "initial_pars" => [0.7, 11.0, 0.0],
#     "s2_shift" => 30.0,
#     "scale_α" => 0.8,
# )

# # cesar model with f2/f2: opposite sign
# settings = Dict(
#     "system" => "compass_ηπ",
#     "pathtodata" => joinpath("data","exp_raw","PLB_shifted"),
#     "fitrange" => [2.4, 3.0],
#     "t2" => -0.2,
#     "tag" => "a2Po-f2f2-PoPo_opposite-sign",
#     "exchanges" => [1,4,5],
#     "initial_pars" => [0.7, -11.0, 0.0],
#     "s2_shift" => 30.0,
#     "scale_α" => 0.8,
# )

# # bottom-Po model opposite sign
# settings = Dict(
#     "system" => "compass_ηπ",
#     "pathtodata" => joinpath("data","exp_raw","PLB_shifted"),
#     "fitrange" => [2.4, 3.0],
#     "t2" => -0.2,
#     "tag" => "a2Po-f2Po-PoPo_opposite-sign_s2shift",
#     "exchanges" => [1,3,5],
#     "initial_pars" => [0.7, -0.7, 0.0 ],
#     "s2_shift" => 30.0,
#     "scale_α" => 0.8,
# )
# four-parameters fit
# settings = Dict(
#     "system" => "compass_ηπ",
#     "pathtodata" => joinpath("data","exp_raw","PLB_shifted"),
#     "fitrange" => [2.4, 3.0],
#     "t2" => -0.2,
#     "tag" => "a2Po-f2Po-a2f2-f2f2_opposite-sign",
#     "exchanges" => [1,2,3,4],
#     "initial_pars" => [0.7, 0.2, -0.7, -0.2],
#     "s2_shift" => 0.0,
#     "scale_α" => 0.8,
# )
# five-parameter fit
# settings = Dict(
#     "system" => "compass_ηπ",
#     "pathtodata" => joinpath("data","exp_raw","PLB_shifted"),
#     "fitrange" => [2.4, 3.0],
#     "t2" => -0.2,
#     "tag" => "a2Po-f2Po-a2f2-f2f2-PoPo_opposite-sign",
#     "exchanges" => [1,2,3,4,5],
#     "initial_pars" => [0.3, 3.7, -0.2, -9.0, 0.0],
#     "s2_shift" => 0.0,
#     "scale_α" => 0.8,
# )
# 
# settings = Dict(
#     "system" => "compass_ηπ",
#     "pathtodata" => joinpath("data","exp_raw","PLB_shifted"),
#     "fitrange" => [2.4, 3.0],
#     "t2" => -0.2,
#     "tag" => "a2Po-f2Po-a2f2-f2f2-Pof2_opposite-sign",
#     "exchanges" => [1,2,3,4,6],
#     "initial_pars" => [0.3, 3.7, -0.2, -9.0, 0.0],
#     "s2_shift" => 0.0,
#     "scale_α" => 0.8,
# )
# #
