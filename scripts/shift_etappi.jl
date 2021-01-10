using DelimitedFiles
using DoubleRegge

const path_from = joinpath("data", "exp_raw", "PLB")
const path_to = joinpath("data", "exp_raw", "PLB_shifted")

function to180(ϕ)
    ϕ >  180 && return ϕ - 360*div(ϕ+180, 360)
    ϕ < -180 && return ϕ - 360*div(ϕ-180, 360)
    return ϕ
end

const regex_phases = r"EtapPi-Ph([1-6]).*\.txt"


function open_shift_save(src, dst)
    d = readdlm(src)
    d[:,2] .= map(ϕ -> iszero(ϕ) ? 0.0 : ϕ+180, d[:,2])
    writedlm(dst, d)
end

for f in readdir(path_from)
    in_f, out_f = joinpath(path_from, f), joinpath(path_to, f)
    rm = match(regex_phases, f)
    if rm === nothing || iseven(Meta.parse(rm[1]))
        cp(in_f, out_f)
        continue
    else
        isodd(Meta.parse(rm[1])) &&
            open_shift_save(in_f, out_f)
    end
end
