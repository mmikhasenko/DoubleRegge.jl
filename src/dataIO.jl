
#
function reshape_compass_format(f_of_LM)
    tmp = [f_of_LM(L,M) for (L,M) in LMs]
    return [getindex.(tmp, i) for i in 1:length(tmp[1])]
end

function x_IδI_ϕδϕ_compass_ηπ(pathtodata)
    # 
    title_tobias(L,M) = "$(L)"*(isodd(L) ? "m" : "p")*"p"*(M==2 ? "M2" : "")
    # 
    get_intesity(L,M) = readdlm(joinpath(pathtodata, "EtaPi-$(title_tobias(L,M)).txt"))
    get_phase(L,M) = (L,M) == (2,1) ? fill(0.0,size(get_intesity(L,M))) :
        readdlm(joinpath(pathtodata, "EtaPi-Ph$(title_tobias(L,M)).txt"))
    # 
    xdata = get_intesity(1,1)[:,1];
    #
    Idata  = reshape_compass_format((L,M)->get_intesity(L,M)[:,2])
    δIdata = reshape_compass_format((L,M)->get_intesity(L,M)[:,3])
    ϕdata  = reshape_compass_format((L,M)->get_phase(L,M)[:,2])
    δϕdata = reshape_compass_format((L,M)->get_phase(L,M)[:,3])
    #
    xdata, Idata, δIdata, ϕdata, δϕdata
end
