

compass_ηπ_LMs = [
    (L = 1, M = 1), (L = 2, M = 1), (L = 2, M = 2),
    (L = 3, M = 1), (L = 4, M = 1),
    (L = 5, M = 1), (L = 6, M = 1)];

#
function reshape_compass_format(f_of_LM, LMs)
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
    Idata  = reshape_compass_format((L,M)->get_intesity(L,M)[:,2], compass_ηπ_LMs)
    δIdata = reshape_compass_format((L,M)->get_intesity(L,M)[:,3], compass_ηπ_LMs)
    ϕdata  = reshape_compass_format((L,M)->get_phase(L,M)[:,2], compass_ηπ_LMs)
    δϕdata = reshape_compass_format((L,M)->get_phase(L,M)[:,3], compass_ηπ_LMs)
    #
    NamedTuple{(:x, :I, :δI, :ϕ, :δϕ)}.(zip(xdata,Idata,δIdata,ϕdata,δϕdata))
end
