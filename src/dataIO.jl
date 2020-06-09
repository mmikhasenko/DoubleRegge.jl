
title_tobias(L,M) = "$(L)"*(isodd(L) ? "m" : "p")*"p"*(M==2 ? "M2" : "")
get_intesity(L,M; pathtodata=joinpath("data", "exp_raw", "PLB_shifted")) = readdlm(joinpath(pathtodata, "EtaPi-$(title_tobias(L,M)).txt"))
get_phase(L,M; pathtodata=joinpath("data", "exp_raw", "PLB_shifted")) = (L,M) == (2,1) ? fill(0.0,size(get_intesity(L,M))) :
    readdlm(joinpath(pathtodata, "EtaPi-Ph$(title_tobias(L,M)).txt"))
