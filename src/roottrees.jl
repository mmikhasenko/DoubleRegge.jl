#
function list_of_vectors(tree_entry;
        varnames=error("list of four-vector names"),
        branchnames=error("list of branch names"),
        components=['X', 'Y', 'Z', 'E'],
        between="",
        before="")
    tpl = NamedTuple{varnames}(
        [1e3 .* SVector([getproperty(tree_entry,Symbol(before*p*between*c)) for c in ['X', 'Y', 'Z', 'E']]...)
            for p in branchnames])
    return tpl
end
#
function broadcast_over_tree(f::Function;
        filename = error("give the name of the root file"),
        treename = error("give the path to the tree, root_directory/tree_name"),
        regex::Regex = r"P_[E,X-Z]$"m,
        Nev::Int = -1) # to search momentum
    #
    !(isfile(filename)) && error("no file $filename")
    file = TFile(filename)
    tree = file[treename]
    filtered_tree = getindex(tree, eachindex(tree)) # , nms[filt]
    f.(filtered_tree)
end

invmasssq(p) = p[4]^2-sum(abs2, p[1:3])

##########################################

invariants(pb,pr,pπ,pη) =
    (s0 = invmasssq(pr+pπ+pη),
     s1 = invmasssq(pη+pπ),
     s2 = invmasssq(pπ+pr),
     t1 = invmasssq(pb-pη),
     t2 = invmasssq(pb-pη-pπ))

# λ(x,y,z) = x^2+y^2+z^2-2x*y-2y*z-2z*x
function cosθη_of_s1t1(vars; m1 = mη, m2 = mπ)
    @unpack s1, t1, t2 = vars
    return (2s1*(t1-m2^2-m1^2)+(s1-t2+mb^2)*(s1+m1^2-m2^2)) / sqrt(λ(s1,t2,mb^2)*λ(s1,m1^2,m2^2))
end

function roty(p,θ)
    c,s = cos(θ),sin(θ)
    p3,p1 = [c -s; s c] * [p[3],p[1]];
    SVector(p1,p[2],p3,p[4])
end
function rotz(p,ϕ)
    c,s = cos(ϕ),sin(ϕ)
    (p1,p2) = [c -s; s c] * [p[1],p[2]]
    SVector(p1,p2,p[3],p[4])
end
function boostz(p,γ)
    βγ = sqrt(γ^2-1)*sign(γ)
    γ = abs(γ)
    (p3,p4) = [γ βγ; βγ γ] * [p[3],p[4]]
    SVector(p[1],p[2],p3,p4)
end
function move_to_rest_frame(p,p0)
    cosθ = p0[3]/sqrt(sum(abs2,p0[1:3]))
    ϕ = atan(p0[2],p0[1])
    γ = p0[4]/sqrt(invmasssq(p0))
    boostz(roty(rotz(p,-ϕ),-acos(cosθ)),-γ)
end
function ϕTY(vectors)
    @unpack pb,pη,pπ,pr = vectors
    pηπ = pη+pπ
    (pb,pη,pπ,pr) = [move_to_rest_frame(p, pηπ) for p in (pb,pη,pπ,pr)]
    z = pb[1:3]; z ./= sqrt(sum(abs2,z))
    xz = -pr[1:3]
    x = xz - z*(xz'*z); x ./= sqrt(sum(abs2,x))
    y = z × x; !(sum(abs2, y) ≈ 1.0) && error("??")
    ϕ = atan(pη[1:3]'*y, pη[1:3]'*x)
    return ϕ
end
