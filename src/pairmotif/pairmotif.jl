import Pkg
Pkg.add("Combinatorics")

using Combinatorics

include("../helpers.jl")

function pairmotif(dna, k, d_thres)
    s1 = dna[1]
    M = []
    motif_positions = Dict()
    for i in 1:length(dna[1]) - k + 1
        x = dna[1][i:i+k-1]
        Cxsi = get_C(x, dna, d_thres)
        r = rand(2:length(dna))
        sr = dna[r]
        Cxsr = Cxsi[r]
        Md = Dict()
        for xp in Cxsr
            pair = (x, xp[1])
            Md[pair] = []
            get_all_md_kmers(x, xp[1], d_thres, ['A','T','G','C'], k, Md[pair])

            Cpxsi = Dict()
            for zi in 2:length(dna)
                if zi != r
                    Cpxsi[zi] = []
                    for z in Cxsi[zi]
                        if filter1(x, xp[1], z[1], d_thres) && filter2(x, xp[1], z[1], d_thres)
                            push!(Cpxsi[zi], [z[1],z[2]])
                        end
                    end
                end
            end
            for y in Md[pair]
                verification_results = verify(y, Cpxsi, r, d_thres)
                if verification_results[1] == true
                    if !(y in M)
                        push!(M, y)
                        verification_results[2][1] = [pair[1],i]
                        verification_results[2][r] = [pair[2],xp[2]]
                        motif_positions[y] = verification_results[2]
                    end
                end
            end
            
        end
    end
    return M, motif_positions
end