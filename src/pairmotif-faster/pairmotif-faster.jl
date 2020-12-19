import Pkg
Pkg.add("Combinatorics")

using Combinatorics

include("../helpers.jl")

# probability and entropy
function init_screen(DNA, k, n, use_prob) # DNA, motif, number of candidate motif, whether use prob or entropy
    l = length(DNA)
    c_A = 0
    c_C = 0
    c_G = 0
    c_T = 0
    nucleo = ["ACGT"]
    for j in 1:l
        if DNA[j] == nucleo[1][1]
            c_A = c_A + 1 
        elseif DNA[j] == nucleo[1][2]
            c_C = c_C + 1
        elseif DNA[j] == nucleo[1][3]
            c_G = c_G + 1
        elseif DNA[j] == nucleo[1][4]
            c_T = c_T + 1
        end
    end
    pro_A = c_A/(c_A + c_C + c_G + c_T)           # probability of the whole DNA seq
    pro_C = c_C/(c_A + c_C + c_G + c_T)
    pro_G = c_G/(c_A + c_C + c_G + c_T)
    pro_T = c_T/(c_A + c_C + c_G + c_T)
    pro = []
    D = Dict()
    for j in 1:(l-k+1)
        if use_prob == true                  # calculate probability
            prob = 1
            corr_motif = DNA[j:j+(k-1)]
            for h in 1:k
                if DNA[j+h-1] == nucleo[1][1]
                    prob = prob*pro_A
                elseif DNA[j+h-1] == nucleo[1][2]
                    prob = prob*pro_C
                elseif DNA[j+h-1] == nucleo[1][3]
                    prob = prob*pro_G
                elseif DNA[j+h-1] == nucleo[1][4]
                    prob = prob*pro_T
                end
            end
            push!(pro, prob)
            if prob in keys(D)
                push!(D[prob], [corr_motif,j])
            else
                D[prob] = [[corr_motif,j]]
            end
        else                                  # calculate entropy
            prob = 0
            corr_motif = DNA[j:j+(k-1)]
            for h in 1:k
                if DNA[j+h-1] == nucleo[1][1]
                    prob = prob - pro_A*log2(pro_A)
                elseif DNA[j+h-1] == nucleo[1][2]
                    prob = prob - pro_C*log2(pro_C)
                elseif DNA[j+h-1] == nucleo[1][3]
                    prob = prob - pro_G*log2(pro_G)
                elseif DNA[j+h-1] == nucleo[1][4]
                    prob = prob - pro_T*log2(pro_T)
                end
            end
            push!(pro, prob)
            if prob in keys(D)
                push!(D[prob], [corr_motif,j])
            else
                D[prob] = [[corr_motif,j]]
            end
        end
    end
    cand_motif = []
    pro = sort(unique(pro), rev=true)
    i_n = 1
    for i in pro
        for j in 1:length(D[i])
            push!(cand_motif, D[i][j])
        end
        i_n += 1
        if n != -1 && i_n == n
            break
        end
    end
    return(cand_motif)  
end

function pairmotif(dna, k, expected_n, use_prob, d_thres)
    s1 = dna[1]
    M = []
    motif_positions = Dict()
    X = init_screen(s1, k, expected_n, use_prob)
    for (x,i) in X
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