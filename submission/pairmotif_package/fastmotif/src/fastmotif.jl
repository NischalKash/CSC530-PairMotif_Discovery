module fastmotif

function d_hamm(x,x1)
    i=1
    count = 0
    while i<=length(x)
        if x[i]!=x1[i]
            count+=1
        end
        i+=1
    end
    return count
end

function get_C(x, dna, d)
    Cxsi = Dict()
    t = length(dna)
    l = length(dna[1])
    k = length(x)
    for i in 2:t
        C = Array[]
        for j in 1:(l-k+1)
            xp = dna[i][j:j+k-1]
            if d_hamm(x,xp) <= 2*d
               push!(C, [xp,j])
            end
        end
        Cxsi[i] = C
    end
    return Cxsi
end

function filter1(x, xp, z, d)
    d1 = 0
    d2 = 0
    for i in 1:length(z)
        if x[i] != z[i]
            d1 += 1
        end
        if xp[i] != z[i]
            d2 += 1
        end
    end
    if d1 > 2*d || d2 > 2*d
       return false
    end
    return true
end

function p_val(x,x1,z,d)
    i = 1
    count1 = 0
    count2 = 0
    while i<=length(x)
        if x[i]==x1[i]
            if z[i]!=x[i] || z[i]!=x1[i]
                count1+=1
            end
        else
            if z[i]!=x[i]  || z[i]!=x1[i]
                count2+=1
            end
        end
        i+=1
    end
    return [count1,count2]
end

function filter2(x,x1,z,d)
    alpha_range = length(x)-d_hamm(x,x1)
    beta_range = d_hamm(x,x1)

    p_values = p_val(x,x1,z,d)

    for i in 0:alpha_range
        for j in 0:beta_range
            k = 2*i+j+d_hamm(x,x1)
            if k<=2*d
                if abs(p_values[1]-i)+abs(p_values[2]-j)<=d
                    return true
                end
            end
        end
    end
    return false
end

function get_all_md_kmers(x, xp, d, ch, k, Md)
    push!(Md, x)
    push!(Md, xp)

    pos = collect(combinations(1:k, d))

    for p in pos
        curr = []
        for p_i in p
            push!(curr, x[p_i])
        end
        rem_neucleo = []
        for c in curr
            nucleo = copy(ch)
            push!(rem_neucleo, deleteat!(nucleo, findfirst(nucleo.==c)))
        end

        prods = Base.Iterators.product(1)
        for d_i in 1:d
            prods = Base.Iterators.product(prods, rem_neucleo[d_i])
        end

        for prod in prods
            for d_i in 1:d
                prod = Tuple(Iterators.flatten(prod))
            end
            prod = prod[2:end]
            seq = x
            i = 1
            for p_i in p
                seq = seq[1:p_i-1] * prod[i] * seq[p_i+1:end]
                if d_hamm(seq, x) <= d && d_hamm(seq, xp) <= d
                    push!(Md, seq)
                end
                i += 1
            end
        end
    end
end

function verify(y, Cp, r, d)
    final_dict = Dict()
    for (i, C) in Cp
        valid = false
        if i != r
            for yi in C
                if d_hamm(y, yi[1]) <= d
                    final_dict[i] = [yi[1],yi[2]]
                    valid = true
                    break
                end
            end
            if !valid
                return [false,[]]
            end
        end
    end
    return [true,final_dict]
end

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

function init_screen(DNA, k, n, use_prob)
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
    pro_A = c_A/(c_A + c_C + c_G + c_T)
    pro_C = c_C/(c_A + c_C + c_G + c_T)
    pro_G = c_G/(c_A + c_C + c_G + c_T)
    pro_T = c_T/(c_A + c_C + c_G + c_T)
    pro = []
    D = Dict()
    for j in 1:(l-k+1)
        if use_prob == true
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
        else
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

function pairmotif_faster(dna, k, expected_n, use_prob, d_thres)
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

function apriori(motifs)
    if length(motifs)==0
        return []
    end
    index_count = 1
    while index_count<=length(motifs[1])
        dict = Dict()
        total_length = length(motifs)
        for i in 1:length(motifs)
            if haskey(dict,motifs[i][1:index_count])
                dict[motifs[i][1:index_count]]+=1
            else
                dict[motifs[i][1:index_count]]=1
            end
        end
        new_list = []
        keys_values = keys(dict)
        condition = false
        for i in keys_values
            consideration_val = dict[i]/total_length
            if consideration_val>1/total_length
                condition = true
                for all in 1:length(motifs)
                    if motifs[all][1:index_count]==i
                        push!(new_list,motifs[all])
                    end
                end
            end
        end
        if condition==false
            return motifs
        end
        motifs = new_list
        index_count+=1
    end
    return motifs
end

function get_score(motifs)
    t = size(motifs)[1]
    k = length(motifs[1])
    symbols = ['A','C','G','T']
    scores = zeros(Number, k)
    for i in 1:k
        counts = zeros(Number, t)
        for j in 1:t
            index = findfirst(x->x==motifs[j][i], symbols)
            counts[index] = counts[index] + 1
        end
        scores[i] = sum(counts[findmax(counts)[2]])
    end
    return sum(scores)
end

end # module