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