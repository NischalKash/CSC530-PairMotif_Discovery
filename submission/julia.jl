import Pkg
Pkg.add("Combinatorics")
Pkg.add("Plots")

using Random
using Combinatorics
using Plots

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

lines = []
open("motif_input1.txt") do f
  while ! eof(f)
     s = readline(f)
     push!(lines,s)
  end
end

motif_length = parse(Int64,lines[1])

t=5
l=100
k=motif_length
d=0

dnas = []
push!(dnas,lines[3])
push!(dnas,lines[5])
push!(dnas,lines[7])
push!(dnas,lines[9])
push!(dnas,lines[11])
(motifs, motif_positions) = pairmotif(dnas, k, d)

apriori_motifs = apriori(motifs)
new_motif_positions = Dict()
for i in apriori_motifs
    new_motif_positions[i] = motif_positions[i]
end

motifs = apriori_motifs
motif_positions = new_motif_positions

strand1 = []
strand2 = []
strand3 = []
strand4 = []
strand5 = []

for i in motif_positions
    for j in i[2]
        if j[1]==1
            push!(strand1,j[2][1])
        end
        if j[1]==2
            push!(strand2,j[2][1])
        end
        if j[1]==3
            push!(strand3,j[2][1])
        end
        if j[1]==4
            push!(strand4,j[2][1])
        end
        if j[1]==5
            push!(strand5,j[2][1])
        end
    end
end

position_strand1 = []
position_strand2 = []
position_strand3 = []
position_strand4 = []
position_strand5 = []

for i in motif_positions
    motif_name = i[1]
    for j in i[2]
        if j[1] == 1
            if !(j[2][2] in position_strand1)
                push!(position_strand1,j[2][2])
            end
        end
        if j[1] == 2
            if !(j[2][2] in position_strand2)
                push!(position_strand2,j[2][2])
            end
        end
        if j[1] == 3
            if !(j[2][2] in position_strand3)
                push!(position_strand3,j[2][2])
            end
        end
        if j[1] == 4
            if !(j[2][2] in position_strand4)
                push!(position_strand4,j[2][2])
            end
        end
        if j[1] == 5
            if !(j[2][2] in position_strand5)
                push!(position_strand5,j[2][2])
            end
        end
    end
end

open("motif_output1.txt","w") do io
   println(io,"Start positions of motifs in each strand are as follows")
   println(io,"In strand1")
   for i in position_strand1
    print(io,"$i,")
   end
    println(io,"")
   println(io,"In strand2")
   for i in position_strand2
    print(io,"$i,")
   end
    println(io,"")
   println(io,"In strand3")
   for i in position_strand3
    print(io,"$i,")
   end
    println(io,"")
   println(io,"In strand4")
   for i in position_strand4
    print(io,"$i,")
   end
    println(io,"")
   println(io,"In strand5")
   for i in position_strand5
    print(io,"$i,")
   end
    println(io,"")
    println(io,"")
   println(io,"Consensus Motifs : ")
   for i in motifs
    println(io,"$i")
   end
   println(io,"")
   println(io,"Motifs:")
   println(io,"S1 : ")
   for i in strand1
        println(io,"$i")
   end
   println(io,"")
   println(io,"S2 : ")
   for i in strand2
        println(io,"$i")
   end
   println(io,"")
   println(io,"S3 : ")
   for i in strand3
        println(io,"$i")
   end
   println(io,"")
   println(io,"S4 : ")
   for i in strand4
        println(io,"$i")
   end
   println(io,"")
   println(io,"S5 : ")
   for i in strand5
        println(io,"$i")
   end
end