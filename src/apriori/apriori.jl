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
