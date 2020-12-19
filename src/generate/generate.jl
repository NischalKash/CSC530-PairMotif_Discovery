
import Pkg

using Random


function generate_motifs(t, k, m, d)
    generated_motifs = Dict()
    neucleotides = ["A","T","G","C"]
    ref_motifs = []
    for i in 1:m
        motif_char = neucleotides[rand(1:length(neucleotides))]
        motif = repeat(motif_char, k)
        deleteat!(neucleotides, findfirst(neucleotides.== motif_char))
        push!(ref_motifs, motif)
        generated_motifs[motif] = []
    end
    for i in 1:t
        for n in 1:m
            motif = ref_motifs[n]
            for j in 1:d
                pos = rand(1:k)
                neucleotide = neucleotides[rand(1:length(neucleotides))]
                motif = motif[1:pos-1] * neucleotide * motif[pos+1:end]
            end
            push!(generated_motifs[ref_motifs[n]], motif)
        end
    end
    return (ref_motifs, generated_motifs)
end

function Arti_DNA(t, l, k, d, gc, m)  
    DNA = []
    position = Dict()
    (ref_motifs, generated_motifs) = generate_motifs(t, k, m, d)
    
    for consensus in ref_motifs
        position[consensus] = Dict()
    end

    for i in 1:t
        m = length(ref_motifs)
        rand_str = randstring("AT", 100-gc)*randstring("GC", gc)
        seq = randstring(rand_str, l)
        push!(DNA, seq)

        positions = collect(1:k:l-k)
        for n in 1:m
            pos = rand(positions)
            motif = generated_motifs[ref_motifs[n]][i]
            DNA[i] = DNA[i][1:pos-1] * motif * DNA[i][pos+k:end]
            deleteat!(positions, findfirst(positions.== pos))
            position[ref_motifs[n]][i] = [motif, pos]
        end
    end
    return(position, ref_motifs, DNA)
end






