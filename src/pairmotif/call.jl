import Pkg
Pkg.add("ArgParse")
Pkg.add("JSON")

using ArgParse
using JSON

include("pairmotif.jl")
include("../apriori/apriori.jl")

function parse_arguments()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--dna-file"
            help = "Relative path of the DNA file"
            arg_type = String
            default = "../data/test/dna.json"
        "--k"
            help = "Length of the motifs"
            arg_type = Int
            default = 8
        "--d-threshold"
            help = "Hamming distance threshold"
            arg_type = Int
            default = 0
        "--output-dir"
            help = "Directory to generate the result in"
            arg_type = String
            default = "../data/test"
    end
    return parse_args(s)
end

function main()
    parsed_args = parse_arguments()
    dna_file = parsed_args["dna-file"]
    k = parsed_args["k"]
    d_thres = parsed_args["d-threshold"]
    output_dir = parsed_args["output-dir"]

    f = open(dna_file)
    lines = read(f, String)
    dna = JSON.parse(lines)

    (motifs, motif_positions) = pairmotif(dna, k, d_thres)
    if !isdir(output_dir)
        mkdir(output_dir)
    end
    apriori_motifs = apriori(motifs)

    new_motif_positions = Dict()
    for i in apriori_motifs
        new_motif_positions[i] = motif_positions[i]
    end

    meta_data = Dict()
    meta_data["motif_positions"] = new_motif_positions
    meta_data["motifs"] = apriori_motifs

    write(output_dir*"/result.json", JSON.json(meta_data))
end

main()