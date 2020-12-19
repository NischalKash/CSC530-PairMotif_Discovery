import Pkg
Pkg.add("ArgParse")
Pkg.add("JSON")

using ArgParse
using JSON

include("pairmotif-faster.jl")


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
        "--n"
            help = "Expected motifs in each strand"
            arg_type = Int
            default = 1
        "--use-prob"
            help = "strategy for initial screening"
            arg_type = Bool
            default = true
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
    n = parsed_args["n"]
    use_prob = Bool(parsed_args["use-prob"])
    d_thres = parsed_args["d-threshold"]
    output_dir = parsed_args["output-dir"]

    f = open(dna_file)
    lines = read(f, String)
    dna = JSON.parse(lines)

    (motifs, motif_positions) = pairmotif(dna, k, n, use_prob, d_thres)
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