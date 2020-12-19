import Pkg
Pkg.add("ArgParse")
Pkg.add("JSON")

using ArgParse
using JSON

include("generate.jl")

function parse_arguments()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--t"
            help = "Number of dna segments"
            arg_type = Int
            default = 4
        "--l"
            help = "Length of each DNA Segments"
            arg_type = Int
            default = 20
        "--k"
            help = "Size of motifs"
            arg_type = Int
            default = 8
        "--d"
            help = "Hamming distance the motifs should be apart"
            arg_type = Int
            default = 0
        "--gc"
            help = "Percentage of GC content in DNA"
            arg_type = Int
            default = 8
        "--m"
            help = "Number of motifs per segment"
            arg_type = Int
            default = 1
        "--output-dir"
            help = "Directory to generate the file in"
            arg_type = String
            default = "../data/test"
    end
    return parse_args(s)
end

function main()
    parsed_args = parse_arguments()
    t = parsed_args["t"]
    l = parsed_args["l"]
    k = parsed_args["k"]
    d = parsed_args["d"]
    gc = parsed_args["gc"]
    m = parsed_args["m"]
    output_dir = parsed_args["output-dir"]
    (positions, motifs, dna) = Arti_DNA(t, l, k, d, gc, m)
    if !isdir(output_dir)
        mkdir(output_dir)
    end
    meta_data = Dict()
    meta_data["motif_positions"] = positions
    meta_data["motifs"] = motifs
    write(output_dir*"/meta.json", JSON.json(meta_data))
    write(output_dir*"/dna.json", JSON.json(dna))
end

main()