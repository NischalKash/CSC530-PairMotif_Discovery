import Pkg
Pkg.add("ArgParse")
Pkg.add("JSON")
Pkg.add("Plots")

using ArgParse
using JSON
using Plots


include("generate/generate.jl")
include("pairmotif/pairmotif.jl")

function run_results(output, implanted, d, t)
    
    average_motifs_in_segments_found = 0
    average_position_found_error = 0
    average_position_error = 0
    average_consensus_identified = 0

    for (consensus, segments) in implanted["motif_positions"]
        for (segment, implants) in segments
            implant =  implants[1]
            implant_pos = implants[2]
            for implanted_consensus in output["motif_positions"]
                identified = output["motif_positions"][consensus][segment][1]
                identified_pos = output["motif_positions"][consensus][segment][2]
                average_position_error += abs(implant_pos - identified_pos)
                if d_hamm(implant, identified) <= d
                    average_motifs_in_segments_found += 1
                    average_position_found_error += abs(implant_pos - identified_pos)
                    break
                end
            end
        end
    end

    for consensus in implanted["motifs"]
        for identified_consensus in output["motifs"]
            if d_hamm(consensus, identified_consensus) <= d
                average_consensus_identified += 1
                deleteat!(output["motifs"], output["motifs"].==identified_consensus)
            end
        end
    end

    total_motifs = t * length(implanted["motifs"])

    average_motifs_in_segments_found = 100*average_motifs_in_segments_found/total_motifs
    average_position_found_error = 100*average_position_found_error/total_motifs
    average_position_error = average_position_error/total_motifs
    average_consensus_identified = 100*average_consensus_identified/length(implanted["motifs"])

    return (average_motifs_in_segments_found, average_position_found_error, average_position_error, average_consensus_identified)
end


function generate_graphs(temp, y_val, x_val, x_label)
    plot(x_val, y_val["times"], label = "Average time taken", ylabel = "Average time taken (s)", xlabel=x_label)
    savefig(temp*"/time-plot.png")

    plot(x_val, y_val["memory"], label = "Average Memory Allocated", ylabel = "Average Memory Allocate (bytes)", xlabel=x_label)
    savefig(temp*"/memory-plot.png")

    plot(x_val, y_val["motifs_in_segments_found_arr"], label = "Percentage of motifs found in each segment", ylabel = "% of motif found", xlabel=x_label)
    savefig(temp*"/motifs-in-segments-found-plot.png")

    plot(x_val, y_val["position_found_error_arr"], label = "Percentage of error in motifs found in each segment where motifs found", ylabel = "% of position error", xlabel=x_label)
    savefig(temp*"/position-found-error-plot.png")

    plot(x_val, y_val["position_error_arr"], label = "Overall Percentage of error in motifs found in each segment", ylabel = "% of position error", xlabel=x_label)
    savefig(temp*"/position-error-plot.png")

    plot(x_val, y_val["consensus_identified_arr"], label = "Percentage of consensus motifs found", ylabel = "% of consensus found", xlabel=x_label)
    savefig(temp*"/consensus-identified-plot.png")
end


function run_evaluations(eval_file, temp_dir)

    f = open(eval_file)
    lines = read(f, String)
    evaluation_options = JSON.parse(lines)
    exp_count = evaluation_options["exp-count"]

    for (eval_proc, params) in evaluation_options["eval-procs"]
        temp = temp_dir * "/" * evaluation_options["algo"] * "/" * eval_proc

        if !isdir(temp)
            mkpath(temp)
        end

        motif_settings = params["generate-settings"]
        if eval_proc == "motif-length"
            k_range = split(motif_settings["k"],":")
            k_range = collect(parse(Int, k_range[1]):parse(Int, k_range[2]):parse(Int, k_range[3]))
            t = motif_settings["t"]
            l = motif_settings["l"]
            d = motif_settings["d"]
            gc = motif_settings["gc"]
            m = motif_settings["m"]
            d_thres = params["algo-settings"]["d_thres"]

            y_val = Dict()
            y_val["times"] = []
            y_val["memory"] = []
            y_val["motifs_in_segments_found_arr"] = []
            y_val["position_found_error_arr"] = []
            y_val["position_error_arr"] = []
            y_val["consensus_identified_arr"] = []

            for k in k_range
                average_time = 0
                average_memory = 0
                motifs_in_segments_found = 0
                position_found_error = 0
                position_error = 0
                consensus_identified = 0
                for c in 1:exp_count
                    (positions, motifs, dna) = Arti_DNA(t, l, k, d, gc, m)
                    meta_data = Dict()
                    meta_data["motif_positions"] = positions
                    meta_data["motifs"] = motifs
                    
                    temp_path = temp * "/" * string(k) * "/" * string(c)
                    if !isdir(temp_path)
                        mkpath(temp_path)
                    end
                    write(temp_path*"/meta.json", JSON.json(meta_data))
                    write(temp_path*"/dna.json", JSON.json(dna))

                    ((motifs, motif_positions), time, bytes) = @timed(pairmotif(dna, k, d_thres))
                    result = Dict()
                    result["motif_positions"] = motif_positions
                    result["motifs"] = motifs
                    result["time"] = time
                    result["memory"] = bytes
                    write(temp_path*"/result.json", JSON.json(result))

                    (average_motifs_in_segments_found, average_position_found_error, average_position_error, average_consensus_identified) = run_results(result, meta_data, d, t)

                    average_time += time
                    average_memory += bytes
                    motifs_in_segments_found += average_motifs_in_segments_found
                    position_found_error += average_position_found_error
                    position_error += average_position_error
                    consensus_identified += average_consensus_identified
                end
                average_time = average_time/exp_count
                average_memory = average_memory/exp_count

                motifs_in_segments_found = motifs_in_segments_found/exp_count
                position_found_error = position_found_error/exp_count
                position_error = position_error/exp_count
                consensus_identified = consensus_identified/exp_count

                y_val["times"] = [y_val["times"]; average_time]
                y_val["memory"] = [y_val["memory"]; average_memory]
                y_val["motifs_in_segments_found_arr"] = [y_val["motifs_in_segments_found_arr"]; motifs_in_segments_found]
                y_val["position_found_error_arr"] = [y_val["position_found_error_arr"]; position_found_error]
                y_val["position_error_arr"] = [y_val["position_error_arr"]; position_error]
                y_val["consensus_identified_arr"] = [y_val["consensus_identified_arr"]; consensus_identified]
            end
            generate_graphs(temp, y_val, k_range, "Length of Motifs")
            

        elseif eval_proc == "num-segments"
            t_range = split(motif_settings["t"],":")
            t_range = collect(parse(Int, t_range[1]):parse(Int, t_range[2]):parse(Int, t_range[3]))
            k = motif_settings["k"]
            l = motif_settings["l"]
            d = motif_settings["d"]
            gc = motif_settings["gc"]
            m = motif_settings["m"]
            d_thres = params["algo-settings"]["d_thres"]

            y_val = Dict()
            y_val["times"] = []
            y_val["memory"] = []
            y_val["motifs_in_segments_found_arr"] = []
            y_val["position_found_error_arr"] = []
            y_val["position_error_arr"] = []
            y_val["consensus_identified_arr"] = []

            for t in t_range
                average_time = 0
                average_memory = 0
                motifs_in_segments_found = 0
                position_found_error = 0
                position_error = 0
                consensus_identified = 0
                for c in 1:exp_count
                    (positions, motifs, dna) = Arti_DNA(t, l, k, d, gc, m)
                    meta_data = Dict()
                    meta_data["motif_positions"] = positions
                    meta_data["motifs"] = motifs
                    
                    temp_path = temp * "/" * string(t) * "/" * string(c)
                    if !isdir(temp_path)
                        mkpath(temp_path)
                    end

                    write(temp_path*"/meta.json", JSON.json(meta_data))
                    write(temp_path*"/dna.json", JSON.json(dna))

                    ((motifs, motif_positions), time, bytes) = @timed(pairmotif(dna, k, d_thres))
                    result = Dict()
                    result["motif_positions"] = motif_positions
                    result["motifs"] = motifs
                    result["time"] = time
                    result["memory"] = bytes
                    write(temp_path*"/result.json", JSON.json(result))

                    (average_motifs_in_segments_found, average_position_found_error, average_position_error, average_consensus_identified) = run_results(result, meta_data, d, t)

                    average_time += time
                    average_memory += bytes
                    motifs_in_segments_found += average_motifs_in_segments_found
                    position_found_error += average_position_found_error
                    position_error += average_position_error
                    consensus_identified += average_consensus_identified
                end
                average_time = average_time/exp_count
                average_memory = average_memory/exp_count

                motifs_in_segments_found = motifs_in_segments_found/exp_count
                position_found_error = position_found_error/exp_count
                position_error = position_error/exp_count
                consensus_identified = consensus_identified/exp_count

                y_val["times"] = [y_val["times"]; average_time]
                y_val["memory"] = [y_val["memory"]; average_memory]
                y_val["motifs_in_segments_found_arr"] = [y_val["motifs_in_segments_found_arr"]; motifs_in_segments_found]
                y_val["position_found_error_arr"] = [y_val["position_found_error_arr"]; position_found_error]
                y_val["position_error_arr"] = [y_val["position_error_arr"]; position_error]
                y_val["consensus_identified_arr"] = [y_val["consensus_identified_arr"]; consensus_identified]
            end
            generate_graphs(temp, y_val, t_range, "Number of segments")


        elseif eval_proc == "dna-length"
            l_range = split(motif_settings["l"],":")
            l_range = collect(parse(Int, l_range[1]):parse(Int, l_range[2]):parse(Int, l_range[3]))
            k = motif_settings["k"]
            t = motif_settings["t"]
            d = motif_settings["d"]
            gc = motif_settings["gc"]
            m = motif_settings["m"]
            d_thres = params["algo-settings"]["d_thres"]

            y_val = Dict()
            y_val["times"] = []
            y_val["memory"] = []
            y_val["motifs_in_segments_found_arr"] = []
            y_val["position_found_error_arr"] = []
            y_val["position_error_arr"] = []
            y_val["consensus_identified_arr"] = []


            for l in l_range
                average_time = 0
                average_memory = 0
                motifs_in_segments_found = 0
                position_found_error = 0
                position_error = 0
                consensus_identified = 0
                for c in 1:exp_count
                    (positions, motifs, dna) = Arti_DNA(t, l, k, d, gc, m)
                    meta_data = Dict()
                    meta_data["motif_positions"] = positions
                    meta_data["motifs"] = motifs
                    
                    temp_path = temp * "/" * string(l) * "/" * string(c)
                    if !isdir(temp_path)
                        mkpath(temp_path)
                    end

                    write(temp_path*"/meta.json", JSON.json(meta_data))
                    write(temp_path*"/dna.json", JSON.json(dna))

                    ((motifs, motif_positions), time, bytes) = @timed(pairmotif(dna, k, d_thres))
                    result = Dict()
                    result["motif_positions"] = motif_positions
                    result["motifs"] = motifs
                    result["time"] = time
                    result["memory"] = bytes
                    write(temp_path*"/result.json", JSON.json(result))

                    (average_motifs_in_segments_found, average_position_found_error, average_position_error, average_consensus_identified) = run_results(result, meta_data, d, t)

                    average_time += time
                    average_memory += bytes
                    motifs_in_segments_found += average_motifs_in_segments_found
                    position_found_error += average_position_found_error
                    position_error += average_position_error
                    consensus_identified += average_consensus_identified
                end

                average_time = average_time/exp_count
                average_memory = average_memory/exp_count

                motifs_in_segments_found = motifs_in_segments_found/exp_count
                position_found_error = position_found_error/exp_count
                position_error = position_error/exp_count
                consensus_identified = consensus_identified/exp_count

                y_val["times"] = [y_val["times"]; average_time]
                y_val["memory"] = [y_val["memory"]; average_memory]
                y_val["motifs_in_segments_found_arr"] = [y_val["motifs_in_segments_found_arr"]; motifs_in_segments_found]
                y_val["position_found_error_arr"] = [y_val["position_found_error_arr"]; position_found_error]
                y_val["position_error_arr"] = [y_val["position_error_arr"]; position_error]
                y_val["consensus_identified_arr"] = [y_val["consensus_identified_arr"]; consensus_identified]
            end
            generate_graphs(temp, y_val, l_range, "DNA Length")
        end
    end
end


function parse_arguments()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--eval_file"
            help = "File containing the evaluation options"
            arg_type = String
            default = "../data/eval.json"
        "--temp-dir"
            help = "Temporary directory to store the generated outputs"
            arg_type = String
            default = "../data/test"
    end
    return parse_args(s)
end

function main()
    parsed_args = parse_arguments()
    eval_file = parsed_args["eval_file"]
    temp_dir = parsed_args["temp-dir"]
    run_evaluations(eval_file, temp_dir)
end

main()
