push!(LOAD_PATH, "src/")
using ArgParse,XAM, BioAlignments, nanoJulia, DataFrames, CSV

# Arguments parser
function parse_commandline()
	s = ArgParseSettings()
	@add_arg_table! s begin
		"--inputbam", "-i"
			help = "The path to BAM file"
			arg_type = String
			required = true
		"--outputdir", "-o"
			help = "Folder to output QC results"
			required = true
	end

	return parse_args(s)
end

function main()
	parsed_args = parse_commandline()
	inputbamfile = parsed_args["inputbam"]
    outputdir = parsed_args["outputdir"]
    bamindexfile = string(inputbamfile, ".bai")
    if !isfile(bamindexfile)
        println("\033[1;31m[Error] Please make sure if the index of input BAM file exists!\033[0m")
        exit(1)
    end

	if isdir(outputdir) 
		println("\033[1;33m[Warning] The output folder has already existed, so it will be renamed with additional suffix \"_bak\"\033[0m")
		new_outputdir = string(outputdir, "_bak")
		mv(outputdir, new_outputdir)
	end
	mkdir(outputdir)

	length2qualityTextFile = joinpath(outputdir, "readlength_vs_readquality.tsv") 

    println("\033[1;32m* Start parsing this BAM file...\033[0m")
	output_table = readBAM(inputbamfile)
	println("\033[1;32m* Writing output file...\033[0m")
	output_table |> CSV.write(length2qualityTextFile)
	totalLength = nanoJulia.totalLen(output_table.length)
	n50 = nanoJulia.calculate_N50(output_table.length, totalLength)
	println("\033[1;32m* Generating summry file and figures...\033[0m")
	generateStatSummary(output_table, totalLength, n50, outputdir)
end

main()