push!(LOAD_PATH, "src/")
using ArgParse, HDF5, DataFrames, CSV, nanoJulia

# Arguments parser
function parse_commandline()
	s = ArgParseSettings()
	@add_arg_table! s begin
		"--recursive", "-r"
			help = "Search recursively through folders for fast5 files"
			action = :store_true
		"--inputdir", "-i"
			help = "The path to directory of fast5 files"
			arg_type = String
			required = true
		"--outputdir", "-o"
			help = "Folder to output QC results"
			required = true
		"--basecall_group", "-b"
			help = "FAST5 group obtain original basecalls (under Analysesgroup). Default: Basecall_1D_000"
			arg_type = String
			default = "Basecall_1D_000"
	end

	return parse_args(s)
end

function main()
	parsed_args = parse_commandline()
	recursive = parsed_args["recursive"]
	inputdir = parsed_args["inputdir"]
	outputdir = parsed_args["outputdir"]
	basecall_group = parsed_args["basecall_group"]

	if isdir(outputdir)
		println("\033[1;31m[Warning] The output folder has already existed, so it will be renamed with additional suffix \"_bak\"\033[0m")
		new_outputdir = string(outputdir, "_bak")
		mv(outputdir, new_outputdir)
	end
	mkdir(outputdir)

	length2qualityTextFile = joinpath(outputdir, "readlength_vs_readquality.tsv") 
	if recursive
		println("\033[1;32m* Search fast5 files recursively.......\033[0m")
		filelist = []
		for (root, dirs, files) in walkdir(inputdir)
			for file in files
				if occursin(r"fast5$|f5$", file)
					push!(filelist, joinpath(root, file))
				end
			end
		end
	else
		println("\033[1;32m* Search fast5 files.......\033[0m")
		f5files = [f5 for f5 in readdir(inputdir) if occursin(r"fast5$|f5$", f5)]
		filelist = joinpath.(inputdir, f5files)
	end
	numFiles = length(filelist)
	println("\033[1m  - There're $numFiles fast5 files required to be parsed!\033[0m")
	println("\033[1;32m* Start parsing fast5 files...\033[0m")
	output_table = mapreduce(f -> readFast5(f, basecall_group, false), vcat, filelist)
	println("\033[1;32m* Writing output file...\033[0m")
	output_table |> CSV.write(length2qualityTextFile)
	totalLength = nanoJulia.totalLen(output_table.length)
	n50 = nanoJulia.calculate_N50(output_table.length, totalLength)
	println("\033[1;32m* Generating summry file and figures...\033[0m")
	generateStatSummary(output_table, totalLength, n50, outputdir)
end

main()
