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

	if !isdir(outputdir) 
		mkdir(outputdir)
	end

	length2qualityTextFile = joinpath(outputdir, "readlength_vs_readquality.tsv") 
	if recursive
		println("Recursive mode launched!")
		filelist = []
		for (root, dirs, files) in walkdir(inputdir)
			for file in files
				if occursin(r"fast5$", file)
					push!(filelist, joinpath(root, file))
				end
			end
		end
	else
		println("Non-recursive mode launched!")
		filelist = joinpath.(inputdir, readdir(inputdir))
	end

	output_table = mapreduce(f -> readFast5(f, basecall_group, false), vcat, filelist)
	output_table |> CSV.write(length2qualityTextFile)
	totalLength = nanoJulia.totalLen(output_table.length)
	n50 = nanoJulia.calculate_N50(output_table.length, totalLength)
	generateStatSummary(output_table, totalLength, n50, outputdir)
end

main()
