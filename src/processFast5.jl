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
	end

	return parse_args(s)
end

function main()
	parsed_args = parse_commandline()
	recursive = parsed_args["recursive"]
	inputdir = parsed_args["inputdir"]
	outputdir = parsed_args["outputdir"]
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

	output_table = mapreduce(f->readFast5(f), vcat, filelist)
	output_table |> CSV.write(length2qualityTextFile)
	totalLength = totalLen(output_table.length)
	n50 = calculate_N50(output_table.length, totalLength)
	generateStatSummary(output_table, totalLength, n50, outputdirs)
end

main()
