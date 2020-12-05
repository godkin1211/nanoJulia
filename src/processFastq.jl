push!(LOAD_PATH, "src/")
using ArgParse,FASTX, nanoJulia, DataFrames, CSV

# Arguments parser
function parse_commandline()
	s = ArgParseSettings()
	@add_arg_table! s begin
		"--recursive", "-r"
			help = "Search recursively through folders for fastq files"
			action = :store_true
		"--inputdir", "-i"
			help = "The path to directory of fastq files"
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

	if isdir(outputdir) 
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
				if occursin(r"fastq$|fq$", file)
					push!(filelist, joinpath(root, file))
				end
			end
		end
	else
		println("\033[1;32m* Search fast5 files.......\033[0m")
		fqfiles = [fq for fq in readdir(inputdir) if occursin(r"fastq$|fq$", fq)]
		filelist = joinpath.(inputdir, readdir(inputdir))
	end

	numFiles = length(filelist)
	println("\033[1m  - There're $numFiles fastq files required to be parsed!\033[0m")
	println("\033[1;32m* Start parsing fastq files...\033[0m")
	output_table = mapreduce(f -> readFastq(f), vcat, filelist)
	println("\033[1;32m* Writing output file...\033[0m")
	output_table |> CSV.write(length2qualityTextFile)
	totalLength = nanoJulia.totalLen(output_table.length)
	n50 = nanoJulia.calculate_N50(output_table.length, totalLength)
	println("\033[1;32m* Generating summry file and figures...\033[0m")
	generateStatSummary(output_table, totalLength, n50, outputdir)
end

main()