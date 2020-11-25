push!(LOAD_PATH, "src/")
using ArgParse, HDF5, nanoJulia

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

	
end

main()