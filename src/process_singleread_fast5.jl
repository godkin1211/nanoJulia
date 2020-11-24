push!(LOAD_PATH, "src/")
using ArgParse, HDF5, DataFrames, nanoJulia

# Arguments parser
function parse_commandline()
	s = ArgParseSettings()
	@add_arg_table! s begin
		"--recursive", "-r"
			help = "Search recursively through folders for SingleRead fast5 files"
			action = :store_true
		"--inputdir", "-i"
			help = "The path to directory of SingleRead fast5 files"
			arg_type = String
			required = true
		"--outputdir", "-o"
			help = "Folder to output QC results"
			required = true
	end

	return parse_args(s)
end

function readFast5(fast5file::String)::DataFrames.DataFrame
	h5open(fast5file) do sf5
		file_version = read(attrs(sf5), "file_version")
		sf5_info = nanoread(sf5, file_version)
		return sf5_info
	end
end

function main()
	parsed_args = parse_commandline()
	recursive = parsed_args["recursive"]
	inputdir = parsed_args["inputdir"]
	outputdir = parsed_args["outputdir"]

	if recursive
		println("Recursive mode launched!")
		filelist = String[]
		idx = 0
		for (root, dirs, files) in walkdir(inputdir)
			for file in files
				if occursin(r"fast5$", file)
					push!(filelist, joinpath(root, file))
					if idx == 4
						break
					end
					idx += 1
				end
			end
		end

		output = mapreduce(f->readFast5(f), vcat, filelist)
		println(output)
	else
		println("Non-recursive mode launched!")
	end
end

main()
