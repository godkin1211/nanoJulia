module nanoJulia
using Statistics, FASTX, DataFrames, Printf, Formatting
export nanoread
# Readinfo types
include("datatype.jl")
include("utilities.jl")
include("record_handler.jl")
include("info_handler.jl")
include("statistics_handler.jl")

# FASTQ File Reader
function nanoread(input::FASTX.FASTQ.Reader)
	records = collect(input);
	remove_bad_read!(records);
	recordnum = length(records);
	readsinfo = Array{Union{FastqInfo, Missing},1}(missing, recordnum);
	Threads.@threads for i in 1:recordnum
		@inbounds readsinfo[i] = get_info(records[i])
	end
	extract_info(readsinfo)
end

end
