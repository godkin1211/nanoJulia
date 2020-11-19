# Define a generic function
function nanoread end

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


# BAM File Reader
function nanoread(input::XAM.BAM.Reader)
    readsinfo = BAMInfo[]
    record = BAM.Record()
    while !eof(input)
        empty!(record)
        read!(input, record)
        push!(readsinfo, get_info(record))
    end
	extract_info(readsinfo)
end