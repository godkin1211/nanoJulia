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


# Miltiple-reads FAST5 File Reader
function nanoread(input::HDF5File, h5ver::String)
    readIDs = names(input)
	readsnum = length(readIDs)
	readsinfo = Array{Union{FastqInfo, Missing},1}(missing, readsnum);
	@inbounds for i in 1:readsnum
		thisread = readIDs[i]
		readrecord = input[thisread]
		readfastq = parse_f5read_record(readrecord, "Basecall_1D_00")
		if FASTQ.seqlen(readfastq) != length(FASTQ.quality(readfastq))
			continue
		end
		readsinfo[i] = get_info(readfastq)
	end
	filter!(e -> e !== missing, readsinfo)
	extract_info(readsinfo)
end


# Single-read FAST5 File Reader
function nanoread(input::HDF5File, h5ver::Float64)
	readfastq = parse_f5read_record(input, "Basecall_1D_000")
	if FASTQ.seqlen(readfastq) != length(FASTQ.quality(readfastq))
		return missing
	else
		readinfo = get_info(readfastq)
		return DataFrame(quality = readinfo.quality, length = readinfo.length)
	end
end