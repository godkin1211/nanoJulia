# Define a generic function
function nanoread end

# SequencingSummary Reader
function nanoread(input::DataFrame)
	return DataFrame(quality=round.(input.mean_qscore_template, digits=1), length=input.sequence_length_template)
end


# Fastq File Reader
function nanoread(input::FASTX.FASTQ.Reader)
	records = collect(input);
	remove_bad_read!(records);
	recordnum = length(records);
	readsinfo = Array{Union{FastqInfo, Missing},1}(missing, recordnum);
	Threads.@threads for i in 1:recordnum
		@inbounds readsinfo[i] = get_info(records[i])
	end
	return extract_info(readsinfo)
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
	return extract_info(readsinfo)
end


# Miltiple-reads FAST5 File Reader
function nanoread(input::HDF5File, h5ver::String, basecallGroup::String)
	dirpath = "Analyses/$basecallGroup/Summary/basecall_1d_template"
    readIDs = names(input)
	readsnum = length(readIDs)
	readsinfo = Array{Union{FastqInfo, Missing},1}(missing, readsnum)
	@inbounds for i in 1:readsnum
		thisread = readIDs[i]
		readrecord = input[thisread]
		summaryInfo = attrs(readrecord[dirpath])
		readsinfo[i] = get_info(summaryInfo)
	end
	filter!(e -> e !== missing, readsinfo)
	return extract_info(readsinfo)
end


# Single-read FAST5 File Reader
function nanoread(input::HDF5.File, h5ver::Float64, basecallGroup::String)
	dirpath = "Analyses/$basecallGroup/Summary/basecall_1d_template"
	summaryInfo = attrs(input[dirpath])
	readinfo = get_info(summaryInfo)
	return DataFrame(quality = readinfo.quality, length = readinfo.length)
end
