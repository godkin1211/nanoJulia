# Remove the reads with the length of sequence which is different from the length of quality
function remove_bad_read!(records::Array{FASTX.FASTQ.Record, 1})
	seq_lengths = FASTQ.seqlen.(records)
	qual_lengths = map(i -> length(FASTQ.quality(i)), records)
	return records[seq_lengths .==  qual_lengths]
end


# Methods for parsing H5Record

function parse_f5read_record end

function parse_f5read_record(record::HDF5Group, signalpath::String)::FASTX.FASTQ.Record
	fastqpath = string("Analyses/$signalpath/BaseCalled_template/Fastq");
	fastq = FASTQ.Record(read(record, fastqpath))
end

function parse_f5read_record(record::HDF5.File, signalpath::String)::FASTX.FASTQ.Record
	fastqpath = string("Analyses/$signalpath/BaseCalled_template/Fastq");
	fastq = FASTQ.Record(read(record, fastqpath))
end


# Methods for Record-parsing
function get_info end

# Parse HDF5.Group
function get_info(record::HDF5.Attributes)
	readqual = read(record, "mean_qscore")
	readlen = read(record, "sequence_length")
	FastqInfo(readqual, readlen)
end

# Parse BAM.Record
function get_info(record::XAM.BAM.Record)
	flag_check = BAM.flag(record) |> Int
	qual_check = BAM.quality(record) |> length
	if flag_check == 4
		return nothing
	end

	if qual_check == 0
		return nothing
	end

	nm = record["NM"] |> Int
	matches, gaps, cgaps = BAM.cigar_rle(record) |> parse_cigar
	readident = calculate_identity(nm, matches, gaps, cgaps)
	readlen = BAM.seqlength(record)
	readqual = BAM.quality(record) |> calculate_mean_readqual
	return BAMInfo(readqual, readlen, readident)
end

# Parse FASTQ.Record
function get_info(record::FASTX.FASTQ.Record)
	readqual = FASTQ.quality(record) |> calculate_mean_readqual
	readlen = FASTQ.seqlen(record)
	FastqInfo(readqual, readlen)
end