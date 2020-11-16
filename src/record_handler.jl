# Remove the reads with the length of sequence which is different from the length of quality
function remove_bad_read!(records::Array{FASTX.FASTQ.Record, 1})
	seq_lengths = FASTQ.seqlen.(records)
	qual_lengths = map(i -> length(FASTQ.quality(i)), records)
	return records[seq_lengths .==  qual_lengths]
end

# Parse FASTQ.Record
function get_info end

function get_info(record::XAM.BAM.Record)
	nm = Int(record["NM"])
	matches, gaps, cgaps = parse_cigar(BAM.cigar_rle(record))
	readident = calculate_identity(nm, matches, gaps, cgaps)
	readlen = BAM.seqlength(record)
	readqual = calculate_mean_readqual(BAM.quality(record))
	BAMInfo(readqual, readlen, readident)
end

function get_info(record::FASTX.FASTQ.Record)
	readqual = calculate_mean_readqual(FASTQ.quality(record))
	readlen = FASTQ.seqlen(record)
	FastqInfo(readqual, readlen)
end