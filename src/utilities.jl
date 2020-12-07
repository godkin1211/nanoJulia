# Calculate GC content
function calculate_GC_content end

function calculate_GC_content(read::FASTX.FASTQ.Record)::Float64
	seq = FASTQ.sequence(read)
	seqlen = FASTQ.seqlen(read)
	nucl_comp = composition(seq)
	gc_content = (nucl_comp[DNA_C] + nucl_comp[DNA_G]) / seqlen |> c->round(c, digits=4)
	return gc_content
end

function calculate_GC_content(read::XAM.BAM.Record)::Float64
	seq = BAM.sequence(read)
	seqlen = BAM.seqlength(read)
	nucl_comp = composition(seq)
	gc_content = (nucl_comp[DNA_C] + nucl_comp[DNA_G]) / seqlen |> c->round(c, digits=4)
	return gc_content
end

# Convert Phred Score into error probabilities
phred2prob(q::UInt8) = 10^(q/(-10))

# Convert error probabilities into Phred Score
prob2phred(p::Float64) = -10log10(p)

# Calculate mean read phred quality score
function calculate_mean_readqual(quals::Array{UInt8, 1})::Float64
	quals .|> phred2prob |> mean |> prob2phred |> q->round(q, digits = 1)
end

# Parse CIGAR string from BAM
function parse_cigar(cigar::Tuple{Array{BioAlignments.Operation,1}, Array{Int64,1}})
	operations = cigar[1];
	op_counts = cigar[2];
	indel_idx = findall(map(i -> i == OP_INSERT || i == OP_DELETE, operations));
	match_idx = findall(map(i -> i == OP_MATCH, operations));
	gaps = sum(op_counts[indel_idx]);
	cgaps = length(indel_idx);
	matches = sum(op_counts[match_idx]);
	(matches, gaps, cgaps)
end

# Calculate read identities from BAM 
# (See 'Gap-compressed identity' from https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity)
# In this equation, 'm' means 'M' operations, 'g' means the sum of 'i' (insertion) and 'd' (deletion) operations and
# 'o' means the number of gap opens. 
calculate_identity(n::Int64, m::Int64, g::Int64, o::Int64) = round((1.0 - (n-g+o)/(m+o)) * 100, digits = 3)


# Calculate sum of total read lengths
totalLen(lengths::Array{Int64,1}) = sum(lengths)

# Calculate N50
function calculate_N50(len_records::Array{Int64,1}, total_length::Int64)::Int64
	recordnum = length(len_records);
	tmpLenSum = 0;
	N50 = 0;
	if recordnum == 1
		return len_records[1]
	else
		half_total_length = total_length/2
		@inbounds for i in 1:recordnum
			tmpLenSum += len_records[i]
			if tmpLenSum >= half_total_length
				N50 = len_records[i]
				break
			end
		end
		return N50
	end
end
