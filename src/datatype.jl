# Readinfo types
# A container for storage of extracted fastq info
struct FastqInfo
	quality::Float64
	length::Int64
	gc::Float64
end

# A container for storage of extracted BAM info
struct BAMInfo
	quality::Float64
	length::Int64
	identity::Float64
	gc::Float64
end