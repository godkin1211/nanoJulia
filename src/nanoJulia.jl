module nanoJulia
using Statistics
export calculate_mean_readqual

phred2prob(q::UInt8) = 10^(q/(-10))

prob2phred(p::Float64) = -10log10(p)

function calculate_mean_readqual(quals::Array{UInt8, 1})::Float64
	quals .|> phred2prob |> mean |> prob2phred |> q->round(q, digits = 1)
end

function remove_bad_read!(records::Array{FASTX.FASTQ.Record, 1})
	seqlengths = FASTQ.seqlen.(records)
	quallengths = map(i -> length(FASTQ.quality(i)), records)
	return records[seqlengths .==  quallengths]
end

end
