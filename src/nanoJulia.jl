module nanoJulia
using Statistics, FASTX, DataFrames, Printf, Formatting, BioAlignments, XAM, Plots, HDF5, CSV, BioSequences
export nanoread, generateStatSummary, plotReadLen2QualScatter, 
       plotReadLen2QualHistogram2D, readFast5, readFastq, readBAM,
       plotReadQual2IdentScatter, plotReadQual2IdentHistogram2D,
       plotReadLenDist, plotSquiggle, readSeqSummary, calculate_GC_content

include("datatype.jl")
include("utilities.jl")
include("record_handler.jl")
include("nanoread.jl")
include("info_handler.jl")
include("statistics_handler.jl")

function readFast5(fast5file::String, basecall_group::String)::DataFrames.DataFrame
	h5open(fast5file) do f5
        file_version = read(attrs(f5), "file_version")
        extracted_info = nanoread(f5, file_version, basecall_group)
        return extracted_info
	end
end


function readFastq(fastqfile::String)::DataFrames.DataFrame
    open(FASTQ.Reader, fastqfile) do fastq
        extracted_info = nanoread(fastq)
        return extracted_info
    end
end

function readBAM(bamfile::String)::DataFrames.DataFrame
    open(BAM.Reader, bamfile) do bam
        extracted_info = nanoread(bam)
        return extracted_info
    end
end


function readSeqSummary(seqsummaryfile::String)::DataFrames.DataFrame
    extracted_info = CSV.File(seqsummaryfile) |> DataFrame |> nanoread
    return extracted_info
end

end