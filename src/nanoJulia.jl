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


"""
    readFast5(f5file, basecall_group)

Read and parse a fast5 file to extract basecalled fastq information. 
The extracted fastq info will be parsed, then the required information
of each read will be extracted, calculated, and exported as a data frame.
The `basecall_group` could be determined simply by using `h5dump` on a 
fast5 file. For example, you can type `h5dump -n single_read.fast5` in 
terminal, and the content of this fast5 will be shown like this:
```
HDF5 "fff07ce7-f378-47ce-a71a-488d9a1a9af8.fast5" {
FILE_CONTENTS {
 group      /
 group      /Analyses
 group      /Analyses/Basecall_1D_000
 group      /Analyses/Basecall_1D_000/Summary
 group      /Analyses/Basecall_1D_000/Summary/basecall_1d_template
 group      /Analyses/Basecall_1D_001
 group      /Analyses/Basecall_1D_001/BaseCalled_template
 dataset    /Analyses/Basecall_1D_001/BaseCalled_template/Fastq
 dataset    /Analyses/Basecall_1D_001/BaseCalled_template/Move
 dataset    /Analyses/Basecall_1D_001/BaseCalled_template/Trace
 group      /Analyses/Basecall_1D_001/Summary
 group      /Analyses/Basecall_1D_001/Summary/basecall_1d_template
 group      /Analyses/Segmentation_000
 group      /Analyses/Segmentation_000/Summary
 group      /Analyses/Segmentation_000/Summary/segmentation
 group      /Analyses/Segmentation_001
 group      /Analyses/Segmentation_001/Summary
 group      /Analyses/Segmentation_001/Summary/segmentation
 group      /Raw
 group      /Raw/Reads
 group      /Raw/Reads/Read_1950
 dataset    /Raw/Reads/Read_1950/Signal
 group      /UniqueGlobalKey
 group      /UniqueGlobalKey/channel_id
 group      /UniqueGlobalKey/context_tags
 group      /UniqueGlobalKey/tracking_id
 }
}
```
You can observe this output, and notice that `Fastq` locates in the directory 
`/Analyses/Basecall_1D_001/BaseCalled_template`. Therefore, the `basecall_group`
sould be set as `Basecall_1D_001` not `Basecall_1D_000`.

# Examples:
1. Read a single-read fast5.
```
julia> readFast5("fff07ce7-f378-47ce-a71a-488d9a1a9af8.fast5", "Basecall_1D_001")
1×3 DataFrame
 Row │ quality  length  gc_content 
     │ Float64  Int64   Float64    
─────┼─────────────────────────────
   1 │ 5.05611    1326      0.4759
```

2. Read a multi-reads fast5.
julia> readFast5("FAM93744_fail_6a4bbe82_1.fast5", "Basecall_1D_000")
4000×3 DataFrame
  Row │ quality  length  gc_content 
      │ Float64  Int64   Float64    
──────┼─────────────────────────────
    1 │     3.7     255      0.4314
    2 │     4.5     331      0.3897
  ⋮   │    ⋮       ⋮         ⋮
 3999 │     6.6     849      0.3604
 4000 │     4.2     164      0.5915
                   3996 rows omitted
"""
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