### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 7dbd82c8-3ae8-11eb-3c06-9f7be3a109a7
begin
	push!(LOAD_PATH, "/Home/godkin/Projects/nanoJulia/src")
	using nanoJulia, FASTX, PlutoUI, DataFrames, CSV, Plots, Statistics, Markdown, InteractiveUtils, Images
end

# ╔═╡ 0ee8b936-3ae7-11eb-2b11-e5a4b8701e09
md"# nanoJulia Tutorial: use fastq-files-processing as an example"

# ╔═╡ 3a5dfc64-3b5c-11eb-3c09-bf19602fa0fb
md"""Copyright: (c) 2020 Health GeneTech Co."""

# ╔═╡ 630678da-3b5c-11eb-2582-d9f4e971900d
md"Author: Michael Nostalgie <nostalgie.chiu@genebook.com.tw>"

# ╔═╡ 5703c6f6-3af2-11eb-1e24-33ae52fe8d26
md"""## 1. Why we do rebuild a new wheel?
### 1.1 You need a new QC tool when shit happens
```
2020-07-12 11:53:42,152 Records in Fastq files should start with '@' character
concurrent.futures.process._RemoteTraceback:
\"\"\"
Traceback (most recent call last):
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/concurrent/futures/process.py", line 175, in _process_worker
    r = call_item.fn(*call_item.args, **call_item.kwargs)
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/concurrent/futures/process.py", line 153, in _process_chunk
    return [fn(*args) for args in chunk]
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/concurrent/futures/process.py", line 153, in <listcomp>
    return [fn(*args) for args in chunk]
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/site-packages/nanoget/extraction_functions.py", line 322, in process_fastq_plain
    data=[res for res in extract_from_fastq(inputfastq) if res],
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/site-packages/nanoget/extraction_functions.py", line 322, in <listcomp>
    data=[res for res in extract_from_fastq(inputfastq) if res],
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/site-packages/nanoget/extraction_functions.py", line 332, in extract_from_fastq
    for rec in SeqIO.parse(fq, "fastq"):
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/site-packages/Bio/SeqIO/QualityIO.py", line 1055, in FastqPhredIterator
    for title_line, seq_string, quality_string in FastqGeneralIterator(handle):
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/site-packages/Bio/SeqIO/QualityIO.py", line 927, in FastqGeneralIterator
    "Records in Fastq files should start with '@' character"
ValueError: Records in Fastq files should start with '@' character
\"\"\"

The above exception was the direct cause of the following exception:

Traceback (most recent call last):
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/site-packages/nanoplot/NanoPlot.py", line 63, in main
    keep_supp=not(args.no_supplementary))
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/site-packages/nanoget/nanoget.py", line 92, in get_input
    dfs=[out for out in executor.map(extraction_function, files)],
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/site-packages/nanoget/nanoget.py", line 92, in <listcomp>
    dfs=[out for out in executor.map(extraction_function, files)],
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/concurrent/futures/process.py", line 366, in _chain_from_iterable_of_lists
    for element in iterable:
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/concurrent/futures/_base.py", line 586, in result_iterator
    yield fs.pop().result()
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/concurrent/futures/_base.py", line 432, in result
    return self.__get_result()
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/concurrent/futures/_base.py", line 384, in __get_result
    raise self._exception
ValueError: Records in Fastq files should start with '@' character
```

Or you probably have saw such type of error

```
If you read this then NanoPlot 1.30.1 has crashed :-(
Please try updating NanoPlot and see if that helps...

If not, please report this issue at https://github.com/wdecoster/NanoPlot/issues
If you could include the log file that would be really helpful.
Thanks!



concurrent.futures.process._RemoteTraceback:
\"\"\"
Traceback (most recent call last):
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/concurrent/futures/process.py", line 175, in _process_worker
    r = call_item.fn(*call_item.args, **call_item.kwargs)
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/concurrent/futures/process.py", line 153, in _process_chunk
    return [fn(*args) for args in chunk]
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/concurrent/futures/process.py", line 153, in <listcomp>
    return [fn(*args) for args in chunk]
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/site-packages/nanoget/extraction_functions.py", line 322, in process_fastq_plain
    data=[res for res in extract_from_fastq(inputfastq) if res],
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/site-packages/nanoget/extraction_functions.py", line 322, in <listcomp>
    data=[res for res in extract_from_fastq(inputfastq) if res],
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/site-packages/nanoget/extraction_functions.py", line 332, in extract_from_fastq
    for rec in SeqIO.parse(fq, "fastq"):
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/site-packages/Bio/SeqIO/QualityIO.py", line 1055, in FastqPhredIterator
    for title_line, seq_string, quality_string in FastqGeneralIterator(handle):
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/site-packages/Bio/SeqIO/QualityIO.py", line 974, in FastqGeneralIterator
    % (title_line, seq_len, len(quality_string))
ValueError: Lengths of sequence and quality values differs for 6bc232d8-2605-44ac-bb65-068712e9d9c0 runid=c31ae3e6906d7dfb719d9ee4903bca77410329f4 read=2020 ch=480 start_time=2020-07-14T09:31:36Z flow_cell_id=FAN38431 protocol_group_id=20200714_RD_WGA10 sample_id=no_sample (411 and 1679).
\"\"\"

The above exception was the direct cause of the following exception:

Traceback (most recent call last):
  File "/home/analysis/nanopore/miniconda3/envs/py3/bin/NanoPlot", line 8, in <module>
    sys.exit(main())
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/site-packages/nanoplot/NanoPlot.py", line 63, in main
    keep_supp=not(args.no_supplementary))
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/site-packages/nanoget/nanoget.py", line 92, in get_input
    dfs=[out for out in executor.map(extraction_function, files)],
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/site-packages/nanoget/nanoget.py", line 92, in <listcomp>
    dfs=[out for out in executor.map(extraction_function, files)],
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/concurrent/futures/process.py", line 366, in _chain_from_iterable_of_lists
    for element in iterable:
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/concurrent/futures/_base.py", line 586, in result_iterator
    yield fs.pop().result()
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/concurrent/futures/_base.py", line 432, in result
    return self.__get_result()
  File "/home/analysis/nanopore/miniconda3/envs/py3/lib/python3.6/concurrent/futures/_base.py", line 384, in __get_result
    raise self._exception
ValueError: Lengths of sequence and quality values differs for 6bc232d8-2605-44ac-bb65-068712e9d9c0 runid=c31ae3e6906d7dfb719d9ee4903bca77410329f4 read=2020 ch=480 start_time=2020-07-14T09:31:36Z flow_cell_id=FAN38431 protocol_group_id=20200714_RD_WGA10 sample_id=no_sample (411 and 1679).
```

- - -

### 1.2 Most of QC tools do not support Fast5-format files
![minionqc](https://www.biorxiv.org/content/biorxiv/early/2019/10/28/818336/F2.large.jpg)

- - -

### 1.3 Why Julia?
  - Write like Python, run like C
  - Easys parallelisation, high performance
  - Static/Dynamics typed
  - Composable
"""

# ╔═╡ 7dcd7d2c-3ae8-11eb-20e0-9d631b456439
md"## 2. Install and load nanoJulia"

# ╔═╡ 5590e294-3ae9-11eb-15bc-7b0ac2a6b03c
md"## 3. Setting input files"

# ╔═╡ 80d4574e-3ae7-11eb-1b6e-27a843ecf98a
begin
	fastq_files_path = "/home/godkin/Projects/test_fastq"
	fastq_files = joinpath.(fastq_files_path,readdir(fastq_files_path))
end

# ╔═╡ 78b7885a-3ae8-11eb-3e57-e7f86097ba73
md"## 4. Read and parse fastq files"

# ╔═╡ 00250c30-3aea-11eb-25c1-fbc18c3b012f
begin
	output_table = mapreduce(f -> readFastq(f), vcat, fastq_files)
	first(output_table, 10)
end

# ╔═╡ f494c894-3aea-11eb-09d2-473d6f516e8d
md"## 5. Write this table into a CSV-format file"

# ╔═╡ df200b08-3aeb-11eb-2b79-61f1aeb0b5d0
begin
	outdir = "/tmp/QC_test"
	if isdir(outdir)
		rm(outdir, force=true, recursive=true)
	end
	mkdir(outdir)
	outputfile = joinpath(outdir, "quantities.csv")
	output_table |> CSV.write(outputfile)
end

# ╔═╡ 042642a6-3aed-11eb-191a-17796e1e69f8
md"## 6. Generating statistics summaries"

# ╔═╡ 213acc4a-3aed-11eb-331d-1b9a6c4cb004
begin
	totalbases = nanoJulia.totalLen(output_table.length)
	n50 = nanoJulia.calculate_N50(output_table.length, totalbases)
	generateStatSummary(output_table, totalbases, n50, outdir)
end

# ╔═╡ 7f0c927c-3aed-11eb-2607-3dd5841ceebb
md"## 7. Plot outputs
### 7.1 Read length vs read quality 2D histogram"

# ╔═╡ 90291f5e-3b55-11eb-02db-8be99e8f78ee
begin
	binNum1 = @bind nbins_1 Slider(20:10:200)
	md"Number of bins: $(binNum1)"
end

# ╔═╡ be5284a0-3aed-11eb-114e-a129fd7542d5
begin
	gr()
	plot(output_table.length,
		output_table.quality,
		seriestype = :histogram2d,
		nbins = nbins_1,
		title = "Read Length vs Read Quality",
		xlabel = "Read length (bp)",
		ylabel = "Phred Score")
end

# ╔═╡ 718873e6-3af0-11eb-3d66-1f88abde9cdc
md"""- - -
### 7.2 Read Length Distribution"""

# ╔═╡ 9b4a1072-3b56-11eb-1e05-f14374fc27b8
begin
	binNum2 = @bind nbins_2 Slider(20:10:200)
	md"Number of bins: $(binNum2)"
end

# ╔═╡ 7ecf441e-3af0-11eb-1f8c-5795a6930048
begin
	gr()
	plot(output_table.length,
		seriestype = :histogram,
		bins = nbins_2,
		minorticks = true,
		xlabel = "Read Length (bp)",
		ylabel = "Counts",
		title = "Read Length Distribution")
	vline!([n50], lw=2, lab="Read N50")
end

# ╔═╡ 92cfb52c-3aee-11eb-23ed-7176e50f1b33
md"""- - -
### 7.3 GC-content distribution"""

# ╔═╡ e9c499c8-3b56-11eb-3e5a-5b887f90e73d
begin
	binNum3 = @bind nbins_3 Slider(20:10:500)
	md"Number of bins: $(binNum3)"
end

# ╔═╡ c85b8216-3aee-11eb-3b00-2f19c6b981c6
begin
	gr()
	gc_contents = output_table.gc_content .* 100
	mean_gc = mean(gc_contents)
	plot(gc_contents,
		seriestype = :histogram,
		bins = nbins_3,
		minorticks = true,
		xlabel = "GC Content (%)",
		ylabel = "Counts",
		title = "GC-content distribution")
	vline!([mean_gc], lw=2, lab="Avg. GC content")
end

# ╔═╡ 6c4282ce-3b5a-11eb-1dfc-a9159b5a16ad
md"## 8. Processing Fast5 files"

# ╔═╡ 6a55b2f0-3cf5-11eb-1218-e16c1df49aab
?readFast5

# ╔═╡ 851efb7e-3b5a-11eb-0c61-49555c0f91bf
md"""```
begin
	fast5_files_path = "/Users/godkin/Projects/test_fast5/fast5_fail"
	fast5_files = joinpath.(fast5_files_path,readdir(fast5_files_path))
	output_table = mapreduce(f -> readFast5(f), vcat, fast5_files)
	first(output_table, 10)
end
```"""

# ╔═╡ dca44304-3b5a-11eb-2b1b-c9b491518fc3
md"## 9. Processing BAM file"

# ╔═╡ 991ceacc-3cf5-11eb-32cd-bf18f3d34e59
?readBAM

# ╔═╡ e8777660-3b5a-11eb-2585-1747d497d802
md"""```
begin
	bam_file_path = "/path/to/your/bamfile.bam"
	output_table = readBAM(bam_file_path)
	first(output_table, 10)
end
```"""

# ╔═╡ 29d19648-3b5b-11eb-038a-3f5216225a12
md"## 10. Processing SequencingSummary file"

# ╔═╡ b220666e-3cf5-11eb-0731-7ba6a69176d8
?readSeqSummary

# ╔═╡ 3baf4808-3b5b-11eb-0332-cd4191c64a26
md"""```
begin
	seqSummary_file_path = "/path/to/your/seqSummary.txt"
	output_table = readSeqSummary(seqSummary_file_path)
	first(output_table, 10)
end
```"""

# ╔═╡ ab40ff42-3b59-11eb-2502-b1a512615492
md"## 11. Speed up these processes"

# ╔═╡ c5301f64-3b59-11eb-1c54-37237f9ee35a
md"""```
# In Julia Console
> using Pkg
> Pkg.add("PackageCompiler")
> Using PackageCompiler, nanoJulia
> create_sysimage(:nanoJulia; sysimage_path="nanoJulia.so")

# In Terminal
$ julia -JnanoJulia.so processFast5.jl -h 
```"""

# ╔═╡ Cell order:
# ╟─0ee8b936-3ae7-11eb-2b11-e5a4b8701e09
# ╟─3a5dfc64-3b5c-11eb-3c09-bf19602fa0fb
# ╟─630678da-3b5c-11eb-2582-d9f4e971900d
# ╟─5703c6f6-3af2-11eb-1e24-33ae52fe8d26
# ╟─7dcd7d2c-3ae8-11eb-20e0-9d631b456439
# ╠═7dbd82c8-3ae8-11eb-3c06-9f7be3a109a7
# ╟─5590e294-3ae9-11eb-15bc-7b0ac2a6b03c
# ╠═80d4574e-3ae7-11eb-1b6e-27a843ecf98a
# ╟─78b7885a-3ae8-11eb-3e57-e7f86097ba73
# ╠═00250c30-3aea-11eb-25c1-fbc18c3b012f
# ╟─f494c894-3aea-11eb-09d2-473d6f516e8d
# ╠═df200b08-3aeb-11eb-2b79-61f1aeb0b5d0
# ╟─042642a6-3aed-11eb-191a-17796e1e69f8
# ╠═213acc4a-3aed-11eb-331d-1b9a6c4cb004
# ╟─7f0c927c-3aed-11eb-2607-3dd5841ceebb
# ╟─90291f5e-3b55-11eb-02db-8be99e8f78ee
# ╠═be5284a0-3aed-11eb-114e-a129fd7542d5
# ╟─718873e6-3af0-11eb-3d66-1f88abde9cdc
# ╟─9b4a1072-3b56-11eb-1e05-f14374fc27b8
# ╠═7ecf441e-3af0-11eb-1f8c-5795a6930048
# ╟─92cfb52c-3aee-11eb-23ed-7176e50f1b33
# ╟─e9c499c8-3b56-11eb-3e5a-5b887f90e73d
# ╠═c85b8216-3aee-11eb-3b00-2f19c6b981c6
# ╟─6c4282ce-3b5a-11eb-1dfc-a9159b5a16ad
# ╠═6a55b2f0-3cf5-11eb-1218-e16c1df49aab
# ╟─851efb7e-3b5a-11eb-0c61-49555c0f91bf
# ╟─dca44304-3b5a-11eb-2b1b-c9b491518fc3
# ╠═991ceacc-3cf5-11eb-32cd-bf18f3d34e59
# ╠═e8777660-3b5a-11eb-2585-1747d497d802
# ╟─29d19648-3b5b-11eb-038a-3f5216225a12
# ╠═b220666e-3cf5-11eb-0731-7ba6a69176d8
# ╟─3baf4808-3b5b-11eb-0332-cd4191c64a26
# ╟─ab40ff42-3b59-11eb-2502-b1a512615492
# ╟─c5301f64-3b59-11eb-1c54-37237f9ee35a
