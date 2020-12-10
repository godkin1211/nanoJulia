### A Pluto.jl notebook ###
# v0.12.16

using Markdown
using InteractiveUtils

# ╔═╡ 7dbd82c8-3ae8-11eb-3c06-9f7be3a109a7
begin
	push!(LOAD_PATH, "/home/godkin/Projects/nanoJulia/src")
	using nanoJulia, FASTX, PlutoUI, DataFrames, CSV, Plots, Statistics, Markdown, InteractiveUtils, Images
end

# ╔═╡ 0ee8b936-3ae7-11eb-2b11-e5a4b8701e09
md"# Tutorial: process Fastq files"

# ╔═╡ 5703c6f6-3af2-11eb-1e24-33ae52fe8d26
md"""## 0. Why we do rebuild a new wheel?
### 0.1 You need a new QC tool when shit happens
![error](https://user-images.githubusercontent.com/59709147/72103466-92ee7c00-3329-11ea-91e2-4717f559328c.png)
### 0.2 Mosts of QC tools do not support Fast5 files
![minionqc](https://www.biorxiv.org/content/biorxiv/early/2019/10/28/818336/F2.large.jpg)
### 0.3 Julia Programming Language"""

# ╔═╡ 7dcd7d2c-3ae8-11eb-20e0-9d631b456439
md"## 1. Load nanoJulia"

# ╔═╡ 5590e294-3ae9-11eb-15bc-7b0ac2a6b03c
md"## 2. Setting input files"

# ╔═╡ 80d4574e-3ae7-11eb-1b6e-27a843ecf98a
begin
	fastq_files_path = "/home/godkin/Projects/test_fastq"
	fastq_files = joinpath.(fastq_files_path,readdir(fastq_files_path))
end

# ╔═╡ 78b7885a-3ae8-11eb-3e57-e7f86097ba73
md"## 3. Read and parse fastq files"

# ╔═╡ 00250c30-3aea-11eb-25c1-fbc18c3b012f
begin
	output_table = mapreduce(f -> readFastq(f), vcat, fastq_files);
	first(output_table, 10)
end

# ╔═╡ f494c894-3aea-11eb-09d2-473d6f516e8d
md"## 4. Write this table into a CSV-format file"

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
md"## 5. Generating statistics summary and plots"

# ╔═╡ 213acc4a-3aed-11eb-331d-1b9a6c4cb004
begin
	totalbases = nanoJulia.totalLen(output_table.length)
	n50 = nanoJulia.calculate_N50(output_table.length, totalbases)
	generateStatSummary(output_table, totalbases, n50, outdir)
end

# ╔═╡ 7f0c927c-3aed-11eb-2607-3dd5841ceebb
md"## 6. Plot outputs
### 6.1 Read length vs read quality 2D histogram"

# ╔═╡ be5284a0-3aed-11eb-114e-a129fd7542d5
begin
	gr()
	plot(output_table.length,
		output_table.quality,
		seriestype = :histogram2d,
		nbins = 50,
		title = "Read Length vs Read Quality",
		legend = false,
		xlabel = "Read length (bp)",
		ylabel = "Phred Score")
end

# ╔═╡ 718873e6-3af0-11eb-3d66-1f88abde9cdc
md"### 6.2 Read Length Distribution"

# ╔═╡ 7ecf441e-3af0-11eb-1f8c-5795a6930048
begin
	gr()
	plot(output_table.length,
		seriestype = :histogram,
		bins = 100,
		minorticks = true,
		xlabel = "Read Length (bp)",
		ylabel = "Counts",
		title = "Read Length Distribution")
	vline!([n50], lw=2, lab="Read N50")
end

# ╔═╡ 92cfb52c-3aee-11eb-23ed-7176e50f1b33
md"### 6.3 GC-content distribution"

# ╔═╡ c85b8216-3aee-11eb-3b00-2f19c6b981c6
begin
	gr()
	gc_contents = output_table.gc_content .* 100
	mean_gc = mean(gc_contents)
	plot(gc_contents,
		seriestype = :histogram,
		bins = 400,
		minorticks = true,
		xlabel = "GC Content (%)",
		ylabel = "Counts",
		title = "GC-content distribution")
	vline!([mean_gc], lw=2, lab="Avg. GC content")
end

# ╔═╡ Cell order:
# ╟─0ee8b936-3ae7-11eb-2b11-e5a4b8701e09
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
# ╠═be5284a0-3aed-11eb-114e-a129fd7542d5
# ╟─718873e6-3af0-11eb-3d66-1f88abde9cdc
# ╠═7ecf441e-3af0-11eb-1f8c-5795a6930048
# ╟─92cfb52c-3aee-11eb-23ed-7176e50f1b33
# ╠═c85b8216-3aee-11eb-3b00-2f19c6b981c6
