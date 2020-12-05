# Subsetting info-table with quality
function getQScoreOverQ(df::DataFrame, qscore::Int64, totallength::Int64) 
	df_subset = df[df.quality .>= qscore, :]
	readNum = size(df_subset, 1)
	baseNum = totalLen(df_subset.length)
	basePerc = round(baseNum * 100 / totallength, digits=2)
	return(readNum, baseNum, "$basePerc%")
end

# Plotting functions
function plotSquiggle(signals::Array{Int16,1})
	gr()
	timepoints = [i for i in 1:length(signals)] ./ 4000
	plot(timepoints, signals, w=3, xlabel = "Time (seconds)", ylabel= "Detected Current (pA)", legend = false)
end

function plotReadQual2IdentScatter(df::DataFrames.DataFrame, outputDir::String)
	gr()
	outputfile = joinpath(outputdir, "read_quality_vs_identity_dot.png")
	plot(df.quality,
		df.identity,
		seriestype = :scatter,
		title = "Read Quality vs Identity",
		xlabel = "Phred Score",
		ylabel = "Identity (%)",
		dpi = 300)
	savefig(outputfile)
end

function plotReadLen2QualScatter(df::DataFrames.DataFrame, N50::Int64, outputDir::String)
	gr()
	outputfile = joinpath(outputDir, "read_length_vs_quality_dot.png")
	if "identity" in names(df)
		plot(df.length,
			df.quality,
			zcolor = df.identity,
			seriestype = :scatter,
			title = "Read Length vs Quality",
			xlabel = "Read Length (bp)",
			ylabel = "Phred Score",
			lab = "Identity (%)",
			dpi = 300)
	else
		plot(df.length, 
			df.quality, 
			seriestype = :scatter, 
			title = "Read Length vs Quality", 
			legend = false, 
			xlabel = "Read Length (bp)", 
			ylabel = "Phred Score", 
			dpi = 300)
	end
	vline!([N50], lw = 2, lab = "N50")
	savefig(outputfile)
end

function plotReadLen2QualHistogram2D(df::DataFrames.DataFrame, outputDir::String)
	gr()
	outputfile = joinpath(outputDir, "read_length_vs_quality_histogram.png")
	histogram2d(df.length, df.quality, nbins = 50)
	savefig(outputfile)
end

function plotReadQual2IdentHistogram2D(df::DataFrames.DataFrame, outputDir::String)
	gr()
	outputfile = joinpath(outputDir, "read_quality_vs_identity_histogram.png")
	histogram2d(df.quality, df.identity, nbins = 50)
	savefig(outputfile)
end

function plotReadLenDist(lengthRecords::Array{Int64,1}, N50::Int64, outputDir::String)
	gr()
	outputfile = joinpath(outputDir, "read_length_distribution_histogram.png")
	plot(lengthRecords, 
		 seriestype=:histogram, 
		 bins=400, 
		 minorticks=true, 
		 legend=false, 
		 xlabel = "Read Length (bp)",
		 ylabel = "Read Counts",
		 dpi = 300)
	vline!([N50], lw=2, lab="N50")
	savefig(outputfile)
end


# Generating statics summary files
function generateStatSummary(df::DataFrames.DataFrame, totalLength::Int64, N50::Int64, outputDir::String)::DataFrames.DataFrame
	stat_summary = describe(df)
	readNum = size(df, 1)
	readsQ7, basesQ7, percQ7 = getQScoreOverQ(df, 7, totalLength)
	readsQ10, basesQ10, percQ10 = getQScoreOverQ(df, 10, totalLength)
	readsQ12, basesQ12, percQ12 = getQScoreOverQ(df, 12, totalLength)
	readsQ15, basesQ15, percQ15 = getQScoreOverQ(df, 15, totalLength)
	output_df = DataFrame(Property = AbstractString[], Value = Number[])
	if ncol(df) == 3
		meanQual, meanLen, meanIdent = tuple(round.(stat_summary[!,:mean], digits=1)...)
		medianQual, medianLen, medianIdent = tuple(stat_summary[!,:median]...)
		push!(output_df, ("Average Identity", meanIdent))
		push!(output_df, ("Median Identity", medianIdent))
		meanIdenttxt = @sprintf "Mean Identity: %24.1f" meanIdent
		medianIdenttxt = @sprintf "Median Identity: %22.1f" medianIdent
	else
		meanQual, meanLen = tuple(round.(stat_summary[!,:mean], digits=1)...)
		medianQual, medianLen = tuple(stat_summary[!,:median]...)
	end
	push!(output_df, ("Mean Read Length", meanLen))
	push!(output_df, ("Median Read Length", medianLen))
	push!(output_df, ("Mean Read Quality", meanQual))
	push!(output_df, ("Median Read Quality", medianQual))
	push!(output_df, ("Read N50", N50))
	push!(output_df, ("Read Numbers", readNum))
	push!(output_df, ("Total Bases", totalLength))

	meanQualtxt = @sprintf "Mean Read Quality: %27.1f" meanQual
	meanLentxt = @sprintf "Mean Read Length: %28.1f" meanLen
	medianQualtxt = @sprintf "Median Read Quality: %25.1f" medianQual
	medianLentxt = @sprintf "Median Read Length: %26.1f" medianLen
	readN50txt = @sprintf "Read N50: %36s" format(N50, commas=true)
	readNumtxt = @sprintf "Read Number: %33s" format(readNum, commas=true)
	totalBasestxt = @sprintf "Total Bases: %33s" format(totalLength, commas=true)
	readsQ7txt = @sprintf "Reads with quality score > 7: %16s" format(readsQ7, commas=true)
	basesQ7txt = @sprintf "Bases with quality score > 7: %16s" format(basesQ7, commas=true)
	percQ7txt = @sprintf "Percentage of bases (Q>7): %19s" percQ7
	readsQ10txt = @sprintf "Reads with quality score > 10: %15s" format(readsQ10, commas=true)
	basesQ10txt = @sprintf "Bases with quality score > 10: %15s" format(basesQ10, commas=true)
	percQ10txt = @sprintf "Percentage of bases (Q>10): %18s" percQ10
	readsQ12txt = @sprintf "Reads with quality score > 12: %15s" format(readsQ12, commas=true)
	basesQ12txt = @sprintf "Bases with quality score > 12: %15s" format(basesQ12, commas=true)
	percQ12txt = @sprintf "Percentage of bases (Q>12): %18s" percQ12
	readsQ15txt = @sprintf "Reads with quality score > 15: %15s" format(readsQ15, commas=true)
	basesQ15txt = @sprintf "Bases with quality score > 15: %15s" format(basesQ15, commas=true)
	percQ15txt = @sprintf "Percentage of bases (Q>15): %18s" percQ15

	outputfile = joinpath(outputDir, "statistics_summary.txt")
	open(outputfile, "w") do io
		if ncol(df) == 3
			write(io, "$meanQualtxt\n$meanLentxt\n$meanIdenttxt\n$medianQualtxt\n$medianLentxt\n$medianIdenttxt\n$readN50txt\n$readNumtxt\n$totalBasestxt\n")
		else
			write(io, "$meanQualtxt\n$meanLentxt\n$medianQualtxt\n$medianLentxt\n$readN50txt\n$readNumtxt\n$totalBasestxt\n")
		end
		write(io, "$readsQ7txt\n$basesQ7txt\n$percQ7txt\n")
		write(io, "$readsQ10txt\n$basesQ10txt\n$percQ10txt\n")
		write(io, "$readsQ12txt\n$basesQ12txt\n$percQ12txt\n")
		write(io, "$readsQ15txt\n$basesQ15txt\n$percQ15txt\n")
	end
	
	plotReadLen2QualScatter(df, N50, outputDir)
	plotReadLen2QualHistogram2D(df, outputDir)
	plotReadLenDist(df.length, N50, outputDir)

	if ncol(df) == 3
		plotReadQual2IdentScatter(df, outputDir)
		plotReadQual2IdentHistogram2D(df, outputDir)
	end
	
	return output_df
end
