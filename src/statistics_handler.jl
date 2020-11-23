# Subsetting info-table with quality
getQScorePart(df::DataFrame, qscore::Int64) = df[df.quality .>= qscore, :]

function plotSquiggle(signals::Array{Int16,1})
	gr()
	timepoints = [i for i in 1:length(signals)] ./ 4000
	plot(timepoints, signals, w=3, xlabel = "Time (seconds)", ylabel= "Detected Current (pA)", legend = false)
end

function plotReadLen2QualScatter(df::DataFrames.DataFrame, N50::Int64)
	gr()
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
	savefig("read_length_vs_quality_dot.png")
end

function plotReadLen2QualHistogram2D(df::DataFrames.DataFrame)
	gr()
	historgam2d(df.length, df.quality, nbins = 50)
	savefig("read_length_vs_quality_histogram.png")
end

function generateStatSummary(df::DataFrames.DataFrame, totalLength::Int64, N50::Int64)::DataFrames.DataFrame
	stat_summary = describe(df)
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
	push!(output_df, ("Total Bases", totalLength))
	meanQualtxt = @sprintf "Mean Read Quality: %20.1f" meanQual
	meanLentxt = @sprintf "Mean Read Length: %21.1f" meanLen
	medianQualtxt = @sprintf "Median Read Quality: %18.1f" medianQual
	medianLentxt = @sprintf "Median Read Length: %19.1f" medianLen
	readN50txt = @sprintf "Read N50: %29s" format(N50, commas=true)
	totalBasestxt = @sprintf "Total Bases: %26s" format(totalLength, commas=true)
	if ncol(df) == 3
		open("statistics_summary.txt", "w") do io
			write(io, "$meanQualtxt\n$meanLentxt\n$meanIdenttxt\n$medianQualtxt\n$medianLentxt\n$medianIdenttxt\n$readN50txt\n$totalBasestxt")
		end
	else
		open("statistics_summary.txt", "w") do io
			write(io, "$meanQualtxt\n$meanLentxt\n$medianQualtxt\n$medianLentxt\n$readN50txt\n$totalBasestxt")
		end
	end
	return output_df
end
