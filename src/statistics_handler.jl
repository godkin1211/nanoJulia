# Subsetting info-table with quality
getQScorePart(df::DataFrame, qscore::Int64) = df[df.quality .>= qscore, :]

totalLen(lengths::Array{Int64,1}) = sum(lengths)

function plotReadLen2Qaul(length_to_quality_df::DataFrames.DataFrame, N50::Int64)
	gr()
	plot(length_to_quality_df.length, 
		length_to_quality_df.quality, 
		seriestype = :scatter, 
		title = "Read Length vs Quality", 
		legend = false, 
		xlabel = "Read Length (bp)", 
		ylabel = "Phred Score", 
		dpi = 300)
	vline!([N50], lw = 2)
	savefig("read_length_vs_quality.png")
end


function generateStatSummary(length_to_quality_df::DataFrames.DataFrame, totalLength::Int64, N50::Int64)
	stat_summary = describe(length_to_quality_df)
	if ncol(length_to_quality_df) == 3
		meanQual, meanLen, meanIdent = tuple(round.(stat_summary[!,:mean], digits=1)...)
		medianQual, medianLen, medianIdent = tuple(stat_summary[!,:median]...)
	else
		meanQual, meanLen = tuple(round.(stat_summary[!,:mean], digits=1)...)
		medianQual, medianLen = tuple(stat_summary[!,:median]...)
	end
	meanQualtxt = @sprintf "Mean Read Quality: %20.1f" meanQual
	meanLentxt = @sprintf "Mean Read Length: %21.1f" meanLen
	medianQualtxt = @sprintf "Median Read Quality: %18.1f" medianQual
	medianLentxt = @sprintf "Median Read Length: %19.1f" medianLen
	readN50txt = @sprintf "Read N50: %29s" format(N50, commas=true)
	totalBasestxt = @sprintf "Total Bases: %26s" format(totalLength, commas=true)
	if ncol(length_to_quality_df) == 3
		meanIdenttxt = @sprintf "Mean Identity: %28.1f" meanIdent
		medianIdenttxt = @sprintf "Median Identity: %20.1f" medianIdent
		open("statistics_summary.txt", "w") do io
			write(io, "$meanQualtxt\n$meanLentxt\n$meanIdenttxt\n$medianQualtxt\n$medianLentxt\n$medianIdenttxt\n$readN50txt\n$totalBasestxt")
		end
	else
		open("statistics_summary.txt", "w") do io
			write(io, "$meanQualtxt\n$meanLentxt\n$medianQualtxt\n$medianLentxt\n$readN50txt\n$totalBasestxt")
		end
	end
end
