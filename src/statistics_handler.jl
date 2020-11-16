# Subsetting info-table with quality
getQScorePart(df::DataFrame, qscore::Int64) = df[df.quality .>= qscore, :]

function generateStatSummary(df::DataFrame)
	stat_summary = describe(df)
	if ncol(df) == 3
		meanQual, meanLen, meanIdent = tuple(round.(stat_summary[!,:mean], digits=1)...)
		medianQual, medianLen, medianIdent = tuple(stat_summary[!,:median]...)
	else
		meanQual, meanLen = tuple(round.(stat_summary[!,:mean], digits=1)...)
		medianQual, medianLen = tuple(stat_summary[!,:median]...)
	end
	totalLength = sum(df.length)
	n50 = calculate_N50(df.length, totalLength)
	meanQualtxt = @sprintf "Mean Read Quality: %20.1f" meanQual
	meanLentxt = @sprintf "Mean Read Length: %21.1f" meanLen
	medianQualtxt = @sprintf "Median Read Quality: %18.1f" medianQual
	medianLentxt = @sprintf "Median Read Length: %19.1f" medianLen
	readN50txt = @sprintf "Read N50: %29d" n50
	totalBasestxt = @sprintf "Total Bases: %26s" format(totalLength, commas=true)
	if ncol(df) == 3
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
