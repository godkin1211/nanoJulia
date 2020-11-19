function extract_info end

# Extract info from FastqInfo
function extract_info(readsinfo::Array{Union{FastqInfo,Missing},1})
	quals = map(e -> e.quality, readsinfo);
	lengths = map(e -> e.length, readsinfo);
	DataFrame(quality = quals, length = lengths)
end

# Extract info from BAMInfo
function extract_info(readsinfo::Array{BAMInfo,1})
	quals = map(e -> e.quality, readsinfo);
	lengths = map(e -> e.length, readsinfo);
	identities = map(e -> round(e.identity, digits=1), readsinfo);
	DataFrame(quality = quals, length = lengths, identity = identities)
end