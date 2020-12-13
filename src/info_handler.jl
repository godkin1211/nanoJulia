function extract_info end

# Extract info from FastqInfo
function extract_info(readsinfo::Array{Union{FastqInfo,Missing},1})
	quals = map(e -> round(e.quality, digits=1), readsinfo)
	lengths = map(e -> e.length, readsinfo)
	gcs = map(e -> e.gc, readsinfo)
	DataFrame(quality = quals, length = lengths, gc_content = gcs)
end

# Extract info from BAMInfo
function extract_info(readsinfo::Array{BAMInfo,1})
	quals = map(e -> round(e.quality, digits=1), readsinfo)
	lengths = map(e -> e.length, readsinfo)
	identities = map(e -> round(e.identity, digits=1), readsinfo)
	gcs = map(e -> e.gc, readsinfo)
	DataFrame(quality = quals, length = lengths, identity = identities, gc_content = gcs)
end