cacheHash = {}

def getValue(key):
	if cacheHash.has_key(key):
		return cacheHash[key]
	return None	

def putValue(key, val):
	cacheHash[key] = val
