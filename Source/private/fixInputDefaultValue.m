function default = fixInputDefaultValue(default)
assert(isscalar(default) && isnumeric(default) && default >= 0, 'KroneckBio:Input:default', 'The default input value must be a nonnegative scalar')
