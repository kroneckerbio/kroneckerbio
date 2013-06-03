function value = fixParameterValue(value)

assert(isnumeric(value) && isscalar(value) && value >= 0, ...
    'KroneckerBio:Parameter:Value', 'Parameter value must be a nonnegative scalar')
