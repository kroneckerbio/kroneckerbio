function value = ParameterValue(value)

assert(isnumeric(value) && isscalar(value) && value >= 0, ...
    'KroneckerBio:FieldValidator:ParameterValue', 'Parameter value must be a nonnegative scalar')
