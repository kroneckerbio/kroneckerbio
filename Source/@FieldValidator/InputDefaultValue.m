function default = InputDefaultValue(default)

assert(isscalar(default) && isnumeric(default) && default >= 0, ...
    'KroneckBio:FieldValidator:InputDefaultValue', 'The default input value must be a nonnegative scalar')
