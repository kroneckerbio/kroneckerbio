function relTol = fixRelTol(relTol)

assert(isnumeric(relTol), 'KroneckerBio:RelTol:Numeric', 'RelTol must be numeric')
assert(numel(relTol) <= 1, 'KroneckerBio:RelTol:Length', 'RelTol cannot be provided as a vector')

if isempty(relTol) || isnan(relTol)
    relTol = 1e-6;
end
