function size = fixCompartmentSize(size)
% Standardize the compartment size as a positive scalar
assert(isscalar(size) && size > 0, 'KroneckerBio:Compartment:Size', 'Compartment size must be a positive scalar')
