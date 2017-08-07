function volume = fixCompartmentSize(volume, dimension)
% Standardize the compartment size as cell array with two columns, where
% the first column is strings and the second column is numeric

if dimension == 0
    assert(isnumeric(volume) && isscalar(volume) && volume == 1, 'KroneckerBio:Compartment:NonUnityZeroCompartment', 'A compartment with dimension 0 must have a size equal to 1.')
end

if isnumeric(volume) && isscalar(volume)
    assert(isscalar(volume) && volume > 0, 'KroneckerBio:Compartment:NonpositiveSize', 'Scalar compartment size must be positive')
    volume = {'', volume};
elseif iscell(volume)
    assert(ndims(volume) == 2 && size(volume,2) == 2, 'KroneckerBio:Compartment:BadCellSize', 'Cell array compartment size must be two-dimensional and have exactly two columns')
    assert(iscellstr(volume(:,1)) && all(cellfun(@(x)(isnumeric(x) && isscalar(x) && x > 0), volume(:,2))), 'KroneckerBio:Compartment:BadCell', 'Rows of cell array compartment size must have a string and a positive scalar')
else
    error('KroneckerBio:Compartment:InvalidSize', 'Not a valid type for a compartment size')
end
