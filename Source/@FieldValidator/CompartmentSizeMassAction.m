function volume = CompartmentSizeMassAction(volume, dimension)
% Standardize the compartment size

if dimension == 0
    assert(isnumeric(volume) && isscalar(volume) && volume == 1, 'KroneckerBio:Compartment:NonUnityZeroCompartment', 'A compartment with dimension 0 must have a size equal to 1.')
end

if isnumeric(volume) && isscalar(volume)
    assert(isscalar(volume) && volume > 0, 'KroneckerBio:Compartment:NonpositiveSize', 'Scalar compartment size must be positive')
elseif iscell(volume)
    assert(ismatrix(volume) && size(volume,2) == 2, 'KroneckerBio:Compartment:BadCellSize', 'Cell array compartment size must be two-dimensional and have exactly two columns')
    assert(iscellstr(volume(:,1)) && all(cellfun(@(x)(isnumeric(x) && isscalar(x) && x > 0), volume(:,2))), 'KroneckerBio:Compartment:BadCell', 'Rows of cell array compartment size must have a string and a positive scalar')
elseif ishandle(volume)
    assert(nargin(volume) == 3, 'KroneckerBio:Compartment:SizeArgumentCount', 'Function handle compartment size must accept three arguments')
else
    error('KroneckerBio:Compartment:InvalidSize', 'Not a valid type for a compartment size')
end