function size = fixCompartmentSizeAnalytic(size, dimension)
% Standard compartment size as a string

if ischar(size)
    warnAboutQuotedFullName(size)
elseif isnumeric(size) && isscalar(size)
    if dimension == 0 && size ~= 1
        error('KroneckerBio:Compartment:NonUnityZeroCompartment', 'A compartment with dimension 0 must have a size equal to 1.')
    end
    size = num2str(size);
else
    error('KroneckerBio:Compartment:size', 'Invalid compartment size')
end
