function size = fixCompartmentSizeAnalytic(size, dimension)

if isnumeric(size) && isscalar(size)
    if dimension == 0 && size ~= 1
        error('KroneckerBio:Compartment:NonUnityZeroCompartment', 'A compartment with dimension 0 must have a size equal to 1.')
    end
else
    error('fixCompartmentSizeAnalytic: invalid compartment size')
end