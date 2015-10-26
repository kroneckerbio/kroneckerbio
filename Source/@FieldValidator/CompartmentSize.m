function size = CompartmentSize(size, dimension)
% Standard compartment size as a string

if ischar(size)
    % Seems legit
elseif isnumeric(size) && isscalar(size)
    if dimension == 0 && size ~= 1
        error('KroneckerBio:FieldValidator:CompartmentSize:NonUnityZeroCompartment', 'A compartment with dimension 0 must have a size equal to 1.')
    end
    size = num2str(size);
else
    error('KroneckerBio:FieldValidator:CompartmentSize:Invalid: Invalid compartment size')
end
