function dimension = CompartmentDimension(dimension)
% Standardize the compartment dimension as 0, 1, 2, or 3

assert(isnumeric(dimension) && isscalar(dimension) && any(dimension == 0:3), ...
    'KroneckerBio:FieldValidator:CompartmentDimension', 'Compartment dimension must be 0, 1, 2, or 3')
