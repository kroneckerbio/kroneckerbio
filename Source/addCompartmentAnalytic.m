function m = addCompartmentAnalytic(m, name, dimension, size)
%AddCompartment Add a compartment to a Model.Analytic
%
%   m = AddCompartment(m, name, dimension, size)
%
%   Compartments hold species and have size. They are purely an
%   organizational feature in analytic models and have no impact on the
%   behavior of the model.
%
%   Inputs
%   m: [ model struct scalar ]
%       The model to which the compartment will be added
%   name: [ string ]
%       A name for the compartment
%   dimension: [ 0 1 2 3 ]
%       The dimensionality of the compartment. Example: the cytoplasm would
%       be 3, the cell membrane would be 2, DNA would be 1, and the
%       centromere would be zero.
%   size: [ positive scalar | string ]
%       The size of the compartment. Example: the volume of the cytoplasm,
%       the surface area of the membrane, or the length of the DNA. The
%       compartment size can either be a constant or a string expression as a
%       function of other components in the model.
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new compartment added.

% Increment counter
nv = m.nv + 1;
m.nv = nv;
m.Compartments = growCompartments(m.Compartments, nv);

% Add item
m.Compartments(nv).Name      = fixCompartmentName(name);
m.Compartments(nv).Dimension = fixCompartmentDimension(dimension);
m.Compartments(nv).Size      = fixCompartmentSizeAnalytic(size, dimension);

m.Ready = false;
