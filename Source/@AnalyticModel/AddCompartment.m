function AddCompartment(this, name, dimension, size)
%AddCompartment Add a compartment to a Model.Analytic
%
%   Compartments hold species and have size. They are purely an
%   organizational feature in analytic models and have no impact on the
%   behavior of the model.
%
%   Inputs
%   this: [ scalar AnalyticModel ]
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

% Increment counter
nv = this.nv + 1;
this.nv = nv;
this.growCompartments;

% Add item
this.Compartments(nv).Name      = FieldValidator.CompartmentName(name);
this.Compartments(nv).Dimension = FieldValidator.CompartmentDimension(dimension);
this.Compartments(nv).Size      = FieldValidator.CompartmentSize(size, dimension);

this.Ready = false;
