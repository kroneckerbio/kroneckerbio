function AddCompartment(this, name, dimension, size)
%AddCompartment Add a compartment to a Model.MassActionAmount
%
%   Compartments hold species and have size. The size affects the rates of
%   bimolecular reactions.
%
%   Inputs
%   m: [ model struct scalar ]
%       The model to which the compartment will be added
%   name: [ string ]
%       A name for the compartment
%   dimension: [ 0 1 2 3 ]
%       The dimensionality of the compartment. Example: the cytoplasm would
%       be 3, the cell membrane would be 2, DNA would be 1, and the
%       centromere would be zero. This is only used to determine which
%       compartment's volume plays a part in the rate of a bimolecular
%       reaction between compartments.
%   size: [ positive scalar | cell array ]
%       The size of the compartment. Example: the volume of the cytoplasm,
%       the surface area of the membrane, or the length of the DNA. The
%       compartment size can either be a constant or a cell array. The
%       cell array must have two columns. Each row is a string and a
%       positive scalar. The string is a regular expression matching the
%       full name of species in the model. The scalar is the amount that
%       each species matched contributes to the size of this compartment.
%       if the string is empty, then this row is a constant contribution to
%       the compartment size.

% Increment counter
nv = this.nv + 1;
this.nv = nv;
this.growCompartments;

% Add item
this.Compartments(nv).Name      = FieldValidator.CompartmentName(name);
this.Compartments(nv).Dimension = FieldValidator.CompartmentDimension(dimension);
this.Compartments(nv).Size      = FieldValidator.CompartmentSizeMassAction(size, dimension);

this.Ready = false;
