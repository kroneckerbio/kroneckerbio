function AddInput(this, name, compartment, default)
%AddInput Add an input species to a Model
%
%   Inputs
%   m: [ scalar Model ]
%       The model to which the input species will be added
%   name: [ string ]
%       A name for the input. This is the name by which reactions will
%       refer to it.
%   compartment: [ string ]
%       The name of the compartment to which it will be added.
%   default: [ nonnegative scalar {0} ]
%       The default value for this input.

if nargin < 4
    default = [];
end
if isempty(default)
    default = 0;
end

% Increment counter
nu = this.nu + 1;
this.nu = nu;
this.growInputs;

% Add item
this.Inputs(nu).Name         = FieldValidator.SpeciesName(name);
this.Inputs(nu).Compartment  = FieldValidator.CompartmentName(compartment);
this.Inputs(nu).DefaultValue = FieldValidator.InputDefaultValue(default);

this.Ready = false;
