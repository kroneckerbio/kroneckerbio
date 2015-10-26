function AddSeed(this, name, value)
%AddSeed Add a seed parameter to a model. The parameter must be a
%nonnegative scalar.
%
%   Seed parameters determine the initial conditions of states that
%   associate with them.
%
%   Inputs
%   this: [ scalar Model ]
%       The model to which the parameter will be added
%   name: [ string ]
%       A name for the parameter. This is the name by which states will
%       refer to it.
%   value: [ nonnegative scalar ]
%       The numeric value of the parameter

% Increment counter
ns = this.ns + 1;
this.ns = ns;
this.growSeeds;

% Add item
this.Seeds(ns).Name  = FieldValidator.ParameterName(name);
this.Seeds(ns).Value = FieldValidator.ParameterValue(value);

this.Ready = false;
