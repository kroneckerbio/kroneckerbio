function m = AddSeed(m, name, value)
%AddSeed Add a seed parameter to a model. The parameter must be a
%nonnegative scalar.
%
%   m = AddSeed(m, name, value)
%
%   Seed parameters determine the initial conditions of states that
%   associate with them.
%
%   Inputs
%   m: [ model struct scalar ]
%       The model to which the parameter will be added
%   name: [ string ]
%       A name for the parameter. This is the name by which states will
%       refer to it.
%   value: [ nonnegative scalar ]
%       The numeric value of the parameter
%
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new parameter added

% Increment counter
ns = m.ns + 1;
m.ns = ns;
m.Seeds = growSeeds(m.Seeds, ns);

% Add item
m.Seeds(ns).Name  = fixParameterName(name);
m.Seeds(ns).Value = fixParameterValue(value);

m.Ready = false;
