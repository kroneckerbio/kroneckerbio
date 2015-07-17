function m = addSeedMassActionAmount(m, name, value)
%AddSeed Add a seed parameter to a KroneckerBio model
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

% (c) 2013 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Increment counter
ns = m.add.ns + 1;
m.add.ns = ns;
m.add.Seeds = growSeeds(m.add.Seeds, ns);

% Add item
m.add.Seeds(ns).Name = fixParameterName(name);
m.add.Seeds(ns).Value = fixParameterValue(value);

m.Ready = false;
