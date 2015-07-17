function m = addSeedAnalytic(m, name, value, id)
%AddSeed Add a seed parameter to a Model.Analytic
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
%   id: [ string {random UUID} ]
%       A unique, valid variable name
%
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new parameter added

% Generate new ID if not supplied
if nargin < 5
    id = [];
end
if isempty(id)
    id = genUID;
end

% Increment counter
ns = m.add.ns + 1;
m.add.ns = ns;
m.add.Seeds = growSeedsAnalytic(m.add.Seeds, ns);

% Add item
m.add.Seeds(ns).Name = fixParameterName(name);
m.add.Seeds(ns).ID = id;
m.add.Seeds(ns).Value = fixParameterValue(value);

m.Ready = false;
