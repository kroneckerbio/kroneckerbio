function m = addParameterAnalytic(m, name, value, id)
%addParameterAnalytic Add a rate parameter to a Model.Analytic
%
%   m = addParameterAnalytic(m, name, value)
%
%   Rate parameters determine the rate of reactions to which they are
%   associated.
%
%   Inputs
%   m: [ model struct scalar ]
%       The model to which the parameter will be added
%   name: [ string ]
%       A name for the parameter. This is the name by which reactions will
%       refer to it.
%   value: [ nonnegative scalar ]
%       The numeric value of the rate parameter.
%   id: [ string {random UUID} ]
%       A unique, valid variable name
%
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new parameter added

% Generate new ID if not supplied
if nargin < 4
    id = [];
end
if isempty(id)
    id = genUID;
end

% Increment counter
nk = m.add.nk + 1;
m.add.nk = nk;
m.add.Parameters = growParametersAnalytic(m.add.Parameters, nk);

% Add item
m.add.Parameters(nk).Name = fixParameterName(name);
m.add.Parameters(nk).ID = id;
m.add.Parameters(nk).Value = fixParameterValue(value);

m.Ready = false;
