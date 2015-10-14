function m = addStatesAsOutputs(m)
%addStatesAsOutputs A quick and dirty helper script that adds one output
%   for each state currently in the model. Note: this
%   function is order-dependent, meaning it only matches species already in the
%   model.
%
%   m = addStatesAsOutputs(m)

% (c) 2015 David R Hagen
% This work is released under the MIT license.

full_names = vec(strcat({m.States(1:m.nx).Compartment}, '.', {m.States(1:m.nx).Name}));

for i = 1:numel(full_names)
    if is(m, 'Model.MassActionAmount')
        m = AddOutput(m, full_names{i});
    elseif is(m, 'Model.Analytic')
        m = AddOutput(m, full_names{i}, ['"' full_names{i} '"']); % quotes around expressions with potentially invalid names
    else
        error('KroneckerBio:AddState:m', 'm must be a model')
    end
end
