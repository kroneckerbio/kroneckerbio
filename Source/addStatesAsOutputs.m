function m = addStatesAsOutputs(m)
%addStatesAsOutputs A quick and dirty helper script that adds one output
%   for each state currently in the model
%
%   m = addStatesAsOutputs(m)

% (c) 2015 David R Hagen
% This work is released under the MIT license.

all_states = [m.States; m.add.States(1:m.add.nx)];

full_names = vec({all_states.Name});

if is(m, 'Model.MassActionAmount')
    full_names_regex = strcat('^', full_names, '$');
    
    for i = 1:numel(full_names)
        m = AddOutput(m, full_names{i}, full_names_regex{i});
    end
elseif is(m, 'Model.MassActionConcentration')
    error('NotImplemented')
elseif is(m, 'Model.Analytic')
    full_names_expression = strcat('"', full_names, '"');

    for i = 1:numel(full_names)
        m = AddOutput(m, full_names{i}, full_names_expression{i});
    end
else
    error('KroneckerBio:AddState:m', 'm must be a model')
end
