function m = addRulesAsOutputs(m)
%addRulesAsOutputs A quick and dirty helper script that adds one output
%   for each rule currently in the model
%
%   m = addRulesAsOutputs(m)

% (c) 2015 David R Hagen
% This work is released under the MIT license.

all_rules = [m.Rules; m.add.Rules(1:m.add.nz)];

full_names = vec({all_rules.Name});

if is(m, 'Model.MassActionAmount')
    error('NotImplemented')
elseif is(m, 'Model.MassActionConcentration')
    error('NotImplemented')
elseif is(m, 'Model.Analytic')
    full_names_expression = strcat('"', full_names, '"');

    for i = 1:numel(full_names)
        m = AddOutput(m, full_names{i}, full_names_expression{i});
    end
else
    error('KroneckerBio:addRulesAsOutputs:m', 'm must be a model')
end
