function m = addRulesAsOutputs(m)
%addRulesAsOutputs A quick and dirty helper script that adds one output
%   for each rule currently in the model
%
%   m = addRulesAsOutputs(m)

% (c) 2015 David R Hagen
% This work is released under the MIT license.

full_names = vec({m.Rules(1:m.nz).Name});

if is(m, 'Model.MassActionAmount')
    error('KroneckerBio:AddRule:MassActionAmount', 'Rules are not implemented for massaction models')
elseif is(m, 'Model.Analytic')
    full_names_expression = strcat('"', full_names, '"');

    for i = 1:numel(full_names)
        m = AddOutput(m, full_names{i}, full_names_expression{i});
    end
else
    error('KroneckerBio:addRulesAsOutputs:m', 'm must be a model')
end
