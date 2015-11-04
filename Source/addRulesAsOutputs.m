function m = addRulesAsOutputs(m)
%addRulesAsOutputs A quick and dirty helper script that adds one output
%   for each rule currently in the model. Note: this function is
%   order-dependent, meaning it only matches rules already in the model.
%
%   m = addRulesAsOutputs(m)

% (c) 2015 David R Hagen
% This work is released under the MIT license.

full_names = vec({m.Rules(1:m.nz).Name});

for i = 1:numel(full_names)
    if is(m, 'Model.MassActionAmount')
        error('KroneckerBio:AddRule:MassActionAmount', 'Rules are not implemented for massaction models')
    elseif is(m, 'Model.Analytic')
        m = AddOutput(m, full_names{i}, ['"' full_names{i} '"']); % quotes around expressions with potentially invalid names
    else
        error('KroneckerBio:AddOutput:m', 'm must be a model')
    end
end
