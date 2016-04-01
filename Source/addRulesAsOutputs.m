function m = addRulesAsOutputs(m)
%addRulesAsOutputs A quick and dirty helper script that adds one output
%   for each rule currently in the model. Note: this function is
%   order-dependent, meaning it only matches rules already in the model.
%
%   m = addRulesAsOutputs(m)

% (c) 2015 David R Hagen
% This work is released under the MIT license.

if is(m, 'Model.MassActionAmount')
    error('KroneckerBio:AddRule:MassActionAmount', 'Rules are not implemented for massaction models')
elseif is(m, 'Model.Analytic')
    for i = 1:m.nz
        name = m.Rules(i).Name;
        expression = quoteIfInvalid(m.Rules(i).Name);
        
        m = AddOutput(m, name, expression);
    end
else
    error('KroneckerBio:AddOutput:m', 'm must be a model')
end

end

function name = quoteIfInvalid(name)
if ~isValidIdentifier(name)
    name = ['"' name '"'];
end
end
