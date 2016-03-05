function m = addInputsAsOutputs(m, include_compartment)
%addInputsAsOutputs A quick and dirty helper script that adds one output
%   for each input currently in the model. Note: this function is
%   order-dependent, meaning it only matches inputs already in the model.
%
%   m = addInputsAsOutputs(m, include_compartment)

% (c) 2015 David R Hagen
% This work is released under the MIT license.

if nargin < 2
    include_compartment = true;
end

if is(m, 'Model.MassActionAmount')
    for i = 1:m.nu
        if include_compartment
            name = [m.Inputs(i).Compartment '.' m.Inputs(i).Name];
        else
            name = m.Inputs(i).Name;
        end
        
        m = AddOutput(m, name, name);
    end
elseif is(m, 'Model.Analytic')
    for i = 1:m.nu
        if include_compartment
            name = [m.Inputs(i).Compartment '.' m.Inputs(i).Name];
            expression = [quoteIfInvalid(m.Inputs(i).Compartment) '.' quoteIfInvalid(m.Inputs(i).Name)];
        else
            name = m.Inputs(i).Name;
            expression = quoteIfInvalid(m.Inputs(i).Name);
        end
        
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
