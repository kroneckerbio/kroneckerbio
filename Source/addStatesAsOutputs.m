function m = addStatesAsOutputs(m, include_compartment)
%addStatesAsOutputs A quick and dirty helper script that adds one output
%   for each state currently in the model. Note: this function is
%   order-dependent, meaning it only matches states already in the model.
%
%   m = addStatesAsOutputs(m, include_compartment)

% (c) 2015 David R Hagen
% This work is released under the MIT license.

if nargin < 2
    include_compartment = true;
end

if is(m, 'Model.MassActionAmount')
    for i = 1:m.nx
        if include_compartment
            name = [m.States(i).Compartment '.' m.States(i).Name];
        else
            name = m.States(i).Name;
        end
        
        m = AddOutput(m, name, name);
    end
elseif is(m, 'Model.Analytic')
    for i = 1:m.nx
        if include_compartment
            name = [m.States(i).Compartment '.' m.States(i).Name];
            expression = [quoteIfInvalid(m.States(i).Compartment) '.' quoteIfInvalid(m.States(i).Name)];
        else
            name = m.States(i).Name;
            expression = quoteIfInvalid(m.States(i).Name);
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
