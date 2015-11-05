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

if include_compartment
    full_names = vec(strcat({m.Inputs(1:m.nu).Compartment}, '.', {m.Inputs(1:m.nu).Name}));
else
    full_names = vec({m.Inputs(1:m.nu).Name});
end

for i = 1:numel(full_names)
    if is(m, 'Model.MassActionAmount')
        m = AddOutput(m, full_names{i});
    elseif is(m, 'Model.Analytic')
        m = AddOutput(m, full_names{i}, ['"' full_names{i} '"']); % quotes around expressions with potentially invalid names
    else
        error('KroneckerBio:AddOutput:m', 'm must be a model')
    end
end