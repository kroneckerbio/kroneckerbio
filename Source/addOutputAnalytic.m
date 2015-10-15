function m = addOutputAnalytic(m, name, expression)
%AddOutput Add an output to a Model.Analytic
%
%   m = AddOutput(m, name, expression)
%
%   Outputs are arbitrary functions of species and rate constants in analytic
%   models.
%
%   Inputs
%   m: [ model struct scalar ]
%       The model to which the output will be added
%   name: [ string ]
%       A name for the output
%   expression: [ string {name} ]
%       Mathematical expression for the output. If blank, expression is set to
%       the state represented by the name, if possible. Species names must be
%       fully qualified with "compartment.species", including the double-quotes,
%       or unqualified with species if the species is globally unique.
%
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new output added.

% Clean up arguments
if nargin < 3
    expression = [];
end

% Set defaults
if isempty(expression)
    expression = name;
end

% Increment counter
ny = m.ny + 1;
m.ny = ny;
m.Outputs = growOutputs(m.Outputs, ny);

% Add item
m.Outputs(ny).Name        = fixOutputName(name);
m.Outputs(ny).Expression  = fixOutputAnalytic(expression);

m.Ready = false;
