function AddOutput(this, name, expression)
%AddOutput Add an output to a  AnalyticModel
%
%   Outputs are arbitrary functions of species and rate constants in analytic
%   models.
%
%   Inputs
%   m: [ scalar AnalyticModel ]
%       The model to which the output will be added
%   name: [ string ]
%       A name for the output
%   expression: [ string {name} ]
%       Mathematical expression for the output. If blank, expression is set to
%       the state represented by the name, if possible. Species names must be
%       fully qualified with "compartment.species", including the double-quotes,
%       or unqualified with species if the species is globally unique.

% Clean up arguments
if nargin < 3
    expression = [];
end

% Set defaults
if isempty(expression)
    expression = name;
end

% Increment counter
ny = this.ny + 1;
this.ny = ny;
this.growOutputs;

% Add item
this.Outputs(ny).Name        = FieldValidator.OutputName(name);
this.Outputs(ny).Expression  = FieldValidator.OutputExpression(expression);

this.Ready = false;
