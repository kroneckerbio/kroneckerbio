function AddOutput(this, name, expression)
%AddOutput Add an output to a MassActionAmountModel
%
%   Outputs are linear combinations of species, with a possibly non-unity
%   coefficient, in massaction models.
%
%   Inputs
%   m: [ scalar MassActionAmountModel ]
%       The model to which the output will be added
%   name: [ string ]
%       A name for the output
%   expression: [ string {name} | cell vector of strings | n x 2 cell matrix of [string, double] pairs ]
%       Names of species comprising the output. May be specified as a single
%       string, a cell vector of species names with implied unity coefficient,
%       or a n x 2 cell matrix where the 1st col contains species names and the
%       2nd col contains the coefficients. If expressions is blank expressions
%       is set to name. Species names must be
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
this.Outputs(ny).Name       = FieldValidator.OutputName(name);
this.Outputs(ny).Expression = FieldValidator.OutputExpressionMassAction(expression);

this.Ready = false;
