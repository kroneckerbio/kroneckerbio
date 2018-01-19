function m = addOutputMassActionAmount(m, name, expression)
%AddOutput Add an output to a Model.MassActionAmount
%
%   m = AddOutput(m, name, expression)
%
%   Outputs are linear combinations of species, with a possibly non-unity
%   coefficient, in massaction models.
%
%   Inputs
%   m: [ model struct scalar ]
%       The model to which the output will be added
%   name: [ string ]
%       A name for the output
%   expression: [ string {name} | cell vector of strings 
%               | n x 2 cell matrix of [string, double] pairs ]
%       Names of species contributing to the output. May be specified as a
%       single string, a cell vector of species names with implied unity
%       coefficient, or a n x 2 cell matrix where the 1st col contains
%       species names and the 2nd col contains the coefficients. If
%       expressions is blank, expressions is set to name. Species names
%       must be fully qualified as compartment.species or unqualified as
%       species if the species is globally unique.
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
m.Outputs(ny).Name       = fixOutputName(name);
m.Outputs(ny).Expression = fixOutputMassAction(expression);

m.Ready = false;
