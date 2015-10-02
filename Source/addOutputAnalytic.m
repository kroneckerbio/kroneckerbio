function m = addOutputAnalytic(m, name, expression, compartment)
%AddOutput Add an output to a Model.Analytic
%
%   m = AddOutput(m, name, expression)
%
%   Outputs are linear combinations of species. In Kronecker, they take
%   advantage of species naming schemes by being defined by regular
%   expressions. Each regular expression of the output has a corresponding
%   value indicating how much each species with a matching substring
%   contributes to this output. The value of the output is the sum of all
%   matches times their values. If a species is matched more than once, the
%   last match is used.
%
%   See regexp for help on how to write Matlab regular expressions.
%
%   Inputs
%   m: [ model struct scalar ]
%       The model to which the output will be added
%   name: [ string ]
%       A name for the output
%   expression: [ string {name} ]
%       Mathematical expression for the output. If blank, expression is set to
%       the state represented by the name, if possible.
%   compartment: [ string | empty {''} ]
%       A compartment name.
%       If empty, then each species referred to in the reactants, products,
%       forward, and reverse must unambiguously refer to a single species
%       in the model, either because it is unique or because it is a
%       species full name.
%       If supplied, then each species in the reactants and products that
%       does not have a compartment will be disambiguated with this
%       compartment. The reactants and products that are disambiguated will
%       also be disambiguated in forward and reverse expressions.
%
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new output added.

% Clean up inputs
if nargin < 4
    compartment = [];
    if nargin < 3
        expression = [];
    end
end

% Set defaults
if isempty(expression)
    expression = name;
end
if isempty(compartment)
    compartment = '';
end

% Increment counter
ny = m.ny + 1;
m.ny = ny;
m.Outputs = growOutputsAnalytic(m.Outputs, ny);

% Add item
m.Outputs(ny).Name        = fixOutputName(name);
m.Outputs(ny).Expression  = fixOutputAnalytic(expression);
m.Outputs(ny).Compartment = compartment;

m.Ready = false;
