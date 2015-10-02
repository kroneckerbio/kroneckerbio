function m = addOutputMassActionAmount(m, name, expressions, compartment)
%AddOutput Add an output to a KroneckerBio model
%
%   m = AddOutput(m, name, expressions)
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
%   expressions: [ string | nonnegative scalar | cell vector of strings
%                 | cell matrix {0} ]
%       The expressions list is a ? by 2 cell matrix where the elements in
%       the first column contain a regular expression and the elements in
%       the second column contain a scalar. If just a string is given, the
%       default scalar is 1. If just a scalar is given or a scalar is
%       paired with an empty string, then that is assumed to be a
%       constitutive amount added to the output value.
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

% (c) 2013 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean-up
if nargin < 4
    compartment = [];
end

if isempty(compartment)
    compartment = '';
end

% Increment counter
ny = m.ny + 1;
m.ny = ny;
m.Outputs = growOutputs(m.Outputs, ny);

% Add item
m.Outputs(ny).Name        = fixOutputName(name);
m.Outputs(ny).Expressions = fixOutputExpressions(expressions);
m.Outputs(ny).Compartment = compartment;

m.Ready = false;
