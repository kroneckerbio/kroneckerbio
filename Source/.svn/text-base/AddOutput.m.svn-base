function m = AddOutput(m, name, expressions, values, units)
%AddOutput Add an output to a KroneckerBio model
%
%   m = AddOutput(m, name, expressions, values, units)
%
%   Outputs are linear combinations of species. In Kronecker, they take
%   advantage of species naming schemes by being defined by regular
%   expressions. Each regular expression of the output has a corresponding
%   value indicating how much each species with a matching substring
%   contributes to this output. The value of the output is the sum of all
%   matches times their values.
%
%   See regexp for help on how to write Matlab regular expressions.
%
%   Inputs
%   m: [ model struct scalar ]
%       The model to which the output will be added
%   name: [ string ]
%       A name for the output
%   expressions: [ string | cell array of strings ]
%       Each species full name will be substring matched against each
%       regular expression provided here. If the string is empty, it means
%       that the output has a constant component.
%   values: [ nonnegative scalar | nonnegative vector numel(expressions)
%             {1} ]
%       Optional
%       Each entry in this vector is associated with one of the
%       expressions. This value tells how much a species will contribute
%       when it matches the corresponding expressions. If a species is
%       matched by multiple expressions, the last expression to match
%       overrides all others.
%   units: [ '' ]
%       This is currently not implemented
%
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new output added.

% (c) 2011 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean-up inputs
if nargin < 5
    units = [];
    if nargin < 4
        values = [];
    end
end

% Increment counter
ny = m.add.ny + 1;
m.add.ny = ny;
m.add.Outputs = growOutputs(m.add.Outputs, ny);

% Add item
m.add.Outputs(ny).Name = fixName(name);
m.add.Outputs(ny).Expressions = fixExpression(expressions);
nExpr = numel(m.add.Outputs(ny).Expressions);
m.add.Outputs(ny).Values = fixOutputValues(nExpr, values, fixUnits(units));

m.Ready = false;