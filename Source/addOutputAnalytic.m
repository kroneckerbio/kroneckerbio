function m = addOutputAnalytic(m, name, expression, id)
%AddOutput Add an output to a Model.Analytic
%
%   m = AddOutput(m, name, expression, id, opts)
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
%   id: [ string {random UUID} ]
%       A unique, valid variable name
%
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new output added.

% Clean up inputs
if nargin < 4
    id = [];
    if nargin < 3
        expression = [];
    end
end

% Set defaults
if isempty(expression)
    expression = name;
end
if isempty(id)
    id = genUID;
end

% Increment counter
ny = m.add.ny + 1;
m.add.ny = ny;
m.add.Outputs = growOutputsAnalytic(m.add.Outputs, ny);

% Add item
m.add.Outputs(ny).Name = fixOutputName(name);
m.add.Outputs(ny).ID = id;
m.add.Outputs(ny).Expression = fixOutputAnalytic(expression);

m.Ready = false;