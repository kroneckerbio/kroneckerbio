function m = addOutputAsRegex(m, name, regex)
%AddOutput Add an output to a Model.MassActionAmount using a regex. Note: this
%   function is order-dependent, meaning it only matches species already in the
%   model.
%
%   m = AddOutput(m, name, expressions)
%
%   Outputs are linear combinations of species, with a possibly non-unity
%   coefficient.
%
%   Inputs
%   m: [ model struct scalar ]
%       The model to which the output will be added
%   name: [ string ]
%       A name for the output
%   expression: [ string regex | cell vector of string regexes 
%               | n x 2 cell matrix of [string regex, double] pairs ]
%       A regex as a string for matching in the model. Matches are
%       performed on the qualified compartment.species names. The output
%       will be the sum of matched species multiplied by the specified
%       coefficient (default 1). May be specified as a single string, a
%       cell vector of strings with implied unity coefficient, or a n x 2
%       cell matrix where the 1st col contains strings and the 2nd col
%       contains the coefficients. No compartment qualification is done.
%       Beware of species whose names are substrings of other species' or
%       compartments' names. Beware: if multiple matches of a species are
%       made, then that species will be counted multiple times.
%
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new output added.

% Clean arguments
regex = fixOutputMassAction(regex);

% Get list of species from model
full_names = vec([strcat({m.States(1:m.nx).Compartment}, '.', {m.States(1:m.nx).Name}), ...
    strcat({m.Inputs(1:m.nu).Compartment}, '.', {m.Inputs(1:m.nu).Name})]);

% Search species for matches
nExpr = size(regex,1);
expression = cell(0,2);
for i = 1:nExpr
    regex_i = regex{i,1};
    coeff_i = regex{i,2};
    match = find(~cellfun(@isempty, regexp(full_names, regex_i, 'once')));
    expression = [expression; full_names(match), num2cell(ones(numel(match),1)*coeff_i)];
end

% Add item
m = addOutputMassActionAmount(m, name, expression);

