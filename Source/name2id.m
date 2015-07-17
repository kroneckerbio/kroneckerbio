function expr = name2id(expr, names, ids, xuvNames)
% Parse arbitrary mathematical expression, replacing user-specified names with
% valid ids
% Inputs:
%   expr [ string ]
%       Input expression containing user-specified component names. Enclose
%       names with invalid characters with double-quotes.
%   names [ cell array of strings ]
%       Cell array of user-specified names. States and inputs appear first.
%       Invalid names should not be surrounded by quotes.
%   ids [ cell array of strings ]
%       Cell array of valid IDs in same order as names
%   xuvNames [ cell array of strings {''} ]
%       Cell array of state and input names, in same order as states and inputs
%       in names and ids. If blank, compartment.names aren't tried.
% Outputs:
%   expr [ string ]
%       Input expression with names replaced by IDs

if nargin < 4
    xuvNames = [];
end

% Make list of compartment.species
nxu = length(xuvNames);
if ~isempty(xuvNames)
    xuFullNames = strcat(xuvNames, '.', names(1:nxu));
end

% Clean up input expression
% Tokenize expression into potentially substitutable parts
parts = regexp(expr, '[\w\.]+|"[^"]+"', 'match');
for i = 1:length(parts)
    
    part = parts{i};
    partNoQuotes = strrep(part, '"', '');
    
    % Match directly
    match = ismember(names, partNoQuotes);
    if sum(match) == 1 % safer than any(match) by deferring to compartment.species for repeated species
        expr = replaceSymbolRegex(expr, part, ids{match});
        continue
    end
    
    % Match compartment.species
    if ~isempty(xuvNames)
        match = ismember(xuFullNames, partNoQuotes);
        if sum(match) == 1
            expr = replaceSymbolRegex(expr, part, ids{match});
            continue
        end
    end
    
    % Fell through - invalid and ignore
    % Includes numbers, extra spaces, and mathematical functions
    % TODO: maybe make the tokenizer/parser smarter to warn about mistyped
    % components - difficult because of mathematical functions and whatnot that
    % look like components
    
end

end
