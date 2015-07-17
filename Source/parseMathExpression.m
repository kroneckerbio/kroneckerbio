function exprOut = parseMathExpression(exprIn, type, m)
% Parse arbitrary mathematical expression in terms of model components for a
% given type of expression, returning expression with user-specified names
% replaced with unique, nice IDs.
% Inputs:
%   exprIn [ string ]
%       Input expression containing user-specified component names. Enclose
%       names with invalid characters with double-quotes.
%   type [ string ]
%       Type of expression being parsed. Allowed values:
%           - 'rule': allows states, inputs, seeds, and rate parameters
%           - 'ic': initial condition, only allows seeds
%           - 'output': allows states, inputs, and rate parameters
%   m [ Model.Analytic struct ]
%       Model the expression is a part of
% Outputs:
%   exprOut [ string ]
%       Output expression with user-specified names replaced by IDs

%% Get lists of allowed component names and corresponding IDs
% Draw from existing and added components
% Make lists of cleaned up model component names
vNames = {m.Compartments.Name};
vIDs   = {m.Compartments.ID}';

if strcmpi(type, 'ic') || strcmpi(type, 'rule')
    sNames = {m.Seeds.Name}';
    sIDs   = {m.Seeds.ID}';
end

if strcmpi(type, 'output') || strcmpi(type, 'rule')
    xuNames        = [{m.States.Name}'; {m.Inputs.Name}'];
    xuCompartments = [{m.States.Compartment}'; {m.Inputs.Compartment}'];
    xuIDs          = [{m.States.ID}'; {m.Inputs.ID}'];
    kNames         = {m.Parameters.Name}';
    kIDs           = {m.Parameters.ID}';
end

%% Clean up input expression
% Tokenize expression into potentially substitutable parts
parts = regexp(exprIn, '[\w\.]+|"[^"]+"', 'match');
exprOut = exprIn;
for i = 1:length(parts)
    
    part = parts{i};
    partNoQuotes = strrep(part, '"', '');
    
    id = getVolume(partNoQuotes);
    if ~isempty(id)
        exprOut = replaceSymbolRegex(exprOut, part, id);
        continue
    end
    
    if strcmpi(type, 'ic') || strcmpi(type, 'rule')
        id = getSeed(partNoQuotes);
        if ~isempty(id)
            exprOut = replaceSymbolRegex(exprOut, part, id);
            continue
        end
    end
    
    if strcmpi(type, 'output') || strcmpi(type, 'rule')
        
        id = getParameter(partNoQuotes);
        if ~isempty(id)
            exprOut = replaceSymbolRegex(exprOut, part, id);
            continue
        end
        
        id = getSpecies(partNoQuotes);
        if ~isempty(id)
            exprOut = replaceSymbolRegex(exprOut, part, id);
            continue
        end
        
    end
    
    % Fell through everything - invalid and ignore
    % Includes numbers, extra spaces, and mathematical functions
    % TODO: maybe make the tokenizer/parser smarter to warn about mistyped
    % components - difficult because of mathematical functions and whatnot that
    % loo like components
    
end

%% Helper functions
% get* functions try to extract the id from the expression, returning the id
% string if found, empty [] otherwise
    function id = getVolume(expr)
        id = [];
        match = ismember(vNames, expr);
        if any(match)
            id = vIDs{match};
        end
    end

    function id = getSeed(expr)
        id = [];
        match = ismember(sNames, expr);
        if any(match)
            id = sIDs{match};
        end
    end

    function id = getParameter(expr)
        id = [];
        match = ismember(kNames, expr);
        if any(match)
            id = kIDs{match};
        end
    end
    
    function id = getSpecies(expr)
        id = [];
        
        % Prepend "compartment." to species, returning if species invalid/not found
        try
            expr = fixSpeciesFullName(expr, [], m); % throws error if invalid/not found species name
        catch
            return
        end
        if isempty(expr)
            return
        end
        
        % Expr is now a valid compartment.species
        fullName = strsplit(expr, '.');
        match = ismember(xuCompartments, fullName{1}) & ismember(xuNames, fullName{2});
        if any(match)
            id = xuIDs{match};
        end
    end

end
