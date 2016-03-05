function expressions = substituteQuotedExpressions(expressions, inputs, outputs, add_output_quotes)
% Replaces instances of inputs with outputs
% If an input is a valid identifier, it may be found in the expression as
%   either a bare identifier or quoted.
% If an input is not a valid identifier, it must be quoted.
% The outputs will only be quoted if both add_quotes is true and the output
%   is not a valid identifier
% Note that fully qualified compartment.species will always be quoted, since it
%   contains a dot, which isn't a valid identifier

if nargin < 4 || isempty(add_output_quotes)
    add_output_quotes = false;
end

if ischar(inputs)
    inputs = {inputs};
end
if ischar(outputs)
    outputs = {outputs};
end

if add_output_quotes
    for i = 1:numel(outputs)
        dot_location = find(outputs{i} == '.', 1);
        if isempty(dot_location)
            % Plain identifier
            if ~isValidIdentifier(outputs{i})
                outputs{i} = ['"' outputs{i} '"'];
            end
        else
            % compartment.species
            compartment = outputs{i}(1:dot_location-1);
            species = outputs{i}(dot_location+1:end);
            if ~isValidIdentifier(compartment)
                compartment = ['"' compartment '"'];
            end
            if ~isValidIdentifier(species)
                species = ['"' species '"'];
            end
            outputs{i} = [compartment '.' species];
        end
    end
end

% Capture all quoted expressions and bare variable names and run them
% through a function that looks up the name. If it is found replace it,
% otherwise leave it the same
temp = @transmute; % Matlab bug can't find local functions
expressions = regexprep(expressions, '((("[^"]*")|(\<[A-Za-z_][A-Za-z0-9_]*\>))(\.(("[^"]*")|(\<[A-Za-z_][A-Za-z0-9_]*\>)))?)', '${temp($1)}');

    function output = transmute(input)
        % Remove quote characters
        stripped_input = strrep(input, '"', '');
        
        index = lookupmember(stripped_input, inputs);
        
        if index ~= 0
            % Replace
            output = outputs{index};
        else
            % Echo
            output = input;
        end
    end
end
