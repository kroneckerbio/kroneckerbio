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
    outputs_need_quotes = ~isValidIdentifier(outputs);
    outputs(outputs_need_quotes) = strcat('"', outputs(outputs_need_quotes), '"');
end

% Capture all quoted expressions and bare variable names and run them
% through a function that looks up the name. If it is found replace it,
% otherwise leave it the same
temp = @transmute; % Matlab bug can't find local functions
expressions = regexprep(expressions, '("[^"]*")|(\<[A-Za-z_][A-Za-z0-9_]*\>)', '${temp($1)}');

    function output = transmute(input)
        if input(1) == '"'
            % Quoted branch
            index = lookupmember(input(2:end-1), inputs);
        else
            % Identifier branch
            index = lookupmember(input, inputs);
        end
        
        if index ~= 0
            % Replace
            output = outputs{index};
        else
            % Echo
            output = input;
        end
    end
end
