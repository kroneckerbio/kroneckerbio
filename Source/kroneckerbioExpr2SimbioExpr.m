function expr = kroneckerbioExpr2SimbioExpr(expr)
% Convert kroneckerbio (double-quoted) name to SimBiology bracketed name.
%   Assumes the input expr is a valid kroneckerbio expression; even if expr is a
%   single term, it must be double-quoted if it contains invalid chars.
% Note, this also takes care of stuff like '"a@".b' -> '[a@].b'

temp = @transmute; % Matlab bug can't find local functions
expr = regexprep(expr, '"(.*?)"', '${temp($1)}');

    function output = transmute(input)
        % Remove quote characters
        input = strrep(input, '"', '');
        
        % Split into compartment.species if applicable
        parts = strsplit(input, '.');
        
        % Surround in brackets if necessary
        for i = 1:length(parts)
            if ~isValidIdentifier(parts{i})
                parts{i} = ['[' parts{i} ']'];
            end
        end
        
        % Combine compartment.species if multiple parts
        output = strjoin(parts, '.');
    end

end

