function expr = replaceSymbolRegex(expr, symbolIn, symbolOut)
% Replace literal 'symbolIn' with 'symbolOut' in 'expr' only when symbolIn occurs as a word, i.e.,
% surrounded by whitespace and mathematical expression or at the beginnning or end of expr.
% Used when substituting names and ids in mathematical expressions.
literal = regexptranslate('escape', symbolIn);

if symbolIn(1) == '"' && symbolIn(end) == '"' % If symbolIn surrounded by quotes
    expr = regexprep(expr, literal, symbolOut);
else  % Not surrounded by quotes
    % Literal contains only valid characters '[\w\.]'
    % Ignore things in quotes by looking for matches between even numbers of quotes
    % Remove quotes placed around valid symbols
    % Repeat as many times as needed for overlapping matches (when 2 instances are separated by a single position)
    regex = ['(?=([^\"]*[\"][^\"]*[\"])*[^\"]*$)(^|\W)(\"?', literal, '\"?)(\W|$)(?=([^\"]*[\"][^\"]*[\"])*[^\"]*$)'];
    while true
        exprOld = expr;
        expr = regexprep(expr, regex, ['$1', symbolOut, '$3']);
        if strcmpi(expr, exprOld) % no more replacements
            break
        end
    end
end
