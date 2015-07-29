function expr = simbioExpr2kroneckerbioExpr(expr)
% Convert SimBio expression to kroneckerbio expression - replaces invalid names
% escaped inside square brackets [] with double-quotes "".
% SimBio only surrounds the component with invalid characters. In the case of
% compartment.species, only [compartment].[species] are bracketed while in
% kroneckerbio, all of "compartment.species" is quoted.

% Regexes that recognize bracketed components
brackets = {'(^|[^w\.])(\[.+?\])([^w\.]|$)', ... % [invalid] or [invalid].[invalid]
    '(^|\W)(\w+\.\[.+?\])(\W|$)', ... % valid.[invalid]
    '(^|\W)(\[.+?\]\.\w+)(\W|$)'}; % [invalid].valid

% Go through different replacements, replacing 1 by 1
for i = 1:length(brackets)
    while true
        tokens = regexp(expr, brackets{i}, 'tokens', 'once');
        if isempty(tokens)
            break
        end
        quoted = ['"' regexprep(tokens{2}, '(\[|\])', '') '"'];
        old = regexptranslate('escape', [tokens{1} tokens{2} tokens{3}]);
        new = [tokens{1} quoted tokens{3}];
        expr = regexprep(expr, old, new, 'once');
    end
end
