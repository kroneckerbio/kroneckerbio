function expr = simbioExpr2kroneckerbioExpr(expr)
% Convert SimBio expression to kroneckerbio expression - replaces invalid names
%   escaped inside square brackets [] with double-quotes "" and double-quotes
%   compartment.species (since this implicitly contains a dot, which is an invalid
%   character)
% SimBio only surrounds the component with invalid characters. In the case of
%   compartment.species, only [compartment].[species] are bracketed while in
%   kroneckerbio, all of "compartment.species" is quoted.

% Regex that matches SimBio species, which are assumed to be qualified compartment.species
% Species names must either be valid Matlab identifiers (can't start with
%   underscore) or bracketed (will be quoted).
% Note: check that SimBio's species naming convention matches Matlab's valid
%   variable naming requirement
% Matches the species in: 'compartment.species + c.s + [c*].[s*] + c.[s*] +
%   [c*].s + 1.2', but not the last term.
match = '\<(([a-zA-Z]\w*|\[[^\]]+\])\.([a-zA-Z]\w*|\[[^\]]+\]))\>';

% Fix all species names
tokens = regexp(expr, match, 'tokens');
for i = 1:length(tokens)
    token = tokens{i}{1};
    quoted = ['"' regexprep(token, '[\[|\]]', '') '"'];
    expr = regexprep(expr, regexptranslate('escape', token), quoted, 'once');
end
