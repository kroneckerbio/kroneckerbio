function warnAboutQuotedFullName(expression)
% Now that issue #70 is implemented, recommend that "compartment.species"
% be switched in favor of compartment.species

for id = regexp(expression, '((("[^"]*")|(\<[A-Za-z_][A-Za-z0-9_]*\>))(\.(("[^"]*")|(\<[A-Za-z_][A-Za-z0-9_]*\>)))?)', 'match')
    id = id{1}; % unpack
    if regexp(id, '^"[^."]*\.[^."]*"$')
        compartment = regexp(id(2:end), '^[^.]*', 'match', 'once');
        if ~isValidIdentifier(compartment)
            compartment = ['"' compartment '"'];
        end
        species = regexp(id(1:end-1), '[^.]*$', 'match', 'once');
        if ~isValidIdentifier(species)
            species = ['"' species '"'];
        end
        warning('KroneckerBio:UnnecessaryQuotedSpecies',  ...
            'Putting quotes around the dot seperating the compartment and species is no longer required. Consider replacing %s with %s', ...
            id, [compartment '.' species])
    end
end
