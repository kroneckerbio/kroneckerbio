function warnAboutQuotedFullName(expression)
% Now that issue #70 is implemented, recommend that "compartment.species"
% be switched in favor of compartment.species

% Extract every legal species in the expression: x1.v1, "x1.v1", "x1"."v1" etc.
identifier = '(\<[A-Za-z_][A-Za-z0-9_]*\>)';
quoted = '("[^"]*")';
name = ['(' quoted '|' identifier ')'];
species = ['(' name '(\.' name ')?)'];

for id = regexp(expression, species, 'match')
    id = id{1}; % unpack
    if ~isempty(regexp(id, ['^"' identifier '"$'], 'once')) || ~isempty(regexp(id, ['^"' identifier '\.' identifier '"$'], 'once'))
        % Totally unecessary quotes "x1" or "v1.x1"
        warning('KroneckerBio:UnnecessaryQuotedSpecies',  ...
            'Unecessary quotes around species found. Consider replacing %s with %s', ...
            id, id(2:end-1))
    elseif regexp(id, ['^' quoted '\.' quoted '$'])
        % Unecessary separation of compartment and species
        warning('KroneckerBio:UnnecessaryQuotedSpecies',  ...
            'Unecessary quotes in species found. Consider replacing %s with %s', ...
            id, ['"' strrep(id, '"', '') '"'])
    end
end
