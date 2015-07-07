function symModel = AddOutputsToSymbolic(symModel, yNames, yExprs, opts)
%AddOutputsToSymbolic Add outputs to symbolic model before converting to
%   symbolic.
% 
%   AddOutputsToSymbolic(symModel, yNames, yExprs)
%
%   Inputs
%   symModel [symbolic model struct]
%       Symbolic model in which to add outputs
%   yNames [ string | cell array of strings {} ]
%       If yNames in nonempty and yExprs is empty, attempts to convert selected
%       expression in states/inputs given in yNames to outputs.
%       If single string, have only 1 output.
%   yExprs [ string | cell array of strings {} ]
%       Mathematical expressions for yNames. Species names with invalid
%       Matlab variable names should be enclosed in double quotes ("").
%       Non-default compartments or compartments for species that appear in
%       multiple compartments should be specified as compartment.species.
%   opts [ struct ]
%       Options struct allowing the following fields:
%       .Verbose: [ integer {0} ]
%           Level of logging detail.
%       .DefaultCompartment: [ string {1st compartment in symModel} ]
%           Name of default compartment for species expressions. If blank, the
%           1st compartment is the default.
%
%   Outputs
%   symModel [ symbolic model struct ]
%       Contains added outputs
%
% Note: make idempotent/give warning, add additional outputs; maybe make default
%   with no yNames or yExprs that turns all states into outputs (or consider 
%   adding a separate function that calls this one)

% Clean up inputs
if nargin < 4
    opts = [];
    if nargin < 3
        yExprs = [];
        if nargin < 2
            yNames = [];
        end
    end
end

% Default options
defaultOpts = [];
defaultOpts.Verbose = 0;
defaultOpts.DefaultCompartment = symModel.vNames{1};

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
if verbose; fprintf('Adding outputs to symbolic model...\n'); end

%% Clean up yNames and yExprs
% Single output string
if ischar(yNames) 
    yNames = {yNames};
end
if ischar(yExprs)
    yExprs = {yExprs};
end

% Make col vectors
[~, nCol] = size(yNames);
if nCol > 1
    yNames = yNames';
end

% If only yNames were provided, then each output expression is directly 
%   written into the name.
if ~isempty(yNames) && isempty(yExprs)
    yExprs = yNames;
end

ny = length(yNames);

%% Make new outputs
delimiter = '[\s+-*/\^\(\)]+'; % whitespace and mathematical operators

y = sym(zeros(ny,1));
yStrings = cell(ny,1); % maintains readable expression and double checking symbolic substitution
for i = 1:ny
    yExpr = yExprs{i};
    
    % Replace all expressions in quotes with cleaned versions
    yExpr = cleanExpr(yExpr);
    
    % Get terms to convert to symbolics
    % Note: captures cruft like blanks, numbers, operators, and invalid
    %   names, which must be checked/removed later
    C = strsplit(yExpr, delimiter, 'DelimiterType', 'RegularExpression');
    
    % Get valid identifiers from tokens: species (x and u) and rate constants (k)
    nC = length(C);
    ids = [];
    cs = [];
    for j = 1:nC
        c = C{j};
        [id, status] = getIDFromName(c, symModel, opts.DefaultCompartment);
        if status >= 0
            ids = [ids; id];
            cs = [cs; sym(c)]; % save good c's for symbolic sub
        end
    end
    
    % Substitute
    yi = sym(yExpr);
    y(i) = subs(yi, cs, ids);
    yStrings{i} = char(yi);
end

%% Write all outputs back to symbolic model
% Append to existing outputs if any exist
symModel.yNames   = [symModel.yNames;   yNames];
symModel.y        = [symModel.y;        y];
symModel.yStrings = [symModel.yStrings; yStrings];

end

%% Helper functions
function [id, status] = getIDFromName(name, symModel, defaultCompartment)
% Returns species or parameter ID from name in symModel
%   If no compartment, either corresponding compartment (unique name) or default
%   compartment (appears in multiple compartments)
% Inputs:
%   name [ string ]
%       Putative model component: parameter or species
%   symModel [ symbolic model struct ]
%       Symbolic model in which to search for component
%   defaultCompartment [ string {1st compartment in symModel} ]
%       Default compartment in which to search for species; ignored for
%       parameters
% Outputs:
%   id [ symbolic var ]
%       Unique model component symbolic variable
%   status [ integer -1,0,1 ] 
%       Code for id type:
%           0: parameter
%           1: species (state or input)
%           -1: error

% Clean up inputs; add default compartment
if nargin < 3
    defaultCompartment = symModel.vNames{1};
end

verbose = false;
defaultId = sym;

% Try parameter lookup
ks = symModel.kNames;
kSyms = symModel.kSyms;
kMask = ismember(ks, name);
if any(kMask)
    id = kSyms(kMask);
    status = 0;
    return
end

% Try species lookup
try
    [compartment, species] = cleanSpeciesName(name, symModel, defaultCompartment);
catch err
    % Rethrow fatal errors
    % TODO: the following error types also match random bits of expressions:
    %   CleanSpeciesName:InvalidNameFormat
    %   CleanSpeciesName:SpeciesNotFoundAny
    %   CleanSpeciesName:SpeciesNotFound
    %   CleanSpeciesName:CompartmentNotFound
    % Fix the parser to warn the user of invalid species and not numbers,
    % mathematical functions, blank spaces, etc.
    switch err.identifier
        case 'CleanSpeciesName:TooManyParts'
            throw(err)
    end
    id = defaultId;
    status = -1;
    return
end
% Valid compartment and species assumed here.

% Get list of species in right compartment
allSpecies = [symModel.xNames; symModel.uNames];
invalid = '[^a-zA-Z0-9_\".]'; % any invalid character
allSpecies = regexprep(allSpecies, invalid, '_');
allCompartments = symModel.vNames;
allSpeciesCompartmentInds = [symModel.vxInd; symModel.vuInd];
compartmentInd = find(ismember(allCompartments, compartment));

speciesMask = allSpeciesCompartmentInds == compartmentInd & ismember(allSpecies, species);

if ~any(speciesMask) % not found; either invalid species or a constant or other control character or function
    if verbose; fprintf('Invalid identifier detected - ignoring.\n'); end
    id = defaultId;
    status = -1;
    return
end

% Get symbolic var of correct species
allSpeciesSyms = [symModel.xSyms; symModel.uSyms];
id = allSpeciesSyms(speciesMask);
status = 1;
end

function [compartment, species] = cleanSpeciesName(name, symModel, defaultCompartment)
% Check and clean up species compartment and name in model. Uses cleaned up
% names in quotes with underscores.
% Inputs:
%   name [ string ]
%       "name" or "compartment.name"
%   symModel [ symbolic model struct ]
%       Symbolic model in which to search for species
%   defaultCompartment [ string {1st compartment in symModel} ]
%       Default compartment in which to search for species
% Outputs:
%   compartment [ string ]
%       Compartment name
%   species [ string ]
%       Species name
% Throws: various errors indicating invalid symbols to ignore and invalid
%   species names to warn user about
% TODO: clean up error handling

% Clean up inputs; add default compartment
if nargin < 3
    defaultCompartment = symModel.vNames{1};
end

allSpecies = [symModel.xNames; symModel.uNames];
allCompartments = symModel.vNames;
allSpeciesCompartmentInds = [symModel.vxInd; symModel.vuInd];

% Cleaned up versions of species names
invalid = '[^a-zA-Z0-9_\".]'; % any invalid character
allSpecies = regexprep(allSpecies, invalid, '_');

% Split into compartment.species components and look up in model
C = strsplit(name, '.');
nC = length(C);
switch nC
    case 1 % "species"
        species = C{1};
        
        % See if species appears in multiple compartments
        speciesMask = ismember(allSpecies, species);
        nSpeciesOcurrences = sum(speciesMask);
        switch nSpeciesOcurrences
            case 0 % not found
                error('CleanSpeciesName:SpeciesNotFoundAny', 'Species %s not found in any compartment', species)
            case 1 % unique species in 1 compartment
                compartment = allCompartments{allSpeciesCompartmentInds(speciesMask)};
            otherwise % > 1, non-unique species in multiple compartments; use default
                % Note: error if species not in default compartment
                compartment = defaultCompartment;
        end
    case 2 % "species.name"
        compartment = C{1};
        species = C{2};
        
        % Verify compartment exists
        compartmentMask = ismember(allCompartments, compartment);
        if ~any(compartmentMask)
            error('CleanSpeciesName:CompartmentNotFound', 'Compartment %s not found', compartment)
        end
        
        % Verify species is in compartment
        compartmentInd = find(compartmentMask);
        allSpeciesInCompartment = allSpecies(allSpeciesCompartmentInds == compartmentInd);
        speciesMask = ismember(allSpeciesInCompartment, species);
        if ~any(speciesMask)
            error('CleanSpeciesName:SpeciesNotFound', 'Species %s not found in compartment %s', species, compartment)
        end
    otherwise % Invalid names with x.x.x format or more
        error('CleanSpeciesName:TooManyParts', 'Invalid name %s', name)
end
end

function expr = cleanExpr(expr)
% Clean up rate form expression, replacing invalid characters in double-quoted
% names with valid names with underscores and no quotes.
% Inputs:
%   expr [ string ]
%       Input expression to be cleaned up
% Outputs:
%   expr [ string ]
%       Cleaned up input expression
inQuotes = '"(.*?)"'; % anything inside double quotes, no nested quotes
invalid = '[^a-zA-Z0-9_\".]'; % any invalid character

qExprs = regexp(expr, inQuotes, 'match');
nQ = length(qExprs);
for i = 1:nQ
    qExpr = qExprs{i};
    cleanedExpr = strrep(regexprep(qExpr, invalid, '_'), '"', '');
    expr = regexprep(expr, regexptranslate('escape', qExpr), cleanedExpr);
end
end
