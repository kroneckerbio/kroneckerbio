function sbml = symbolic2sbml(symbolic, opts)
%SYMBOLIC2SBML Convert intermediate Model.Symbolic to libSBML Matlab SBML model
%struct
% Inputs:
%   symbolic [ Model.Symbolic struct ]
%   opts [ options struct ]
%       Options struct allowing the following fields:
%       .Verbose [ nonnegative integer {0} ]
%       .Validate [ true | {false} ]
%       .SBML_level [ positive integer {2} ]
%       .SBML_version [ positive integer {1} ]
%       .UseConcentrations [ true | {false} ]
% Outputs:
%   sbml [ libSBML Matlab model struct ]
%
% NOTE: Passing isSBML_Model check in libSBML matlab binding can still lead to
%   models segfaulting on export

%% Options
% Resolve missing inputs
if nargin < 2
    opts = [];
end

% Options for displaying progress
opts_.Verbose = 0;
opts_.Validate = false;
opts_.SBML_level = 2;
opts_.SBML_version = 1;
opts_.UseConcentrations = false; % true = species use initial concs.; false = species use initial amounts; Note/TODO: doesn't do anything now

opts = mergestruct(opts_, opts);

verbose = logical(opts.Verbose);

% Sanity check
if ~strcmpi(symbolic.Type, 'Model.Symbolic')
    error('symbolic2sbml: input model is not a Model.Symbolic')
end

if verbose; fprintf('Converting symbolic model to SBML...'); end

%% Base struct with fields needed for all SBML sub-struct fields
baseStruct = struct('typecode', '', 'metaid', '', 'notes', '', 'annotation', '', ...
    'name', '', 'id', '', 'SBML_level', opts.SBML_level, 'SBML_version', opts.SBML_version);
emptyStruct = struct('typecode', {}, 'metaid', {}, 'notes', {}, 'annotation', {}, ...
    'name', {}, 'id', {}, 'level', {}, 'version', {});

%% Main SBML struct
% TODO: annotate with version of kroneckerbio: release and commit
modelName = symbolic.Name;
modelID   = genUID;

sbml = baseStruct;
sbml.typecode = 'SBML_MODEL';
sbml.notes    = 'Converted from kroneckerbio model';
sbml.id       = modelID;
sbml.name     = modelName;

namespaces = [];
namespaces.prefix = '';
namespaces.uri = 'http://www.sbml.org/sbml/level2'; % fix this for now
sbml.namespaces = namespaces;

%% Misc fields
% kroneckerbio doesn't use these
functionDefinition = emptyStruct;
functionDefinition().math = [];
sbml.functionDefinition = functionDefinition;

unitDefinition = emptyStruct;
unitDefinition().unit = [];
sbml.unitDefinition = unitDefinition;

sbml.time_symbol = '';
sbml.delay_symbol = '';

%% Compartments
vStruct = baseStruct;
vStruct.typecode = 'SBML_COMPARTMENT';
vStruct.units = '';
vStruct.outside = '';
vStruct.constant = 1;
vStruct.isSetSize = 1;
vStruct.isSetVolume = 1;

nv = symbolic.nv;
vEntries = [];
vIDs   = symbolic.vIDs;
vNames = symbolic.vNames;
dv     = symbolic.dv;
v      = symbolic.v;
for i = 1:nv
    vEntry = vStruct;
    vEntry.id = vIDs{i};
    vEntry.name = vNames{i};
    vEntry.spatialDimensions = dv(i);
    vEntry.size = v(i);
    
    vEntries = [vEntries, vEntry];
end

sbml.compartment = vEntries;

%% Species
% kroneckerbio states + inputs
% TODO: difference between boundary condition and constant - all inputs assumed
%   constant for now
xuStruct = baseStruct;
xuStruct.typecode = 'SBML_SPECIES';
xuStruct.isSetCharge = 0;
xuStruct.charge = 0;
xuStruct.substanceUnits = '';
xuStruct.spatialSizeUnits = '';

nx = symbolic.nx;
nu = symbolic.nu;

xuEntries = [];

% Add states
xIDs    = symbolic.xIDs;
xNames  = symbolic.xNames;
xvNames = symbolic.xvNames;
x       = symbolic.x; % string number of expression of seeds
sIDs   = symbolic.sIDs;
sSyms  = sym(sIDs);
sNames = symbolic.sNames;
s      = symbolic.s;
for i = 1:nx
    xuEntry = xuStruct;
    xuEntry.id = xIDs{i};
    xuEntry.name = xNames{i};
    xuEntry.compartment = vIDs{ismember(xvNames{i}, vNames)};
    
    % Get state initial condition
    %   Invalid names should be quoted
    %   Symbolic sub seeds for values
    %   Shouldn't be anything but seeds in expression
    % TODO:Initial conditions that are more complicated should get initial assignment rules
    x0 = x{i};
    x0sym = sym(name2id(x0, sNames, sIDs));
    x0val = eval(subs(x0sym, sSyms, s));
    
    if opts.UseConcentrations
        xuEntry.hasOnlySubstanceUnits = 0;
        xuEntry.isSetInitialConcentration = 1;
        xuEntry.isSetInitialAmount = 0;
        xuEntry.initialConcentration = x0val;
        xuEntry.initialAmount = NaN;
    else % UseAmounts
        xuEntry.hasOnlySubstanceUnits = 1;
        xuEntry.isSetInitialConcentration = 0;
        xuEntry.isSetInitialAmount = 1;
        xuEntry.initialConcentration = NaN;
        xuEntry.initialAmount = x0val;
    end
    
    xuEntry.boundaryCondition = 0;
    xuEntry.constant = 0;
    
    xuEntries = [xuEntries, xuEntry];
end

% Add inputs
uIDs    = symbolic.uIDs;
uNames  = symbolic.uNames;
uvNames = symbolic.uvNames;
u       = symbolic.u; % doubles of default values
for i = 1:nu
    xuEntry = xuStruct;
    xuEntry.id = uIDs{i};
    xuEntry.name = uNames{i};
    xuEntry.compartment = vIDs{ismember(uvNames{i}, vNames)};
    
    if opts.UseConcentrations
        xuEntry.hasOnlySubstanceUnits = 0;
        xuEntry.isSetInitialConcentration = 1;
        xuEntry.isSetInitialAmount = 0;
        xuEntry.initialConcentration = u(i);
        xuEntry.initialAmount = NaN;
    else % UseAmounts
        xuEntry.hasOnlySubstanceUnits = 1;
        xuEntry.isSetInitialConcentration = 0;
        xuEntry.isSetInitialAmount = 1;
        xuEntry.initialConcentration = NaN;
        xuEntry.initialAmount = u(i);
    end
    
    xuEntry.boundaryCondition = 0;
    xuEntry.constant = 1;
    
    xuEntries = [xuEntries, xuEntry];
end

sbml.species = xuEntries;

%% Parameters
% All parameters are called k for sbml
% Assume all parameters are constant
kStruct = baseStruct;
kStruct.typecode = 'SBML_PARAMETER';
kStruct.units = '';
kStruct.constant = 1;
kStruct.isSetValue = 1;

kStructEmpty = emptyStruct; % For unused empty parameter struct in reactions
kStructEmpty().id = [];
kStructEmpty().name = [];
kStructEmpty().value = [];
kStructEmpty().units = [];
kStructEmpty().constant = [];
kStructEmpty().isSetValue = [];

nk = symbolic.nk;
ns = symbolic.ns;

kEntries = [];

kIDs   = symbolic.kIDs;
kNames = symbolic.kNames;
k      = symbolic.k;
for i = 1:nk
    kEntry = kStruct;
    kEntry.id = kIDs{i};
    kEntry.name = kNames{i};
    kEntry.value = k(i);
    
    kEntries = [kEntries, kEntry];
end

% Seeds are added but may not refer to anything if initial assignment rule isn't added
sIDs   = symbolic.sIDs;
sNames = symbolic.sNames;
s      = symbolic.s;
for i = 1:ns
    kEntry = kStruct;
    kEntry.id = sIDs{i};
    kEntry.name = sNames{i};
    kEntry.value = s(i);
    
    kEntries = [kEntries, kEntry];
end

sbml.parameter = kEntries;

%% Rules
% TODO: implement
zStruct = emptyStruct;
zStruct().formula = [];
zStruct().variable = [];
zStruct().species = [];
zStruct().compartment = [];
zStruct().units = [];

sbml.rule = zStruct;

%% Reactions
% Rates contain rules implicitly
% Make base reactions sub-structs
rStruct = baseStruct;
rStruct.typecode = 'SBML_REACTION';
rStruct.reversible = 1; % SBML Level 2 default; positive/negative fluxes determine this
rStruct.isSetFast = 0; % kronecker has no concept of fast reactions
rStruct.fast = 0;

reactantStruct = baseStruct;
reactantStruct.typecode = 'SBML_SPECIES_REFERENCE';
reactantStruct.denominator = 1;
reactantStruct.stoichiometryMath = '';

productStruct = reactantStruct;

% Modifier struct missing some usual fields and different names for others (is it going to be deprecated?)
modifierStruct = []; % implied in kronecker rate expressions but explictly required in SBML
modifierStruct.typecode = 'SBML_MODIFIER_SPECIES_REFERENCE';
modifierStruct.metaid = '';
modifierStruct.notes = '';
modifierStruct.annotation = '';
modifierStruct.level = opts.SBML_level;
modifierStruct.version = opts.SBML_version;

emptyModifierStruct = struct('typecode', {}, 'metaid', {}, 'notes', {}, 'annotation', {}, 'species', {}, 'level', {}, 'version', {});

kineticLawStruct = baseStruct;
kineticLawStruct.typecode = 'SBML_KINETIC_LAW';
kineticLawStruct.parameter = kStructEmpty; % local parameters not used in kronecker
kineticLawStruct.timeUnits = '';
kineticLawStruct.substanceUnits = '';

allReactions = [];



% Get names/ids for symbolic substitutions
xuNames  = [xNames; uNames];
xuvNames = [xvNames; uvNames];
xuIDs    = [xIDs; uIDs];
allNames = [xuNames; vNames; kNames; sNames];
allIDs   = [xuIDs; vIDs; kIDs; sIDs];

nr = symbolic.nr;
rIDs   = symbolic.rIDs;
rNames = symbolic.rNames;
r      = symbolic.r;
for i = 1:nr
    
    reaction = rStruct;
    
    reaction.id = rIDs{i};
    reaction.name = rNames{i};
    
    %% Reactants and products
    % Convert reactant and product names to IDs
    reactants = vec(r{i,1});
    products = vec(r{i,2});
    nReac = length(reactants);
    nProd = length(products);
    for j = 1:nReac
        reactants{j} = name2id(quoteInvalid(reactants{j}), xuNames, xuIDs, xuvNames);
    end
    for j = 1:nProd
        products{j} = name2id(quoteInvalid(products{j}), xuNames, xuIDs, xuvNames);
    end
    
    % Combine reactants and products according to stoichiometries for
    % multi-consumption and production
    [reactants, ~, ic] = unique(reactants);
    Sr = histc(ic, 1:numel(reactants));
    [products, ~, ic] = unique(products);
    Sp = histc(ic, 1:numel(products));
    
    nReac = length(reactants);
    nProd = length(products);
    
    allReactants = [];
    for j = 1:nReac
        reactant = reactantStruct;
        reactant.species = reactants{j};
        reactant.stoichiometry = Sr(j);
        allReactants = [allReactants, reactant];
    end
    
    allProducts = [];
    for j = 1:nProd
        product = productStruct;
        product.species = products{j};
        product.stoichiometry = Sp(j);
        allProducts = [allProducts, product];
    end
    
    reaction.reactant = allReactants;
    reaction.product = allProducts;
    
    %% Kinetic Law
    % Convert names to IDs in rate expression
    rate = r{i,3};
    rate = name2id(rate, allNames, allIDs, xuvNames);
    
    kineticLaw = kineticLawStruct;
    kineticLaw.formula = rate;
    kineticLaw.math = rate; % should these be the same?
    
    reaction.kineticLaw = kineticLaw;
    
    %% Modifiers
    % Need to identify species present in kinetic law but not reactants or products
    vars = symvar(rate);
    modifierIDs = setdiff(vars, unique([vIDs; kIDs; sIDs; reactants; products])); % exclude compartment, param, and participating species
    nModifiers = length(modifierIDs);
    
    allModifiers = [];
    if nModifiers == 0
        allModifiers = emptyModifierStruct;
    else
        for j = 1:nModifiers
            modifier = modifierStruct;
            modifier.species = modifierIDs{j};
            allModifiers = [allModifiers, modifier];
        end
    end
    reaction.modifier = allModifiers;
    
    allReactions = [allReactions, reaction];
end

sbml.reaction = allReactions;

%% Events
% Events not implemented in models. kroneckerbio uses them in experiments but
% they aren't exposed to the symbolic model intermediate
eventStruct = emptyStruct;
eventStruct().trigger = [];
eventStruct().delay = [];
eventStruct().timeUnits = [];
eventStruct().eventAssignment = [];

sbml.event = eventStruct;

if verbose; fprintf('done.\n'); end

%% Optional validation
% This doesn't seem to be very thorough...
if opts.Validate
    if verbose; fprintf('Validating SBML model struct...'); end
    [valid, message] = isSBML_Model(sbml);
    if logical(valid)
        fprintf('Model appears to be valid...')
    else
        warning('Validation failed: %s', message)
    end
    if verbose; fprintf('done.\n'); end
end

end

function expr = quoteInvalid(expr)
% Wrap expression (single identifier) in quotes if it contains invalid
% characters ('\W')
if regexp(expr, '\W') % wrap in quotes if invalid chars present in state name
    expr = ['"', expr, '"'];
end
end