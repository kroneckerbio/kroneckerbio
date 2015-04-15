function sbmlModel = symbolic2Sbml(symModel, opts)
%SYMBOLIC2SBML Reverse of sbml2Symbolic
% Output to libSBML format
% Note: Use generated IDs (UUIDs or sym names) for safety for some fields.
% Note: requires Java for UUID generation (TODO: remove this dependence)
% NOTE: passes isSBML_Model check in libSBML matlab binding but segfaults when
%   using OutputSBML

%% Options
% Resolve missing inputs
if nargin < 2
    opts = [];
end

% Options for displaying progress
defaultOpts.Verbose = 0;
defaultOpts.Validate = false;
defaultOpts.UseNames = false; % true uses Names as IDs - potentially dangerous
defaultOpts.SBML_level = 2;
defaultOpts.SBML_version = 1;
defaultOpts.UseConcentrations = true; % species use initial concs.; false = species use initial amounts

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

if verbose; fprintf('Converting symbolic model to SBML.\n'); end

%% Base struct with fields needed for all SBML sub-struct fields
baseStruct = struct('typecode', '', 'metaid', '', 'notes', '', 'annotation', '', ...
    'name', '', 'id', '', 'SBML_level', opts.SBML_level, 'SBML_version', opts.SBML_version);
emptyStruct = struct('typecode', {}, 'metaid', {}, 'notes', {}, 'annotation', {}, ...
    'name', {}, 'id', {}, 'level', {}, 'version', {});

%% IDs of all compartments (v), species (x, u), parameters (k, s, q)
% For comparisons when finding reaction modifier species
compartmentIds = arrayfun(@char, symModel.vSyms, 'UniformOutput', false);
speciesIDs = [arrayfun(@char, symModel.xSyms, 'UniformOutput', false); ...
    arrayfun(@char, symModel.uSyms, 'UniformOutput', false)];
parameterIDs = [arrayfun(@char, symModel.kSyms, 'UniformOutput', false); ...
    arrayfun(@char, symModel.sSyms, 'UniformOutput', false); ...
    arrayfun(@char, symModel.qSyms, 'UniformOutput', false)];

%% Main SBML struct
if opts.UseNames
    modelID = symModel.Name;
else
    modelID = genUUID;
end
modelName = symModel.Name;

sbmlModel = baseStruct;
sbmlModel.typecode = 'SBML_MODEL';
sbmlModel.notes    = 'Converted from kroneckerbio symbolic model.';
sbmlModel.id       = modelID;
sbmlModel.name     = modelName;

%% Misc fields
% kroneckerbio doesn't use these
functionDefinition = emptyStruct;
functionDefinition().math = [];
sbmlModel.functionDefinition = functionDefinition;

unitDefinition = emptyStruct;
unitDefinition().unit = [];
sbmlModel.unitDefinition = unitDefinition;

sbmlModel.time_symbol = '';
sbmlModel.delay_symbol = '';

%% Compartments
compartmentStruct = baseStruct;
compartmentStruct.typecode = 'SBML_COMPARTMENT';
compartmentStruct.units = [];
compartmentStruct.outside = [];
compartmentStruct.constant = 1;
compartmentStruct.isSetSize = 1;
compartmentStruct.isSetVolume = 1;

nv = symModel.nv;
allCompartments = [];
compartmentIDs = cell(nv,1); % for species compartment lookup
for i = 1:nv
    compartment = compartmentStruct;
    compartment.id = char(symModel.vSyms(i));
    compartment.name = symModel.vNames{i};
    compartment.spatialDimensions = symModel.dv(i);
    compartment.size = symModel.v(i);
    
    compartmentIDs{i} = compartment.id;
    
    allCompartments = [allCompartments, compartment];
end

sbmlModel.compartment = allCompartments;

%% Species
% kroneckerbio states + inputs
% TODO: difference between boundary condition and constant - all inputs assumed
%   constant for now
speciesStruct = baseStruct;
speciesStruct.typecode = 'SBML_SPECIES';
speciesStruct.isSetCharge = 0;
speciesStruct.charge = 0;
speciesStruct.substanceUnits = [];
speciesStruct.spatialSizeUnits = [];

nx = symModel.nx;
nu = symModel.nu;

allSpecies = [];

% Add states
for i = 1:nx
    species = speciesStruct;
    species.id = char(symModel.xSyms(i));
    species.name = symModel.xNames{i};
    species.compartment = compartmentIDs{symModel.vxInd(i)};
    
    if opts.UseConcentrations
        species.hasOnlySubstanceUnits = 0;
        species.isSetInitialConcentration = 1;
        species.isSetInitialAmount = 0;
        species.initialConcentration = double(symModel.x0(i));
        species.initialAmount = NaN;
    else % UseAmounts
        species.hasOnlySubstanceUnits = 1;
        species.isSetInitialConcentration = 0;
        species.isSetInitialAmount = 1;
        species.initialConcentration = NaN;
        species.initialAmount = double(symModel.x0(i));
    end
    
    species.boundaryCondition = 0;
    species.constant = 0;
    
    allSpecies = [allSpecies, species];
end

% Add inputs
for i = 1:nu
    species = speciesStruct;
    species.id = char(symModel.uSyms(i));
    species.name = symModel.uNames{i};
    species.compartment = compartmentIDs{symModel.vuInd(i)};
    
    if opts.UseConcentrations
        species.hasOnlySubstanceUnits = 0;
        species.isSetInitialConcentration = 1;
        species.isSetInitialAmount = 0;
        species.initialConcentration = double(symModel.u(i));
        species.initialAmount = NaN;
    else % UseAmounts
        species.hasOnlySubstanceUnits = 1;
        species.isSetInitialConcentration = 0;
        species.isSetInitialAmount = 1;
        species.initialConcentration = NaN;
        species.initialAmount = double(symModel.u(i));
    end
    
    species.boundaryCondition = 0;
    species.constant = 1;
    
    allSpecies = [allSpecies, species];
end

sbmlModel.species = allSpecies;

%% Parameters
% Note: regular parameters k are tested for now; s and q can be added but aren't
%   implemented yet in other code
% Assume all parameters are constant
parameterStruct = baseStruct;
parameterStruct.typecode = 'SBML_PARAMETER';
parameterStruct.units = [];
parameterStruct.constant = 1;
parameterStruct.isSetValue = 1;

emptyParameterStruct = parameterStruct; % For unused empty parameter struct in reactions
emptyParameterStruct().id = [];
emptyParameterStruct().name = [];
emptyParameterStruct().value = [];

nk = symModel.nk;
ns = symModel.ns;
nq = symModel.nq;

allParameters = [];

for i = 1:nk
    parameter = parameterStruct;
    parameter.id = char(symModel.kSyms(i));
    parameter.name = symModel.kNames{i};
    parameter.value = symModel.k(i);
    
    allParameters = [allParameters, parameter];
end

for i = 1:ns
    parameter = parameterStruct;
    parameter.id = char(symModel.sSyms(i));
    parameter.name = symModel.sNames{i};
    parameter.value = symModel.s(i);
    
    allParameters = [allParameters, parameter];
end

for i = 1:nq
    parameter = parameterStruct;
    parameter.id = char(symModel.qSyms(i));
    parameter.name = symModel.qNames{i};
    parameter.value = symModel.q(i);
    
    allParameters = [allParameters, parameter];
end

sbmlModel.parameter = allParameters;

%% Rules
% Rules are stripped from the model when converting to kroneckerbio formats
% Reaction rates contain the rules
ruleStruct = emptyStruct;
ruleStruct().formula = [];
ruleStruct().variable = [];
ruleStruct().species = [];
ruleStruct().compartment = [];
ruleStruct().units = [];

sbmlModel.rule = ruleStruct;

%% Reactions
% Rates contain rules implicitly
% Make base reactions sub-structs
reactionStruct = baseStruct;
reactionStruct.typecode = 'SBML_REACTION';
reactionStruct.reversible = 1; % SBML Level 2 default; positive/negative fluxes determine this
reactionStruct.isSetFast = 0; % kronecker has no concept of fast reactions
reactionStruct.fast = 0;

reactantStruct = baseStruct;
reactantStruct.typecode = 'SBML_SPECIES_REFERENCE';
reactantStruct.denominator = 1;
reactantStruct.stoichiometryMath = [];

productStruct = reactantStruct;

% Modifier struct missing some usual fields and different names for others (is it going to be deprecated?)
modifierStruct = []; % implied in kronecker rate expressions but explictly required in SBML
modifierStruct.typecode = 'SBML_MODIFIER_SPECIES_REFERENCE';
modifierStruct.metaid = '';
modifierStruct.notes = '';
modifierStruct.annotation = '';
modifierStruct.level = opts.SBML_level;
modifierStruct.version = opts.SBML_version;

kineticLawStruct = baseStruct;
kineticLawStruct.typecode = 'SBML_KINETIC_LAW';
kineticLawStruct.parameter = emptyParameterStruct; % local parameters not used in kronecker
kineticLawStruct.timeUnits = '';
kineticLawStruct.substanceUnits = '';

nr = symModel.nr;
allReactions = [];
for i = 1:nr
    reaction = reactionStruct;
    
    reaction.id = genUUID;
    reaction.name = symModel.rNames{i};
    
    %% Get reactants and products from states and inputs
    allReactants = [];
    allProducts = [];
    reactantIDs = {};
    productIDs = {};
    
    %% Get appropriate column of states stoichiometry matrix
    % Note that reactant stoichiometries are reported as positive in reactant struct
    species = symModel.S(:,i);
    
    % Get reactants from states negative row indices
    reactantIdxs    = find(species < 0);
    reactantStoichs = species(reactantIdxs);
    nReactants = length(reactantIdxs);
    for j = 1:nReactants
        reactantIdx = reactantIdxs(j);
        reactant = reactantStruct;
        reactant.species = char(symModel.xSyms(reactantIdx));
        if opts.UseConcentrations
            reactant.stoichiometry = -reactantStoichs(j)*symModel.v(symModel.vxInd(reactantIdx));
        else % UseAmounts
            reactant.stoichiometry = -reactantStoichs(j);
        end
        allReactants = [allReactants, reactant];
        reactantIDs = [reactantIDs; reactant.species];
    end
    
    % Get products from states positive row indices
    productIdxs    = find(species < 0);
    productStoichs = species(productIdxs);
    nProducts = length(productIdxs);
    for j = 1:nProducts
        productIdx = productIdxs(j);
        product = productStruct;
        product.species = char(symModel.xSyms(productIdx));
        if opts.UseConcentrations
            product.stoichiometry = productStoichs(j)*symModel.v(symModel.vxInd(productIdx));
        else % UseAmounts
            product.stoichiometry = productStoichs(j);
        end
        allProducts = [allProducts, product];
        productIDs = [productIDs; product.species];
    end
    
    %% Get appropriate column of inputs stoichiometry matrix
    species = symModel.Su(:,i);
    
    % Get reactants from inputs negative row indices
    reactantIdxs    = find(species < 0);
    reactantStoichs = species(reactantIdxs);
    nReactants = length(reactantIdxs);
    for j = 1:nReactants
        reactantIdx = reactantIdxs(j);
        reactant = reactantStruct;
        reactant.species = char(symModel.xSyms(reactantIdx));
        if opts.UseConcentrations
            reactant.stoichiometry = -reactantStoichs(j)*symModel.v(symModel.vuInd(reactantIdx));
        else % UseAmounts
            reactant.stoichiometry = -reactantStoichs(j);
        end
        allReactants = [allReactants, reactant];
        reactantIDs = [reactantIDs; reactant.species];
    end
    
    % Get products from inputs positive row indices
    productIdxs    = find(species < 0);
    productStoichs = species(productIdxs);
    nProducts = length(productIdxs);
    for j = 1:nProducts
        productIdx = productIdxs(j);
        product = productStruct;
        product.species = char(symModel.xSyms(productIdx));
        if opts.UseConcentrations
            product.stoichiometry = productStoichs(j)*symModel.v(symModel.vuInd(productIdx));
        else % UseAmounts
            product.stoichiometry = productStoichs(j);
        end
        allProducts = [allProducts, product];
        productIDs = [productIDs; product.species];
    end
    
    reaction.reactant = allReactants;
    reaction.product = allProducts;
    reactantAndProductIDs = [reactantIDs; productIDs];
    
    %% Kinetic Law
    kineticLaw = kineticLawStruct;
    kineticLaw.formula = char(symModel.r(i));
    kineticLaw.math = kineticLaw.formula; % should these be the same?
    
    reaction.kineticLaw = kineticLaw;
    
    %% Modifiers
    % Need to identify species present in kinetic law but not reactants or products
    vars = symvar(symModel.r(i));
    varIDs = arrayfun(@char,vars, 'UniformOutput', false)';
    modifierIDs = setdiff(varIDs, [compartmentIds; parameterIDs; reactantAndProductIDs]); % exclude compartment, param, and participating species
    
    nModifiers = length(modifierIDs);
    allModifiers = [];
    for j = 1:nModifiers
        modifier = modifierStruct;
        modifier.species = modifierIDs{j};
        allModifiers = [allModifiers, modifier];
    end
    reaction.modifier = allModifiers;
    
    allReactions = [allReactions, reaction];
end

sbmlModel.reaction = allReactions;

%% Events
% Events not implemented in models. kronecker uses them in experiments?
eventStruct = emptyStruct;
eventStruct().trigger = [];
eventStruct().delay = [];
eventStruct().timeUnits = [];
eventStruct().eventAssignment = [];

sbmlModel.event = eventStruct;

if verbose; fprintf('done.\n'); end

end


%% Helper functions
function uuid = genUUID
    uuid = strrep(char(java.util.UUID.randomUUID), '-', '_');
end