function m = LoadModelSbmlAnalytic(model, opts)
%LoadModelSbmlAnalytic Import SBML model and covert to kroneckerbio analytic
%model. Modify the model and add outputs after calling this.
%
%   m = LoadModelSbmlAnalytic(model, yNames, yMembers, yValues, opts)
%
%   Inputs
%   model: [ Simbiology model scalar ]
%       A Simbiology model object
%   opts: [ options struct scalar {} ]
%       .Verbose [ logical scalar {false} ]
%       	Print progress to command window
%       .Validate [ logical scalar {false} ]
%           Whether to use libSBML's model validation tool
%       .UseNames [logical scalar {false} ]
%           Whether to convert SBML IDs to Names and autogenerate new IDs
%           Use this when the supplied SBML model uses "nice" names as IDs
%
%   Outputs
%   m: [ Kronecker model struct scalar ]
%       An analytic pseudo-Kronecker model
%
%   Notes: 
%   - Inputs are any species that have "constant" or "boundaryCondition" set.
%   - Seeds are generated from any parameter that appears in an InitialAssignment.
%
%   Limitations:
%   Not all Simbiology features are compatible with this converter. This
%   function ignores events, some rules, and functions in the model.

%% Clean up inputs
if nargin < 2
    opts = [];
end

% Default options
opts_.Verbose = 0;
opts_.Validate = false;
opts_.UseNames = false;

opts = mergestruct(opts_, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

%% Call libSBML to import SBML model
if verbose; fprintf('Convert SBML model using libSBML...'); end

sbml = TranslateSBML(model, double(opts.Validate), opts.Verbose);

if verbose; fprintf('done.\n'); end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Part 1: Extracting the Model Variables %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose; fprintf('Extracting model components...'); end

%% Model name and initialization
if opts.UseNames
    name = sbml.id;
else
    name = sbml.name;
end

m = InitializeModelAnalytic(name);

%% Compartments
nv = length(sbml.compartment);
for i = 1:nv
    
    compartment = sbml.compartment(i);
    
    if opts.UseNames
        name = compartment.id;
        id = genUID;
    else
        name = compartment.name;
        id = compartment.id;
    end
    
    dv = compartment.spatialDimensions;
    
    if compartment.isSetSize
        vNames = compartment.size;
    else
        warning('Warning:LoadModelSbmlAnalytic:CompartmentSizeNotSet: Compartment size not set, setting default size = 1.')
        vNames = 1;
    end
    
    m = AddCompartment(m, name, dv, vNames, id);
    
end

vNames = {m.add.Compartments.Name}';
vIDs   = {m.add.Compartments.ID}';
v      = [m.add.Compartments.Size]';

%% Species
nxu = length(sbml.species);
xuIDs   = cell(nxu,1);
xuNames = cell(nxu,1);
xuCompartments = cell(nxu,1);
xuSubstanceUnits = false(nxu,1);
xuvInd = zeros(nxu,1);
sIDs = cell(0,1);
sNames = cell(0,1);
for i = 1:nxu
    
    species = sbml.species(i);
    
    if opts.UseNames
        name = species.id;
        id = genUID;
    else
        name = species.name;
        id = species.id;
    end
    
    % Get species compartment
    xuvID = species.compartment;
    if opts.UseNames
        xuvID = vIDs(ismember(xuvID, vNames));
    end
    [~, vInd] = ismember(xuvID, vIDs);
    xuvName   = vNames{vInd};
    xuv       = v(vInd);
    
    % Get initial amount, converting initial concentrations if necessary
    if species.isSetInitialAmount
        xu0 = species.initialAmount;
    elseif species.isSetInitialConcentration
        xu0 = species.initialConcentration;
    else
        warning('LoadModelSbmlAnalytic:InitialConcentrationNotSet: Initial species conc. not set for %s, setting default conc. = 0.', xuIDs{i})
        xu0 = 0;
    end
    
    % Species substance units in amount/true or conc./false
    xuSubstanceUnits(i) = logical(species.hasOnlySubstanceUnits);
    
    % Species is input/conc. doesn't change due to reactions, etc.
    isu = species.boundaryCondition || species.constant;
    if isu
        m = AddInput(m, name, xuvName, xu0, id);
    else % x
        % Make a seed for each state's initial condition
        sName = [name '_0']; % hopefully unique
        sExpr = sName;
        if regexp(sName, '\W') % wrap in quotes if invalid chars present in state name
            sExpr = ['"', sExpr, '"'];
        end
        sID = genUID;
        m = AddSeed(m, sName, xu0, sID);
        m = AddState(m, name, xuvName, sExpr, id);
        sIDs   = [sIDs;   sID];
        sNames = [sNames; sName];
    end
    
    xuIDs{i} = id;
    xuNames{i} = name;
    xuCompartments{i} = xuvName;
    xuvInd(i) = vInd;
end

%% Parameters
% Global
nk  = length(sbml.parameter);
for i = 1:nk
    
    parameter = sbml.parameter(i);
    
    if opts.UseNames
        name = parameter.id;
        id = genUID;
    else
        name = parameter.name;
        id = parameter.id;
    end
    
    value = parameter.value;
    
    m = AddParameter(m, name, value, id);
    
end

% Local to reactions
% TODO: prefix reaction-local parameters and modify reaction rate in case of
% local vs global name clashes
nr = length(sbml.reaction);
for i = 1:nr
    
    reaction = sbml.reaction(i);
    
    kineticLaw = reaction.kineticLaw; % Will be empty if no kinetic law parameters exist
    if ~isempty(kineticLaw)
        parameters = kineticLaw.parameter; % Will only fetch parameters unique to this kinetic law
        nkl = length(parameters);
        for j = 1:nkl
            parameter = parameters(j);
            
            if opts.UseNames
                name = parameter.id;
                id = genUID;
            else
                name = parameter.name;
                id = parameter.id;
            end
            
            value = parameter.value;
            
            m = AddParameter(m, name, value, id);
        end
    end
    
end

kNames = {m.add.Parameters.Name}';
kIDs   = {m.add.Parameters.ID}';

%% Reactions, local parameters
% Note: ignores reversible flag - not sure where the reverse rate is specified
nr = length(sbml.reaction);
for i = 1:nr
    
    reaction = sbml.reaction(i);
    
    if opts.UseNames
        name = reaction.id;
        id = genUID;
    else
        name = reaction.name;
        id = reaction.id;
    end
    
    % Get reactant and product names
    [r1,r2] = processSpecies(reaction.reactant);
    [p1,p2] = processSpecies(reaction.product);
    
    % Reaction rate - keep as IDs
    %   UUID-like IDs won't have any problems when IDs are attempted to be
    %   subbed in here again when the model is parsed
    kineticLaw = reaction.kineticLaw;
    rate = kineticLaw.math;
    
    m = AddReaction(m, name, r1, r2, p1, p2, rate, [], [], id);
    
end

    function [r1Name, r2Name] = processSpecies(species)
        % Get reactant and product species by ID in a reaction
        % Note: only accepts up to 2 reactants and products
        % TODO: make this as flexible as SBML, with stoichiometries of arbitrary
        % numbers of species in m.Reactions.Stoichiometry or as another field of
        % m.Reactions.Reactants/Products
        nSpecies = numel(species);
        switch nSpecies % rxn order
            case 0
                r1Name = [];
                r2Name = [];
            case 1
                r1ID = species(1).species;
                if opts.UseNames
                    r1ID = xuIDs(ismember(r1ID, xuNames));
                end
                [~, r1Ind] = ismember(r1ID, xuIDs);
                r1Name = strcat(xuCompartments{r1Ind}, '.', xuNames{r1Ind});
                switch species(1).stoichiometry;
                    case 1
                        r2Name = [];
                    case 2
                        r2Name = r1Name;
                    otherwise
                        error('LoadModelSbmlAnalytic:processSpecies: reaction has species with stoichiometry ~= 1 or 2')
                end
            case 2
                r1ID = species(1).species;
                if opts.UseNames
                    r1ID = xuIDs(ismember(r1ID, xuNames));
                end
                [~, r1Ind] = ismember(r1ID, xuIDs);
                r1Name = strcat(xuCompartments{r1Ind}, '.', xuNames{r1Ind});
                
                r2ID = species(2).species;
                if opts.UseNames
                    r2ID = xuIDs(ismember(r2ID, xuNames));
                end
                [~, r2Ind] = ismember(r2ID, xuIDs);
                r2Name = strcat(xuCompartments{r2Ind}, '.', xuNames{r2Ind});
            otherwise
                error('LoadModelSbmlAnalytic:processSpecies: reaction has > 2 reactants or products, exiting')
        end
    end

%% Switch from old names to new IDs if requested
if opts.UseNames
    % Assemble mapping
    ids = [vIDs; xuIDs; kIDs; sIDs];
    ids(cellfun('isempty', ids)) = [];
    ids = sym(ids);
    names = [vNames; xuNames; kNames; sNames];
    names(cellfun('isempty', names)) = [];
    names = sym(names);
    
    % Get symbolic expressions for rates
    rates = {m.add.Reactions.Rate}';
    rates(cellfun('isempty', rates)) = [];
    ratesSym = sym(rates);
    
    % Perform substitutions
    ratesSym = subs(ratesSym, names, ids);
    
    % Convert rate expressions back to strings
    for i = 1:length(rates)
        m.add.Reactions(i).Rate = char(ratesSym(i));
    end
end

if verbose; fprintf('done.\n'); end

%% Rules
% Repeated assignment
if isfield(sbml, 'rule') && ~isempty(sbml.rule)
    for i = 1:length(sbml.rule)
        rule = sbml.rule(i);
        name = rule.name; % doesn't map to anything?
        id = genUID; % libSBML Matlab loader doesn't generate rule IDs
        target = rule.variable;
        expression = rule.formula;
        m = AddRuleAnalytic(m, name, target, expression, 'repeated assignment', id);
    end
end

% Initial assignment
if isfield(sbml, 'initialAssignment') && ~isempty(sbml.initialAssignment)
    for i = 1:length(sbml.rule)
        rule = sbml.initialAssignment(i);
        name = rule.symbol; % no name field
        id = genUID; % libSBML Matlab loader doesn't generate rule IDs
        target = rule.symbol;
        expression = rule.math;
        m = AddRuleAnalytic(m, name, target, expression, 'initial assignment', id);
    end
end

%% Done
if verbose; fprintf('SBML model loaded.\n'); end
end
