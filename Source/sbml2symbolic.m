function symbolic = sbml2symbolic(sbml, opts)
% Inputs:
%   sbml [ libSBML Matlab converter SBML model struct ]
%   opts [ options struct ]
% Outputs:
%   symbolic [ Model.Symbolic struct ]

if nargin < 2
    opts = [];
end

% Default options
opts_.Verbose = 0;
opts_.UseNames = false;

opts = mergestruct(opts_, opts);

verbose = logical(opts.Verbose);

%% Extract model components
if verbose; fprintf('Extracting model components...'); end

%% Compartments
nv = length(sbml.compartment);
vNames = cell(nv,1);
vIDs   = cell(nv,1);
v      = zeros(nv,1);
dv     = zeros(nv,1);
for i = 1:nv
    
    compartment = sbml.compartment(i);
    
    if opts.UseNames
        name = compartment.id;
        id = genUID;
    else
        name = compartment.name;
        id = compartment.id;
    end
    vNames{i} = name;
    vIDs{i} = id;
    
    if compartment.isSetSize
        v(i) = compartment.size;
    else
        warning('Warning:LoadModelSbmlAnalytic:CompartmentSizeNotSet: Compartment size not set, setting default size = 1.')
        v(i) = 1;
    end
    
    dv(i) = compartment.spatialDimensions;
    
end

%% Species, turned into states and inputs below
nxu = length(sbml.species);
xuNames          = cell(nxu,1);
xuIDs            = cell(nxu,1);
xu0              = zeros(nxu,1);
xuvNames         = cell(nxu,1);
xuSubstanceUnits = false(nxu,1);
isu              = false(nxu,1);
for i = 1:nxu
    
    species = sbml.species(i);
    
    if opts.UseNames
        name = species.id;
        id = genUID;
    else
        name = species.name;
        id = species.id;
    end
    xuNames{i} = name;
    xuIDs{i} = id;
    
    % Get initial amount, converting initial concentrations if necessary
    if species.isSetInitialAmount
        xu0i = species.initialAmount;
    elseif species.isSetInitialConcentration
        xu0i = species.initialConcentration;
    else
        warning('LoadModelSbmlAnalytic:InitialConcentrationNotSet: Initial species conc. not set for %s, setting default conc. = 0.', xuIDs{i})
        xu0i = 0;
    end
    xu0(i) = xu0i;
    
    % Get species compartment
    xuvID = species.compartment;
    if opts.UseNames
        xuvID = vIDs(ismember(xuvID, vNames));
    end
    [~, vInd] = ismember(xuvID, vIDs);
    xuvNames{i} = vNames{vInd};
    
    % Species substance units in amount/true or conc./false
    xuSubstanceUnits(i) = logical(species.hasOnlySubstanceUnits);
    
    % Species is input/conc. doesn't change due to reactions, etc.
    isu(i) = species.boundaryCondition || species.constant;
    
end

%% States and seeds
% Make a seed for each state, equaling initial condition
nx = sum(~isu);
xNames  = xuNames(~isu);
xIDs    = xuIDs(~isu);
xvNames = xuvNames(~isu);
ns = nx;
sNames = strcat(xNames, '_0');
sNamesQuoted = sNames;
for i = 1:ns
    if regexp(sNamesQuoted{i}, '\W') % wrap in quotes if invalid chars present in state name
        sNamesQuoted{i} = ['"', sNamesQuoted{i}, '"'];
    end
end
sIDs   = cellfun(@(x)genUID, sNames, 'UniformOutput', false);
x      = sNamesQuoted;
s      = xu0(~isu);

%% Inputs
nu = sum(isu);
uNames  = xuNames(isu);
uIDs    = xuIDs(isu);
uvNames = xuvNames(isu);
u       = xu0(isu);

%% Parameters
% Global
nk  = length(sbml.parameter);
kNames = cell(nk,1);
kIDs   = cell(nk,1);
k      = zeros(nk,1);
for i = 1:nk
    
    parameter = sbml.parameter(i);
    
    if opts.UseNames
        name = parameter.id;
        id = genUID;
    else
        name = parameter.name;
        id = parameter.id;
    end
    kNames{i} = name;
    kIDs{i} = id;
    
    k(i) = parameter.value;
    
end

% Reaction-local
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
            
            nk     = nk + 1;
            kNames = [kNames; name];
            kIDs   = [kIDs; id];
            k      = [k; value];
        end
    end
    
end

%% Reactions
rNames = cell(nr,1);
rIDs   = cell(nr,1);
r      = cell(nr,3); % [reactant cells, product cells, rate expression]
for i = 1:nr
    
    reaction = sbml.reaction(i);
    
    if opts.UseNames
        name = reaction.id;
        id = genUID;
    else
        name = reaction.name;
        id = reaction.id;
    end
    rNames{i} = name;
    rIDs{i} = id;
    
    % Get reactant and product names
    reactants = getSpeciesNames(reaction.reactant);
    products = getSpeciesNames(reaction.product);
    
    % Reaction rate - keep as IDs
    %   UUID-like IDs won't have any problems when IDs are attempted to be
    %   subbed in here again when the model is parsed
    % Formula uses "power" which is recognized by Matlab's symbolic toolbox
    %   Math uses "pow" which is not
    rate = reaction.kineticLaw.formula;
    
    r(i,:) = {reactants, products, rate};
    
end

if verbose; fprintf('done.\n'); end

    function speciesNames = getSpeciesNames(species)
        % Get reactant and product species' names from their IDs. Needed as
        % kroneckerbio analytic model parsing uses species names in reactants,
        % not IDs. Species with stoichiometries > 1 are repeated.
        % Inputs:
        %   species [ SBML reaction.reactant|product struct ]
        %       libSBML Matlab loader reaction.* struct
        % Outputs:
        %   speciesNames [ 1 x nSpecies cell array of strings ]
        %       Cell array of compartment.name strings with species repeated
        %       according to stoichiometry.
        nSpecies = numel(species);
        speciesNames = cell(1,0);
        for iSpecies = 1:nSpecies
            xuID = species(iSpecies).species;
            stoich = species(iSpecies).stoichiometry;
            if opts.UseNames
                xuID = xuIDs(ismember(xuID, xuNames));
            end
            [~, xuInd] = ismember(xuID, xuIDs);
            xuName = strcat(xuvNames{xuInd}, '.', xuNames{xuInd});
            speciesNames = [speciesNames, repmat({xuName},1,stoich)];
        end
    end

%% Switch from old names to new IDs if requested
if opts.UseNames
    if verbose; fprintf('Converting IDs to Names...'); end
    
    % Assemble mapping
    ids = [vIDs; xuIDs; kIDs; sIDs];
    names = [vNames; xuNames; kNames; sNames];
    
    ids = sym(ids);
    names = sym(names); % safe because opts.UseNames assumes names are valid IDs
    
    % Get symbolic expressions for rates
    rates = r(:,3);
    rates = sym(rates);
    
    % Perform substitutions
    rates = subs(rates, names, ids);
    
    % Convert rate expressions back to strings
    for i = 1:length(rates)
        r{i,3} = char(rates(i));
    end
    if verbose; fprintf('done.\n'); end
end

%% Rules
nz     = 0;
zNames = cell(nz,1);
zIDs   = cell(nz,1);
z      = cell(nz,3); % [target, expression, type]

% Repeated assignment
if isfield(sbml, 'rule') && ~isempty(sbml.rule)
    for i = 1:length(sbml.rule)
        rule = sbml.rule(i);
        name = rule.name; % doesn't map to anything?
        id = genUID; % libSBML Matlab loader doesn't generate rule IDs
        target = rule.variable;
        expression = rule.formula;
        
        nz     = nz + 1;
        zNames = [zNames; name];
        zIDs   = [zIDs; id];
        z      = [z; {target, expression, 'repeated assignment'}];
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
        
        nz     = nz + 1;
        zNames = [zNames; name];
        zIDs   = [zIDs; id];
        z      = [z; {target, expression, 'initial assignment'}];
    end
end

%% Assemble symbolic model
if verbose; fprintf('Assembling symbolic model...'); end
if opts.UseNames
    name = sbml.id;
else
    name = sbml.name;
end

symbolic = [];

symbolic.Type = 'Model.Symbolic';
symbolic.Name = name;

symbolic.nv     = nv;
symbolic.vNames = vNames;
symbolic.vIDs   = vIDs;
symbolic.v      = v; % doubles
symbolic.dv     = dv;

symbolic.nx      = nx;
symbolic.xNames  = xNames;
symbolic.xIDs    = xIDs;
symbolic.xvNames = xvNames;
symbolic.x       = x; % cell array of strings

symbolic.nu      = nu;
symbolic.uNames  = uNames;
symbolic.uIDs    = uIDs;
symbolic.uvNames = uvNames;
symbolic.u       = u; % doubles

symbolic.ns     = ns;
symbolic.sNames = sNames;
symbolic.sIDs   = sIDs;
symbolic.s      = s; % doubles

symbolic.nk     = nk;
symbolic.kNames = kNames;
symbolic.kIDs   = kIDs;
symbolic.k      = k; % doubles

symbolic.nr     = nr;
symbolic.rNames = rNames;
symbolic.rIDs   = rIDs;
symbolic.r      = r; % cell matrix

symbolic.nz     = nz;
symbolic.zNames = zNames;
symbolic.zIDs   = zIDs;
symbolic.z      = z; % cell matrix

if verbose; fprintf('done.\n'); end

end