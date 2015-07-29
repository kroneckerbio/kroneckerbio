function symbolic = simbio2symbolic(simbio, opts)
% Inputs:
%   sbml [ simbio model object ]
%   opts [ options struct ]
%       Options struct allowing the following fields:
%       .Verbose
%       .OutputsAsRules [ true | {false} ]
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
nv     = length(simbio.compartment);
vNames = cell(nv,1);
vIDs   = cell(nv,1);
v      = zeros(nv,1);
dv     = zeros(nv,1);
for i = 1:nv
    
    compartment = simbio.Compartments(i);
    
    vNames{i} = compartment.Name;
    vIDs{i}   = genUID;
    v(i)      = compartment.Capacity;
    dv(i)     = 3; % default; TODO: from units if specified
    
end

%% Species, turned into states and inputs below
nxu       = length(simbio.species);
xuNames   = cell(nxu,1);
xuIDs     = cell(nxu,1);
xu0       = zeros(nxu,1);
xuvNames  = cell(nxu,1);
isu       = false(nxu,1);
for i = 1:nxu
    
    species = simbio.Species(i);
    
    xuNames{i}  = species.Name;
    xuIDs{i}    = genUID;
    xu0(i)      = species.InitialAmount;
    xuvNames{i} = vNames{ismember(vNames, species.Parent.Name)};
    isu(i)      = species.BoundaryCondition || species.ConstantAmount;
    
end

%% States and seeds
% Make a seed for each state, equaling initial condition
nx           = sum(~isu);
xNames       = xuNames(~isu);
xIDs         = xuIDs(~isu);
xvNames      = xuvNames(~isu);
ns           = nx;
sNames       = strcat(xvNames, '_', xNames, '_0');
sNamesQuoted = sNames; % initial expression requires quotes around invalid names
for i = 1:ns
    if regexp(sNames{i}, '\W')
        sNamesQuoted{i} = ['"', sNames{i}, '"'];
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
nk     = length(simbio.Parameters);
kNames = cell(nk,1);
kIDs   = cell(nk,1);
k      = zeros(nk,1);
for i = 1:nk
    
    parameter = simbio.Parameters(i);
    
    kNames{i} = parameter.Name;
    kIDs{i}   = genUID;
    k(i)      = parameter.Value;
    
end

% Reaction-local
nr = length(simbio.Reactions);
for i = 1:nr
    
    reaction = simbio.Reactions(i);
    
    kineticLaw = reaction.kineticLaw; % Will be empty if no kinetic law parameters exist
    if ~isempty(kineticLaw)
        parameters = kineticLaw.Parameters; % Will only fetch parameters unique to this kinetic law
        nkl = length(parameters);
        for j = 1:nkl
            parameter = parameters(j);
            
            name = parameter.Name;
            id = genUID;
            value = parameter.Value;
            
            nk     = nk + 1;
            kNames = [kNames; name];
            kIDs   = [kIDs; id];
            k      = [k; value];
            
        end
    end
    
end

%% Reactions
% Note/TODO: amounts/concentrations and their effect on stoichiometry
%   For now, just using stoichiometry directly (assumed to be integers)
rNames = cell(nr,1);
rIDs   = cell(nr,1);
r      = cell(nr,3); % [reactant cells, product cells, rate expression]
for i = 1:nr
    
    reaction = simbio.Reactions(i);
    
    rNames{i} = reaction.Name;
    rIDs{i} = genUID;
    
    % Get reactant and product names
    %   Prepend compartment names
    Sindex = 1; % keep track of stoichiometry to determine how many species to add
    assert(all(mod(reaction.Stoichiometry, 1) == 0), 'simbio2symbolic:invalidStoich', 'Stoichiometries of reaction %s not all integers', rNames{i})
    
    nReactants = length(reaction.Reactants);
    reactants = cell(1,0);
    for j = 1:nReactants
        S = abs(reaction.Stoichiometry(Sindex));
        reactant = [reaction.Reactants(j).Parent.Name '.' reaction.Reactants(j).Name];
        reactants = [reactants, repmat({reactant}, [1,S])];
        Sindex = Sindex + 1;
    end
    
    nProducts = length(reaction.Products);
    products = cell(1,0);
    for j = 1:nProducts
        S = abs(reaction.Stoichiometry(Sindex));
        product = [reaction.Products(j).Parent.Name '.' reaction.Products(j).Name];
        products = [products, repmat({product}, [1,S])];
        Sindex = Sindex + 1;
    end
    assert(Sindex == length(reaction.Stoichiometry)+1, 'simbio2symbolic:missingReactantsOrProducts' ,'Missing reactants or products') % didn't get them all somehow
    
    % Reaction rates are strings containing names
    %   SimBio wraps invalid names with square brackets [] - replace with quotes
    %   If name is invalid, brackets surround only the invalid part, not
    %       the full name, i.e. cell.[C:C]. Change this to surround the full 
    %       name as kroneckerbio expects
    rate = simbioExpr2kroneckerbioExpr(reaction.ReactionRate);
    
    r(i,:) = {reactants, products, rate};
    
end

%% Rules
% Skip rules that aren't supported and give a warning
% Currently only supports repeatedAssignment and 
nz     = length(simbio.Rules);
zNames = cell(nz,1);
zIDs   = cell(nz,1);
z      = cell(nz,3); % [target, expression, type]
for i = 1:nz
    
    rule = simbio.Rules(i);
    
    zNames{i} = rule.name;
    zIDs{i} = genUID;
    
    % Special handling of different rule types
    type = rule.RuleType;
    switch type
        case 'repeatedAssignment'
            % pass
        case 'initialAssignment'
            % pass
        otherwise
            warning('simbio2symbolic:invalidRuleType', 'Rule %s with type %s not recognized. Ignoring', zNames{i}, type)
            continue
    end
    
    % Get rule expressions
    splits = regexp(rule.Rule, '=', 'split');
    assert(numel(splits) == 2, 'simbio2symbolic: rule %s had an unparsible expression', zNames{i})
    target = cleanSimBioExpr(strtrim(splits{1}));
    expression = cleanSimBioExpr(strtrim(splits{2}));
    
    z(i,:) = {target, expression, type};
end

% Clean up rules
zValid = cellfun(@isempty, z(:,3));
nz     = sum(zValid);
zNames = zNames(zValid);
zIDs   = zIDs(zValid);
z      = z(zValid,:);

%% Assemble symbolic model
if verbose; fprintf('Assembling symbolic model...'); end

symbolic = [];

symbolic.Type = 'Model.Symbolic';
symbolic.Name = simbio.name;

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
