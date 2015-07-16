function symbolic = analytic2symbolic(m, opts)
% Inputs:
%   m [ Model.Analytic struct ]
%   opts [ options struct ]
%       Options struct allowing the following fields:
%       .Verbose
%       .Finalized [ {true} | false ]
%           Whether to only use model components from a finalized model. True
%           means m.* components are used. False means m.* and m.add.*
%           components are used.
% Outputs:
%   symbolic [ Model.Symbolic struct ]

if nargin < 2
    opts = [];
end

% Default options
opts_.Verbose = 0;
opts_.Finalized = false;

opts = mergestruct(opts_, opts);

verbose = logical(opts.Verbose);

%% Copy m.add.* fields in non-finalized models
if ~opts.Finalized
    if verbose; fprintf('Copying nonfinalized componenets...'); end
    m.Compartments = combineComponents(m.Compartments, m.add.Compartments, opts.Verbose);
    m.Parameters = combineComponents(m.Parameters, m.add.Parameters, opts.Verbose);
    m.Seeds = combineComponents(m.Seeds, m.add.Seeds, opts.Verbose);
    m.States = combineComponents(m.States, m.add.States, opts.Verbose);
    m.Inputs = combineComponents(m.Inputs, m.add.Inputs, opts.Verbose);
    m.Reactions = combineComponents(m.Reactions, m.add.Reactions, opts.Verbose);
    m.Outputs = combineComponents(m.Outputs, m.add.Outputs, opts.Verbose);
    m.Rules = combineComponents(m.Rules, m.add.Rules, opts.Verbose);
    m.nv = length(m.Compartments);
    m.nk = length(m.Parameters);
    m.ns = length(m.Seeds);
    m.nx = length(m.States);
    m.nu = length(m.Inputs);
    m.nr = length(m.Reactions);
    m.ny = length(m.Outputs);
    m.nz = length(m.Rules);
    if verbose; fprintf('done.\n'); end
end


%% Extract model components
if verbose; fprintf('Extracting model components...'); end

% Compartments
nv     = m.nv;
vNames = {m.Compartments.Name}';
vIDs   = {m.Compartments.ID}';
v      = [m.Compartments.Size]';
dv     = [m.Compartments.Dimension]';

% States
nx      = m.nx;
xNames  = {m.States.Name}';
xIDs    = {m.States.ID}';
xvNames = {m.States.Compartment}';
x       = {m.States.InitialValue}';

% Inputs
nu      = m.nu;
uNames  = {m.Inputs.Name}';
uIDs    = {m.Inputs.ID}';
uvNames = {m.Inputs.Compartment}';
u       = [m.Inputs.DefaultValue]';

% Seeds
ns      = m.ns;
sNames  = {m.Seeds.Name}';
sIDs    = {m.Seeds.ID}';
s       = [m.Seeds.Value]';

% Parameters
nk      = m.nk;
kNames  = {m.Parameters.Name}';
kIDs    = {m.Parameters.ID}';
k       = [m.Parameters.Value]';

% Reactions
nr      = m.nr;
rNames  = {m.Reactions.Name}';
rIDs    = {m.Reactions.ID}';
r = cell(nr,3);
for i = 1:nr
    reaction = m.Reactions(i);
    
    % Replace cell arrays of empty double vector (no product or reactant) with
    % empty spot
    reactants = reaction.Reactants;
    if iscell(reactants) && isempty(reactants{1})
        reactants = [];
    end
    products = reaction.Products;
    if iscell(products) && isempty(products{1})
        products = [];
    end
    r(i,:) = {reactants, products, reaction.Rate};
end

% Rules
nz      = m.nz;
zNames  = {m.Rules.Name}';
zIDs    = {m.Rules.ID}';
z = cell(nz,3);
for i = 1:nz
    rule = m.Rules(i);
    
    % Check rule type validity
    switch rule.Type
        case 'repeated assignment'
            % pass
        case 'initial assignment'
            % pass
        otherwise
            warning('analytic2symbolic: rule %s of type %s not recognized. Ignoring.', rule.Name, rule.Type)
    end
    
    z(i,:) = {rule.Target, rule.Expression, rule.Type};
end
% Remove invalid, empty rule spots
zInvalid = cellfun(@isempty,z(:,3));
nz = nz - sum(zInvalid);
zNames(zInvalid) = [];
zIDs(zInvalid)   = [];
z(zInvalid,:)    = []; 

if verbose; fprintf('done.\n'); end

%% Assemble symbolic model
if verbose; fprintf('Assembling symbolic model...'); end

symbolic = [];

symbolic.Type = 'Model.Symbolic';
symbolic.Name = m.Name;

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