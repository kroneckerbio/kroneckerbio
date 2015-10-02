function symbolic = analytic2symbolic(m, opts)
% Inputs:
%   m [ Model.Analytic struct ]
%   opts [ options struct ]
%       Options struct allowing the following fields:
%       .Verbose [ nonnegative integer ]
%           The level of debug information verbosity, with higher values giving
%           more output.
%
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

%% Extract model components
if verbose; fprintf('Extracting model components...'); end

%% Compartments
nv     = m.nv;
vNames = {m.Compartments.Name}';
vIDs   = {m.Compartments.ID}';
v      = [m.Compartments.Size]';
dv     = [m.Compartments.Dimension]';

%% States
nx      = m.nx;
xNames  = {m.States.Name}';
xIDs    = {m.States.ID}';
xvNames = {m.States.Compartment}';
x       = {m.States.InitialValue}';

%% Inputs
nu      = m.nu;
uNames  = {m.Inputs.Name}';
uIDs    = {m.Inputs.ID}';
uvNames = {m.Inputs.Compartment}';
u       = [m.Inputs.DefaultValue]';
% Handle blank inputs
if nu == 0
    uNames  = cell(0,1);
    uIDs    = cell(0,1);
    uvNames = cell(0,1);
    u       = zeros(0,1);
end

%% Seeds
ns      = m.ns;
sNames  = {m.Seeds.Name}';
sIDs    = {m.Seeds.ID}';
s       = [m.Seeds.Value]';
% Handle blank seeds
if ns == 0
    sNames  = cell(0,1);
    sIDs    = cell(0,1);
    s       = zeros(0,1);
end

%% Parameters
nk      = m.nk;
kNames  = {m.Parameters.Name}';
kIDs    = {m.Parameters.ID}';
k       = [m.Parameters.Value]';

%% Reactions
nr      = m.nr;
rNames  = {m.Reactions.Name}';
rIDs    = {m.Reactions.ID}';
r = cell(nr,3);
for i = 1:nr
    reaction = m.Reactions(i);
    
    % Replace cell arrays of empty double vector (no product or reactant) with
    % 0x1 cell array
    reactants = reaction.Reactants;
    if iscell(reactants) && isempty(reactants{1})
        reactants = cell(1,0);
    end
    products = reaction.Products;
    if iscell(products) && isempty(products{1})
        products = cell(1,0);
    end
    r(i,:) = {reactants, products, reaction.Rate};
end

%% Rules
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
% Handle blank rules
if nz == 0
    zNames = cell(0,1);
    zIDs   = cell(0,1);
end

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