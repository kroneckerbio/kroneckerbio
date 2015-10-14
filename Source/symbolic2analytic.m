function m = symbolic2analytic(symbolic, opts)
% Inputs:
%   symbolic [ Model.Symbolic struct ]
%   opts [ options struct ]
% Outputs:
%   m [ Model.Analytic struct ]

if nargin < 2
    opts = [];
end

% Default options
opts_.Verbose = 0;

opts = mergestruct(opts_, opts);

verbose = logical(opts.Verbose);

% Sanity check
if ~strcmpi(symbolic.Type, 'Model.Symbolic')
    error('symbolic2analytic: input model is not a Model.Symbolic')
end

%% Add symbolic components to model
if verbose; fprintf('Adding components to analytic model...'); end

% Initialize model
m = InitializeModelAnalytic(symbolic.Name);

% Compartments
nv     = symbolic.nv;
vNames = symbolic.vNames;
v      = symbolic.v; % doubles
dv     = symbolic.dv;
for i = 1:nv
    m = AddCompartment(m, vNames{i}, dv(i), v(i));
end

% States
nx      = symbolic.nx;
xNames  = symbolic.xNames;
xvNames = symbolic.xvNames;
x       = symbolic.x; % cell array of strings
for i = 1:nx
    m = AddState(m, xNames{i}, xvNames{i}, x{i});
end

% Inputs
nu      = symbolic.nu;
uNames  = symbolic.uNames;
uvNames = symbolic.uvNames;
u       = symbolic.u; % doubles
for i = 1:nu
    m = AddInput(m, uNames{i}, uvNames{i}, u(i));
end

% Seeds
ns      = symbolic.ns;
sNames  = symbolic.sNames;
s       = symbolic.s; % doubles
for i = 1:ns
    m = AddSeed(m, sNames{i}, s(i));
end

% Parameters
nk     = symbolic.nk;
kNames = symbolic.kNames;
k      = symbolic.k; % doubles
for i = 1:nk
    m = AddParameter(m, kNames{i}, k(i));
end

% Reactions
nr     = symbolic.nr;
rNames = symbolic.rNames;
r      = symbolic.r; % cell matrix
for i = 1:nr
    m = AddReaction(m, rNames{i}, r{i,1}, r{i,2}, r{i,3}, []);
end

% Rules
nz     = symbolic.nz;
zNames = symbolic.zNames;
z      = symbolic.z; % cell matrix
for i = 1:nz
    m = AddRuleAnalytic(m, zNames{i}, z{i,1}, z{i,2}, z{i,3});
end

if verbose; fprintf('done.\n'); end
