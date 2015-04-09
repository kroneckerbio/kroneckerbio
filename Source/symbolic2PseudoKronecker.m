function m = symbolic2PseudoKronecker(SymModel, opts)
%symbolic2PseudoKronecker converts a symbolic model into a pseudo-kronecker
%   model, which interacts with the Kronecker Bio toolbox much like a
%   Kronecker model.
% 
%   m = symbolic2PseudoKronecker(SymModel, opts)
% 
%   Inputs
%   SymModel: [ symbolic model scalar ]
%       A symbolic model
%   opts: [ options struct scalar ]
%       Optional
%       .Order [ 0 | 1 | {2} | 3 ]
%           Determines how deep the derivatives should be taken with
%           respect to x and p. Each level increases the cost
%           exponentially, but increases the number of Kronecker Bio
%           functions that can be run on the model.
%       .VolumeToParameter [ true | {false} ]
%           Should the volume of the compartments be converted to
%           parameters or converted to immutable constants
%       .Verbose [ nonnegative integer scalar {1} ]
%           Bigger number displays more progress information
%
%   Outputs
%   m: [ psuedo-kronecker model scalar ]
%       The useable form of the model

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Clean up inputs
if nargin < 2
    opts = [];
end

% Default options
defaultOpts.Order             = 2;
defaultOpts.VolumeToParameter = false;
defaultOpts.Verbose           = 0;

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

order = opts.Order;

%% Extract symbolic values
if isfield(SymModel, 'Name')
    name = SymModel.Name;
else
    name = '';
end

kSyms   = SymModel.kSyms;
k       = SymModel.k;
nk      = numel(k);

if isfield(SymModel, 'kNames')
    kNames = SymModel.kNames;
else
    kNames = cell(nk,1);
    for ik = 1:nk
        kNames{ik} = char(kSyms(ik));
    end
end

sSyms   = SymModel.sSyms;
s       = SymModel.s;
ns      = numel(s);

if isfield(SymModel, 'sNames')
    sNames = SymModel.sNames;
else
    sNames = cell(ns,1);
    for is = 1:ns
        sNames{is} = char(sSyms(is));
    end
end

uSyms   = SymModel.uSyms;
u       = SymModel.u;
nu      = numel(u);

if isfield(SymModel, 'uNames')
    uNames = SymModel.uNames;
else
    uNames = cell(nu,1);
    for iu = 1:nu
        uNames{iu} = char(uSyms(iu));
    end
end

xSyms   = SymModel.xSyms;
x0      = SymModel.x0;
nx      = numel(x0);

if isfield(SymModel, 'xNames')
    xNames = SymModel.xNames;
else
    xNames = cell(nx,1);
    for ix = 1:nx
        xNames{ix} = char(xSyms(ix));
    end
end

if isfield(SymModel, 'v')
    vSyms   = SymModel.vSyms;
    vNames  = SymModel.vNames;
    dv      = SymModel.dv;
    v       = SymModel.v;
    nv      = numel(v);
    vuInd   = SymModel.vuInd;
    vxInd   = SymModel.vxInd;
else
    vSyms   = sym('co1x');
    vNames  = {'v'};
    dv      = 3;
    v       = 1;
    nv      = 1;
    vuInd   = ones(nu,1);
    vxInd   = ones(nx,1);
end

if isfield(SymModel, 'f')
    f = SymModel.f;
    if isfield(SymModel, 'r')
        r = SymModel.r;
        S = SymModel.S;
        nr = numel(r);
        
        if isfield(SymModel, 'rNames')
            rNames = SymModel.rNames;
        else
            rNames = repmat({''}, [nr,1]);
        end
    else
        r = sym(zeros(0,1));
        S = zeros(nx,0);
        nr = 0;
        rNames = cell(0,1);
    end
else
    % If there is no f, then r and S must be supplied
    r = SymModel.r;
    S = SymModel.S;
    nr = numel(r);
    
    if isfield(SymModel, 'rNames')
        rNames = SymModel.rNames;
    else
        rNames = repmat({''}, [nr,1]);
    end

    f = S*r;
end

y       = SymModel.y;
ny      = numel(y);

if isfield(SymModel, 'yNames')
    yNames = SymModel.yNames;
elseif isfield(SymModel, 'ySyms')
    yNames = cell(ny,1);
    for iy = 1:ny
        yNames{iy} = char(SymModel.ySyms(iy));
    end
else
    % A reasonable name for each y is not strictly necessary
    yNames = cell(ny,1);
    for iy = 1:ny
        yNames{iy} = char(y(iy));
    end
end

% Convert compartment volumes to parameters or constants
if opts.VolumeToParameter
    nk = nk + nv;
    kNames = [kNames; vNames];
    k = [k; v];
    kSyms = [kSyms; vSyms];
else% ~VolumeToParameter
    r = subs(r, vSyms, v);
    f = subs(f, vSyms, v);
end

% Sanitize symbols
old_symbols = [kSyms; sSyms; xSyms; uSyms];
new_symbols = sym(zeros(nk+ns+nx+nu,1));
for i = 1:nk+ns+nx+nu
    new_symbols(i) = sym(sprintf('sy%dx', i));
end

kSyms = subs(kSyms, old_symbols, new_symbols);
sSyms = subs(sSyms, old_symbols, new_symbols);
xSyms = subs(xSyms, old_symbols, new_symbols);
uSyms = subs(uSyms, old_symbols, new_symbols);
f = subs(f, old_symbols, new_symbols);
r = subs(r, old_symbols, new_symbols);
x0 = subs(x0, old_symbols, new_symbols);
y =  subs(y, old_symbols, new_symbols);

% String representations that will be useful
vStrs = cell(nv,1);
for iv = 1:nv
    vStrs{iv} = char(vSyms(iv));
end

kStrs = cell(nk,1);
for ik = 1:nk
    kStrs{ik} = char(kSyms(ik));
end

sStrs = cell(ns,1);
for is = 1:ns
    sStrs{is} = char(sSyms(is));
end

uStrs = cell(nu,1);
uNamesFull = cell(nu,1);
for iu = 1:nu
    uStrs{iu} = char(uSyms(iu));
    uNamesFull{iu} = [vStrs{vuInd(iu)} '.' uStrs{iu}];
end

xStrs = cell(nx,1);
xNamesFull = cell(nx,1);
for ix = 1:nx
    xStrs{ix} = char(xSyms(ix));
    xNamesFull{ix} = [vStrs{vxInd(ix)} '.' xStrs{ix}];
end

% Generate derivatives of desired order
if order >= 1
    % Gradient of r with respect to x
    drdx = jacobian(r, xSyms);
    
    % Gradient of r with respect to u
    drdu = jacobian(r, uSyms);
    
    % Gradient of r with respect to k
    drdk = jacobian(r, kSyms);
    
    % Gradient of f with respect to x
    dfdx = S*drdx;
    
    % Gradient of f with respect to u
    dfdu = S*drdu;
    
    % Gradient of f with respect to k
    dfdk = S*drdk;
    
    % Gradient of y with respect to x
    dydx = jacobian(y, xSyms);
    
    % Gradient of y with respect to u
    dydu = jacobian(y, uSyms);
else
    drdx = '';
    drdk = '';
    dfdx = '';
    dfdk = '';
end

if order >= 2
    % Gradient of drdx with respect to x
    d2rdx2 = jacobian(vec(drdx), xSyms);
    
    % Gradient of drdu with respect to u
    d2rdu2 = jacobian(vec(drdu), uSyms);
    
    % Gradient of drdu with respect to x
    d2rdxdu = jacobian(vec(drdu), xSyms);
    
    % Gradient of drdx with respect to u
    d2rdudx = jacobian(vec(drdx), uSyms);
    
    % Gradient of drdk with respect to k
    d2rdk2 = jacobian(vec(drdk), kSyms);
    
    % Gradient of drdx with respect to k
    d2rdkdx = jacobian(vec(drdx), kSyms);
    
    % Gradient of drdu with respect to k
    d2rdkdu = jacobian(vec(drdu), kSyms);
    
    % Gradient of drdk with respect to x
    d2rdxdk = jacobian(vec(drdk), xSyms);
    
    % Gradient of drdk with respect to u
    d2rdudk = jacobian(vec(drdk), uSyms);
    
    % Gradient of dfdx with respect to x
    d2fdx2 = S*reshape(d2rdx2, nr,nx*nx);
    d2fdx2 = reshape(d2fdx2, nx*nx,nx);
    
    % Gradient of dfdu with respect to u
    d2fdu2 = S*reshape(d2rdu2, nr,nu*nu);
    d2fdu2 = reshape(d2fdu2, nx*nu,nu);
    
    % Gradient of dfdu with respect to x
    d2fdxdu = S*reshape(d2rdxdu, nr,nu*nx);
    d2fdxdu = reshape(d2fdxdu, nx*nu,nx);
    
    % Gradient of dfdx with respect to u
    d2fdudx = S*reshape(d2rdudx, nr,nx*nu);
    d2fdudx = reshape(d2fdudx, nx*nx,nu);
    
    % Gradient of dfdk with respect to k
    d2fdk2 = S*reshape(d2rdk2, nr,nk*nk);
    d2fdk2 = reshape(d2fdk2, nx*nk,nk);
    
    % Gradient of dfdx with respect to k
    d2fdkdx = S*reshape(d2rdkdx, nr,nx*nk);
    d2fdkdx = reshape(d2fdkdx, nx*nx,nk);
    
    % Gradient of dfdu with respect to k
    d2fdkdu = S*reshape(d2rdkdu, nr,nu*nk);
    d2fdkdu = reshape(d2fdkdu, nx*nu,nk);
    
    % Gradient of dfdk with respect to x
    d2fdxdk = S*reshape(d2rdxdk, nr,nk*nx);
    d2fdxdk = reshape(d2fdxdk, nx*nk,nx);

    % Gradient of dfdk with respect to u
    d2fdudk = S*reshape(d2rdudk, nr,nk*nu);
    d2fdudk = reshape(d2fdudk, nx*nk,nu);
    
    % Gradient of dydx with respect to x
    d2ydx2 = jacobian(vec(dydx), xSyms);
    
    % Gradient of dydu with respect to u
    d2ydu2 = jacobian(vec(dydu), uSyms);
    
    % Gradient of dydu with respect to x
    d2ydxdu = jacobian(vec(dydu), xSyms);
    
    % Gradient of dydx with respect to u
    d2ydudx = jacobian(vec(dydx), uSyms);
else
    d2rdx2  = '';
    d2rdu2  = '';
    d2rdxdu = '';
    d2rdudx = '';
    d2rdk2  = '';
    d2rdkdx = '';
    d2rdkdu = '';
    d2rdxdk = '';
    d2rdudk = '';
    d2fdx2  = '';
    d2fdu2  = '';
    d2fdxdu = '';
    d2fdudx = '';
    d2fdk2  = '';
    d2fdkdx = '';
    d2fdkdu = '';
    d2fdxdk = '';
    d2fdudk = '';
end

if order >= 3
    %Gradient of d2rdx2 with respect to x
    d3rdx3 = jacobian(vec(d2rdx2), xSyms);
    
    %Gradient of d2rdx2 with respect to k
    d3rdkdx2 = jacobian(vec(d2rdx2), kSyms);
    
    %Gradient of d2fdx2 with respect to x
    d3fdx3 = S*reshape(d3rdx3, nr,nx*nx*nx);
    d3fdx3 = reshape(d3fdx3, xn*nx*nx,nx);
    
    %Gradient of d2fdx2 with respect to k
    d3fdkdx2 = S*reshape(d3rdkdx2, nr,nx*nx*nk);
    d3fdkdx2 = reshape(d3fdkdx2, nx*nx*nx,nk);
else
    d3rdx3   = '';
    d3rdkdx2 = '';
    d3fdx3   = '';
    d3fdkdx2 = '';
end

%% Extract seed information
dx0ds = double(jacobian(x0, sSyms));
x0c = double(x0 - dx0ds*sSyms);

initial_values = cell(nx,1);
for ix = 1:nx
    % Add constant value first
    if x0c(ix)
        values_i = {'', x0c(ix)};
    else
        values_i = cell(0,2);
    end
    
    % Append seed parameters
    nonzero_seed_indexes = logical(dx0ds(ix,:));
    values_i = [values_i; 
                vec(sStrs(nonzero_seed_indexes)), vec(num2cell(dx0ds(ix,nonzero_seed_indexes)))];
    
    % Store values
    initial_values{ix} = values_i;
end

%% Replace symbolic names with systematic
if verbose; fprintf('Converting symbolics to strings...\n'); end
% Convert the symbolics into strings
f        = symbolic2string('f', nx);
dfdx     = symbolic2string('dfdx', nx,nx);
dfdu     = symbolic2string('dfdu', nx,nu);
dfdk     = symbolic2string('dfdk', nx,nk);
d2fdx2   = symbolic2string('d2fdx2', nx,nx,nx);
d2fdu2   = symbolic2string('d2fdu2', nx,nu,nu);
d2fdk2   = symbolic2string('d2fdk2', nx,nk,nk);
d2fdudx  = symbolic2string('d2fdudx', nx,nx,nu);
d2fdxdu  = symbolic2string('d2fdxdu', nx,nu,nx);
d2fdk2   = symbolic2string('d2fdk2', nx,nk,nk);
d2fdkdx  = symbolic2string('d2fdkdx', nx,nx,nk);
d2fdxdk  = symbolic2string('d2fdxdk', nx,nk,nx);
d3fdx3   = symbolic2string('d3fdx3', nx,nx,nx,nx);
d3fdkdx2 = symbolic2string('d3fdkdx2', nx,nx,nx,nk);

r        = symbolic2string('r', nr);
drdx     = symbolic2string('drdx', nr,nx);
drdu     = symbolic2string('drdu', nr,nu);
drdk     = symbolic2string('drdk', nr,nk);
d2rdx2   = symbolic2string('d2rdx2', nr,nx,nx);
d2rdu2   = symbolic2string('d2rdu2', nr,nu,nu);
d2rdxdu  = symbolic2string('d2rdxdu', nr,nu,nx);
d2rdudx  = symbolic2string('d2rdudx', nr,nx,nu);
d2rdk2   = symbolic2string('d2rdk2', nr,nk,nk);
d2rdkdx  = symbolic2string('d2rdkdx', nr,nx,nk);
d2rdkdu  = symbolic2string('d2rdkdu', nr,nu,nk);
d2rdxdk  = symbolic2string('d2rdxdk', nr,nk,nx);
d2rdudk  = symbolic2string('d2rdudk', nr,nk,nu);

y        = symbolic2string('y', ny);
dydx     = symbolic2string('dydx', ny,nx);
dydu     = symbolic2string('dydu', ny,nu);
d2ydx2   = symbolic2string('d2ydx2', ny,nx,nx);
d2ydu2   = symbolic2string('d2ydu2', ny,nu,nu);
d2ydxdu  = symbolic2string('d2ydxdu', ny,nu,nx);
d2ydudx  = symbolic2string('d2ydudx', ny,nx,nu);

% Replace species names with vector index names
if verbose; fprintf('Replacing names with vectorized variables...\n'); end
if verbose; fprintf('   states...\n');end
for i = 1:nx
    sys_string = sprintf('x(%d)', i);
    f        = regexprep(f, xStrs{i}, sys_string, 0);
    dfdx     = regexprep(dfdx, xStrs{i}, sys_string, 0);
    dfdu     = regexprep(dfdu, xStrs{i}, sys_string, 0);
    dfdk     = regexprep(dfdk, xStrs{i}, sys_string, 0);
    d2fdx2   = regexprep(d2fdx2, xStrs{i}, sys_string, 0);
    d2fdu2   = regexprep(d2fdu2, xStrs{i}, sys_string, 0);
    d2fdk2   = regexprep(d2fdk2, xStrs{i}, sys_string, 0);
    d2fdudx  = regexprep(d2fdudx, xStrs{i}, sys_string, 0);
    d2fdxdu  = regexprep(d2fdxdu, xStrs{i}, sys_string, 0);
    d2fdkdx  = regexprep(d2fdkdx, xStrs{i}, sys_string, 0);
    d2fdxdk  = regexprep(d2fdxdk, xStrs{i}, sys_string, 0);
    d3fdx3   = regexprep(d3fdx3, xStrs{i}, sys_string, 0);
    d3fdkdx2 = regexprep(d3fdkdx2, xStrs{i}, sys_string, 0);
    r        = regexprep(r, xStrs{i}, sys_string, 0);
    drdx     = regexprep(drdx, xStrs{i}, sys_string, 0);
    drdu     = regexprep(drdu, xStrs{i}, sys_string, 0);
    drdk     = regexprep(drdk, xStrs{i}, sys_string, 0);
    d2rdx2   = regexprep(d2rdx2, xStrs{i}, sys_string, 0);
    d2rdu2   = regexprep(d2rdu2, xStrs{i}, sys_string, 0);
    d2rdxdu  = regexprep(d2rdxdu, xStrs{i}, sys_string, 0);
    d2rdudx  = regexprep(d2rdudx, xStrs{i}, sys_string, 0);
    d2rdk2   = regexprep(d2rdk2, xStrs{i}, sys_string, 0);
    d2rdkdx  = regexprep(d2rdkdx, xStrs{i}, sys_string, 0);
    d2rdkdu  = regexprep(d2rdkdu, xStrs{i}, sys_string, 0);
    d2rdxdk  = regexprep(d2rdxdk, xStrs{i}, sys_string, 0);
    d2rdudk  = regexprep(d2rdudk, xStrs{i}, sys_string, 0);
    y        = regexprep(y, xStrs{i}, sys_string, 0);
    dydx     = regexprep(dydx, xStrs{i}, sys_string, 0);
    dydu     = regexprep(dydu, xStrs{i}, sys_string, 0);
    d2ydx2   = regexprep(d2ydx2, xStrs{i}, sys_string, 0);
    d2ydu2   = regexprep(d2ydu2, xStrs{i}, sys_string, 0);
    d2ydxdu  = regexprep(d2ydxdu, xStrs{i}, sys_string, 0);
    d2ydudx  = regexprep(d2ydudx, xStrs{i}, sys_string, 0);
end

% Replace input namus with vector index names
if verbose; fprintf('   input...\n'); end
for i = 1:nu
    sys_string = sprintf('u(%d)', i);
    f        = regexprep(f, uStrs{i}, sys_string, 0);
    dfdx     = regexprep(dfdx, uStrs{i}, sys_string, 0);
    dfdu     = regexprep(dfdu, uStrs{i}, sys_string, 0);
    dfdk     = regexprep(dfdk, uStrs{i}, sys_string, 0);
    d2fdx2   = regexprep(d2fdx2, uStrs{i}, sys_string, 0);
    d2fdu2   = regexprep(d2fdu2, uStrs{i}, sys_string, 0);
    d2fdk2   = regexprep(d2fdk2, uStrs{i}, sys_string, 0);
    d2fdudx  = regexprep(d2fdudx, uStrs{i}, sys_string, 0);
    d2fdxdu  = regexprep(d2fdxdu, uStrs{i}, sys_string, 0);
    d2fdkdx  = regexprep(d2fdkdx, uStrs{i}, sys_string, 0);
    d2fdxdk  = regexprep(d2fdxdk, uStrs{i}, sys_string, 0);
    d3fdx3   = regexprep(d3fdx3, uStrs{i}, sys_string, 0);
    d3fdkdx2 = regexprep(d3fdkdx2, uStrs{i}, sys_string, 0);
    r        = regexprep(r, uStrs{i}, sys_string, 0);
    drdx     = regexprep(drdx, uStrs{i}, sys_string, 0);
    drdu     = regexprep(drdu, uStrs{i}, sys_string, 0);
    drdk     = regexprep(drdk, uStrs{i}, sys_string, 0);
    d2rdx2   = regexprep(d2rdx2, uStrs{i}, sys_string, 0);
    d2rdu2   = regexprep(d2rdu2, uStrs{i}, sys_string, 0);
    d2rdxdu  = regexprep(d2rdxdu, uStrs{i}, sys_string, 0);
    d2rdudx  = regexprep(d2rdudx, uStrs{i}, sys_string, 0);
    d2rdk2   = regexprep(d2rdk2, uStrs{i}, sys_string, 0);
    d2rdkdx  = regexprep(d2rdkdx, uStrs{i}, sys_string, 0);
    d2rdkdu  = regexprep(d2rdkdu, uStrs{i}, sys_string, 0);
    d2rdxdk  = regexprep(d2rdxdk, uStrs{i}, sys_string, 0);
    d2rdudk  = regexprep(d2rdudk, uStrs{i}, sys_string, 0);
    y        = regexprep(y, uStrs{i}, sys_string, 0);
    dydx     = regexprep(dydx, uStrs{i}, sys_string, 0);
    dydu     = regexprep(dydu, uStrs{i}, sys_string, 0);
    d2ydx2   = regexprep(d2ydx2, uStrs{i}, sys_string, 0);
    d2ydu2   = regexprep(d2ydu2, uStrs{i}, sys_string, 0);
    d2ydxdu  = regexprep(d2ydxdu, uStrs{i}, sys_string, 0);
    d2ydudx  = regexprep(d2ydudx, uStrs{i}, sys_string, 0);
end

% Replace parameters with vector index names
if verbose; fprintf('   kinetic parameters...\n');end
for i = 1:nk
    sys_string = sprintf('k(%d)', i);
    f        = regexprep(f, kStrs{i}, sys_string, 0);
    dfdx     = regexprep(dfdx, kStrs{i}, sys_string, 0);
    dfdu     = regexprep(dfdu, kStrs{i}, sys_string, 0);
    dfdk     = regexprep(dfdk, kStrs{i}, sys_string, 0);
    d2fdx2   = regexprep(d2fdx2, kStrs{i}, sys_string, 0);
    d2fdu2   = regexprep(d2fdu2, kStrs{i}, sys_string, 0);
    d2fdk2   = regexprep(d2fdk2, kStrs{i}, sys_string, 0);
    d2fdudx  = regexprep(d2fdudx, kStrs{i}, sys_string, 0);
    d2fdxdu  = regexprep(d2fdxdu, kStrs{i}, sys_string, 0);
    d2fdkdx  = regexprep(d2fdkdx, kStrs{i}, sys_string, 0);
    d2fdxdk  = regexprep(d2fdxdk, kStrs{i}, sys_string, 0);
    d3fdx3   = regexprep(d3fdx3, kStrs{i}, sys_string, 0);
    d3fdkdx2 = regexprep(d3fdkdx2, kStrs{i}, sys_string, 0);
    r        = regexprep(r, kStrs{i}, sys_string, 0);
    drdx     = regexprep(drdx, kStrs{i}, sys_string, 0);
    drdu     = regexprep(drdu, kStrs{i}, sys_string, 0);
    drdk     = regexprep(drdk, kStrs{i}, sys_string, 0);
    d2rdx2   = regexprep(d2rdx2, kStrs{i}, sys_string, 0);
    d2rdu2   = regexprep(d2rdu2, kStrs{i}, sys_string, 0);
    d2rdxdu  = regexprep(d2rdxdu, kStrs{i}, sys_string, 0);
    d2rdudx  = regexprep(d2rdudx, kStrs{i}, sys_string, 0);
    d2rdk2   = regexprep(d2rdk2, kStrs{i}, sys_string, 0);
    d2rdkdx  = regexprep(d2rdkdx, kStrs{i}, sys_string, 0);
    d2rdkdu  = regexprep(d2rdkdu, kStrs{i}, sys_string, 0);
    d2rdxdk  = regexprep(d2rdxdk, kStrs{i}, sys_string, 0);
    d2rdudk  = regexprep(d2rdudk, kStrs{i}, sys_string, 0);
end

if verbose; fprintf('   ...done.\n'); end

%% Convert strings into function handles
if verbose; fprintf('Converting expressions into handles...'); end

% Clear unnecessary variables from scope
clear SymModel vSyms kSyms sSyms qSyms xuSyms xSyms uSyms x0 ...
    vStrs kStrs sStrs qStrs xuStrs xStrs uStrs ...
    xNamesFull uNamesFull vxNames vuNames ...
    C1Entries C1Values C2Entries C2Values cEntries cValues ...
    nC1Entries nC2Entries ncEntries defaultOpts opts ...
    addlength currentLength iExpr ind match nAdd nExpr uAddInd sys_string ...
    i iv ik is iq iu ix ir iy
    
f = eval(['@(t,x,u,k) [' f ']']);

dfdx = regexprep(dfdx, '[', ''); %remove extra "[" from front of lines
dfdx = regexprep(dfdx, ']', ''); %and "]"
dfdx = strtrim(dfdx);            %trim excess new lines from the end
dfdx = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' dfdx '])))']);

dfdu = regexprep(dfdu, '[', '');
dfdu = regexprep(dfdu, ']', '');
dfdu = strtrim(dfdu);
dfdu = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' dfdu '])))']);

dfdk = regexprep(dfdk, '[', '');
dfdk = regexprep(dfdk, ']', '');
dfdk = strtrim(dfdk);
dfdk = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' dfdk '])))']);

d2fdx2 = regexprep(d2fdx2, '[', '');
d2fdx2 = regexprep(d2fdx2, ']', '');
d2fdx2 = strtrim(d2fdx2);
d2fdx2 = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2fdx2 '])))']);

d2fdu2 = regexprep(d2fdu2, '[', '');
d2fdu2 = regexprep(d2fdu2, ']', '');
d2fdu2 = strtrim(d2fdu2);
d2fdu2 = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2fdu2 '])))']);

d2fdk2 = regexprep(d2fdk2, '[', '');
d2fdk2 = regexprep(d2fdk2, ']', '');
d2fdk2 = strtrim(d2fdk2);
d2fdk2 = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2fdk2 '])))']);

d2fdudx = regexprep(d2fdudx, '[', '');
d2fdudx = regexprep(d2fdudx, ']', '');
d2fdudx = strtrim(d2fdudx);
d2fdudx = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2fdudx '])))']);

d2fdxdu = regexprep(d2fdxdu, '[', '');
d2fdxdu = regexprep(d2fdxdu, ']', '');
d2fdxdu = strtrim(d2fdxdu);
d2fdxdu = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2fdxdu '])))']);

d2fdkdx = regexprep(d2fdkdx, '[', '');
d2fdkdx = regexprep(d2fdkdx, ']', '');
d2fdkdx = strtrim(d2fdkdx);
d2fdkdx = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2fdkdx '])))']);

d2fdxdk = regexprep(d2fdxdk, '[', '');
d2fdxdk = regexprep(d2fdxdk, ']', '');
d2fdxdk = strtrim(d2fdxdk);
d2fdxdk = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2fdxdk '])))']);

d3fdx3 = regexprep(d3fdx3, '[', '');
d3fdx3 = regexprep(d3fdx3, ']', '');
d3fdx3 = strtrim(d3fdx3);
d3fdx3 = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d3fdx3 '])))']);

d3fdkdx2 = regexprep(d3fdkdx2, '[', '');
d3fdkdx2 = regexprep(d3fdkdx2, ']', '');
d3fdkdx2 = strtrim(d3fdkdx2);
d3fdkdx2 = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d3fdkdx2 '])))']);

r = eval(['@(t,x,u,k) [' r ']']);

% drdx = regexprep(drdx, '[', ''); %remove extra "[" from front of lines
% drdx = regexprep(drdx, ']', ''); %and "]"
% drdx = strtrim(drdx);            %trim excess new lines from the end
drdx = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' drdx '])))']);

% drdu = regexprep(drdu, '[', '');
% drdu = regexprep(drdu, ']', '');
% drdu = strtrim(drdu);
drdu = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' drdu '])))']);

% drdk = regexprep(drdk, '[', '');
% drdk = regexprep(drdk, ']', '');
% drdk = strtrim(drdk);
drdk = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' drdk '])))']);

% d2rdx2 = regexprep(d2rdx2, '[', '');
% d2rdx2 = regexprep(d2rdx2, ']', '');
% d2rdx2 = strtrim(d2rdx2);
d2rdx2 = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2rdx2 '])))']);

% d2rdu2 = regexprep(d2rdu2, '[', '');
% d2rdu2 = regexprep(d2rdu2, ']', '');
% d2rdu2 = strtrim(d2rdu2);
d2rdu2 = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2rdu2 '])))']);

% d2rdxdu = regexprep(d2rdxdu, '[', '');
% d2rdxdu = regexprep(d2rdxdu, ']', '');
% d2rdxdu = strtrim(d2rdxdu);
d2rdxdu = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2rdxdu '])))']);

% d2rdudx = regexprep(d2rdudx, '[', '');
% d2rdudx = regexprep(d2rdudx, ']', '');
% d2rdudx = strtrim(d2rdudx);
d2rdudx = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2rdudx '])))']);

% d2rdk2 = regexprep(d2rdk2, '[', '');
% d2rdk2 = regexprep(d2rdk2, ']', '');
% d2rdk2 = strtrim(d2rdk2);
d2rdk2 = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2rdk2 '])))']);

% d2rdkdx = regexprep(d2rdkdx, '[', '');
% d2rdkdx = regexprep(d2rdkdx, ']', '');
% d2rdkdx = strtrim(d2rdkdx);
d2rdkdx = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2rdkdx '])))']);

% d2rdkdu = regexprep(d2rdkdu, '[', '');
% d2rdkdu = regexprep(d2rdkdu, ']', '');
% d2rdkdu = strtrim(d2rdkdu);
d2rdkdu = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2rdkdu '])))']);

% d2rdxdk = regexprep(d2rdxdk, '[', '');
% d2rdxdk = regexprep(d2rdxdk, ']', '');
% d2rdxdk = strtrim(d2rdxdk);
d2rdxdk = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2rdxdk '])))']);

% d2rdudk = regexprep(d2rdudk, '[', '');
% d2rdudk = regexprep(d2rdudk, ']', '');
% d2rdudk = strtrim(d2rdudk);
d2rdudk = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2rdudk '])))']);

y = eval(['@(t,x,u) [' y ']']);

dydx = regexprep(dydx, '[', ''); %remove extra "[" from front of lines
dydx = regexprep(dydx, ']', ''); %and "]"
dydx = strtrim(dydx);            %trim excess new lines from the end
dydx = eval(['@(t,x,u) inf2big(nan2zero(sparse([' dydx '])))']);

dydu = regexprep(dydu, '[', '');
dydu = regexprep(dydu, ']', '');
dydu = strtrim(dydu);
dydu = eval(['@(t,x,u) inf2big(nan2zero(sparse([' dydu '])))']);

d2ydx2 = regexprep(d2ydx2, '[', '');
d2ydx2 = regexprep(d2ydx2, ']', '');
d2ydx2 = strtrim(d2ydx2);
d2ydx2 = eval(['@(t,x,u) inf2big(nan2zero(sparse([' d2ydx2 '])))']);

d2ydu2 = regexprep(d2ydu2, '[', '');
d2ydu2 = regexprep(d2ydu2, ']', '');
d2ydu2 = strtrim(d2ydu2);
d2ydu2 = eval(['@(t,x,u) inf2big(nan2zero(sparse([' d2ydu2 '])))']);

d2ydxdu = regexprep(d2ydxdu, '[', '');
d2ydxdu = regexprep(d2ydxdu, ']', '');
d2ydxdu = strtrim(d2ydxdu);
d2ydxdu = eval(['@(t,x,u) inf2big(nan2zero(sparse([' d2ydxdu '])))']);

d2ydudx = regexprep(d2ydudx, '[', '');
d2ydudx = regexprep(d2ydudx, ']', '');
d2ydudx = strtrim(d2ydudx);
d2ydudx = eval(['@(t,x,u) inf2big(nan2zero(sparse([' d2ydudx '])))']);

if verbose; fprintf('done.\n'); end

%% Replace rate constants with their values
if verbose; fprintf('Evaluating handles...'); end

m = InitializeModel();

m.Type = 'Model.AnalyticReactions';
m.Name = name;

m.Compartments = struct('Name', vNames, 'Dimension', num2cell(dv));
m.Parameters   = struct('Name', kNames, 'Value', num2cell(k));
m.Seeds        = struct('Name', sNames, 'Value', num2cell(s));
m.Inputs       = struct('Name', uNames, 'Compartment', vNames(vuInd));
m.States       = struct('Name', xNames, 'Compartment', vNames(vxInd), 'InitialValue', initial_values);
m.Reactions    = struct('Name', rNames);
m.Outputs      = struct('Name', yNames);

m.nv = nv;
m.nk = nk;
m.ns = ns;
m.nu = nu;
m.nx = nx;
m.nr = nr;
m.ny = ny;

m.dv = dv;
m.k  = k;
m.s  = s;
m.dx0ds = dx0ds;
m.x0c = x0c;

m.u = u;

m.vxInd = vxInd;
m.vuInd = vuInd;
%rOrder
%krInd

%As
%Bs
%Cs

clear vNames kNames sNames uNames xNames rNames yNames

m.f         = @(t,x,u)f(t,x,u,k);

if order >= 1
    m.dfdx      = @(t,x,u)dfdx(t,x,u,k);
    m.dfdk      = @(t,x,u)dfdk(t,x,u,k);
    m.dfdu      = @(t,x,u)dfdu(t,x,u,k);
end

if order >= 2
    m.d2fdx2    = @(t,x,u)d2fdx2(t,x,u,k);
    m.d2fdu2    = @(t,x,u)d2fdu2(t,x,u,k);
    m.d2fdk2    = @(t,x,u)d2fdk2(t,x,u,k);
    m.d2fdudx   = @(t,x,u)d2fdudx(t,x,u,k);
    m.d2fdxdu   = @(t,x,u)d2fdxdu(t,x,u,k);
    m.d2fdkdx   = @(t,x,u)d2fdkdx(t,x,u,k);
    m.d2fdkdu   = @(t,x,u)d2fdkdu(t,x,u,k);
    m.d2fdxdk   = @(t,x,u)d2fdxdk(t,x,u,k);
    m.d2fdudk   = @(t,x,u)d2fdudk(t,x,u,k);
end

if order >= 3
    m.d3fdx3    = @(t,x,u)d3fdx3(t,x,u,k);
    m.d3fdkdx2  = @(t,x,u)d3fdkdx2(t,x,u,k);
end

m.S = S;
m.r = @(t,x,u)r(t,x,u,k);

if order >= 1
    m.drdx      = @(t,x,u)drdx(t,x,u,k);
    m.drdk      = @(t,x,u)drdk(t,x,u,k);
    m.drdu      = @(t,x,u)drdu(t,x,u,k);
end

if order >= 2
    m.d2rdx2    = @(t,x,u)d2rdx2(t,x,u,k);
    m.d2rdu2    = @(t,x,u)d2rdu2(t,x,u,k);
    m.d2rdk2    = @(t,x,u)d2rdk2(t,x,u,k);
    m.d2rdudx   = @(t,x,u)d2rdudx(t,x,u,k);
    m.d2rdxdu   = @(t,x,u)d2rdxdu(t,x,u,k);
    m.d2rdkdx   = @(t,x,u)d2rdkdx(t,x,u,k);
    m.d2rdkdu   = @(t,x,u)d2rdkdu(t,x,u,k);
    m.d2rdxdk   = @(t,x,u)d2rdxdk(t,x,u,k);
    m.d2rdudk   = @(t,x,u)d2rdudk(t,x,u,k);
end

m.y = @(t,x,u)vectorize_y(y,t,x,u);

if order >= 1
    m.dydx      = @(t,x,u)dydx(t,x,u);
    m.dydu      = @(t,x,u)dydu(t,x,u);
end

if order >= 2
    m.d2ydx2    = @(t,x,u)d2ydx2(t,x,u);
    m.d2ydu2    = @(t,x,u)d2ydu2(t,x,u);
    m.d2ydudx   = @(t,x,u)d2ydudx(t,x,u);
    m.d2ydxdu   = @(t,x,u)d2ydxdu(t,x,u);
end

m.Ready  = true;
%m.add
m.Update = @update;

if verbose; fprintf('done.\n'); end

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Update function %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function varargout = update(newk)
        % Apply changes
        k = newk;
        m.k  = k;
        
        % Update function handles
        m.f             = @(t,x,u)f(t,x,u,k);
        m.r             = @(t,x,u)r(t,x,u,k);
        m.y             = @(t,x,u)vectorize_y(y,t,x,u);

        if order >= 1
            m.dfdx      = @(t,x,u)dfdx(t,x,u,k);
            m.dfdk      = @(t,x,u)dfdk(t,x,u,k);
            m.dfdu      = @(t,x,u)dfdu(t,x,u,k);
            m.drdx      = @(t,x,u)drdx(t,x,u,k);
            m.drdk      = @(t,x,u)drdk(t,x,u,k);
            m.drdu      = @(t,x,u)drdu(t,x,u,k);
            m.dydx      = @(t,x,u)dydx(t,x,u);
            m.dydu      = @(t,x,u)dydu(t,x,u);
        end
        
        if order >= 2
            m.d2fdx2    = @(t,x,u)d2fdx2(t,x,u,k);
            m.d2fdu2    = @(t,x,u)d2fdu2(t,x,u,k);
            m.d2fdk2    = @(t,x,u)d2fdk2(t,x,u,k);
            m.d2fdudx   = @(t,x,u)d2fdudx(t,x,u,k);
            m.d2fdxdu   = @(t,x,u)d2fdxdu(t,x,u,k);
            m.d2fdkdx   = @(t,x,u)d2fdkdx(t,x,u,k);
            m.d2fdkdu   = @(t,x,u)d2fdkdu(t,x,u,k);
            m.d2fdxdk   = @(t,x,u)d2fdxdk(t,x,u,k);
            m.d2fdudk   = @(t,x,u)d2fdudk(t,x,u,k);
            m.d2rdx2    = @(t,x,u)d2rdx2(t,x,u,k);
            m.d2rdu2    = @(t,x,u)d2rdu2(t,x,u,k);
            m.d2rdk2    = @(t,x,u)d2rdk2(t,x,u,k);
            m.d2rdudx   = @(t,x,u)d2rdudx(t,x,u,k);
            m.d2rdxdu   = @(t,x,u)d2rdxdu(t,x,u,k);
            m.d2rdkdx   = @(t,x,u)d2rdkdx(t,x,u,k);
            m.d2rdkdu   = @(t,x,u)d2rdkdu(t,x,u,k);
            m.d2rdxdk   = @(t,x,u)d2rdxdk(t,x,u,k);
            m.d2rdudk   = @(t,x,u)d2rdudk(t,x,u,k);
            m.d2ydx2    = @(t,x,u)d2ydx2(t,x,u);
            m.d2ydu2    = @(t,x,u)d2ydu2(t,x,u);
            m.d2ydudx   = @(t,x,u)d2ydudx(t,x,u);
            m.d2ydxdu   = @(t,x,u)d2ydxdu(t,x,u);
        end
        
        if order >= 3
            m.d3fdx3    = @(t,x,u)d3fdx3(t,x,u,k);
            m.d3fdkdx2  = @(t,x,u)d3fdkdx2(t,x,u,k);
        end
        
        m.Update = @update;
        
        varargout{1} = m;
    end

%% %%%%%%%%%%%%%%%%%%%%%%%
%%%% Helper functions %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
    function string_rep = symbolic2string(variable_name, varargin)
        dimensions = [varargin{:}];
        if numel(dimensions) == 1
            dimensions = [dimensions, 1];
        end
        
        if any(dimensions == 0)
            string_rep = ['zeros([' num2str(dimensions) '])'];
        else
            string_rep = evalc(['disp(' variable_name ')']);
        end
    end

    function val = vectorize_y(y, t, x, u)
        nt = numel(t);
        val = zeros(ny,nt);
        for it = 1:nt
            val(:,it) = y(t(it), x(:,it), u(:,it));
        end
    end

end