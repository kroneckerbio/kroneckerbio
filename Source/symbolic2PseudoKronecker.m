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
%       .NewJacobianMethod [ true | {false} ]
%           Determines whether to use the more memory-efficient Jacobian
%           method or the original method. For smaller models ( < 100
%           reactions or species), there is little difference between the
%           methods.
%       .UseMEX [ true | {false} ]
%           Set to true to write MEX functions in C that calculate f, r,
%           and their derivatives with respect to x, u, and k. For larger
%           models, this will speed up most of Kronecker's functionality.
%           Following building the model, the .mex files will need to be
%           compiled using the compileMEXFunctions function. Currently
%           requires gcc, a C compiler, to be configured to be used by
%           MATLAB. Following compilation, the .mex files must be in the
%           current path for the model to work.
%       .MEXDirectory [ string {['mexfuns_' dateandtimestring]} ]
%           If .UseMEX is set to TRUE, the string provided here sets the
%           directory to which the .mex files will be written. If .UseMEX
%           is FALSE, this option is ignored.
%
%   Outputs
%   m: [ psuedo-kronecker model scalar ]
%       The useable form of the model

% (c) 2015 David Flowers, David R Hagen, & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Clean up inputs
if nargin < 2
    opts = [];
end

% Set up default directory name, based on the clock
thistime = num2cell(clock);
thistime = cellfun(@num2str,thistime,'UniformOutput',false);
thistime = [thistime{:}];
defaultMEXdirectory = ['mexfuns_' regexprep(thistime,'\.','_') filesep];

% Default options
defaultOpts.Order             = 2;
defaultOpts.VolumeToParameter = false;
defaultOpts.Verbose           = 0;
defaultOpts.MultipleyMembers  = false;
defaultOpts.NewJacobianMethod = false;
defaultOpts.UseMEX            = false;
defaultOpts.MEXDirectory      = defaultMEXdirectory;

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

order = opts.Order;

if opts.UseMEX && exist(opts.MEXDirectory,'dir') ~= 7
    mkdir(opts.MEXDirectory);
end

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
    dv      = zeros(3,1);
    v       = zeros(1,1);
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
    uNamesFull{iu} = [vNames{vuInd(iu)} '.' uNames{iu}];
end

xStrs = cell(nx,1);
xNamesFull = cell(nx,1);
for ix = 1:nx
    xStrs{ix} = char(xSyms(ix));
    xNamesFull{ix} = [vNames{vxInd(ix)} '.' xNames{ix}];
end

% Convert compartment volumes to parameters or constants
if opts.VolumeToParameter
    nk = nk + nv;
    kNames = [kNames; vNames];
    k = [k; v];
    kSyms = [kSyms; vSyms];
else% ~VolumeToParameter
    r = subs(r, vSyms, v, 0);
end

f = S*r;

q = SymModel.q;
nq = SymModel.nq;
qStrs = cell(nq,1);
qSyms = SymModel.qSyms;

% Standardize names of input parameters
for iq = 1:nq
    newqSym = sym(['q' num2str(iq) 'x']);
    u = subs(u, qSyms(iq), newqSym);
    qSyms(iq) = newqSym;
    qStrs{iq} = char(qSyms(iq));
end

%% Determine the species and parameters in each reaction

rstates = cell(nr,1);
rinputs = cell(nr,1);
rparams = cell(nr,1);
% For each reaction...
if verbose; fprintf('Determining what terms are present in each reaction rate...\n'); end
% Substitute 1's for parameters and inputs to find states
statesubs = subs(r,[kSyms;uSyms],ones(size([kSyms;uSyms])));
% Substitute 1's for parameters and states to find inputs
inputsubs = subs(r,[kSyms;xSyms],ones(size([kSyms;xSyms])));
% Substitute 1's for species to find parameters
paramssubs =  subs(r,[xSyms;uSyms],ones(size([xSyms;uSyms])));
ndigitstoprint = 4;
strtoprint = ['Reaction %' num2str(ndigitstoprint) 'g'];
lenstr = length(strtoprint) - 3 + ndigitstoprint;
fprintf(repmat(' ',1,lenstr))
for ri = 1:nr
    if verbose; fprintf([repmat('\b',1,lenstr) strtoprint],ri); end
    rstates{ri} = symvar(statesubs(ri));
    rinputs{ri} = symvar(inputsubs(ri));
    rparams{ri} = symvar(paramssubs(ri));
end
if verbose; fprintf('\n'); end

%% Find nonzero derivative terms

% Convert the species and parameter indicators into logical matrices
rstatesi = false(nr,nx);
rinputsi = false(nr,nu);
rparamsi = false(nr,nk);
if verbose; fprintf('Determining approximate sparsity patterns for derivative matrices...\n'); end
ndigitstoprint = 4;
strtoprint = ['Reaction %' num2str(ndigitstoprint) 'g'];
lenstr = length(strtoprint) - 3 + ndigitstoprint;
fprintf(repmat(' ',1,lenstr))
for ri = 1:nr
    if verbose; fprintf([repmat('\b',1,lenstr) strtoprint],ri); end
    for xii = 1:length(rstates{ri})
        rstatesi(ri,:) = rstatesi(ri,:) + logical(rstates{ri}(xii) == xSyms)';
    end
    for ui = 1:length(rinputs{ri})
        rinputsi(ri,:) = rinputsi(ri,:) + logical(rinputs{ri}(ui) == uSyms)';
    end
    for ki = 1:length(rparams{ri})
        rparamsi(ri,:) = rparamsi(ri,:) + logical(rparams{ri}(ki) == kSyms)';
    end
end

uqi = false(nu,nq);
for ui = 1:nu
    qsinui = symvar(u(ui));
    for qi = 1:length(qsinui)
        uqi(ui,:) = uqi(ui,:) + logical(qsinui(qi) == qSyms)';
    end
end

yxi = false(ny,nx);
yui = false(ny,nu);
for yi = 1:ny
    xsandusinyi = symvar(y(yi));
    for ii = 1:length(xsandusinyi)
        yxi(yi,:) = yxi(yi,:) + logical(xsandusinyi(ii) == xSyms)';
        yui(yi,:) = yui(yi,:) + logical(xsandusinyi(ii) == uSyms)';
    end
end
% The above logical arrays are the nonzero elements of the first derivative
% matrix for r, u, and y. Record these in a struct, and construct a
% function that can calculate whether higher-order derivatives contain a
% given parameter, input, or state.
nz.r.x = rstatesi;
nz.r.u = rinputsi;
nz.r.k = rparamsi;
nz.f.x = logical(abs(S)*rstatesi); % Take the absolute value of S so that there are no accidental cancellations between positive and negative terms
nz.f.u = logical(abs(S)*rinputsi);
nz.f.k = logical(abs(S)*rparamsi);
nz.u.q = uqi;
nz.y.x = yxi;
nz.y.u = yui;

% Set up sizes struct
sizes.f = nx;
sizes.r = nr;
sizes.k = nk;
sizes.x = nx;
sizes.u = nu;
sizes.q = nq;
sizes.y = ny;

    function nzout = getNonZeroEntries(num,dens)
        % Inputs:
        %   num: "numerator" of derivative, either 'r', 'f', 'u', or 'y'
        %   dens: "denominator" terms of derivative, any combination of 'x',
        %   'u', and 'k' in a cell array for 'r' and 'f' numerators, 'q'
        %   for 'u' numerators, and 'x' and 'u' for 'y' numerators
        % Outputs:
        %   nzout: a logical array of size n(r, f, u, or y)-by-n(x or u or
        %   k)-by-n(x or u or k)-by-..., where nzout(i,j,k,...) is true if
        %   d(r or f)/d(x or u or k)(x or u or k)... is potentially
        %   nonzero, based on what terms are present in each of the
        %   reaction or state first derivative terms. Note that the first term
        %   appearing in the "denominator" string of the derivative is the
        %   last dimension, since it is the derivative taken last. I.E.,
        %   dr/dxdk has dimensions nr-by-nk-by-nx.
        %
        % Example:
        %   nzout = getNonZeroEntries('f',{'x','k'}) would return nzout,
        %   the logical array of size nf-by-nx-by-nk indicating the
        %   potentially nonzero entries of df/dkdx.
        if ischar(dens)
            dens = {dens};
        end
        ndens = numel(dens);
        singletondims = 3:(ndens+1);
        nzout = true(sizes.(num),1); % Start with this value
        for nzi = 1:ndens
            nztoadd = nz.(num).(dens{nzi});
            if length(singletondims) <= 1
                permutedims = [1 singletondims(1:nzi-1) 2];
            else
                permutedims = [1 singletondims(1:nzi-1) 2 singletondims(nzi:end)];
            end
            permutednztoadd = permute(nztoadd,permutedims);
            nzout = bsxfun(@and,nzout,permutednztoadd);
        end
    end

% if order > 1
%     nz_rxx = nz_addorder(nz_rx,nz_rx);
%     nz_rkx = nz_addorder(nz_rx,nz_rk);
%     nz_rxk = nz_addorder(nz_rk,nz_rx);
%     nz_rkk = nz_addorder(nz_rk,nz_rk);
%     nz_ruu = nz_addorder(nz_ru,nz_ru);
%     nz_rux = nz_addorder(nz_rx,nz_ru);
%     nz_rxu = nz_addorder(nz_ru,nz_rx);
%     nz_ruk = nz_addorder(nz_rk,nz_ru);
%     nz_rku = nz_addorder(nz_ru,nz_rk);
% end
% if order > 2
%     error('This function does not support orders greater than 2 at this time.')
% end
% 
%     function nzout = nz_addorder(nz1,nz2)
%         nz1dimsizes = cell(ndims(nz1),1);
%         nz2dimsizes = cell(ndims(nz2),1);
%         [nz1dimsizes{:}] = size(nz1);
%         [nz2dimsizes{:}] = size(nz2);
%         permutedims1 = [nz1dimsizes{:} 1];
%         nzout = bsxfun(@and,nz1,permute(nz2,
%     end
fprintf('\n')

%% Generate derivatives of desired order

% Determine which method to use to calculate Jacobians
if opts.NewJacobianMethod
    jacobianfun = @calcDerivative;
    reshapefun = @reshapeDerivative;
else
    jacobianfun = @(x, syms) jacobian(vec(x), syms);
    reshapefun = @reshape_extraarguments;
end

if order >= 1
    % Gradient of u with respect to q
    if verbose; fprintf('Calculating dudq...'); end
    dudq = jacobianfun(u, qSyms);
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of r with respect to x
    if verbose; fprintf('Calculating drdx...'); end
    drdx = jacobianfun(r, xSyms);
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of r with respect to u
    if verbose; fprintf('Calculating drdu...'); end
    drdu = jacobianfun(r, uSyms);
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of r with respect to k
    if verbose; fprintf('Calculating drdk...'); end
    drdk = jacobianfun(r, kSyms);
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of f with respect to x
    if verbose; fprintf('Calculating dfdx...'); end
    dfdx = S*drdx;
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of f with respect to u
    if verbose; fprintf('Calculating dfdu...'); end
    dfdu = S*drdu;
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of f with respect to k
    if verbose; fprintf('Calculating dfdk...'); end
    dfdk = S*drdk;
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of y with respect to x
    if verbose; fprintf('Calculating dydx...'); end
    dydx = jacobianfun(y, xSyms);
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of y with respect to u
    if verbose; fprintf('Calculating dydu...'); end
    dydu = jacobian(y, uSyms);
    if verbose; fprintf('Done.\n'); end

else
    drdx = '';
    drdk = '';
    dfdx = '';
    dfdk = '';
    dydx = '';
    dydu = '';
end

if order >= 2
    % Gradient of drdx with respect to x
    if verbose; fprintf('Calculating d2rdx2...'); end
    d2rdx2 = jacobianfun(drdx, xSyms);
    if verbose; fprintf('\n'); end
    
    % Gradient of drdu with respect to u
    if verbose; fprintf('Calculating d2rdu2...'); end
    d2rdu2 = jacobianfun(drdu, uSyms);
    if verbose; fprintf('\n'); end
    
    % Gradient of drdu with respect to x
    if verbose; fprintf('Calculating d2rdxdu...'); end
    d2rdxdu = jacobianfun(drdu, xSyms);
    if verbose; fprintf('\n'); end
    
    % Gradient of drdx with respect to u
    if verbose; fprintf('Calculating d2rdudx...'); end
    d2rdudx = jacobianfun(drdx, uSyms);
    if verbose; fprintf('\n'); end
    
    % Gradient of drdk with respect to k
    if verbose; fprintf('Calculating d2rdk2...'); end
    d2rdk2 = jacobianfun(drdk, kSyms);
    if verbose; fprintf('\n'); end
    
    % Gradient of drdx with respect to k
    if verbose; fprintf('Calculating d2rdkdx...'); end
    d2rdkdx = jacobianfun(drdx, kSyms);
    if verbose; fprintf('\n'); end
    
    % Gradient of drdu with respect to k
    if verbose; fprintf('Calculating d2rdkdu...'); end
    d2rdkdu = jacobianfun(drdu, kSyms);
    if verbose; fprintf('\n'); end
    
    % Gradient of drdk with respect to x
    if verbose; fprintf('Calculating d2rdxdk...'); end
    d2rdxdk = jacobianfun(drdk, xSyms);
    if verbose; fprintf('\n'); end
    
    % Gradient of drdk with respect to u
    if verbose; fprintf('Calculating d2rdudk...'); end
    d2rdudk = jacobianfun(drdk, uSyms);
    if verbose; fprintf('\n'); end
    
    % Gradient of dfdx with respect to x
    if verbose; fprintf('Calculating d2fdx2...'); end
    d2fdx2 = S*reshapefun(d2rdx2, [nr nx*nx], 'r', {'x' 'x'});
    d2fdx2 = reshapefun(d2fdx2, [nx*nx nx], 'f', {'x' 'x'});
    if verbose; fprintf('\n'); end
    
    % Gradient of dfdu with respect to u
    if verbose; fprintf('Calculating d2fdu2...'); end
    d2fdu2 = S*reshapefun(d2rdu2, [nr,nu*nu], 'r', {'u' 'u'});
    d2fdu2 = reshapefun(d2fdu2, [nx*nu,nu], 'f', {'u' 'u'});
    if verbose; fprintf('\n'); end
    
    % Gradient of dfdu with respect to x
    if verbose; fprintf('Calculating d2fdxdu...'); end
    d2fdxdu = S*reshapefun(d2rdxdu, [nr,nu*nx], 'r', {'u' 'x'});
    d2fdxdu = reshapefun(d2fdxdu, [nx*nu,nx], 'f', {'u' 'x'});
    if verbose; fprintf('\n'); end
    
    % Gradient of dfdx with respect to u
    if verbose; fprintf('Calculating d2fdudx...'); end
    d2fdudx = S*reshapefun(d2rdudx, [nr,nx*nu], 'r', {'x' 'u'});
    d2fdudx = reshapefun(d2fdudx, [nx*nx,nu], 'f', {'x' 'u'});
    if verbose; fprintf('\n'); end
    
    % Gradient of dfdk with respect to k
    if verbose; fprintf('Calculating d2fdk2...'); end
    d2fdk2 = S*reshapefun(d2rdk2, [nr,nk*nk], 'r', {'k' 'k'});
    d2fdk2 = reshapefun(d2fdk2, [nx*nk,nk], 'f', {'k' 'k'});
    if verbose; fprintf('\n'); end
    
    % Gradient of dfdx with respect to k
    if verbose; fprintf('Calculating d2fdkdx...'); end
    d2fdkdx = S*reshapefun(d2rdkdx, [nr,nx*nk], 'r', {'x' 'k'});
    d2fdkdx = reshapefun(d2fdkdx, [nx*nx,nk], 'f', {'x' 'k'});
    if verbose; fprintf('\n'); end
    
    % Gradient of dfdu with respect to k
    if verbose; fprintf('Calculating d2fdkdu...'); end
    d2fdkdu = S*reshapefun(d2rdkdu, [nr,nu*nk], 'r',{'u' 'k'});
    d2fdkdu = reshapefun(d2fdkdu, [nx*nu,nk], 'f',{'u' 'k'});
    if verbose; fprintf('\n'); end
    
    % Gradient of dfdk with respect to x
    if verbose; fprintf('Calculating d2fdxdk...'); end
    d2fdxdk = S*reshapefun(d2rdxdk, [nr,nk*nx], 'r',{'k' 'x'});
    d2fdxdk = reshapefun(d2fdxdk, [nx*nk,nx], 'f',{'k' 'x'});
    if verbose; fprintf('\n'); end
    
    % Gradient of dfdk with respect to u
    if verbose; fprintf('Calculating d2fdudk...'); end
    d2fdudk = S*reshapefun(d2rdudk, [nr,nk*nu], 'r',{'k' 'u'});
    d2fdudk = reshapefun(d2fdudk, [nx*nk,nu], 'f',{'k' 'u'});
    if verbose; fprintf('\n'); end
    
    % Output second derivatives
    if verbose; fprintf('Calculating d2ydx2...'); end
    d2ydx2 = jacobianfun(dydx, xSyms);
    if verbose; fprintf('\n'); end
    
    if verbose; fprintf('Calculating d2ydu2...'); end
    d2ydu2 = jacobianfun(dydu, uSyms);
    if verbose; fprintf('\n'); end
    
    if verbose; fprintf('Calculating d2ydxdu...'); end
    d2ydxdu = jacobianfun(dydu, xSyms);
    if verbose; fprintf('\n'); end
    
    if verbose; fprintf('Calculating d2ydudx...'); end
    d2ydudx = jacobianfun(dydx, uSyms);
    if verbose; fprintf('\n'); end
    
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
    d2ydx2 = '';
    d2ydu2 = '';
    d2ydxdu = '';
    d2ydudx = '';
end

if order >= 3
    %Gradient of d2rdx2 with respect to x
    d3rdx3 = jacobianfun(vec(d2rdx2), xSyms);
    
    %Gradient of d2rdx2 with respect to k
    d3rdkdx2 = jacobianfun(vec(d2rdx2), kSyms);
    
    %Gradient of d2fdx2 with respect to x
    d3fdx3 = S*reshapefun(d3rdx3, [nr,nx*nx*nx], 'r',{'x' 'x' 'x'});
    d3fdx3 = reshapefun(d3fdx3, [xn*nx*nx,nx], 'f',{'x' 'x' 'x'});
    
    %Gradient of d2fdx2 with respect to k
    d3fdkdx2 = S*reshapefun(d3rdkdx2, [nr,nx*nx*nk], 'r',{'x' 'x' 'k'});
    d3fdkdx2 = reshapefun(d3fdkdx2, [nx*nx*nx,nk], 'f',{'x' 'x' 'k'});
else
    d3rdx3   = '';
    d3rdkdx2 = '';
    d3fdx3   = '';
    d3fdkdx2 = '';
end

%% Extract seed information
dx0ds = double(jacobianfun(x0, sSyms));
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

if opts.UseMEX
    symbolic2stringmethod = 'mex';
else
    symbolic2stringmethod = 'efficient';
end

if verbose; fprintf('Converting symbolics to functions...\n'); end
% Convert the symbolics into strings
u        = symbolic2function('u', 'u', {});
dudq     = symbolic2function('dudq', 'u', 'q');

f        = symbolic2function('f', 'f', {});
r        = symbolic2function('r', 'r', {});
y        = symbolic2function('y', 'y', {});

if order >= 1
    dfdx     = symbolic2function('dfdx', 'f', 'x');
    dfdu     = symbolic2function('dfdu', 'f', 'u');
    dfdk     = symbolic2function('dfdk', 'f', 'k');

    drdx     = symbolic2function('drdx', 'r', 'x');
    drdu     = symbolic2function('drdu', 'r', 'u');
    drdk     = symbolic2function('drdk', 'r', 'k');
    
    dydx     = symbolic2function('dydx', 'y', 'x');
    dydu     = symbolic2function('dydu', 'y', 'u');
end

if order >= 2
    d2fdx2   = symbolic2function('d2fdx2', 'f', {'x' 'x'});
    d2fdu2   = symbolic2function('d2fdu2', 'f', {'u' 'u'});
    d2fdk2   = symbolic2function('d2fdk2', 'f', {'k' 'k'});
    d2fdudx  = symbolic2function('d2fdudx', 'f', {'x' 'u'});
    d2fdxdu  = symbolic2function('d2fdxdu', 'f', {'u' 'x'});
    d2fdkdx  = symbolic2function('d2fdkdx', 'f', {'x' 'k'});
    d2fdxdk  = symbolic2function('d2fdxdk', 'f', {'k' 'x'});

    d2rdx2   = symbolic2function('d2rdx2', 'r', {'x' 'x'});
    d2rdu2   = symbolic2function('d2rdu2', 'r', {'u' 'u'});
    d2rdxdu  = symbolic2function('d2rdxdu', 'r', {'u' 'x'});
    d2rdudx  = symbolic2function('d2rdudx', 'r', {'x' 'u'});
    d2rdk2   = symbolic2function('d2rdk2', 'r', {'k' 'k'});
    d2rdkdx  = symbolic2function('d2rdkdx', 'r', {'x' 'k'});
    d2rdkdu  = symbolic2function('d2rdkdu', 'r', {'u' 'k'});
    d2rdxdk  = symbolic2function('d2rdxdk', 'r', {'k' 'x'});
    d2rdudk  = symbolic2function('d2rdudk', 'r', {'k' 'u'});
    
    d2ydx2   = symbolic2function('d2ydx2',  'y', {'x' 'x'});
    d2ydu2   = symbolic2function('d2ydu2',  'y', {'u' 'u'});
    d2ydudx  = symbolic2function('d2ydudx', 'y', {'x' 'u'});
    d2ydxdu  = symbolic2function('d2ydxdu', 'y', {'u' 'x'});
end

if order >= 3
    d3fdx3   = symbolic2function('d3fdx3', 'f', {'x' 'x' 'x'});
    d3fdkdx2 = symbolic2function('d3fdkdx2', 'f', {'x' 'x' 'k'});
end

clear symbolic2function % clears the persistent save directory variable from memory

if verbose; fprintf('   ...done.\n'); end

%% Set up model structure

% Clear unnecessary variables from scope
clear SymModel vSyms kSyms sSyms qSyms xuSyms xSyms uSyms x0 ...
    vStrs kStrs sStrs qStrs xuStrs xStrs uStrs ...
    xNamesFull uNamesFull vxNames vuNames ...
    C1Entries C1Values C2Entries C2Values cEntries cValues ...
    nC1Entries nC2Entries ncEntries defaultOpts opts ...
    addlength currentLength iExpr ind match nAdd nExpr uAddInd sys_string ...
    i iv ik is iq iu ix ir iy ...
    rstates rinputs rparams statesubs inputsubs paramssubs rstatesi rinputsi rparamsi...
    uqi qsinui nz sizes jacobianfun reshapefun symbolic2stringmethod...
    thistime defaultMEXdirectory

m = InitializeModel();

m.Type = 'Model.AnalyticReactions';
m.Name = name;

m.Compartments = struct('Name', vNames, 'Dimension', num2cell(dv));
m.Parameters   = struct('Name', kNames, 'Value', num2cell(k));
m.Seeds        = struct('Name', sNames, 'Value', num2cell(s));
m.Inputs       = struct('Name', uNames, 'Compartment', vNames(vuInd));
m.States       = struct('Name', xNames, 'Compartment', vNames(vxInd), 'InitialValue', initial_values);
m.Reactions    = struct('Name', rNames);
% Old code for when y values were not functions
% if isempty(yMembers)
%     Expressions = [];
% else
%     Expressions = mat2cell([vertcat(yMembers{:}), num2cell(vertcat(yValues{:}))], cellfun(@length,yMembers), 2);
% end
% m.Outputs      = struct('Name', yNames, 'Expressions', Expressions );
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

m.u    = setfun_u(u,q);
% Not sure if we still need the following values
% m.q    = q;
% m.dudq = dudq;
% m.nqu  = zeros(nu,1);

m.vxInd = vxInd;
m.vuInd = vuInd;
%rOrder
%krInd

%As
%Bs
%Cs

clear vNames kNames sNames uNames xNames rNames yNames yMembers yValues

m.f         = setfun_rf(f,k);

if order >= 1
    m.dfdx      = setfun_rf(dfdx,k);
    m.dfdk      = setfun_rf(dfdk,k);
    m.dfdu      = setfun_rf(dfdu,k);
end

if order >= 2
    m.d2fdx2    = setfun_rf(d2fdx2,k);
    m.d2fdu2    = setfun_rf(d2fdu2,k);
    m.d2fdk2    = setfun_rf(d2fdk2,k);
    m.d2fdudx   = setfun_rf(d2fdudx,k);
    m.d2fdxdu   = setfun_rf(d2fdxdu,k);
    m.d2fdkdx   = setfun_rf(d2fdkdx,k);
    m.d2fdkdu   = setfun_rf(d2fdkdu,k);
    m.d2fdxdk   = setfun_rf(d2fdxdk,k);
    m.d2fdudk   = setfun_rf(d2fdudk,k);
end

if order >= 3
    m.d3fdx3    = setfun_rf(d3fdx3,k);
    m.d3fdkdx2  = setfun_rf(d3fdkdx2,k);
end

m.S = S;
m.r = setfun_rf(r,k);

if order >= 1
    m.drdx      = setfun_rf(drdx,k);
    m.drdk      = setfun_rf(drdk,k);
    m.drdu      = setfun_rf(drdu,k);
end

if order >= 2
    m.d2rdx2    = setfun_rf(d2rdx2,k);
    m.d2rdu2    = setfun_rf(d2rdu2,k);
    m.d2rdk2    = setfun_rf(d2rdk2,k);
    m.d2rdudx   = setfun_rf(d2rdudx,k);
    m.d2rdxdu   = setfun_rf(d2rdxdu,k);
    m.d2rdkdx   = setfun_rf(d2rdkdx,k);
    m.d2rdkdu   = setfun_rf(d2rdkdu,k);
    m.d2rdxdk   = setfun_rf(d2rdxdk,k);
    m.d2rdudk   = setfun_rf(d2rdudk,k);
end

m.y = setfun_y(y,true,ny);

if order >= 1
    m.dydx      = setfun_y(dydx,false,ny);
    m.dydu      = setfun_y(dydu,false,ny);
end

if order >= 2
    m.d2ydx2    = setfun_y(d2ydx2,false,ny);
    m.d2ydu2    = setfun_y(d2ydu2,false,ny);
    m.d2ydudx   = setfun_y(d2ydudx,false,ny);
    m.d2ydxdu   = setfun_y(d2ydxdu,false,ny);
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
        m.f             = setfun_rf(f,k);
        m.r             = setfun_rf(r,k);
        m.y             = setfun_y(y,true,ny);
        
        if order >= 1
            m.dfdx      = setfun_rf(dfdx,k);
            m.dfdk      = setfun_rf(dfdk,k);
            m.dfdu      = setfun_rf(dfdu,k);
            m.drdx      = setfun_rf(drdx,k);
            m.drdk      = setfun_rf(drdk,k);
            m.drdu      = setfun_rf(drdu,k);
            m.dydx      = setfun_y(dydx,false,ny);
            m.dydu      = setfun_y(dydu,false,ny);
        end
        
        if order >= 2
            m.d2fdx2    = setfun_rf(d2fdx2,k);
            m.d2fdu2    = setfun_rf(d2fdu2,k);
            m.d2fdk2    = setfun_rf(d2fdk2,k);
            m.d2fdudx   = setfun_rf(d2fdudx,k);
            m.d2fdxdu   = setfun_rf(d2fdxdu,k);
            m.d2fdkdx   = setfun_rf(d2fdkdx,k);
            m.d2fdkdu   = setfun_rf(d2fdkdu,k);
            m.d2fdxdk   = setfun_rf(d2fdxdk,k);
            m.d2fdudk   = setfun_rf(d2fdudk,k);
            m.d2rdx2    = setfun_rf(d2rdx2,k);
            m.d2rdu2    = setfun_rf(d2rdu2,k);
            m.d2rdk2    = setfun_rf(d2rdk2,k);
            m.d2rdudx   = setfun_rf(d2rdudx,k);
            m.d2rdxdu   = setfun_rf(d2rdxdu,k);
            m.d2rdkdx   = setfun_rf(d2rdkdx,k);
            m.d2rdkdu   = setfun_rf(d2rdkdu,k);
            m.d2rdxdk   = setfun_rf(d2rdxdk,k);
            m.d2rdudk   = setfun_rf(d2rdudk,k);
            m.d2ydx2    = setfun_y(d2ydx2,false,ny);
            m.d2ydu2    = setfun_y(d2ydu2,false,ny);
            m.d2ydxdu   = setfun_y(d2ydxdu,false,ny);
            m.d2ydudx   = setfun_y(d2ydudx,false,ny);
        end
        
        if order >= 3
            m.d3fdx3    = setfun_rf(d3fdx3,k);
            m.d3fdkdx2  = setfun_rf(d3fdkdx2,k);
        end
        
        m.Update = @update;
        
        varargout{1} = m;
    end

%% %%%%%%%%%%%%%%%%%%%%%%%
%%%% Helper functions %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

    function fun = symbolic2function(variable_name, num, dens)
        
        if ischar(dens)
            dens = {dens};
        end
        
        switch symbolic2stringmethod
            case 'efficient'
                string_rep = symbolic2string_new(variable_name, num, dens);
                string_rep = convertSymStrtoIndexedStr(string_rep, num);
                fun = string2fun(string_rep, num, dens);
            case 'original'
                string_rep = symbolic2string_old(variable_name, num, dens);
                string_rep = convertSymStrtoIndexedStr(string_rep, num);
                fun = string2fun(string_rep, num, dens);
            case 'mex'
                
                % Don't create MEX files for u and dudq. Use the
                % 'efficient' method instead
                if any(strcmp(variable_name,{'u';'dudq'}))
                    string_rep = symbolic2string_new(variable_name, num, dens);
                    string_rep = convertSymStrtoIndexedStr(string_rep, num);
                    fun = string2fun(string_rep, num, dens);
                    return
                end
                
                % Prepare inputs to C code generation function
                dsym = eval(variable_name);
                nzlogical = getNonZeroEntries(num, dens);
                nzi = cell(ndims(nzlogical),1);
                nziind = find(nzlogical);
                [nzi{:}] = ind2sub(size(nzlogical),nziind);
                nzi = [nzi{:}];
                nzsizes = [sizes.(num) cellfun(@(den)sizes.(den),dens(:)')];
                
                % Generate mex C code
                getMexReadyCode(dsym,nzi,nzsizes,xSyms,uSyms,kSyms,variable_name,opts.MEXDirectory)
             
                fun = str2func([variable_name 'fun']);
        end
        
    end
    
    function string_rep = symbolic2string_old(variable_name, num, dens)
        
        densizes = cellfun(@(den)sizes.(den),dens);
        dimensions = [sizes.(num) densizes(:)'];
        
        if numel(dimensions) == 1
            dimensions = [dimensions, 1];
        end
        
        if any(dimensions == 0)
            string_rep = ['zeros([' num2str(dimensions) '])'];
        else
            string_rep = evalc(['disp(' variable_name ')']);
        end
        
        string_rep = regexprep(string_rep, '[', '');
        string_rep = regexprep(string_rep, ']', '');
        string_rep = strtrim(string_rep);
    end

    function string_rep = symbolic2string_new(variable_name, num, dens)
        
        % If provided an empty input, return an empty output
        if isempty(variable_name)
            string_rep = '';
            return
        end
        
        % For f, r, u, and y, just use the old method, since they are not sparse
        if any(strcmp(variable_name,{'f';'r';'u';'y'}))
            string_rep = symbolic2string_old(variable_name, num, dens);
            return
        end
        
        % Evaluate variable name to get symbolic array
        variable = eval(variable_name);
        
        % Get nonzero derivative logical matrix
        nzlogical = getNonZeroEntries(num, dens);
        
        % Reshape nzlogical to the same size as that of the provided value
        varsizes = size(variable);
        nzlogical = reshape(nzlogical,varsizes);
        
        % Find symbolic derivatives for nonzero elements. This outputs an
        % oddly-formatted string that needs to be stripped of some extra characters.
        str_unformatted = char(variable(nzlogical(:)));
        
        % Get string elements between square brackets, not including square
        % brackets
        nnzelements = sum(nzlogical(:));
        if nnzelements == 1
            % Special case for one element: the odd formatting isn't used, so just
            % copy the result
            strelements = {str_unformatted};
        else
            [~,~,~,strelements] = regexp(str_unformatted,'(?<=\[)[^\[\]]+(?=\])');
        end
        
        assert(length(strelements) == nnzelements, ['The number of expressions in d' num 'd' dens{:} ' does not match the number of nonzero elements expected for it.'])
        
        strelements = strjoin(strelements,',');
        
        % Convert subscripts into strings
        [isubscripts,jsubscripts] = find(nzlogical);
        isubscriptstrs = strtrim(cellstr(num2str(isubscripts)));
        jsubscriptstrs = strtrim(cellstr(num2str(jsubscripts)));
        isubstring = strjoin(isubscriptstrs,',');
        jsubstring = strjoin(jsubscriptstrs,',');
        
        % Write sparse initialization string, which will fit in the
        % following expression:
        % '@(t,x,u,k) inf2big(nan2zero(sparse([' STRING_REP '])))'
        % To fit this expression, note the square brackets are left off the
        % left and right sides of string_rep.
        string_rep = sprintf([isubstring '],[' jsubstring '],[' strelements '],[' num2str(varsizes(1)) '],[' num2str(varsizes(2))]);
        
    end

    function string_rep = convertSymStrtoIndexedStr(string_rep, num)
        if any(strcmp(num, {'f' 'r' 'y'}))
            % Replace expressions for states, inputs, and parameters with
            % indexed calls to x, u, and k, respectively
            xindexstring = sprintf('x(%d)\n', 1:nx);
            xindexstring = textscan(xindexstring,'%s','Delimiter','\n');
            xindexstring = xindexstring{1};
            uindexstring = sprintf('u(%d)\n', 1:nu);
            uindexstring = textscan(uindexstring,'%s','Delimiter','\n');
            uindexstring = uindexstring{1};
            kindexstring = sprintf('k(%d)\n', 1:nk);
            kindexstring = textscan(kindexstring,'%s','Delimiter','\n');
            kindexstring = kindexstring{1};
            string_rep = regexprep(string_rep, [xStrs; uStrs; kStrs], [xindexstring; uindexstring; kindexstring], 0);
        elseif strcmp(num, 'u')
            % Replace expressions for input parameters with indexed calls to
            % q
            qindexstring = sprintf('q(%d)\n', 1:nq);
            qindexstring = textscan(qindexstring,'%s','Delimiter','\n');
            qindexstring = qindexstring{1};
            string_rep = regexprep(string_rep, qStrs, qindexstring, 0);
        end
    end

    function d2rdkdx_ = calcDerivative(drdx,kSyms)
        if size(drdx,2) == 1
            d2rdkdx_ = jacobian(drdx,kSyms);
        elseif size(drdx,2) > 1
            d2rdkdx_ = calcHigherOrderDerivative(drdx,kSyms);
        end
    end
            

    function d2rdkdx_ = calcHigherOrderDerivative(drdx, kSyms)
        % Note that in the input argument names here, I use the
        % prototypical example of taking the derivative of drdx with
        % respect to k. This means "r" stands for the variable whose
        % derivative is being taken, "x" stands for the first derivative
        % variable, and "k" stands for the second derivative variable.
        nk_ = length(kSyms);
        nr_ = size(drdx,1);
        nxu_ = size(drdx,2);
        d2rdkdxtemp = cell(nxu_,1);
        if verbose; fprintf(repmat(' ',1,25)); end
        for xi = 1:nxu_
            if verbose; fprintf([repmat('\b',1,25) '%25s'],sprintf('%3g%% complete', 50*xi/nxu_)); end
            d2rdkdxtemp{xi} = jacobian(drdx(:,xi), kSyms);
        end
        d2rdkdx_ = d2rdkdxtemp{1};
        for xi = 1:nxu_-1
            if verbose; fprintf([repmat('\b',1,25) '%25s'],sprintf('%3g%% complete', 50+50*xi/nxu_)); end
            d2rdkdx_ = feval(symengine,'linalg::stackMatrix',d2rdkdx_,d2rdkdxtemp{xi+1});
        end
        
    end
    
    function matout = reshape_extraarguments(mat, newsize, ~, ~)
        matout = reshape(mat,newsize);
    end

    function matout = reshapeDerivative(mat, newsize, num, dens)
        oldsize = size(mat);
        nzlogical = getNonZeroEntries(num,dens);
        newsubscripts = cell(length(newsize),1);
        nzindices = find(nzlogical);
        oldsubscripts = cell(length(oldsize),1);
        newsubscripts = cell(length(newsize),1);
        [oldsubscripts{:}] = ind2sub(oldsize,nzindices);
        [newsubscripts{:}] = ind2sub(newsize,nzindices);
        filename = which('initializematrix.mu');
        read(symengine, filename);
        matout = feval(symengine,'initializematrix',mat,[oldsubscripts{:}],newsize(1),newsize(2),[newsubscripts{:}]);
    end

end

function fun = string2fun(string_rep, num, dens)
% Note that string2fun is a subfunction instead of a nested function to
% prevent the anonymous functions created here from saving copies of
% the primary function workspace variables.

if any(strcmp(num, {'f' 'r'}))
    
    % Set up the function handle by evaluating the string
    if isempty(dens)
        fun = eval(['@(t,x,u,k) [' string_rep ']']);
    else
        fun = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' string_rep '])))']);
    end
    
elseif strcmp(num, 'u')
    
    % Set up the function handle by evaluating the string
    if isempty(dens) % u
        fun = eval(['@(t,q) repmat([' string_rep '],1,numel(t))']);
    else % dudq
        fun = eval(['@(t,q) inf2big(nan2zero(sparse([' string_rep '])))']);
    end
    
elseif strcmp(num, 'y')
    
    if isempty(dens)
        fun = eval(['@(t,x,u) [' string_rep ']']);
    else
        fun = eval(['@(t,x,u) inf2big(nan2zero(sparse([' string_rep '])))']);
    end
end

% Convert function to and from function handle to ensure that
% MATLAB recognizes and stores the workspace for the anonymous
% functions. Without this, saving and loading a model will
% cause MATLAB to lose the function.
fun = func2str(fun);
fun = str2func(fun);

end

function fun = setfun_rf(basefun, k)
    fun = @(t,x,u)basefun(t,x,u,k);
end

function fun = setfun_u(basefun, q)
    fun = @(t)basefun(t,q);
end

function fun = setfun_y(basefun,is0order,ny)
    if is0order
        fun = @(t,x,u) vectorize_y(basefun, t, x, u, ny);
    else
        fun = @(t,x,u) basefun(t,x,u);
    end
end

function val = vectorize_y(y, t, x, u, ny)
    nt = numel(t);
    val = zeros(ny,nt);
    for it = 1:nt
        val(:,it) = y(t(it), x(:,it), u(:,it));
    end
end