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
%       .NewJacobianMethod [ {true} | false ]
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
defaultOpts.NewJacobianMethod = true;
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
        Ssym = initSsym(S);
        nr = numel(r);
        
        if isfield(SymModel, 'rNames')
            rNames = SymModel.rNames;
        else
            rNames = repmat({''}, [nr,1]);
        end
    else
        r = sym(zeros(0,1));
        S = zeros(nx,0);
        Ssym = initSsym(S);
        nr = 0;
        rNames = cell(0,1);
    end
else
    % If there is no f, then r and S must be supplied
    r = SymModel.r;
    S = SymModel.S;
    Ssym = initSsym(S);
    nr = numel(r);
    
    if isfield(SymModel, 'rNames')
        rNames = SymModel.rNames;
    else
        rNames = repmat({''}, [nr,1]);
    end
    
    f = Ssym*r;
end

    function Ssym = initSsym(S)
        Sind = find(S ~= 0);
        [Si,Sj] = ind2sub(size(S),Sind);
        Ssize = size(S);
        Ssym = initializeMatrixMupad(Si,Sj,S(Sind),Ssize(1),Ssize(2));
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
new_symbols = sym('sy%dx',[nk+ns+nx+nu,1]);
% new_symbols = sym(zeros(nk+ns+nx+nu,1));
% for i = 1:nk+ns+nx+nu
%     new_symbols(i) = sym(sprintf('sy%dx', i));
% end

kSyms = new_symbols(1:nk);
sSyms = new_symbols(nk+1:nk+ns);
xSyms = new_symbols(nk+ns+1:nk+ns+nx);
uSyms = new_symbols(nk+ns+nx+1:nk+ns+nx+nu);
f = fastsubs(f, old_symbols, new_symbols);
r = fastsubs(r, old_symbols, new_symbols);
x0 = fastsubs(x0, old_symbols, new_symbols);
y =  fastsubs(y, old_symbols, new_symbols);

% String representations that will be useful
vStrs = fastchar(vSyms);
% vStrs = cell(nv,1);
% for iv = 1:nv
%     vStrs{iv} = char(vSyms(iv));
% end

kStrs = fastchar(kSyms);
% kStrs = cell(nk,1);
% for ik = 1:nk
%     kStrs{ik} = char(kSyms(ik));
% end

sStrs = fastchar(sSyms);
% sStrs = cell(ns,1);
% for is = 1:ns
%     sStrs{is} = char(sSyms(is));
% end

uStrs = fastchar(uSyms);
uNamesFull = cell(nu,1);
for iu = 1:nu
    uNamesFull{iu} = [vNames{vuInd(iu)} '.' uNames{iu}];
end

xStrs = fastchar(xSyms);
xNamesFull = cell(nx,1);
for ix = 1:nx
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

%f = Ssym*r;

q = SymModel.q;
nq = SymModel.nq;
qStrs = cell(nq,1);
qSyms = SymModel.qSyms;

% Standardize names of input parameters
newqSyms = sym('q%dx',[nq,1]);
u = fastsubs(u, qSyms, newqSyms);
qSyms = newqSyms; clear newqSyms
qStrs = fastchar(qSyms);
% 
% for iq = 1:nq
%     u = subs(u, qSyms(iq), newqSym(iq));
%     qSyms(iq) = newqSym;
%     qStrs{iq} = char(qSyms(iq));
% end

%% Determine the species and parameters in each reaction

rstr = fastchar(r);
ystr = fastchar(y);
ustr = fastchar(u);
x0str = fastchar(x0);

rhasx = exprhasvar(rstr,xStrs,nr,nx);
rhasu = exprhasvar(rstr,uStrs,nr,nu);
rhask = exprhasvar(rstr,kStrs,nr,nk);

uhasq = exprhasvar(ustr,qStrs,nu,nq);

yhasx = exprhasvar(ystr,xStrs,ny,nx);
yhasu = exprhasvar(ystr,uStrs,ny,nu);

x0hass = exprhasvar(x0str,sStrs,nx,ns);

    function out = exprhasvar(exprStrs,varStrs,nexprs,nvars)
        
        if isempty(exprStrs) || isempty(varStrs)
            out = false(nexprs,nvars);
        else
        
            % Determine where the variable's substrings appear in the r
            % expressions. Cell array varpositionsinr element i,j contains an
            % array indicating at what positions varStr(j) appears in r(i)
            varpositionsinexpr = cellfun(@strfind,...
                repmat(exprStrs,1,nvars),...
                repmat(varStrs(:)',nexprs,1),...
                'UniformOutput',false...
                );
            
            % Determine which elements are empty. If varpositionsinr{i,j} is
            % not empty, then reaction i contains variable j.
            out = ~cellfun(@isempty,varpositionsinexpr);
            
        end
        
    end

% Set up sizes struct
sizes.f = nx;
sizes.r = nr;
sizes.k = nk;
sizes.x = nx;
sizes.u = nu;
sizes.q = nq;
sizes.y = ny;

% The above logical arrays are the nonzero elements of the first derivative
% matrix for r, u, and y. Record these in a map, and construct a
% function that can retrieve which elements in a derivative matrix are 
% nonzero or, in the absence of this information, calculate whether 
% higher-order derivatives can contain a given parameter, input, or state 
% based on the lower derivatives.
nz = containers.Map;
nz('r') = true(nr,1);
nz('f') = true(nx,1);
nz('y') = true(ny,1);
nz('u') = true(nu,1);
nz('rx') = rhasx;
nz('ru') = rhasu;
nz('rk') = rhask;
nz('fx') = logical(abs(S)*rhasx); % Take the absolute value of S so that there are no accidental cancellations between positive and negative terms
nz('fu') = logical(abs(S)*rhasu);
nz('fk') = logical(abs(S)*rhask);
nz('uq') = uhasq;
nz('yx') = yhasx;
nz('yu') = yhasu;
nz('xs') = x0hass;


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
        
        % Standardize denominators as a cell array
        if ischar(dens)
            dens = {dens};
        end
        
        % Find the highest-order derivative we have nonzero entry data for
        nzkey = [num dens{:}];
        nzout = [];
        nzfound = false;
        while ~nzfound
            if isKey(nz,nzkey)
                nzout = nz(nzkey);
                nzfound = true;
            else
                nzkey = nzkey(1:end-1);
            end
            if isempty(nzkey)
                error('No data found for this derivative. This should never happen and needs to be debugged.')
            end
        end
        
        ndens = length(dens);
        ndensknown = length(nzkey) - 1;
        denstoestimate = (ndensknown+1):ndens;
        singletondims = 3:(ndens+1);
        for nzi = denstoestimate(:)'
            nzaddkey = [num dens{nzi}];
            nztoadd = nz(nzaddkey);
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

%% Analyze common subexpressions in reaction rates
% Reaction rates will typically consist of a few rate forms repeated many
% times. I identify repeated rate forms here and only take the symbolic
% derivatives of each rate form once.
%
% This isn't finished yet.

% Add indicators of variable type (parameter, state, input) to
% standardized names by replacing y with x, u, or k
% xStrs_types = strrep(xStrs,'y','x');
% uStrs_types = strrep(uStrs,'y','u');
% kStrs_types = strrep(kStrs,'y','k');
% 
% % Replace variable names with indicator names
% rstr_types = regexprep(rstr,[xStrs;uStrs;kStrs],[xStrs_types;uStrs_types;kStrs_types]);
% 
% % Determine the order of appearance of each variable type in each
% % expression
% vartypes = {'x','u','k'};
% rstr_forms = rstr_types;
% for vi = 1:length(vartypes)
%     vartype = vartypes{vi};
%     rstr_forms = getRateForms_onevartype(rstr_forms,vartype);
% end
% 
%     function rstr_forms = getRateForms_onevartype(rstr_types,vartype)
%         rstr_varStrs = regexp(rstr_types,['s' vartype '\d+x'],'match');
%         lengths = cellfun(@length,rstr_varStrs);
%         rstr_repvarStrs = strsplit(strtrim(sprintf(['s' vartype '%dx\n'],1:max(lengths))))';
%         rstr_forms = cell(sizes.r,1);
%         for ri = 1:sizes.r
%             rstr_forms{ri} = regexprep(rstr_types{ri},rstr_varStrs{ri},rstr_repvarStrs(1:lengths(ri)));
%         end
%     end

% Find rates with same form 

%% Generate derivatives of desired order

% Determine which method to use to calculate Jacobians
if opts.NewJacobianMethod
    jacobianfun = @calcDerivative;
    reshapefun = @reshapeDerivative;
else
    jacobianfun = @(x, syms, num, dens) jacobian(vec(x), syms);
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
    dfdx = Ssym*drdx;
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of f with respect to u
    if verbose; fprintf('Calculating dfdu...'); end
    dfdu = Ssym*drdu;
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of f with respect to k
    if verbose; fprintf('Calculating dfdk...'); end
    dfdk = Ssym*drdk;
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of y with respect to x
    if verbose; fprintf('Calculating dydx...'); end
    dydx = jacobianfun(y, xSyms);
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of y with respect to u
    if verbose; fprintf('Calculating dydu...'); end
    dydu = jacobianfun(y, uSyms);
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
    d2rdx2 = jacobianfun(drdx, xSyms, 'r', {'x','x'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of drdu with respect to u
    if verbose; fprintf('Calculating d2rdu2...'); end
    d2rdu2 = jacobianfun(drdu, uSyms, 'r', {'u','u'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of drdu with respect to x
    if verbose; fprintf('Calculating d2rdxdu...'); end
    d2rdxdu = jacobianfun(drdu, xSyms, 'r', {'u','x'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of drdx with respect to u
    if verbose; fprintf('Calculating d2rdudx...'); end
    d2rdudx = jacobianfun(drdx, uSyms, 'r', {'x','u'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of drdk with respect to k
    if verbose; fprintf('Calculating d2rdk2...'); end
    d2rdk2 = jacobianfun(drdk, kSyms, 'r', {'k','k'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of drdx with respect to k
    if verbose; fprintf('Calculating d2rdkdx...'); end
    d2rdkdx = jacobianfun(drdx, kSyms, 'r', {'x','k'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of drdu with respect to k
    if verbose; fprintf('Calculating d2rdkdu...'); end
    d2rdkdu = jacobianfun(drdu, kSyms, 'r', {'u','k'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of drdk with respect to x
    if verbose; fprintf('Calculating d2rdxdk...'); end
    d2rdxdk = jacobianfun(drdk, xSyms, 'r', {'k','x'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of drdk with respect to u
    if verbose; fprintf('Calculating d2rdudk...'); end
    d2rdudk = jacobianfun(drdk, uSyms, 'r', {'k','u'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of dfdx with respect to x
    if verbose; fprintf('Calculating d2fdx2...'); end
    d2fdx2 = Ssym*reshapefun(d2rdx2, [nr nx*nx], 'r', {'x' 'x'});
    d2fdx2 = reshapefun(d2fdx2, [nx*nx nx], 'f', {'x' 'x'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of dfdu with respect to u
    if verbose; fprintf('Calculating d2fdu2...'); end
    d2fdu2 = Ssym*reshapefun(d2rdu2, [nr,nu*nu], 'r', {'u' 'u'});
    d2fdu2 = reshapefun(d2fdu2, [nx*nu,nu], 'f', {'u' 'u'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of dfdu with respect to x
    if verbose; fprintf('Calculating d2fdxdu...'); end
    d2fdxdu = Ssym*reshapefun(d2rdxdu, [nr,nu*nx], 'r', {'u' 'x'});
    d2fdxdu = reshapefun(d2fdxdu, [nx*nu,nx], 'f', {'u' 'x'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of dfdx with respect to u
    if verbose; fprintf('Calculating d2fdudx...'); end
    d2fdudx = Ssym*reshapefun(d2rdudx, [nr,nx*nu], 'r', {'x' 'u'});
    d2fdudx = reshapefun(d2fdudx, [nx*nx,nu], 'f', {'x' 'u'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of dfdk with respect to k
    if verbose; fprintf('Calculating d2fdk2...'); end
    d2fdk2 = Ssym*reshapefun(d2rdk2, [nr,nk*nk], 'r', {'k' 'k'});
    d2fdk2 = reshapefun(d2fdk2, [nx*nk,nk], 'f', {'k' 'k'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of dfdx with respect to k
    if verbose; fprintf('Calculating d2fdkdx...'); end
    d2fdkdx = Ssym*reshapefun(d2rdkdx, [nr,nx*nk], 'r', {'x' 'k'});
    d2fdkdx = reshapefun(d2fdkdx, [nx*nx,nk], 'f', {'x' 'k'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of dfdu with respect to k
    if verbose; fprintf('Calculating d2fdkdu...'); end
    d2fdkdu = Ssym*reshapefun(d2rdkdu, [nr,nu*nk], 'r',{'u' 'k'});
    d2fdkdu = reshapefun(d2fdkdu, [nx*nu,nk], 'f',{'u' 'k'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of dfdk with respect to x
    if verbose; fprintf('Calculating d2fdxdk...'); end
    d2fdxdk = Ssym*reshapefun(d2rdxdk, [nr,nk*nx], 'r',{'k' 'x'});
    d2fdxdk = reshapefun(d2fdxdk, [nx*nk,nx], 'f',{'k' 'x'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of dfdk with respect to u
    if verbose; fprintf('Calculating d2fdudk...'); end
    d2fdudk = Ssym*reshapefun(d2rdudk, [nr,nk*nu], 'r',{'k' 'u'});
    d2fdudk = reshapefun(d2fdudk, [nx*nk,nu], 'f',{'k' 'u'});
    if verbose; fprintf('Done.\n'); end
    
    % Output's second derivatives
    if verbose; fprintf('Calculating d2ydx2...'); end
    d2ydx2 = jacobianfun(dydx, xSyms, 'y', {'x','x'});
    if verbose; fprintf('Done.\n'); end
    
    if verbose; fprintf('Calculating d2ydu2...'); end
    d2ydu2 = jacobianfun(dydu, uSyms, 'y', {'u','u'});
    if verbose; fprintf('Done.\n'); end
    
    if verbose; fprintf('Calculating d2ydxdu...'); end
    d2ydxdu = jacobianfun(dydu, xSyms, 'y', {'u','x'});
    if verbose; fprintf('Done.\n'); end
    
    if verbose; fprintf('Calculating d2ydudx...'); end
    d2ydudx = jacobianfun(dydx, uSyms, 'y', {'x','u'});
    if verbose; fprintf('Done.\n'); end
    
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
    d3rdx3 = jacobianfun(vec(d2rdx2), xSyms, 'r', {'x','x','x'});
    
    %Gradient of d2rdx2 with respect to k
    d3rdkdx2 = jacobianfun(vec(d2rdx2), kSyms, 'r', {'x','x','k'});
    
    %Gradient of d2fdx2 with respect to x
    d3fdx3 = Ssym*reshapefun(d3rdx3, [nr,nx*nx*nx], 'r',{'x' 'x' 'x'});
    d3fdx3 = reshapefun(d3fdx3, [xn*nx*nx,nx], 'f',{'x' 'x' 'x'});
    
    %Gradient of d2fdx2 with respect to k
    d3fdkdx2 = Ssym*reshapefun(d3rdkdx2, [nr,nx*nx*nk], 'r',{'x' 'x' 'k'});
    d3fdkdx2 = reshapefun(d3fdkdx2, [nx*nx*nx,nk], 'f',{'x' 'x' 'k'});
else
    d3rdx3   = '';
    d3rdkdx2 = '';
    d3fdx3   = '';
    d3fdkdx2 = '';
end

%% Extract seed information

% Take derivative of initial condition expressions with respect to seed
% parameters
dx0ds = jacobianfun(x0, sSyms);

% Convert nonzero terms to doubles (we currently expect only linear
% combinations of seeds with added constant factors for each x0)
dx0ds_logical = getNonZeroEntries('x','s');
dx0ds_index = find(dx0ds_logical(:));
dx0ds_nz = double(dx0ds(dx0ds_index));

% Initialize sparse matrix with these double values
[dx0ds_i,dx0ds_j] = find(dx0ds_logical);
dx0ds = sparse(dx0ds_i,dx0ds_j,dx0ds_nz,nx,ns);

% Calculate constant terms of ICs by setting seeds to 0
x0c = double(fastsubs(x0,sSyms,zeros(size(sSyms))));

% Store information about seeds and ICs
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
f        = symbolic2function('f', 'f', {});
r        = symbolic2function('r', 'r', {});
y        = symbolic2function('y', 'y', {});

if order >= 1
    dudq     = symbolic2function('dudq', 'u', 'q');
    
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

m.y = setfun_y(y,true,ny,nk);

if order >= 1
    m.dydx      = setfun_y(dydx,false,ny,nk);
    m.dydu      = setfun_y(dydu,false,ny,nk);
end

if order >= 2
    m.d2ydx2    = setfun_y(d2ydx2,false,ny,nk);
    m.d2ydu2    = setfun_y(d2ydu2,false,ny,nk);
    m.d2ydudx   = setfun_y(d2ydudx,false,ny,nk);
    m.d2ydxdu   = setfun_y(d2ydxdu,false,ny,nk);
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
        m.y             = setfun_y(y,true,ny,nk);
        
        if order >= 1
            m.dfdx      = setfun_rf(dfdx,k);
            m.dfdk      = setfun_rf(dfdk,k);
            m.dfdu      = setfun_rf(dfdu,k);
            m.drdx      = setfun_rf(drdx,k);
            m.drdk      = setfun_rf(drdk,k);
            m.drdu      = setfun_rf(drdu,k);
            m.dydx      = setfun_y(dydx,false,ny,nk);
            m.dydu      = setfun_y(dydu,false,ny,nk);
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
            m.d2ydx2    = setfun_y(d2ydx2,false,ny,nk);
            m.d2ydu2    = setfun_y(d2ydu2,false,ny,nk);
            m.d2ydxdu   = setfun_y(d2ydxdu,false,ny,nk);
            m.d2ydudx   = setfun_y(d2ydudx,false,ny,nk);
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
        
        if verbose; fprintf([variable_name '...']); end
        
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
                    if verbose; fprintf('Done.\n'); end
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
        
        if verbose; fprintf('Done.\n'); end
        
    end
    
    function string_rep = symbolic2string_old(variable_name, num, dens)
        
        densizes = cellfun(@(den)sizes.(den),dens);
        dimensions = [sizes.(num) densizes(:)'];
        
        if numel(dimensions) == 1
            dimensions = [dimensions, 1];
        end
        
        if any(dimensions == 0)
            string_rep = ['zeros([' num2str(dimensions,['%d' repmat(',%d',1,length(dimensions)-1)]) '])'];
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
        nzindices = find(nzlogical(:));
        strelements = fastchar(variable(nzindices)); % Don't replace indices with logicals here. The conversion of the logical to sym takes too long.
        
%         % Get string elements between square brackets, not including square
%         % brackets
%         nnzelements = sum(nzlogical(:));
%         if nnzelements == 1
%             % Special case for one element: the odd formatting isn't used, so just
%             % copy the result
%             strelements = {str_unformatted};
%         else
%             [~,~,~,strelements] = regexp(str_unformatted,'(?<=\[)[^\[\]]+(?=\])');
%         end
%         
%         assert(length(strelements) == nnzelements, ['The number of expressions in d' num 'd' dens{:} ' does not match the number of nonzero elements expected for it.'])
%         
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
            if nx == 0
                xindexstring = {};
            else
                xindexstring = sprintf('x(%d)\n', 1:nx);
                xindexstring = textscan(xindexstring,'%s','Delimiter','\n');
                xindexstring = xindexstring{1};
            end
            if nu == 0
                uindexstring = {};
            else
                uindexstring = sprintf('u(%d)\n', 1:nu);
                uindexstring = textscan(uindexstring,'%s','Delimiter','\n');
                uindexstring = uindexstring{1};
            end
            if nk == 0
                kindexstring = {};
            else
                kindexstring = sprintf('k(%d)\n', 1:nk);
                kindexstring = textscan(kindexstring,'%s','Delimiter','\n');
                kindexstring = kindexstring{1};
            end
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

    function d2rdkdx_ = calcDerivative(drdx,kSyms,num,dens)
        if size(drdx,2) <= 1
            d2rdkdx_ = jacobian(drdx,kSyms);
        elseif size(drdx,2) > 1
            d2rdkdx_ = calcHigherOrderDerivative(drdx,kSyms,num,dens);
        end
    end
            

    function d2rdkdx_ = calcHigherOrderDerivative(drdx_, kSyms_, num, dens)
        % Note that in the input argument names here, I use the
        % prototypical example of taking the derivative of drdx with
        % respect to k. This means "r" stands for the variable whose
        % derivative is being taken, "x" stands for the first derivative
        % variable, and "k" stands for the second derivative variable.
        %
        % num and dens should be for the output derivative matrix.
        
        nk_ = length(kSyms_);
        nr_ = size(drdx_,1);
        nxu_ = size(drdx_,2);
        
        % Standardize dens as a cell array
        if ischar(dens)
            dens = {dens};
        end
        
        % Find the entries of d2rdkdx that might be nonzero
        nze = getNonZeroEntries(num,dens);
        nze = reshape(nze,nr_*nxu_,nk_);
        [nzterms,nzdens] = find(nze);
        
        % Take derivatives of the possibly nonzero entries
        nzders = diff_vectorized(drdx_(nzterms),kSyms_(nzdens));
        
        % Of the supposedly nonzero derivatives, find the ones that are
        % actually nonzero, and only keep those
        %iszero = cellfun(@(nzders)logical(nzders == 0),nzders);
        iszero = logical(nzders == 0);
        
        if all(iszero) || isempty(iszero)
            d2rdkdx_ = initializeMatrixMupad([],[],[],nr_*nxu_, nk_);
            return
        end
        
        % Update the nz map with the new information
        nzeiszero = sub2ind([nr_*nxu_,nk_],nzterms(iszero),nzdens(iszero));
        nze(nzeiszero) = false;
        nze = reshape(nze,nr_,nxu_,nk_);
        nzkey = [num dens{:}];
        nz(nzkey) = nze;
        % Also get information about f derivative, if numerator is r
        if strcmp(num,'r')
            nzkey_f = strrep(nzkey,'r','f');
            % Check that we don't already know more about f than we can
            % gather from the current derivative
            if ~isKey(nz,nzkey_f)
                nztemp = logical(abs(S)*reshape(nze,[nr_,nxu_*nk_]));
                nz(nzkey_f) = reshape(nztemp,[nx,nxu_,nk_]);
            end
        end
        
        % Remove zero terms
        nzders(iszero) = [];
        nzterms(iszero) = [];
        nzdens(iszero) = [];
        
        d2rdkdx_ = initializeMatrixMupad(nzterms, nzdens, nzders, nr_*nxu_, nk_);
        
    end
    
    function d2rdkdx_ = calcHigherOrderDerivative_old(drdx, kSyms)
        
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
        
        % Get logical matrix indicating non-zero entries
        nzlogical = getNonZeroEntries(num,dens);
        
        % Convert logicals into linear indices
        nzindices = find(nzlogical);
        
        % Convert linear indices into subscripts in new reshaped matrix
        [nzsub_i,nzsub_j] = ind2sub(newsize,nzindices);
        
        % Get nonzero terms
        nzterms = mat(nzindices);
        
        % Initialize matrix 
        matout = initializeMatrixMupad(nzsub_i,nzsub_j,nzterms,newsize(1),newsize(2));
        
    end

    function out = fastsubs(f, old, new)
        filename = which('fastsubs.mu');
        read(symengine, filename);
        out = feval(symengine,'fastsubs',f,old,new);
    end

    function fstr = fastchar(f)
        filename = which('fastchar.mu');
        read(symengine, filename);
        fstr = feval(symengine,'fastchar',f);
        fstr = eval(fstr);
    end

    function out = getsymelements(s,indices)
        filename = which('getsymelements.mu');
        read(symengine, filename);
        out = feval(symengine,'getsymelements',s(:),indices);
    end

    function matout = initializeMatrixMupad(i,j,s,m,n)
        % More efficient method for initializing sparse symbolic matrices.
        % Utilizes a MuPAD function that generates a MuPAD table,
        % then uses the table to initialize the matrix.
        filename = which('sparse_mupad.mu');
        read(symengine, filename);
        matout = feval(symengine,'sparse_mupad',i, j, s, m, n);
    end

    function D = diff_vectorized(nums, dens)
        filename = which('diff_vectorized.mu');
        read(symengine, filename);
        D = feval(symengine,'diff_vectorized',nums,dens);
    end

    function nzstring = generateNonZeroElementCCode(mat,num,dens)
        nze = getNonZeroEntries(num,dens);
        [i,j] = find(reshape(nze,size(mat)));
        filename = which('generateNonZeroElementCCode.mu');
        read(symengine, filename);
        nzstring = char(feval(symengine,'generateNonZeroElementCCode',i,j,mat,xSyms,uSyms,kSyms));
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
        fun = eval(['@(t,x,u,k) [' string_rep ']']);
    else
        fun = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' string_rep '])))']);
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

function fun = setfun_y(basefun,is0order,ny,nk)
% nk is needed to initialize a dummy k vector, which is needed to avoid
% having to treat the y functions differently in the MEX generating code.
% This could be done better. I will leave it like this to make it easier to
% add k values to outputs later, if desired.
    if is0order
        fun = @(t,x,u) vectorize_y(basefun, t, x, u, ny, nk);
    else
        fun = @(t,x,u) basefun(t,x,u,ones(nk,1)); % Currently, output dependence on k is not supported
    end
end

function val = vectorize_y(y, t, x, u, ny, nk)
    nt = numel(t);
    val = zeros(ny,nt);
    for it = 1:nt
        val(:,it) = y(t(it), x(:,it), u(:,it), ones(nk,1));
    end
end