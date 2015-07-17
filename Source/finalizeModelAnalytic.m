function m = finalizeModelAnalytic(m, opts)
%finalizeModelAnalytic prepares a constructed Model.Analytic for use with the
%rest of kroneckerbio. Generates symbolic expressions and requires the symbolic
%toolbox.
% 
%   m = finalizeModelAnalytic(m, opts)
% 
%   Inputs
%   m: [ Model.Analytic ]
%       Analytic model with components to add in the m.add.* fields
%   opts: [ options struct ]
%       Options struct with the following optional fields:
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
%       .EvaluateExternalFunctions [ logical scalar {false} ]
%           Determines whether to evaluate calls to external functions in
%           the reaction rate expressions. The external functions are
%           evaluated with symbolic input arguments to obtain a symbolic
%           version of the output that can be differentiated. If set to
%           false, derivatives of external function calls cannot be
%           computed. Note: required to evaluate exponents written usen the pow
%           function. Try setting this to true if pow functions aren't
%           recognized in final symbolic expressions and function handles.
%
%   Outputs
%   m: [ Model.Analytic ]
%       The useable form of the model

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
opts_.Order                     = 2;
opts_.VolumeToParameter         = false;
opts_.Verbose                   = 0;
opts_.UseMEX                    = false;
opts_.MEXDirectory              = defaultMEXdirectory;
opts_.EvaluateExternalFunctions = false; % needed for calls to power() and other functions

opts = mergestruct(opts_, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

order = opts.Order;

if opts.UseMEX && exist(opts.MEXDirectory,'dir') ~= 7
    mkdir(opts.MEXDirectory);
end

%% Extract symbolic values from IDs
% Compartments
m.Compartments = combineComponents(m.Compartments, m.add.Compartments, opts.Verbose);
vStrs  = {m.Compartments.ID}';
vSyms  = sym(vStrs);
vNames = {m.Compartments.Name}';
dv     = [m.Compartments.Dimension]';
v      = [m.Compartments.Size]';
nv     = numel(v);

% Parameters
m.Parameters = combineComponents(m.Parameters, m.add.Parameters, opts.Verbose);
kStrs  = {m.Parameters.ID}';
kSyms  = sym(kStrs);
kNames = {m.Parameters.Name}';
k      = [m.Parameters.Value]';
nk     = numel(k);

% Seeds
m.Seeds = combineComponents(m.Seeds, m.add.Seeds, opts.Verbose);
sStrs  = {m.Seeds.ID}';
sSyms  = sym(sStrs);
sNames = {m.Seeds.Name}';
s      = [m.Seeds.Value]';
ns     = numel(s);

% States
m.States = combineComponents(m.States, m.add.States, opts.Verbose);
xStrs   = {m.States.ID}';
xSyms   = sym(xStrs);
xNames  = {m.States.Name}';
xvNames = {m.States.Compartment}';
x0      = {m.States.InitialValue}';
nx      = numel(x0);

% Inputs
m.Inputs = combineComponents(m.Inputs, m.add.Inputs, opts.Verbose);
uStrs   = {m.Inputs.ID}';
uSyms   = sym(uStrs);
uNames  = {m.Inputs.Name}';
uvNames = {m.Inputs.Compartment}';
u       = [m.Inputs.DefaultValue]';
nu      = numel(u);

% Species compartments
[~, vxInd] = ismember(xvNames, vNames);
[~, vuInd] = ismember(uvNames, vNames);

% Reactions
rprovided = true;
m.Reactions = combineComponents(m.Reactions, m.add.Reactions, opts.Verbose);
rNames = {m.Reactions.Name}';
r = {m.Reactions.Rate}';
nr = numel(r);

% Outputs
m.Outputs = combineComponents(m.Outputs, m.add.Outputs, opts.Verbose);
yStrs  = {m.Outputs.ID}';
ySyms  = sym(yStrs);
yNames = {m.Outputs.Name}';
y      = {m.Outputs.Expression}';
ny     = numel(y);

% Initialize outputs properly if no outputs specified in model
if ny == 0
    yStrs  = cell(0,1);
    ySyms  = sym(yStrs);
    yNames = cell(0,1);
    y      = cell(0,1);
end

% Rules
m.Rules = combineComponents(m.Rules, m.add.Rules, opts.Verbose);
zStrs    = {m.Rules.ID}';
zSyms    = sym(zStrs);
zNames   = {m.Rules.Name}';
z        = {m.Rules.Expression}';
zTargets = {m.Rules.Target}';
zTypes   = {m.Rules.Type}';
nRules   = numel(z); % nz clashes with another variable below

% Useful combinations of names and IDs
outputNames = [xNames; uNames; vNames; kNames];
outputIDs   = [xStrs; uStrs; vStrs; kStrs];
allNames    = [outputNames; sNames];
allIDs      = [outputIDs; sStrs];
xuvNames = [xvNames; uvNames];

%% Clean up m.add.* fields
m.add.Compartments = growCompartmentsAnalytic;
m.add.Parameters   = growParametersAnalytic;
m.add.Seeds        = growSeedsAnalytic;
m.add.Inputs       = growInputsAnalytic;
m.add.States       = growStatesAnalytic;
m.add.Reactions    = growReactionsAnalytic;
m.add.Outputs      = growOutputsAnalytic;
m.add.Rules        = growRulesAnalytic;
m.add.nv = 0;
m.add.nk = 0;
m.add.ns = 0;
m.add.nu = 0;
m.add.nx = 0;
m.add.nr = 0;
m.add.nz = 0;

%% Process expressions and turn strings into symbolics
% Initial conditions - arbitrary expressions of seeds
for i = 1:nx
    x0{i} = name2id(x0{i}, sNames, sStrs);
end

% Reaction rates - arbitrary expressions of anything
for i = 1:nr
    r{i} = name2id(r{i}, allNames, allIDs, xuvNames);
end

% Outputs - arbitrary expressions of states and parameters
for i = 1:ny
    y{i} = name2id(y{i}, outputNames, outputIDs, xuvNames);
end

% Convert to symbolics
x0 = sym(x0);
r = sym(r);
y = sym(y);

%% Perform rule substitutions
% Note/TODO: zTargets substitution probably too permissive - lets anything be replaced, but backend won't handle this properly
for i = 1:nRules
    zTargets{i} = name2id(zTargets{i}, allNames, allIDs, xuvNames);
    z{i} = name2id(z{i}, allNames, allIDs, xuvNames);
end
zTargets = sym(zTargets);
z = sym(z);

% Repeated assignment subs in reaction rates
%   Makes n repeated assignment rules passes, O(n^2) to ensure all subs
zRAInds = strcmpi(zTypes, 'repeated assignment');
zTargetsRA = zTargets(zRAInds);
zRA        = z(zRAInds);
for i = 1:length(zTargetsRA)
    r = subs(r, zTargetsRA, zRA);
end

% Repeated assignment subs in initial assignment rules
zIAInds = strcmpi(zTypes, 'initial assignment');
zTargetsIA = zTargets(zIAInds);
zIA        = z(zIAInds);
for i = 1:length(zTargetsRA)
    zTargetsIA = subs(zTargetsIA, zTargetsRA, zRA);
end

% Initial assignment and repeated assignment subs in initial conditions
%   Note: kronecker currently separates k (params valid at all times) from s
%   (seeds valid only at initial time). If a model has initial assignment rules
%   that depend on k, the initial condition will respect substituting k.
%   However, if k changes, the entire model needs to be rebuilt as during
%   fitting, and the constant initial condition will no longer be valid.
%   Workaround: make 2 parameters, 1 k and 1 s, and set them equal with linear
%   equality constraints when fitting.
%   Note/TODO: inputs can't be assigned with initial assignment rules yet. SBML
%   models with states set to constant that have initial assignment rules won't
%   convert correctly for now. Should be a quick change to assign input default
%   values, though.
%   TODO: consider merging k and s
[~, IAInd] = ismember(char(zTargetsIA), xStrs);
subTarget = [vSyms; kSyms; xSyms; uSyms];
subExpr   = [sym(v); sym(k); x0; sym(u)];
for i = 1:length(zTargetsIA);
    x0(IAInd(i)) = subs(zIA(i), subTarget, subExpr);
end

%% Evaluate external functions
if opts.EvaluateExternalFunctions
    r = evaluateExternalFunctions(r, [vStrs; xStrs; uStrs; kStrs]);
end

%% Process stoichiometry and rate forms/RHS's
% Make stoichiometry matrix nSpecies x nReactions
xuFullNames = strcat([xvNames; uvNames], '.', [xNames; uNames]);
SEntries = zeros(0,3);
for i = 1:nr
    for j = 1:length(m.Reactions(i).Reactants)
        if ~isempty(m.Reactions(i).Reactants{j})
            [~, ind] = ismember(m.Reactions(i).Reactants{j}, xuFullNames);
            SEntries = [SEntries; ind, i, -1];
        end
    end
    for j = 1:length(m.Reactions(i).Products)
        if ~isempty(m.Reactions(i).Products{j})
            [~, ind] = ismember(m.Reactions(i).Products{j}, xuFullNames);
            SEntries = [SEntries; ind, i, 1];
        end
    end
end
S = sparse(SEntries(:,1), SEntries(:,2), SEntries(:,3), nx+nu, nr);
Su = S(nx+1:end,:);
S = S(1:nx,:);

StoichiometricMatrix = S;

StoichiometricSym = initSsym(StoichiometricMatrix);
f = StoichiometricSym*r;

    function StoichiometricSym = initSsym(StoichiometricMatrix)
        Sind = find(StoichiometricMatrix(:) ~= 0);
        [Si,Sj] = ind2sub(size(StoichiometricMatrix),Sind);
        Ssize = size(StoichiometricMatrix);
        Ss = StoichiometricMatrix(Sind);
        StoichiometricSym = initializeMatrixMupad(Si,Sj,Ss,Ssize(1),Ssize(2));
    end

%% Convert compartment volumes to parameters or constants
if opts.VolumeToParameter
    nk = nk + nv;
    kNames = [kNames; vNames];
    k = [k; v];
    kSyms = [kSyms; vSyms];
else
    r = subs(r, vSyms, v);
    f = subs(f, vSyms, v);
    y = subs(y, vSyms, v);
end

%% Sanity checks
% TODO: make sure initial conditions are only functions of seeds

% Check original symbols for reserved MuPAD names
% Note: shouldn't be needed since symbols are usually autogenerated UUIDs
old_symbols_strs = fastchar([kSyms; sSyms; xSyms; uSyms]);
reservednames = {'pi','PI','eulergamma','EULER','catalan','CATALAN'};
for rni = 1:length(reservednames)
    isreservedname = strcmp(old_symbols_strs,reservednames{rni});
    if any(isreservedname)
        warning(['A symbolic variable in the model was named ' ...
            reservednames{rni} ', which is a reserved variable name. ' ...
            'This can cause undesired behavior in rate expressions. ' ...
            'It is recommended to change the symbolic variable''s name.'])
    end
end

%% Determine the species and parameters in each reaction

rstr = fastchar(r);
ystr = fastchar(y);

x0str = fastchar(x0);

rhasx = exprhasvar(rstr,xStrs,nr,nx);
rhasu = exprhasvar(rstr,uStrs,nr,nu);
rhask = exprhasvar(rstr,kStrs,nr,nk);

if ~rprovided
    fstr = fastchar(f);
    fhasx = exprhasvar(fstr,xStrs,nx,nx);
    fhasu = exprhasvar(fstr,uStrs,nx,nu);
    fhask = exprhasvar(fstr,kStrs,nx,nk);
end

yhasx = exprhasvar(ystr,xStrs,ny,nx);
yhasu = exprhasvar(ystr,uStrs,ny,nu);
yhask = exprhasvar(ystr,kStrs,ny,nk);

x0hass = exprhasvar(x0str,sStrs,nx,ns);

    function out = exprhasvar(exprStrs, varStrs, nExprs, nVars)
        
        if isempty(exprStrs) || isempty(varStrs)
            out = false(nExprs, nVars);
        else
        
            % Determine where the variable's substrings appear in the r
            % expressions. Cell array varPositionsInExpr element i,j contains an
            % array indicating at what positions varStr(j) appears in r(i)
            varPositionsInExpr = cellfun(@strfind,...
                repmat(exprStrs, 1, nVars),...
                repmat(varStrs(:)', nExprs, 1),...
                'UniformOutput', false...
                );
            
            % Determine which elements are empty. If varPositionsInExpr{i,j} is
            % not empty, then reaction i contains variable j.
            out = ~cellfun(@isempty, varPositionsInExpr);
            
        end
        
    end

% Set up sizes struct
sizes.f = nx;
sizes.r = nr;
sizes.k = nk;
sizes.x = nx;
sizes.u = nu;
sizes.s = ns;
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
nz('rx') = rhasx;
nz('ru') = rhasu;
nz('rk') = rhask;
if rprovided
    nz('fx') = logical(abs(StoichiometricMatrix)*rhasx); % Take the absolute value of S so that there are no accidental cancellations between positive and negative terms
    nz('fu') = logical(abs(StoichiometricMatrix)*rhasu);
    nz('fk') = logical(abs(StoichiometricMatrix)*rhask);
else
    nz('fx') = fhasx;
    nz('fu') = fhasu;
    nz('fk') = fhask;
end
nz('yx') = yhasx;
nz('yu') = yhasu;
nz('yk') = yhask;
nz('xs') = x0hass;

    function nzout = getNonZeroEntries(num, dens)
        % Inputs:
        %   num: "numerator" of derivative, either 'r', 'f', 'u', or 'y'
        %   dens: "denominator" terms of derivative, any combination of 'x',
        %       'u', and 'k' in a cell array for 'r' and 'f' numerators, 'q'
        %       for 'u' numerators, and 'x', 'u', and 'k' for 'y' numerators
        % Outputs:
        %   nzout: a logical array of size n(r, f, u, or y)-by-n(x or u or
        %       k)-by-n(x or u or k)-by-..., where nzout(i,j,k,...) is true if
        %       d(r or f)/d(x or u or k)(x or u or k)... is potentially
        %       nonzero, based on what terms are present in each of the
        %       reaction or state first derivative terms. Note that the first term
        %       appearing in the "denominator" string of the derivative is the
        %       last dimension, since it is the derivative taken last. I.E.,
        %       dr/dxdk has dimensions nr-by-nk-by-nx.
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

if verbose; fprintf('\n'); end

%% Generate derivatives of desired order
if order >= 1

    % Gradient of r with respect to x
    if verbose; fprintf('Calculating drdx...'); end
    drdx = calcDerivative(r, xSyms, 'r', {'x'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of r with respect to u
    if verbose; fprintf('Calculating drdu...'); end
    drdu = calcDerivative(r, uSyms, 'r', {'u'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of r with respect to k
    if verbose; fprintf('Calculating drdk...'); end
    drdk = calcDerivative(r, kSyms, 'r', {'k'});
    if verbose; fprintf('Done.\n'); end
    
    if rprovided
        
        % Gradient of f with respect to x
        if verbose; fprintf('Calculating dfdx...'); end
        dfdx = StoichiometricSym*drdx;
        if verbose; fprintf('Done.\n'); end

        % Gradient of f with respect to u
        if verbose; fprintf('Calculating dfdu...'); end
        dfdu = StoichiometricSym*drdu;
        if verbose; fprintf('Done.\n'); end

        % Gradient of f with respect to k
        if verbose; fprintf('Calculating dfdk...'); end
        dfdk = StoichiometricSym*drdk;
        if verbose; fprintf('Done.\n'); end
    
    else
        
        % Gradient of f with respect to x
        if verbose; fprintf('Calculating dfdx...'); end
        dfdx = calcDerivative(f, xSyms, 'f', {'x'});
        if verbose; fprintf('Done.\n'); end

        % Gradient of f with respect to u
        if verbose; fprintf('Calculating dfdu...'); end
        dfdu = calcDerivative(f, uSyms, 'f', {'u'});
        if verbose; fprintf('Done.\n'); end

        % Gradient of f with respect to k
        if verbose; fprintf('Calculating dfdk...'); end
        dfdk = calcDerivative(f, kSyms, 'f', {'k'});
        if verbose; fprintf('Done.\n'); end
        
    end
        
    % Gradient of y with respect to x
    if verbose; fprintf('Calculating dydx...'); end
    dydx = calcDerivative(y, xSyms, 'y', {'x'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of y with respect to u
    if verbose; fprintf('Calculating dydu...'); end
    dydu = calcDerivative(y, uSyms, 'y', {'u'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of y with respect to k
    if verbose; fprintf('Calculating dydk...'); end
    dydk = calcDerivative(y, kSyms, 'y', {'k'});
    if verbose; fprintf('Done.\n'); end

    % Gradient of x0 with respect to s
    if verbose; fprintf('Calculating dx0ds...'); end
    dx0ds = calcDerivative(x0, sSyms, 'x', {'s'});
    if verbose; fprintf('Done.\n'); end

else
    drdx = '';
    drdk = '';
    dfdx = '';
    dfdk = '';
    dydx = '';
    dydu = '';
    dydk = '';
    dx0ds = '';
end

if order >= 2
    % Gradient of drdx with respect to x
    if verbose; fprintf('Calculating d2rdx2...'); end
    d2rdx2 = calcDerivative(drdx, xSyms, 'r', {'x','x'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of drdu with respect to u
    if verbose; fprintf('Calculating d2rdu2...'); end
    d2rdu2 = calcDerivative(drdu, uSyms, 'r', {'u','u'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of drdu with respect to x
    if verbose; fprintf('Calculating d2rdxdu...'); end
    d2rdxdu = calcDerivative(drdu, xSyms, 'r', {'u','x'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of drdx with respect to u
    if verbose; fprintf('Calculating d2rdudx...'); end
    d2rdudx = calcDerivative(drdx, uSyms, 'r', {'x','u'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of drdk with respect to k
    if verbose; fprintf('Calculating d2rdk2...'); end
    d2rdk2 = calcDerivative(drdk, kSyms, 'r', {'k','k'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of drdx with respect to k
    if verbose; fprintf('Calculating d2rdkdx...'); end
    d2rdkdx = calcDerivative(drdx, kSyms, 'r', {'x','k'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of drdu with respect to k
    if verbose; fprintf('Calculating d2rdkdu...'); end
    d2rdkdu = calcDerivative(drdu, kSyms, 'r', {'u','k'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of drdk with respect to x
    if verbose; fprintf('Calculating d2rdxdk...'); end
    d2rdxdk = calcDerivative(drdk, xSyms, 'r', {'k','x'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of drdk with respect to u
    if verbose; fprintf('Calculating d2rdudk...'); end
    d2rdudk = calcDerivative(drdk, uSyms, 'r', {'k','u'});
    if verbose; fprintf('Done.\n'); end
    
    if rprovided
        
        % Gradient of dfdx with respect to x
        if verbose; fprintf('Calculating d2fdx2...'); end
        d2fdx2 = StoichiometricSym*reshapeDerivative(d2rdx2, [nr nx*nx], 'r', {'x' 'x'});
        d2fdx2 = reshapeDerivative(d2fdx2, [nx*nx nx], 'f', {'x' 'x'});
        if verbose; fprintf('Done.\n'); end

        % Gradient of dfdu with respect to u
        if verbose; fprintf('Calculating d2fdu2...'); end
        d2fdu2 = StoichiometricSym*reshapeDerivative(d2rdu2, [nr,nu*nu], 'r', {'u' 'u'});
        d2fdu2 = reshapeDerivative(d2fdu2, [nx*nu,nu], 'f', {'u' 'u'});
        if verbose; fprintf('Done.\n'); end

        % Gradient of dfdu with respect to x
        if verbose; fprintf('Calculating d2fdxdu...'); end
        d2fdxdu = StoichiometricSym*reshapeDerivative(d2rdxdu, [nr,nu*nx], 'r', {'u' 'x'});
        d2fdxdu = reshapeDerivative(d2fdxdu, [nx*nu,nx], 'f', {'u' 'x'});
        if verbose; fprintf('Done.\n'); end

        % Gradient of dfdx with respect to u
        if verbose; fprintf('Calculating d2fdudx...'); end
        d2fdudx = StoichiometricSym*reshapeDerivative(d2rdudx, [nr,nx*nu], 'r', {'x' 'u'});
        d2fdudx = reshapeDerivative(d2fdudx, [nx*nx,nu], 'f', {'x' 'u'});
        if verbose; fprintf('Done.\n'); end

        % Gradient of dfdk with respect to k
        if verbose; fprintf('Calculating d2fdk2...'); end
        d2fdk2 = StoichiometricSym*reshapeDerivative(d2rdk2, [nr,nk*nk], 'r', {'k' 'k'});
        d2fdk2 = reshapeDerivative(d2fdk2, [nx*nk,nk], 'f', {'k' 'k'});
        if verbose; fprintf('Done.\n'); end

        % Gradient of dfdx with respect to k
        if verbose; fprintf('Calculating d2fdkdx...'); end
        d2fdkdx = StoichiometricSym*reshapeDerivative(d2rdkdx, [nr,nx*nk], 'r', {'x' 'k'});
        d2fdkdx = reshapeDerivative(d2fdkdx, [nx*nx,nk], 'f', {'x' 'k'});
        if verbose; fprintf('Done.\n'); end

        % Gradient of dfdu with respect to k
        if verbose; fprintf('Calculating d2fdkdu...'); end
        d2fdkdu = StoichiometricSym*reshapeDerivative(d2rdkdu, [nr,nu*nk], 'r',{'u' 'k'});
        d2fdkdu = reshapeDerivative(d2fdkdu, [nx*nu,nk], 'f',{'u' 'k'});
        if verbose; fprintf('Done.\n'); end

        % Gradient of dfdk with respect to x
        if verbose; fprintf('Calculating d2fdxdk...'); end
        d2fdxdk = StoichiometricSym*reshapeDerivative(d2rdxdk, [nr,nk*nx], 'r',{'k' 'x'});
        d2fdxdk = reshapeDerivative(d2fdxdk, [nx*nk,nx], 'f',{'k' 'x'});
        if verbose; fprintf('Done.\n'); end

        % Gradient of dfdk with respect to u
        if verbose; fprintf('Calculating d2fdudk...'); end
        d2fdudk = StoichiometricSym*reshapeDerivative(d2rdudk, [nr,nk*nu], 'r',{'k' 'u'});
        d2fdudk = reshapeDerivative(d2fdudk, [nx*nk,nu], 'f',{'k' 'u'});
        if verbose; fprintf('Done.\n'); end
        
    else % not used
        
        % Gradient of dfdx with respect to x
        if verbose; fprintf('Calculating d2fdx2...'); end
        d2fdx2 = calcDerivative(dfdx, xSyms, 'f', {'x','x'});
        if verbose; fprintf('Done.\n'); end

        % Gradient of dfdu with respect to u
        if verbose; fprintf('Calculating d2fdu2...'); end
        d2fdu2 = calcDerivative(dfdu, uSyms, 'f', {'u','u'});
        if verbose; fprintf('Done.\n'); end

        % Gradient of dfdu with respect to x
        if verbose; fprintf('Calculating d2fdxdu...'); end
        d2fdxdu = calcDerivative(dfdu, xSyms, 'f', {'u','x'});
        if verbose; fprintf('Done.\n'); end

        % Gradient of dfdx with respect to u
        if verbose; fprintf('Calculating d2fdudx...'); end
        d2fdudx = calcDerivative(dfdx, uSyms, 'f', {'x','u'});
        if verbose; fprintf('Done.\n'); end

        % Gradient of dfdk with respect to k
        if verbose; fprintf('Calculating d2fdk2...'); end
        d2fdk2 = calcDerivative(dfdk, kSyms, 'f', {'k','k'});
        if verbose; fprintf('Done.\n'); end

        % Gradient of dfdx with respect to k
        if verbose; fprintf('Calculating d2fdkdx...'); end
        d2fdkdx = calcDerivative(dfdx, kSyms, 'f', {'x','k'});
        if verbose; fprintf('Done.\n'); end

        % Gradient of dfdu with respect to k
        if verbose; fprintf('Calculating d2fdkdu...'); end
        d2fdkdu = calcDerivative(dfdu, kSyms, 'f', {'u','k'});
        if verbose; fprintf('Done.\n'); end

        % Gradient of dfdk with respect to x
        if verbose; fprintf('Calculating d2fdxdk...'); end
        d2fdxdk = calcDerivative(dfdk, xSyms, 'f', {'k','x'});
        if verbose; fprintf('Done.\n'); end

        % Gradient of dfdk with respect to u
        if verbose; fprintf('Calculating d2fdudk...'); end
        d2fdudk = calcDerivative(dfdk, uSyms, 'f', {'k','u'});
        if verbose; fprintf('Done.\n'); end
        
    end
    
    % Output's second derivatives
    if verbose; fprintf('Calculating d2ydx2...'); end
    d2ydx2 = calcDerivative(dydx, xSyms, 'y', {'x','x'});
    if verbose; fprintf('Done.\n'); end
    
    if verbose; fprintf('Calculating d2ydu2...'); end
    d2ydu2 = calcDerivative(dydu, uSyms, 'y', {'u','u'});
    if verbose; fprintf('Done.\n'); end
    
    if verbose; fprintf('Calculating d2ydk2...'); end
    d2ydk2 = calcDerivative(dydk, kSyms, 'y', {'k','k'});
    if verbose; fprintf('Done.\n'); end
    
    if verbose; fprintf('Calculating d2ydudx...'); end
    d2ydudx = calcDerivative(dydx, uSyms, 'y', {'x','u'});
    if verbose; fprintf('Done.\n'); end
    
    if verbose; fprintf('Calculating d2ydkdx...'); end
    d2ydkdx = calcDerivative(dydx, kSyms, 'y', {'x','k'});
    if verbose; fprintf('Done.\n'); end
    
    if verbose; fprintf('Calculating d2ydxdu...'); end
    d2ydxdu = calcDerivative(dydu, xSyms, 'y', {'u','x'});
    if verbose; fprintf('Done.\n'); end
    
    if verbose; fprintf('Calculating d2ydkdu...'); end
    d2ydkdu = calcDerivative(dydu, kSyms, 'y', {'u','k'});
    if verbose; fprintf('Done.\n'); end
    
    if verbose; fprintf('Calculating d2ydxdk...'); end
    d2ydxdk = calcDerivative(dydk, xSyms, 'y', {'k','x'});
    if verbose; fprintf('Done.\n'); end
    
    if verbose; fprintf('Calculating d2ydudk...'); end
    d2ydudk = calcDerivative(dydk, uSyms, 'y', {'k','u'});
    if verbose; fprintf('Done.\n'); end

    if verbose; fprintf('Calculating d2x0ds2...'); end
    d2x0ds2 = calcDerivative(dx0ds, sSyms, 'x', {'s','s'});
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
    d2ydk2 = '';
    d2ydudx = '';
    d2ydkdx = '';
    d2ydxdu = '';
    d2ydkdu = '';
    d2ydxdk = '';
    d2ydudk = '';
    d2x0ds2 = '';
end

if order >= 3
    %Gradient of d2rdx2 with respect to x
    d3rdx3 = calcDerivative(d2rdx2, xSyms, 'r', {'x','x','x'});
    
    %Gradient of d2rdx2 with respect to k
    d3rdkdx2 = calcDerivative(d2rdx2, kSyms, 'r', {'x','x','k'});
    
    %Gradient of d2fdx2 with respect to x
    d3fdx3 = StoichiometricSym*reshapeDerivative(d3rdx3, [nr,nx*nx*nx], 'r',{'x' 'x' 'x'});
    d3fdx3 = reshapeDerivative(d3fdx3, [xn*nx*nx,nx], 'f',{'x' 'x' 'x'});
    
    %Gradient of d2fdx2 with respect to k
    d3fdkdx2 = StoichiometricSym*reshapeDerivative(d3rdkdx2, [nr,nx*nx*nk], 'r',{'x' 'x' 'k'});
    d3fdkdx2 = reshapeDerivative(d3fdkdx2, [nx*nx*nx,nk], 'f',{'x' 'x' 'k'});
    
    % Gradient of d2x0ds2 with respect to s
    d3x0ds3 = calcDerivative(d2x0ds2, sSyms, 'x', {'s','s','s'});
else
    d3rdx3   = '';
    d3rdkdx2 = '';
    d3fdx3   = '';
    d3fdkdx2 = '';
    d3x0ds3  = '';
end

%% Extract seed information
% initial_values = fastchar(fastsubs(x0,sSyms,sNames)); % TODO: what's this for?

%% Convert symbolic expressions to function handles

if opts.UseMEX
    symbolic2stringmethod = 'mex';
else
    symbolic2stringmethod = 'efficient';
end

if verbose; fprintf('Converting symbolics to functions...\n'); end
% for i = 1:length(f)
%     i
%     test1 = symbolic2function(f(i), 'f', {});
% end
f        = symbolic2function(f, 'f', {});
r        = symbolic2function(r, 'r', {});
y        = symbolic2function(y, 'y', {});
x0       = symbolic2function(x0, 'x', {});

if order >= 1
    dfdx     = symbolic2function(dfdx, 'f', 'x');
    dfdu     = symbolic2function(dfdu, 'f', 'u');
    dfdk     = symbolic2function(dfdk, 'f', 'k');

    drdx     = symbolic2function(drdx, 'r', 'x');
    drdu     = symbolic2function(drdu, 'r', 'u');
    drdk     = symbolic2function(drdk, 'r', 'k');
    
    dydx     = symbolic2function(dydx, 'y', 'x');
    dydu     = symbolic2function(dydu, 'y', 'u');
    dydk     = symbolic2function(dydk, 'y', 'k');
    
    dx0ds    = symbolic2function(dx0ds, 'x', 's');
end

if order >= 2
    d2fdx2   = symbolic2function(d2fdx2, 'f', {'x' 'x'});
    d2fdu2   = symbolic2function(d2fdu2, 'f', {'u' 'u'});
    d2fdk2   = symbolic2function(d2fdk2, 'f', {'k' 'k'});
    d2fdudx  = symbolic2function(d2fdudx, 'f', {'x' 'u'});
    d2fdxdu  = symbolic2function(d2fdxdu, 'f', {'u' 'x'});
    d2fdkdx  = symbolic2function(d2fdkdx, 'f', {'x' 'k'});
    d2fdxdk  = symbolic2function(d2fdxdk, 'f', {'k' 'x'});
    d2fdkdu  = symbolic2function(d2fdkdu, 'f', {'u' 'k'});
    d2fdudk  = symbolic2function(d2fdudk, 'f', {'k' 'u'});

    d2rdx2   = symbolic2function(d2rdx2, 'r', {'x' 'x'});
    d2rdu2   = symbolic2function(d2rdu2, 'r', {'u' 'u'});
    d2rdxdu  = symbolic2function(d2rdxdu, 'r', {'u' 'x'});
    d2rdudx  = symbolic2function(d2rdudx, 'r', {'x' 'u'});
    d2rdk2   = symbolic2function(d2rdk2, 'r', {'k' 'k'});
    d2rdkdx  = symbolic2function(d2rdkdx, 'r', {'x' 'k'});
    d2rdkdu  = symbolic2function(d2rdkdu, 'r', {'u' 'k'});
    d2rdxdk  = symbolic2function(d2rdxdk, 'r', {'k' 'x'});
    d2rdudk  = symbolic2function(d2rdudk, 'r', {'k' 'u'});
    
    d2ydx2   = symbolic2function(d2ydx2,  'y', {'x' 'x'});
    d2ydu2   = symbolic2function(d2ydu2,  'y', {'u' 'u'});
    d2ydk2   = symbolic2function(d2ydk2,  'y', {'k' 'k'});
    d2ydudx  = symbolic2function(d2ydudx, 'y', {'x' 'u'});
    d2ydkdx  = symbolic2function(d2ydkdx, 'y', {'x' 'k'});
    d2ydxdu  = symbolic2function(d2ydxdu, 'y', {'u' 'x'});
    d2ydkdu  = symbolic2function(d2ydkdu, 'y', {'u' 'k'});
    d2ydxdk  = symbolic2function(d2ydxdk, 'y', {'k' 'x'});
    d2ydudk  = symbolic2function(d2ydudk, 'y', {'k' 'u'});
    
    d2x0ds2  = symbolic2function(d2x0ds2, 'x', {'s' 's'});
end

if order >= 3
    d3fdx3   = symbolic2function(d3fdx3, 'f', {'x' 'x' 'x'});
    d3fdkdx2 = symbolic2function(d3fdkdx2, 'f', {'x' 'x' 'k'});
    
    d3x0ds3  = symbolic2function(d3x0ds3, 'x', {'s' 's' 's'});
end

if verbose; fprintf('   ...done.\n'); end

%% Set up model structure

% Clear unnecessary variables from scope
clear SymModel vSyms kSyms sSyms qSyms xuSyms xSyms uSyms ...
    vStrs kStrs sStrs qStrs xuStrs xStrs uStrs ...
    xNamesFull uNamesFull vxNames vuNames ...
    C1Entries C1Values C2Entries C2Values cEntries cValues ...
    nC1Entries nC2Entries ncEntries defaultOpts opts ...
    addlength currentLength iExpr ind match nAdd nExpr uAddInd sys_string ...
    i iv ik is iq iu ix ir iy ...
    rstates rinputs rparams statesubs inputsubs paramssubs rstatesi rinputsi rparamsi...
    uqi qsinui nz sizes symbolic2stringmethod...
    thistime defaultMEXdirectory

m.nv = nv;
m.nk = nk;
m.ns = ns;
m.nu = nu;
m.nx = nx;
m.nr = nr;
m.ny = ny;
m.nz = nRules;

m.dv = dv;
m.k  = k;
m.s  = s;
m.u  = u;

m.vxInd = vxInd;
m.vuInd = vuInd;

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

m.S = StoichiometricMatrix;
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

m.y = setfun_y(y,true,k,ny);

if order >= 1
    m.dydx      = setfun_y(dydx,false,k,ny);
    m.dydu      = setfun_y(dydu,false,k,ny);
    m.dydk      = setfun_y(dydk,false,k,ny);
end

if order >= 2
    m.d2ydx2    = setfun_y(d2ydx2,false,k,ny);
    m.d2ydu2    = setfun_y(d2ydu2,false,k,ny);
    m.d2ydk2    = setfun_y(d2ydk2,false,k,ny);
    m.d2ydudx   = setfun_y(d2ydudx,false,k,ny);
    m.d2ydkdx   = setfun_y(d2ydkdx,false,k,ny);
    m.d2ydxdu   = setfun_y(d2ydxdu,false,k,ny);
    m.d2ydkdu   = setfun_y(d2ydkdu,false,k,ny);
    m.d2ydxdk   = setfun_y(d2ydxdk,false,k,ny);
    m.d2ydudk   = setfun_y(d2ydudk,false,k,ny);
end

m.x0            = x0;

if order >= 1
    m.dx0ds     = dx0ds;
end

if order >= 2
    m.d2x0ds2   = d2x0ds2;
end

if order >= 3
    m.d3x0ds3   = d3x0ds3;
end

m.Ready  = true;
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
        
        kAssign = num2cell(k);
        [m.Parameters.Value] = kAssign{:};
        
        m.k  = k;
        
        % Update function handles
        m.f             = setfun_rf(f,k);
        m.r             = setfun_rf(r,k);
        m.y             = setfun_y(y,true,k,ny);
        
        if order >= 1
            m.dfdx      = setfun_rf(dfdx,k);
            m.dfdk      = setfun_rf(dfdk,k);
            m.dfdu      = setfun_rf(dfdu,k);
            m.drdx      = setfun_rf(drdx,k);
            m.drdk      = setfun_rf(drdk,k);
            m.drdu      = setfun_rf(drdu,k);
            m.dydx      = setfun_y(dydx,false,k,ny);
            m.dydu      = setfun_y(dydu,false,k,ny);
            m.dydk      = setfun_y(dydk,false,k,ny);
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
            m.d2ydx2    = setfun_y(d2ydx2,false,k,ny);
            m.d2ydu2    = setfun_y(d2ydu2,false,k,ny);
            m.d2ydk2    = setfun_y(d2ydk2,false,k,ny);
            m.d2ydudx   = setfun_y(d2ydudx,false,k,ny);
            m.d2ydkdx   = setfun_y(d2ydkdx,false,k,ny);
            m.d2ydxdu   = setfun_y(d2ydxdu,false,k,ny);
            m.d2ydkdu   = setfun_y(d2ydkdu,false,k,ny);
            m.d2ydxdk   = setfun_y(d2ydxdk,false,k,ny);
            m.d2ydudk   = setfun_y(d2ydudk,false,k,ny);
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

    function fun = symbolic2function(dsym, num, dens)
        
        % Standardize independent variable names as a cell array
        if ischar(dens)
            dens = {dens};
        end
        
        % Get string of variable name for this function (i.e., 'dfdx')
        if ~isempty(dens)
            order_ = length(dens);
            if order_ >= 2
                orderstr = int2str(length(dens));
            else
                orderstr = '';
            end
            [denterms,~,denoccurenceindex] = unique(dens,'stable');
            dencounts = histcounts(denoccurenceindex,0.5:length(denterms)+0.5);
            dencounts = strtrim(cellstr(int2str(dencounts(:))));
            dencounts(strcmp(dencounts,'1')) = {''};
            denstrs = strcat('d',denterms(:),dencounts(:));
            variable_name = ['d' orderstr num denstrs{:}];
        else
            variable_name = num;
        end
        
        if verbose; fprintf([variable_name '...']); end
        
        switch symbolic2stringmethod
            case 'efficient'
                
                string_rep = symbolic2string(dsym, num, dens);
                fun = string2fun(string_rep, num, dens);
                
            case 'mex'

                % Don't create MEX functions for x0 functions, because they
                % aren't called very many times
                if strcmp(num,'x')
                    string_rep = symbolic2string(dsym, num, dens);
                    fun = string2fun(string_rep, num, dens);
                    return
                end
                
                % Prepare inputs to C code generation function
                nzlogical = getNonZeroEntries(num, dens);
                nzi = cell(ndims(nzlogical),1);
                nziind = find(nzlogical);
                [nzi{:}] = ind2sub(size(nzlogical),nziind(:));
                nzi = [nzi{:}];
                nzsizes = [sizes.(num) cellfun(@(den)sizes.(den),dens(:)')];
                
                % Generate mex C code
                getMexReadyCode(dsym,nzi,nzsizes,xSyms,uSyms,kSyms,variable_name,opts.MEXDirectory)
             
                fun = str2func([variable_name 'fun']);
        end
        
        if verbose; fprintf('Done.\n'); end
        
    end

    function string_rep = symbolic2string(dsym, num, dens)
        
        % If provided an empty input, return an empty output
        if isempty(dsym)
            dsymsize = size(dsym);
            dsymsize = strtrim(cellstr(int2str(dsymsize(:))));
            if isempty(dens)
                string_rep = ['zeros(' dsymsize{:} ')'];
                return
            else
                string_rep = ['[],[],[],' dsymsize{1} ',' dsymsize{2}];
                return
            end
        end
        
        if isempty(dens)
            
            % For zero-order derivatives, don't create a sparse matrix
            strelements = fastchar(dsym);
            string_rep = strjoin(row(strelements),';');
            
        else
    
            % Get nonzero derivative logical matrix
            nzlogical = getNonZeroEntries(num, dens);
            
            % Reshape nzlogical to the same size as that of the provided value
            varsizes = size(dsym);
            nzlogical = reshape(nzlogical,varsizes);
            
            % Convert nonzero terms to strings
            nzindices = find(nzlogical(:));
            strelements = fastchar(dsym(nzindices)); %#ok % Don't replace indices with logicals here. The conversion of the logical to sym takes too long.
            
            strelements = strjoin(row(strelements),',');
            
            % Convert subscripts into strings
            [isubscripts,jsubscripts] = find(nzlogical);
            isubscriptstrs = strtrim(cellstr(num2str(isubscripts)));
            jsubscriptstrs = strtrim(cellstr(num2str(jsubscripts)));
            isubstring = strjoin(row(isubscriptstrs),',');
            jsubstring = strjoin(row(jsubscriptstrs),',');
            
            % Write sparse initialization string, which will fit in the
            % following expression:
            % '@(t,x,u,k) sparse(' STRING_REP ')'
            string_rep = sprintf(['[' isubstring '],[' jsubstring '],[' strelements '],[' num2str(varsizes(1)) '],[' num2str(varsizes(2)) ']']);
            
        end
        
        % Replace symbolic names with indexed references to vectors
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
        if ns == 0
            sindexstring = {};
        else
            sindexstring = sprintf('s(%d)\n', 1:ns);
            sindexstring = textscan(sindexstring,'%s','Delimiter','\n');
            sindexstring = sindexstring{1};
        end
        string_rep = regexprep(string_rep, [xStrs; uStrs; kStrs; sStrs], [xindexstring; uindexstring; kindexstring; sindexstring], 0);
        
    end         

    function d2ydx2dx1_ = calcDerivative(dydx1_, x2Syms_, num, dens)
        % y = dependent variable
        % x1 = 1st derivative variable
        % x2 = 2nd derivative variable
        % num and dens should be for the output derivative matrix.
        
        ny_  = size(dydx1_,1);
        nx1_ = size(dydx1_,2);
        nx2_ = length(x2Syms_);
        
        % If any dimensions are zero, return an empty derivative
        if any([nx1_ nx2_ ny_] == 0, 2)
            d2ydx2dx1_ = initializeMatrixMupad([], [], [], ny_*nx1_, nx2_);
            return
        end
        
        % Standardize dens as a cell array
        if ischar(dens)
            dens = {dens};
        end
        
        % Find the entries of d2ydx2dx1_ that might be nonzero
        nze = getNonZeroEntries(num, dens);
        nze = reshape(nze, ny_*nx1_, nx2_);
        [nzterms, nzdens] = find(nze);
        
        % If there aren't any nonzero terms, return an all-zero derivative
        if isempty(nzterms)
            d2ydx2dx1_ = initializeMatrixMupad([], [], [], ny_*nx1_, nx2_);
            return
        end
        
        % Take derivatives of the possibly nonzero entries
        nzders = diff_vectorized(vec(dydx1_(nzterms)), vec(x2Syms_(nzdens)));
        
        % Of the supposedly nonzero derivatives, find the ones that are
        % actually nonzero, and only keep those
        iszero = logical(nzders == 0);
        nzeiszero = sub2ind([ny_*nx1_, nx2_], nzterms(iszero), nzdens(iszero));
        nze(nzeiszero) = false;
        
        % Get sizes of dependent (numerator) and independent (denominator)
        % variables in the derivative
        numsize = sizes.(num);
        densizes = zeros(1, length(dens));
        for di = 1:length(dens)
            densizes(di) = sizes.(dens{di});
        end
        
        % Reshape nonzero entries to n-dimensional matrix
        nze = reshape(nze, [numsize,densizes]);
        
        % Update the non-zero map to account for newly discovered zero
        % terms
        nzkey = [num dens{:}];
        nz(nzkey) = nze;
        
        % Also get information about f derivative, if numerator is r
        if strcmp(num,'r')
            nzkey_f = strrep(nzkey,'r','f');
            % Only update the key for f if no information about this
            % particular derivative was previously stored in nz
            if ~isKey(nz, nzkey_f)
                nztemp = full(logical(abs(StoichiometricMatrix)*reshape(nze, [ny_, nx1_*nx2_])));
                nz(nzkey_f) = reshape(nztemp, [nx, nx1_, nx2_]);
            end
        end
        
        % Remove zero terms
        nzders(iszero) = [];
        nzterms(iszero) = [];
        nzdens(iszero) = [];
        
        d2ydx2dx1_ = initializeMatrixMupad(nzterms, nzdens, nzders, ny_*nx1_, nx2_);
        
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

end

function fun = string2fun(string_rep, num, dens)
% Note that string2fun is a subfunction instead of a nested function to
% prevent the anonymous functions created here from saving copies of
% the primary function workspace variables.

% If the dependent variable is x0, the input argument is s. Otherwise,
% t,x,u,k.
if strcmp(num,'x')
    inputargstr = 's';
else
    inputargstr = 't,x,u,k';
end

% If this isn't a derivative (no independent variables specified), make a
% full array. Otherwise make a sparse array.
if isempty(dens)
    sparsestr = {'[' ']'};
else
    sparsestr = {'sparse(' ')'};
end

% Set up the function handle
fun = eval(['@(' inputargstr ') ' sparsestr{1} string_rep sparsestr{2}]);

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

function fun = setfun_y(basefun, is0order, k, ny)
    if is0order
        fun = @(t,x,u) vectorize_y(basefun, t, x, u, k, ny);
    else
        fun = @(t,x,u) basefun(t, x, u, k);
    end
end

function val = vectorize_y(y, t, x, u, k, ny)
    nt = numel(t);
    val = zeros(ny,nt);
    if isempty(x)
        x = zeros(0,nt);
    end
    if isempty(u)
        u = zeros(0,nt);
    end
    if isempty(k)
        k = 0; % doesn't change in time
    end
    for it = 1:nt
        val(:,it) = y(t(it), x(:,it), u(:,it), k);
    end
end

function rOut = evaluateExternalFunctions(rIn, ids)
% Evaluate symbolic functions/pull in functions defined in path
%   Necessary for "power" and other MathML function translation
% Initialize symbolic variables
syms(ids{:});

% Evaluate the expressions to remove function calls
rOut = eval(rIn);

% Clear the symbolic variables
clear(ids{:})
end