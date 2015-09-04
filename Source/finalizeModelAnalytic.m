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
default_opts.Order                     = 2;
default_opts.Verbose                   = 0;
default_opts.UseMEX                    = false;
default_opts.MEXDirectory              = defaultMEXdirectory;
default_opts.EvaluateExternalFunctions = true; % needed for calls to power() and other functions

opts = mergestruct(default_opts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

order = opts.Order;

if opts.UseMEX && exist(opts.MEXDirectory,'dir') ~= 7
    mkdir(opts.MEXDirectory);
end

%% Trim m.add.* components to only those actually added
m.add.Compartments = m.add.Compartments(1:m.add.nv);
m.add.Seeds = m.add.Seeds(1:m.add.ns);
m.add.Parameters = m.add.Parameters(1:m.add.nk);
m.add.Inputs = m.add.Inputs(1:m.add.nu);
m.add.States = m.add.States(1:m.add.nx);
m.add.Reactions = m.add.Reactions(1:m.add.nr);
m.add.Rules = m.add.Rules(1:m.add.nz);
m.add.Outputs = m.add.Outputs(1:m.add.ny);

%% Check for uniqueness of names and IDs
% Blank IDs are ignored as they will be automatically generated later
all_existing_ids = vec({m.Compartments.ID, m.Seeds.ID, m.Parameters.ID, m.Inputs.ID, m.States.ID, m.Reactions.ID, m.Rules.ID, m.Outputs.ID});
all_existing_ids = check_unique_name_id('compartment', {m.add.Compartments().Name}, {m.add.Compartments.ID}, {m.Compartments.Name}, all_existing_ids);
all_existing_ids = check_unique_name_id('seed', {m.add.Seeds.Name}, {m.add.Seeds.ID}, {m.Seeds.Name}, all_existing_ids);
all_existing_ids = check_unique_name_id('parameter', {m.add.Parameters.Name}, {m.add.Parameters.ID}, {m.Parameters.Name}, all_existing_ids);
all_existing_ids = check_unique_name_id('input', {m.add.Inputs.Name}, {m.add.Inputs.ID}, {m.Inputs.Name, m.States.Name}, all_existing_ids);
all_existing_ids = check_unique_name_id('state', {m.add.States.Name}, {m.add.States.ID}, {m.Inputs.Name, m.States.Name}, all_existing_ids);
all_existing_ids = check_unique_name_id('reaction', {}, {m.add.Reactions.ID}, {}, all_existing_ids); % Reactions don't need names
all_existing_ids = check_unique_name_id('rule', {m.add.Rules.Name}, {m.add.Rules.ID}, {m.Rules.Name}, all_existing_ids);
all_existing_ids = check_unique_name_id('output', {m.add.Outputs.Name}, {m.add.Outputs.ID}, {m.Outputs.Name}, all_existing_ids);

%% Add in m.add.* components
m.Compartments = [m.Compartments; m.add.Compartments];
m.Seeds = [m.Seeds; m.add.Seeds];
m.Parameters = [m.Parameters; m.add.Parameters];
m.Inputs = [m.Inputs; m.add.Inputs];
m.States = [m.States; m.add.States];
m.Reactions = [m.Reactions; m.add.Reactions];
m.Rules = [m.Rules; m.add.Rules];
m.Outputs = [m.Outputs; m.add.Outputs];

m.add.Compartments = growCompartmentsAnalytic;
m.add.Parameters = growParametersAnalytic;
m.add.Seeds = growSeedsAnalytic;
m.add.Inputs = growInputsAnalytic;
m.add.States = growStatesAnalytic;
m.add.Reactions = growReactionsAnalytic;
m.add.Rules = growRulesAnalytic;
m.add.Outputs = growOutputsAnalytic;

m.add.nv = 0;
m.add.nk = 0;
m.add.ns = 0;
m.add.nu = 0;
m.add.nx = 0;
m.add.nz = 0;
m.add.nr = 0;
m.add.ny = 0;

%% Extract values
nv = numel(m.Compartments);
nk = numel(m.Parameters);
ns = numel(m.Seeds);
nu = numel(m.Inputs);
nx = numel(m.States);
nr = numel(m.Reactions);
nz = numel(m.Rules);
ny = numel(m.Outputs);

% String representations of symbolics. Some may be empty and will be assigned later
v_strs = vec({m.Compartments.ID});
k_strs = vec({m.Parameters.ID});
s_strs = vec({m.Seeds.ID});
u_strs = vec({m.Inputs.ID});
x_strs = vec({m.States.ID});
z_strs = vec({m.Rules.ID});
y_strs = vec({m.Outputs.ID});

v_names = vec({m.Compartments.Name});
k_names = vec({m.Parameters.Name});
s_names = vec({m.Seeds.Name});
u_names = vec({m.Inputs.Name});
x_names = vec({m.States.Name});
z_names = vec({m.Rules.Name});
y_names = vec({m.Outputs.Name});

% Some of these will be converted to symbolics later after all default symbolics have been added
dv = vec([m.Compartments.Dimension]);
v = vec({m.Compartments.Size});

s = vec([m.Seeds.Value]);

k = vec([m.Parameters.Value]);

vu_names = vec({m.Inputs.Compartment});
[~, vuInd] = ismember(vu_names, v_names);
u = vec([m.Inputs.DefaultValue]);

vx_names = vec({m.States.Compartment});
[~, vxInd] = ismember(vx_names, v_names);
x0 = vec({m.States.InitialValue});

r = vec({m.Reactions.Rate});

z = {m.Rules.Expression}';

y = vec({m.Outputs.Expression});

%% Assign missing IDs
[v_strs, all_existing_ids] = update_unique_ids(v_names, v_strs, all_existing_ids);
if ~isempty(v_strs); [m.Compartments.ID] = deal(v_strs{:}); end % Matlab bug will not deal to empty stucts
[k_strs, all_existing_ids] = update_unique_ids(k_names, k_strs, all_existing_ids);
if ~isempty(k_strs); [m.Parameters.ID] = deal(k_strs{:}); end
[s_strs, all_existing_ids] = update_unique_ids(s_names, s_strs, all_existing_ids);
if ~isempty(s_strs); [m.Seeds.ID] = deal(s_strs{:}); end
[u_strs, all_existing_ids] = update_unique_ids(u_names, u_strs, all_existing_ids);
if ~isempty(u_strs); [m.Inputs.ID] = deal(u_strs{:}); end
[x_strs, all_existing_ids] = update_unique_ids(x_names, x_strs, all_existing_ids);
if ~isempty(x_strs); [m.States.ID] = deal(x_strs{:}); end
[y_strs, all_existing_ids] = update_unique_ids(y_names, y_strs, all_existing_ids);
if ~isempty(y_strs); [m.Outputs.ID] = deal(y_strs{:}); end

%% Assemble symbolic representations of model components
v_syms = sym(v_strs);
k_syms = sym(k_strs);
s_syms = sym(s_strs);
u_syms = sym(u_strs);
x_syms = sym(x_strs);
z_syms = sym(z_strs);
y_syms = sym(y_strs);

%% Build map of names to IDs
all_names = [v_names; s_names; k_names; u_names; x_names; z_names; y_names];
all_ids = [v_strs; s_strs; k_strs; u_strs; x_strs; z_strs; y_strs];

% Remove all names containing a double quote
bad_names = cellfun(@(str)ismember('"',str), all_names);
all_names = all_names(~bad_names);
all_ids = all_ids(~bad_names);

%% Replace names with symbolics and convert to symbolics
v = replace_names(v, all_names, all_ids);
x0 = replace_names(x0, all_names, all_ids);
z = replace_names(z, all_names, all_ids);
r = replace_names(r, all_names, all_ids);
y = replace_names(y, all_names, all_ids);

v = sym(v);
x0 = sym(x0);
z = sym(z);
r = sym(r);
y = sym(y);

%% Standardize all symbolic names
% Changing all the symbols to a standard format ensures that string replacement works later
% TODO: check the unlikely case that someone is calling a function sy#x
old_symbols = [s_syms; k_syms; u_syms; x_syms];
new_symbols = sym(zeros(nk+ns+nx+nu,1));
for i = 1:ns+nk+nu+nx
    new_symbols(i) = sym(sprintf('sy%dx', i));
end

v_syms = fastsubs(v_syms, old_symbols, new_symbols);
k_syms = fastsubs(k_syms, old_symbols, new_symbols);
s_syms = fastsubs(s_syms, old_symbols, new_symbols);
u_syms = fastsubs(u_syms, old_symbols, new_symbols);
x_syms = fastsubs(x_syms, old_symbols, new_symbols);
z_syms = fastsubs(z_syms, old_symbols, new_symbols);
y_syms = fastsubs(y_syms, old_symbols, new_symbols);

v = fastsubs(v, old_symbols, new_symbols);
x0 = fastsubs(x0, old_symbols, new_symbols);
z = fastsubs(z, old_symbols, new_symbols);
r = fastsubs(r, old_symbols, new_symbols);
y = fastsubs(y, old_symbols, new_symbols);

v_strs = fastchar(v_syms);
k_strs = fastchar(k_syms);
s_strs = fastchar(s_syms);
u_strs = fastchar(u_syms);
x_strs = fastchar(x_syms);
z_strs = fastchar(z_syms);
y_strs = fastchar(y_syms);

%% Substitute in expressions
% Everything that is substitutable
substitutable_ids = [v_syms; z_syms];
substitutable_exps = [v; z];

% Need to substitute n times to ensure that all substitutions are complete
n_subs = numel(substitutable_ids);
for i = 1:n_subs
    substitutable_exps = subs(substitutable_exps, substitutable_ids, substitutable_exps);
end

% Extract or replace
v = substitutable_exps(1:nv,1);
x0 = fastsubs(x0, substitutable_ids, substitutable_exps);
z = substitutable_exps(nv+(1:nz),1);
r = fastsubs(r, substitutable_ids, substitutable_exps);
y = fastsubs(y, substitutable_ids, substitutable_exps);

%% Evaluate external functions
if opts.EvaluateExternalFunctions
    x0 = evaluate_external_functions(x0, [s_strs; x_strs; u_strs; k_strs]);
    r = evaluate_external_functions(r, [s_strs; x_strs; u_strs; k_strs]);
    y = evaluate_external_functions(y, [s_strs; x_strs; u_strs; k_strs]);
end

%% Process stoichiometry and rate forms/RHS's
% Make stoichiometry matrix nSpecies x nReactions
xu_names = [x_names; u_names];
S_entries = zeros(0,3);
for i = 1:nr
    for j = 1:numel(m.Reactions(i).Reactants)
        ind = lookup(m.Reactions(i).Reactants{j}, xu_names);
        S_entries = [S_entries; ind, i, -1];
    end
    for j = 1:numel(m.Reactions(i).Products)
        ind = lookup(m.Reactions(i).Products{j}, xu_names);
        S_entries = [S_entries; ind, i, 1];
    end
end
Sxu = sparse(S_entries(:,1), S_entries(:,2), S_entries(:,3), nx+nu, nr);
Su = Sxu(nx+1:end,:);
S = Sxu(1:nx,:);

% Convert stoichiometrix matrix to symbolic so that it doesn't have to be
% converted later for each multiplication
[i_S, j_S, val_S] = find(S);
size_S = size(S);
S_sym = initializeMatrixMupad(i_S, j_S, val_S, size_S(1), size_S(2));

%% Construct ODE system
f = S_sym*r;

%% Sanity checks
% TODO: make sure initial conditions are only functions of seeds

% Check original symbols for reserved MuPAD names
old_symbols_strs = fastchar([k_syms; s_syms; x_syms; u_syms]);
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

%% Determine the variables in each expression
x0str = fastchar(x0);
rstr = fastchar(r);
ystr = fastchar(y);

x0hass = expression_has_variable(x0str, s_strs);
x0hask = expression_has_variable(x0str, k_strs);

rhasx = expression_has_variable(rstr, x_strs);
rhasu = expression_has_variable(rstr, u_strs);
rhask = expression_has_variable(rstr, k_strs);

yhasx = expression_has_variable(ystr, x_strs);
yhasu = expression_has_variable(ystr, u_strs);
yhask = expression_has_variable(ystr, k_strs);

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
nonzero_map = containers.Map;
nonzero_map('r') = true(nr,1);
nonzero_map('f') = true(nx,1);
nonzero_map('y') = true(ny,1);
nonzero_map('rx') = rhasx;
nonzero_map('ru') = rhasu;
nonzero_map('rk') = rhask;
nonzero_map('fx') = logical(full(abs(S)*rhasx)); % Take the absolute value of S so that there are no accidental cancellations between positive and negative terms
nonzero_map('fu') = logical(full(abs(S)*rhasu));
nonzero_map('fk') = logical(full(abs(S)*rhask));
nonzero_map('yx') = yhasx;
nonzero_map('yu') = yhasu;
nonzero_map('yk') = yhask;
nonzero_map('xs') = x0hass;
nonzero_map('xk') = x0hask;

if verbose; fprintf('\n'); end

%% Generate derivatives of desired order
if order >= 1

    % Gradient of r with respect to x
    if verbose; fprintf('Calculating drdx...'); end
    drdx = calculate_derivative(r, x_syms, 'r', {'x'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of r with respect to u
    if verbose; fprintf('Calculating drdu...'); end
    drdu = calculate_derivative(r, u_syms, 'r', {'u'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of r with respect to k
    if verbose; fprintf('Calculating drdk...'); end
    drdk = calculate_derivative(r, k_syms, 'r', {'k'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of f with respect to x
    if verbose; fprintf('Calculating dfdx...'); end
    dfdx = S_sym*drdx;
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of f with respect to u
    if verbose; fprintf('Calculating dfdu...'); end
    dfdu = S_sym*drdu;
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of f with respect to k
    if verbose; fprintf('Calculating dfdk...'); end
    dfdk = S_sym*drdk;
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of y with respect to x
    if verbose; fprintf('Calculating dydx...'); end
    dydx = calculate_derivative(y, x_syms, 'y', {'x'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of y with respect to u
    if verbose; fprintf('Calculating dydu...'); end
    dydu = calculate_derivative(y, u_syms, 'y', {'u'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of y with respect to k
    if verbose; fprintf('Calculating dydk...'); end
    dydk = calculate_derivative(y, k_syms, 'y', {'k'});
    if verbose; fprintf('Done.\n'); end

    % Gradient of x0 with respect to s
    if verbose; fprintf('Calculating dx0ds...'); end
    dx0ds = calculate_derivative(x0, s_syms, 'x', {'s'});
    if verbose; fprintf('Done.\n'); end

    % Gradient of x0 with respect to k
    if verbose; fprintf('Calculating dx0k...'); end
    dx0dk = calculate_derivative(x0, k_syms, 'x', {'k'});
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
    dx0k = '';
end

if order >= 2
    % Gradient of drdx with respect to x
    if verbose; fprintf('Calculating d2rdx2...'); end
    d2rdx2 = calculate_derivative(drdx, x_syms, 'r', {'x','x'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of drdu with respect to u
    if verbose; fprintf('Calculating d2rdu2...'); end
    d2rdu2 = calculate_derivative(drdu, u_syms, 'r', {'u','u'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of drdu with respect to x
    if verbose; fprintf('Calculating d2rdxdu...'); end
    d2rdxdu = calculate_derivative(drdu, x_syms, 'r', {'u','x'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of drdx with respect to u
    if verbose; fprintf('Calculating d2rdudx...'); end
    d2rdudx = calculate_derivative(drdx, u_syms, 'r', {'x','u'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of drdk with respect to k
    if verbose; fprintf('Calculating d2rdk2...'); end
    d2rdk2 = calculate_derivative(drdk, k_syms, 'r', {'k','k'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of drdx with respect to k
    if verbose; fprintf('Calculating d2rdkdx...'); end
    d2rdkdx = calculate_derivative(drdx, k_syms, 'r', {'x','k'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of drdu with respect to k
    if verbose; fprintf('Calculating d2rdkdu...'); end
    d2rdkdu = calculate_derivative(drdu, k_syms, 'r', {'u','k'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of drdk with respect to x
    if verbose; fprintf('Calculating d2rdxdk...'); end
    d2rdxdk = calculate_derivative(drdk, x_syms, 'r', {'k','x'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of drdk with respect to u
    if verbose; fprintf('Calculating d2rdudk...'); end
    d2rdudk = calculate_derivative(drdk, u_syms, 'r', {'k','u'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of dfdx with respect to x
    if verbose; fprintf('Calculating d2fdx2...'); end
    d2fdx2 = S_sym*reshape_derivative(d2rdx2, [nr nx*nx], 'r', {'x' 'x'});
    d2fdx2 = reshape_derivative(d2fdx2, [nx*nx nx], 'f', {'x' 'x'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of dfdu with respect to u
    if verbose; fprintf('Calculating d2fdu2...'); end
    d2fdu2 = S_sym*reshape_derivative(d2rdu2, [nr,nu*nu], 'r', {'u' 'u'});
    d2fdu2 = reshape_derivative(d2fdu2, [nx*nu,nu], 'f', {'u' 'u'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of dfdu with respect to x
    if verbose; fprintf('Calculating d2fdxdu...'); end
    d2fdxdu = S_sym*reshape_derivative(d2rdxdu, [nr,nu*nx], 'r', {'u' 'x'});
    d2fdxdu = reshape_derivative(d2fdxdu, [nx*nu,nx], 'f', {'u' 'x'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of dfdx with respect to u
    if verbose; fprintf('Calculating d2fdudx...'); end
    d2fdudx = S_sym*reshape_derivative(d2rdudx, [nr,nx*nu], 'r', {'x' 'u'});
    d2fdudx = reshape_derivative(d2fdudx, [nx*nx,nu], 'f', {'x' 'u'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of dfdk with respect to k
    if verbose; fprintf('Calculating d2fdk2...'); end
    d2fdk2 = S_sym*reshape_derivative(d2rdk2, [nr,nk*nk], 'r', {'k' 'k'});
    d2fdk2 = reshape_derivative(d2fdk2, [nx*nk,nk], 'f', {'k' 'k'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of dfdx with respect to k
    if verbose; fprintf('Calculating d2fdkdx...'); end
    d2fdkdx = S_sym*reshape_derivative(d2rdkdx, [nr,nx*nk], 'r', {'x' 'k'});
    d2fdkdx = reshape_derivative(d2fdkdx, [nx*nx,nk], 'f', {'x' 'k'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of dfdu with respect to k
    if verbose; fprintf('Calculating d2fdkdu...'); end
    d2fdkdu = S_sym*reshape_derivative(d2rdkdu, [nr,nu*nk], 'r',{'u' 'k'});
    d2fdkdu = reshape_derivative(d2fdkdu, [nx*nu,nk], 'f',{'u' 'k'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of dfdk with respect to x
    if verbose; fprintf('Calculating d2fdxdk...'); end
    d2fdxdk = S_sym*reshape_derivative(d2rdxdk, [nr,nk*nx], 'r',{'k' 'x'});
    d2fdxdk = reshape_derivative(d2fdxdk, [nx*nk,nx], 'f',{'k' 'x'});
    if verbose; fprintf('Done.\n'); end
    
    % Gradient of dfdk with respect to u
    if verbose; fprintf('Calculating d2fdudk...'); end
    d2fdudk = S_sym*reshape_derivative(d2rdudk, [nr,nk*nu], 'r',{'k' 'u'});
    d2fdudk = reshape_derivative(d2fdudk, [nx*nk,nu], 'f',{'k' 'u'});
    if verbose; fprintf('Done.\n'); end
    
    % Output's second derivatives
    if verbose; fprintf('Calculating d2ydx2...'); end
    d2ydx2 = calculate_derivative(dydx, x_syms, 'y', {'x','x'});
    if verbose; fprintf('Done.\n'); end
    
    if verbose; fprintf('Calculating d2ydu2...'); end
    d2ydu2 = calculate_derivative(dydu, u_syms, 'y', {'u','u'});
    if verbose; fprintf('Done.\n'); end
    
    if verbose; fprintf('Calculating d2ydk2...'); end
    d2ydk2 = calculate_derivative(dydk, k_syms, 'y', {'k','k'});
    if verbose; fprintf('Done.\n'); end
    
    if verbose; fprintf('Calculating d2ydudx...'); end
    d2ydudx = calculate_derivative(dydx, u_syms, 'y', {'x','u'});
    if verbose; fprintf('Done.\n'); end
    
    if verbose; fprintf('Calculating d2ydkdx...'); end
    d2ydkdx = calculate_derivative(dydx, k_syms, 'y', {'x','k'});
    if verbose; fprintf('Done.\n'); end
    
    if verbose; fprintf('Calculating d2ydxdu...'); end
    d2ydxdu = calculate_derivative(dydu, x_syms, 'y', {'u','x'});
    if verbose; fprintf('Done.\n'); end
    
    if verbose; fprintf('Calculating d2ydkdu...'); end
    d2ydkdu = calculate_derivative(dydu, k_syms, 'y', {'u','k'});
    if verbose; fprintf('Done.\n'); end
    
    if verbose; fprintf('Calculating d2ydxdk...'); end
    d2ydxdk = calculate_derivative(dydk, x_syms, 'y', {'k','x'});
    if verbose; fprintf('Done.\n'); end
    
    if verbose; fprintf('Calculating d2ydudk...'); end
    d2ydudk = calculate_derivative(dydk, u_syms, 'y', {'k','u'});
    if verbose; fprintf('Done.\n'); end

    if verbose; fprintf('Calculating d2x0ds2...'); end
    d2x0ds2 = calculate_derivative(dx0ds, s_syms, 'x', {'s','s'});
    if verbose; fprintf('Done.\n'); end
    
    if verbose; fprintf('Calculating d2x0dk2...'); end
    d2x0dk2 = calculate_derivative(dx0dk, k_syms, 'x', {'k','k'});
    if verbose; fprintf('Done.\n'); end
    
    if verbose; fprintf('Calculating d2x0dkds...'); end
    d2x0dkds = calculate_derivative(dx0ds, k_syms, 'x', {'s','k'});
    if verbose; fprintf('Done.\n'); end
    
    if verbose; fprintf('Calculating d2x0dsdk...'); end
    d2x0dsdk = calculate_derivative(dx0dk, s_syms, 'x', {'k','s'});
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
    d2x0dk2 = '';
    d2x0dkds = '';
    d2x0dsdk = '';
end

if order >= 3
    %Gradient of d2rdx2 with respect to x
    d3rdx3 = calculate_derivative(d2rdx2, x_syms, 'r', {'x','x','x'});
    
    %Gradient of d2rdx2 with respect to k
    d3rdkdx2 = calculate_derivative(d2rdx2, k_syms, 'r', {'x','x','k'});
    
    %Gradient of d2fdx2 with respect to x
    d3fdx3 = S_sym*reshape_derivative(d3rdx3, [nr,nx*nx*nx], 'r',{'x' 'x' 'x'});
    d3fdx3 = reshape_derivative(d3fdx3, [xn*nx*nx,nx], 'f',{'x' 'x' 'x'});
    
    %Gradient of d2fdx2 with respect to k
    d3fdkdx2 = S_sym*reshape_derivative(d3rdkdx2, [nr,nx*nx*nk], 'r',{'x' 'x' 'k'});
    d3fdkdx2 = reshape_derivative(d3fdkdx2, [nx*nx*nx,nk], 'f',{'x' 'x' 'k'});
    
    % Gradient of d2x0ds2 with respect to s
    d3x0ds3 = calculate_derivative(d2x0ds2, s_syms, 'x', {'s','s','s'});
else
    d3rdx3   = '';
    d3rdkdx2 = '';
    d3fdx3   = '';
    d3fdkdx2 = '';
    d3x0ds3  = '';
end

%% Convert symbolic expressions to function handles

if opts.UseMEX
    symbolic2stringmethod = 'mex';
else
    symbolic2stringmethod = 'efficient';
end

if verbose; fprintf('Converting symbolics to functions...\n'); end
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
    dx0dk    = symbolic2function(dx0dk, 'x', 'k');
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
    d2x0dk2  = symbolic2function(d2x0dk2, 'x', {'k' 'k'});
    d2x0dkds = symbolic2function(d2x0dkds, 'x', {'s' 'k'});
    d2x0dsdk = symbolic2function(d2x0dsdk, 'x', {'k' 's'});
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
    uqi qsinui nonzero_map sizes symbolic2stringmethod...
    thistime defaultMEXdirectory

m.nv = nv;
m.nk = nk;
m.ns = ns;
m.nu = nu;
m.nx = nx;
m.nr = nr;
m.ny = ny;
m.nz = nz;

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

m.x0            = setfun_x0(x0,k);

if order >= 1
    m.dx0ds     = setfun_x0(dx0ds,k);
    m.dx0dk     = setfun_x0(dx0dk,k);
end

if order >= 2
    m.d2x0ds2   = setfun_x0(d2x0ds2,k);
    m.d2x0dk2   = setfun_x0(d2x0dk2,k);
    m.d2x0dkds  = setfun_x0(d2x0dkds,k);
    m.d2x0dsdk  = setfun_x0(d2x0dsdk,k);
end

if order >= 3
    m.d3x0ds3   = setfun_x0(d3x0ds3);
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
        
        k_assign = num2cell(k);
        [m.Parameters.Value] = k_assign{:};
        
        m.k  = k;
        
        % Update function handles
        m.x0            = setfun_x0(x0,k);
        m.f             = setfun_rf(f,k);
        m.r             = setfun_rf(r,k);
        m.y             = setfun_y(y,true,k,ny);
        
        if order >= 1
            m.dx0ds     = setfun_x0(dx0ds,k);
            m.dx0dk     = setfun_x0(dx0dk,k);
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
            m.d2x0ds2   = setfun_x0(d2x0ds2,k);
            m.d2x0dk2   = setfun_x0(d2x0dk2,k);
            m.d2x0dkds  = setfun_x0(d2x0dkds,k);
            m.d2x0dsdk  = setfun_x0(d2x0dsdk,k);
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
            order_ = numel(dens);
            if order_ >= 2
                orderstr = int2str(numel(dens));
            else
                orderstr = '';
            end
            [denterms, ~, denoccurenceindex] = unique(dens, 'stable');
            dencounts = histcounts(denoccurenceindex, 0.5:numel(denterms)+0.5);
            dencounts = strtrim(cellstr(int2str(dencounts(:))));
            dencounts(strcmp(dencounts, '1')) = {''};
            denstrs = strcat('d', denterms(:), dencounts(:));
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
                [nzi{:}] = ind2sub(size(nzlogical), nziind(:));
                nzi = [nzi{:}];
                nzsizes = [sizes.(num) cellfun(@(den)sizes.(den), row(dens))];
                
                % Generate mex C code
                getMexReadyCode(dsym, nzi, nzsizes, x_syms, u_syms, k_syms, variable_name, opts.MEXDirectory)
             
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
                string_rep = ['zeros(' strjoin(row(dsymsize), ',') ')'];
                return
            else
                string_rep = ['[],[],[],' dsymsize{1} ',' dsymsize{2}];
                return
            end
        end
        
        if isempty(dens)
            
            % For zero-order derivatives, don't create a sparse matrix
            strelements = fastchar(dsym);
            string_rep = strjoin(row(strelements), ';');
            
        else
    
            % Get nonzero derivative logical matrix
            nzlogical = getNonZeroEntries(num, dens);
            
            % Reshape nzlogical to the same size as that of the provided value
            varsizes = size(dsym);
            nzlogical = reshape(nzlogical, varsizes);
            
            % Convert nonzero terms to strings
            nzindices = find(nzlogical(:));
            strelements = fastchar(dsym(nzindices)); %#ok % Don't replace indices with logicals here. The conversion of the logical to sym takes too long.
            
            strelements = strjoin(row(strelements), ',');
            
            % Convert subscripts into strings
            [isubscripts,jsubscripts] = find(nzlogical);
            isubscriptstrs = strtrim(cellstr(num2str(isubscripts)));
            jsubscriptstrs = strtrim(cellstr(num2str(jsubscripts)));
            isubstring = strjoin(row(isubscriptstrs), ',');
            jsubstring = strjoin(row(jsubscriptstrs), ',');
            
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
            xindexstring = textscan(xindexstring, '%s', 'Delimiter', '\n');
            xindexstring = xindexstring{1};
        end
        if nu == 0
            uindexstring = {};
        else
            uindexstring = sprintf('u(%d)\n', 1:nu);
            uindexstring = textscan(uindexstring, '%s', 'Delimiter', '\n');
            uindexstring = uindexstring{1};
        end
        if nk == 0
            kindexstring = {};
        else
            kindexstring = sprintf('k(%d)\n', 1:nk);
            kindexstring = textscan(kindexstring, '%s', 'Delimiter', '\n');
            kindexstring = kindexstring{1};
        end
        if ns == 0
            sindexstring = {};
        else
            sindexstring = sprintf('s(%d)\n', 1:ns);
            sindexstring = textscan(sindexstring, '%s', 'Delimiter', '\n');
            sindexstring = sindexstring{1};
        end
        string_rep = regexprep(string_rep, [x_strs; u_strs; k_strs; s_strs], [xindexstring; uindexstring; kindexstring; sindexstring], 0);
        
    end         

    function d2ydx2dx1_ = calculate_derivative(dydx1_, x2Syms_, num, dens)
        % y = dependent variable
        % x1 = 1st derivative variable
        % x2 = 2nd derivative variable
        % num and dens should be for the output derivative matrix.
        
        ny_  = size(dydx1_,1);
        nx1_ = size(dydx1_,2);
        nx2_ = numel(x2Syms_);
        
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
        densizes = zeros(1, numel(dens));
        for di = 1:numel(dens)
            densizes(di) = sizes.(dens{di});
        end
        
        % Reshape nonzero entries to n-dimensional matrix
        nze = reshape(nze, [numsize,densizes]);
        
        % Update the non-zero map to account for newly discovered zero
        % terms
        nzkey = [num dens{:}];
        nonzero_map(nzkey) = nze;
        
        % Also get information about f derivative, if numerator is r
        if strcmp(num,'r')
            nzkey_f = strrep(nzkey,'r','f');
            % Only update the key for f if no information about this
            % particular derivative was previously stored in nz
            if ~isKey(nonzero_map, nzkey_f)
                nztemp = full(logical(abs(S)*reshape(nze, [ny_, nx1_*nx2_])));
                nonzero_map(nzkey_f) = reshape(nztemp, [nx, nx1_, nx2_]);
            end
        end
        
        % Remove zero terms
        nzders(iszero) = [];
        nzterms(iszero) = [];
        nzdens(iszero) = [];
        
        d2ydx2dx1_ = initializeMatrixMupad(nzterms, nzdens, nzders, ny_*nx1_, nx2_);
        
    end

    function matout = reshape_derivative(mat, newsize, num, dens)
        
        % Get logical matrix indicating non-zero entries
        nzlogical = getNonZeroEntries(num, dens);
        
        % Convert logicals into linear indices
        nzindices = find(nzlogical);
        
        % Convert linear indices into subscripts in new reshaped matrix
        [nzsub_i,nzsub_j] = ind2sub(newsize, nzindices);
        
        % Get nonzero terms
        nzterms = mat(nzindices);
        
        % Initialize matrix 
        matout = initializeMatrixMupad(nzsub_i, nzsub_j, nzterms, newsize(1), newsize(2));
        
    end

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
            if isKey(nonzero_map,nzkey)
                nzout = nonzero_map(nzkey);
                nzfound = true;
            else
                nzkey = nzkey(1:end-1);
            end
            if isempty(nzkey)
                error('No data found for this derivative. This should never happen and needs to be debugged.')
            end
        end
        
        ndens = numel(dens);
        ndensknown = numel(nzkey) - 1;
        denstoestimate = (ndensknown+1):ndens;
        singletondims = 3:(ndens+1);
        for nzi = denstoestimate(:)'
            nzaddkey = [num dens{nzi}];
            nztoadd = nonzero_map(nzaddkey);
            if numel(singletondims) <= 1
                permutedims = [1, singletondims(1:nzi-1), 2];
            else
                permutedims = [1, singletondims(1:nzi-1), 2, singletondims(nzi:end)];
            end
            permutednztoadd = permute(nztoadd,permutedims);
            nzout = bsxfun(@and,nzout,permutednztoadd);
        end
    end
end

function out = expression_has_variable(expr_strs, var_strs)
% Determine which expressions have which variables. Return a logical matrix
% n_expr by n_var with true whenever the ith expression contains the jth
% variable.

n_exprs = numel(expr_strs);
n_vars = numel(var_strs);

if isempty(expr_strs) || isempty(var_strs)
    out = false(n_exprs, n_vars);
else
    
    % Determine where the variable's substrings appear in the r
    % expressions. Cell array varPositionsInExpr element i,j contains an
    % array indicating at what positions var_strs(j) appears in expr_strs(i)
    varPositionsInExpr = cellfun(@strfind,...
        repmat(expr_strs, 1, n_vars),...
        repmat(row(var_strs), n_exprs, 1),...
        'UniformOutput', false...
        );
    
    % Determine which elements are empty. If varPositionsInExpr{i,j} is
    % not empty, then reaction i contains variable j.
    out = ~cellfun(@isempty, varPositionsInExpr);
end
end

function fun = string2fun(string_rep, num, dens)
% Note that string2fun is a subfunction instead of a nested function to
% prevent the anonymous functions created here from saving copies of
% the primary function workspace variables.

% If the dependent variable is x0, the input argument is s,k. Otherwise,
% t,x,u,k.
if strcmp(num,'x')
    inputargstr = 's,k';
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

function fun = setfun_x0(basefun, k)
fun = @(s)basefun(s,k);
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

function rOut = evaluate_external_functions(rIn, ids)
if ~isempty(rIn)
    % Evaluate symbolic functions/pull in functions defined in path
    %   Necessary for "power" and other MathML function translation
    % Initialize symbolic variables
    syms(ids{:});
    
    % Evaluate the expressions to remove function calls
    rOut = eval(rIn);
    
    % Clear the symbolic variables
    clear(ids{:})
else
    % eval does not preserve size if rIn is empty
    rOut = rIn;
end
end

function name = replace_special_chars(name)
name = regexprep(name, '[^a-zA-Z0-9]', '_');
end

function all_ids = check_unique_name_id(type, new_names, new_ids, old_names, all_ids)
new_names = vec(new_names);
new_ids = vec(new_ids);
old_names = vec(old_names);
all_ids = vec(all_ids);

n_new = numel(new_names);

for i = 1:n_new
    % Check if item by this name already exists
    if ~isempty(new_names)
        name = new_names{i};
        if any(strcmp(name, [old_names; new_names(1:i-1)]))
            error('KroneckerBio:FinalizeModel:RepeatName', ['There is already a ' type ' with the name %s (new ' type ' #%i)'], name, i)
        end
    end
    
    % Check if ID has already been used
    id = new_ids{i};
    if ~isempty(id) && any(strcmp(id, [all_ids; new_ids(1:i-1)]))
        error('KroneckerBio:FinalizeModel:RepeatID', ['There is already a component with the id %s (new ' type ' #%i)'], id, i)
    end
end

all_ids = [all_ids; new_ids(~strcmp('', new_ids))];
end

function [ids, all_ids] = update_unique_ids(names, ids, all_ids)
n = numel(names);

for i = 1:n
    if isempty(ids{i})
        % The standard ID is the name with all symbols converted to underscores
        root_id = replace_special_chars(names{i});
        id = root_id;
        
        % If the ID already exists, append a different number until it is unique
        index = 1;
        while ismember(id, all_ids)
            id = [root_id '_' num2str(index)];
            index = index + 1;
        end
        
        % Apply the new ID
        all_ids = [all_ids; {id}];
        ids{i} = id;
    end
end
end

function strs = replace_names(strs, names, ids)
% Replace all quoted names in strs with their respective IDs

if ischar(strs)
    strs = {strs};
end
n = numel(strs);

% Find all places where there are names embedded in the string expression
[all_starts, all_ends, all_matches] = regexp(strs, '"[^"]*"', 'start', 'end', 'match');

for i_str = 1:n
    starts = all_starts{i_str};
    ends = all_ends{i_str};
    matches = cellfun(@(str){str(2:end-1)}, all_matches{i_str}); % Strip the quotes
    inds_matching = lookup(matches, names);
    
    assert(all(inds_matching), 'KroneckerBio:finalizeModelAnalytic:UnknownName', 'In an expression (%s) a name was not found among the names of the model', strs{i_str})
    
    n_matches = numel(matches);
    for i_match = n_matches:-1:1
        % Keep the string before the quote and after the quote and replace the middle with the ID
        strs{i_str} = [strs{i_str}(1:starts(i_match)-1), ids{inds_matching(i_match)}, strs{i_str}(ends(i_match)+1:end)];
    end
end
end
