function [m, expectedexprs] = analytic_model_syms(opts)
% Builds a simple analytic model and symbolically calculates expected
% components. Note: this function is slow due to the FinalizeModel.
%
% A -> B -> C + DD
%
% opts (optional)
%   .Inputs [ cell array of strings ]
%       Species to be declared inputs at constant concentrations.
%   .Outputs [ cell array of strings [{'A+B'};{'C+DD'}] ]
%       Expressions of species to be set as outputs.

if nargin < 1
    opts = [];
end
if isempty(opts)
    opts = struct;
end

% Default inputs and outputs
if ~isfield(opts,'Inputs')
    opts.Inputs = {'A';'C'};
end
if ~isfield(opts,'Outputs')
    opts.Outputs = {'A+B';'C+DD'};
end

% Standardized empty inputs and outputs as empty cell arrays
if isempty(opts.Inputs)
    opts.Inputs = {};
end
if isempty(opts.Outputs)
    opts.Outputs = {};
end

% Silence the symbolic warning for certain Matlab versions
if verLessThan('matlab', '9.3') && ~verLessThan('matlab', '9.0')
    state = warning('off', 'symbolic:sym:sym:DeprecateExpressions');
    finished = onCleanup(@() warning(state));
end

%% Initialize models
sm.Type       = 'Model.SymbolicReactions';
sm.Name       = 'Symbolic Model Test';

m = InitializeModelAnalytic('Symbolic Model Test');

%% Compartments
compartmenttable = {
    'solution'  1
    };

nv = size(compartmenttable,1);
v = cell2mat(compartmenttable(:,2));
vNames = compartmenttable(:,1);

sm.nv         = nv;
sm.vSyms      = sym(vNames);
sm.vNames     = vNames;
sm.dv         = zeros(nv,1) + 3;
sm.v          = v;

m = AddCompartment(m, vNames{1}, 3, v(1));

%% Parameters
parametertable = {
    'kf'        1
    'kr'        1
    'koff'      1
    'kon'       1
    };

nk = size(parametertable,1);
k = cell2mat(parametertable(:,2));
kNames = parametertable(:,1);

sm.nk         = nk;
sm.kSyms      = sym(kNames);
sm.kNames     = kNames;
sm.k          = k;

for i = 1:nk
    m = AddParameter(m, kNames{i}, k(i));
end

%% Species
% By default, A and C are inputs (values as below) and B and DD are states 
%   (with ICs as seeds, squaring the below)
speciestable = {
    'A'     1 
    'B'     3       
    'C'     2
    'DD'    4
    };

isu = ismember(speciestable(:,1),opts.Inputs);
nx = sum(~isu);
nu = sum(isu);
xNames = speciestable(~isu,1);
uNames = speciestable(isu,1);
xuNames = speciestable(:,1);
u = cell2mat(speciestable(isu,2));

sm.nx         = nx;
sm.xSyms      = sym(xNames);
sm.xNames     = xNames;
sm.vxInd      = ones(nx,1);

sm.nu         = nu;
sm.uSyms      = sym(uNames);
sm.uNames     = uNames;
sm.vuInd      = ones(nu,1);
sm.u          = u;

% Seeds
ns = nx;
s = cell2mat(speciestable(~isu,2));
sNames = strcat(xNames,'_0');

sm.sSyms      = sym(sNames);
sm.sNames     = sNames;
sm.s          = s;
sm.x0         = sym(sNames).^2;

sm.ns         = ns;

for i = 1:nx+nu
    xu0 = speciestable{i,2};
    if isu(i) % input
        m = AddInput(m, xuNames{i}, 'solution', xu0);
    else % state
        m = AddSeed(m, [xuNames{i} '_0'], xu0);
        m = AddState(m, xuNames{i}, 'solution', [xuNames{i} '_0^2']);
    end
end

%% Reactions
reactiontable = {
    {'A'}       {'B'}       'kf*A'  'A -> B'
    {'B'}       {'A'}       'kr*B'  'B -> A'
    {'B'}       {'C' 'DD'}   'pi*koff*B' 'B -> C + DD' % pi included in rate form to test behavior of predefined constants 
    {'C' 'DD'}   {'B'}       'kon*C*DD' 'C + DD -> B'
    };

nr = size(reactiontable,1);
if ~verLessThan('matlab', '9.3')
    r = str2sym(reactiontable(:,3));
else
    r = sym(reactiontable(:,3));
end

% Get indices of reactants and products in each reaction
getspeciesindices = @(x)intersect(x,xuNames,'stable');
[~,~,reactants] = cellfun(getspeciesindices,reactiontable(:,1),'UniformOutput',false);
[~,~,products] = cellfun(getspeciesindices,reactiontable(:,2),'UniformOutput',false);

% Get linear indices of nonzero entries in the stoichiometric matrix
getSindices = @(x,ri) sub2ind([nx+nu,nr],x,repmat(ri,size(x))); % index linearization function
reactantindices = cellfun(getSindices,reactants,num2cell(1:nr)','UniformOutput',false);
productindices = cellfun(getSindices,products,num2cell(1:nr)','UniformOutput',false);

rNames = reactiontable(:,4);

% Create stoichiometric matrix
S = zeros(nx+nu,nr);
S(vertcat(reactantindices{:})) = -1;
S(vertcat(productindices{:})) = 1;

sm.nr         = nr;
sm.rNames     = rNames;
sm.r          = r;
sm.S          = S(~isu,:);
sm.Su         = S(isu,:);

for i = 1:nr
    m = AddReaction(m, reactiontable{i,4}, reactiontable{i,1}, reactiontable{i,2}, reactiontable{i,3});
end

%% Outputs
outputtable = opts.Outputs(:);
yNames = outputtable(:,1);
yStrings = outputtable(:,1);
sm.yNames = yNames;
sm.yStrings = yStrings;
if ~verLessThan('matlab', '9.3')
    sm.y = str2sym(yStrings);
else
    sm.y = sym(yStrings);
end
y = sm.y;

ny = length(y);
for i = 1:ny
    m = AddOutput(m, yNames{i}, yStrings{i});
end

%% Finalize test model
m = FinalizeModel(m);

%% Get expected values for expressions in model, if requested
if nargout > 1
    
    norder = 2;

    f = S(~isu,:)*r;
    x0 = sm.x0;
    
    % Set some random test values for comparing constructed model and symbolic
    % solutions
    x = 2*sm.s;
    u = 3*sm.u;
    k = sm.k;

    %%% Get expected values for f, r, and y functions %%%
    % Names of symbolic variables to be substituted for values
    subvars = [sm.xSyms;sm.uSyms;sm.kSyms];

    subfun = @(expr) double(subs(expr,subvars,[x;u;k])); % sym-to-double function
    exprvals = cell(norder+1,1);
    exprnames = cell(norder+1,1);
    exprs = {f;r;y}; % Symbolic expressions
    exprvars = {'f';'r';'y'}; % short expressions representing which derivative is which, i.e., 'fxu' is d2f/dudx
    for oi = 0:norder
        if oi > 0
            % Take symbolic derivatives and update names to
            % reflect derivatives
            [exprs,exprvars] = getDerivatives(exprs,exprvars);
        end
        % Substitute in values for symbolics to generate values
        exprvals{oi+1} = cellfun(subfun,exprs,'UniformOutput',false);
        exprnames{oi+1} = cellfun(@getDerivativeName,exprvars,'UniformOutput',false);
    end
    %%%
    
    %%% Get expected values for x0 functions %%%
    subvars_x0 = sm.sSyms;
    
    subfun_x0 = @(expr) double(subs(expr,subvars_x0,s)); % sym-to-double function
    exprvals_x0 = cell(norder+1,1);
    exprnames_x0 = cell(norder+1,1);
    exprs_x0 = {x0};
    exprvars_x0 = {'x'};
    for oi = 0:norder
        if oi > 0
            % Take symbolic derivatives and update names to
            % reflect derivatives
            [exprs_x0,exprvars_x0] = getDerivatives(exprs_x0,exprvars_x0,true);
        end
        % Substitute in values for symbolics to generate values
        exprvals_x0{oi+1} = cellfun(subfun_x0,exprs_x0,'UniformOutput',false);
        exprnames_x0{oi+1} = cellfun(@getDerivativeName,exprvars_x0,'UniformOutput',false);
    end
    %%%
    
    % Concatenate f, r, y, and x0 function values and names
    exprvals = vertcat(exprvals{:},exprvals_x0{:});
    exprnames = vertcat(exprnames{:},exprnames_x0{:});

    expectedexprs = cell2struct(exprvals,exprnames,1);

    expectedexprs.x = x;
    expectedexprs.u = u;
    expectedexprs.k = k;
    expectedexprs.s = s;

end

    function [dexprs, dexprvars] = getDerivatives(exprs, exprvars, isx0)
        % Takes a list of symbolic expressions exprs and takes the
        % derivatives of them wrt x, k, and u (or s if isx0 is true),
        % returning them in a 3*length(exprs) cell array dexprs. Also takes
        % a list of short names exprvars (i.e., fx) and updates them to
        % reflect the names of the derivatives taken in dexprvars.
        if nargin < 3
            isx0 = false;
        end
        
        getder = @(expr,vars) reshape(jacobian(expr(:),vars),numel(expr),numel(vars));
        if isx0
            dexprs = cell(length(exprs),1);
            dexprvars = cell(size(dexprs));
            for ei = 1:length(exprs)
                dexprs{ei} = getder(exprs{ei},sm.sSyms);
                dexprvars{ei} = [exprvars{ei} 's'];
            end
        else
            dexprs = cell(3*length(exprs),1);
            dexprvars = cell(size(dexprs));
            for ei = 1:length(exprs)
                dexprs{3*(ei-1)+1} = getder(exprs{ei},sm.xSyms);
                dexprvars{3*(ei-1)+1} = [exprvars{ei} 'x'];
                dexprs{3*(ei-1)+2} = getder(exprs{ei},sm.kSyms);
                dexprvars{3*(ei-1)+2} = [exprvars{ei} 'k'];
                dexprs{3*(ei-1)+3} = getder(exprs{ei},sm.uSyms);
                dexprvars{3*(ei-1)+3} = [exprvars{ei} 'u'];
            end
        end
    end

    function exprname = getDerivativeName(exprvar)
        % Converts short representation of derivative name to full name,
        % i.e., fxu -> d2fdudx
        dnum = length(exprvar)-1;
        if dnum == 0
            exprname = exprvar(1);
            exprname = regexprep(exprname,'x','x0'); % Replace x with x0 (I use x as a single-character representation of x0)
            return
        elseif dnum == 1
            numstr = ['d' exprvar(1)];
        else
            numstr = ['d' num2str(dnum) exprvar(1)];
        end
        dens = exprvar(end:-1:2);
        if dnum == 2 && strcmp(dens(1),dens(2))
            denstr = ['d' dens(1) '2'];
        elseif dnum == 2
            denstr = ['d' dens(1) 'd' dens(2)];
        elseif dnum == 1
            denstr = ['d' dens(1)];
        else
            error('Higher order deriviatives not currently supported.')
        end
        numstr = regexprep(numstr,'x','x0'); % Replace x with x0 (I use x as a single-character representation of x0)
        exprname = [numstr denstr];
    end

end

