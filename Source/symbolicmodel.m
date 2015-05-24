function [ m, expectedexprs ] = symbolicmodel(opts)
% m = symbolicmodel(opts)
% Builds a simple symbolic model for tests
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

m.Type       = 'Model.SymbolicReactions';
m.Name       = 'Symbolic Model Test';

% Compartments
compartmenttable = {
    'solution'  1
    };

nv = size(compartmenttable,1);
v = cell2mat(compartmenttable(:,2));
vNames = compartmenttable(:,1);

m.nv         = nv;
m.vSyms      = sym(vNames);
m.vNames     = vNames;
m.dv         = zeros(nv,1) + 3;
m.v          = v;

% Parameters
parametertable = {
    'kf'        1
    'kr'        1
    'koff'      1
    'kon'       1
    };

nk = size(parametertable,1);
k = cell2mat(parametertable(:,2));
kNames = parametertable(:,1);

m.nk         = nk;
m.kSyms      = sym(kNames);
m.kNames     = kNames;
m.k          = k;

% Species
speciestable = {
    'A'     1 
    'B'     3       
    'C'     2
    'DD'     4
    };

isu = ismember(speciestable(:,1),opts.Inputs);
nx = sum(~isu);
nu = sum(isu);
xNames = speciestable(~isu,1);
uNames = speciestable(isu,1);
xuNames = speciestable(:,1);
u = cell2mat(speciestable(isu,2));

m.nx         = nx;
m.xSyms      = sym(xNames);
m.xNames     = xNames;
m.vxInd      = ones(nx,1);

m.nu         = nu;
m.uSyms      = sym(uNames);
m.uNames     = uNames;
m.vuInd      = ones(nu,1);
m.u          = u;

% Seeds
ns = nx;
s = cell2mat(speciestable(~isu,2));
sNames = strcat(xNames,'_0');

m.sSyms      = sym(sNames);
m.sNames     = sNames;
m.s          = s;
m.x0         = sym(sNames);

m.ns         = ns;

% Reactions
reactiontable = {
    {'A'}       {'B'}       'kf*A'  'A -> B'
    {'B'}       {'A'}       'kr*B'  'B -> A'
    {'B'}       {'C' 'DD'}   'koff*B' 'B -> C + DD'
    {'C' 'DD'}   {'B'}       'kon*C*DD' 'C + DD -> B'
    };

nr = size(reactiontable,1);
r = sym(reactiontable(:,3));

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

m.nr         = nr;
m.rNames     = rNames;
m.r          = r;
m.S          = S(~isu,:);
m.Su         = S(isu,:);

% Outputs
outputtable = opts.Outputs(:);
y = sym(outputtable(:,1));
yNames = outputtable(:,1);

m.y          = y;
m.yNames     = yNames;

% Get expected values for expressions in model, if requested
if nargout > 1
    
norder = 2;
exprvals = cell(norder+1,1);
exprnames = cell(norder+1,1);

f = S(~isu,:)*r;

x = 2*m.s;
u = 3*m.u;
k = m.k;

% Names of symbolic variables to be substituted for values
subvars = [m.xSyms;m.uSyms;m.kSyms];

subfun = @(expr) double(subs(expr,subvars,[x;u;k])); % sym-to-double function
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

exprvals = vertcat(exprvals{:});
exprnames = vertcat(exprnames{:});

% Remove derivatives of y wrt k
isykexpr = ismember(exprnames,{'dydk';'d2ydk2';'d2ydkdx';'d2ydkdu';'d2ydxdk';'d2ydudk'});
exprvals(isykexpr) = [];
exprnames(isykexpr) = [];

expectedexprs = cell2struct(exprvals,exprnames,1);

expectedexprs.x = x;
expectedexprs.u = u;
expectedexprs.k = k;

end

    function [dexprs,dexprvars] = getDerivatives(exprs,exprvars)
        % Takes a list of symbolic expressions exprs and takes the
        % derivatives of them wrt x, k, and u, returning them in a
        % 3*length(exprs) cell array dexprs. Also takes a list of short names
        % exprvars (i.e., fx) and updates them to reflect the names of the
        % derivatives taken in dexprvars.
        dexprs = cell(3*length(exprs),1);
        dexprvars = cell(size(dexprs));
        getder = @(expr,vars) reshape(jacobian(expr(:),vars),numel(expr),numel(vars));
        for ei = 1:length(exprs)
            dexprs{3*(ei-1)+1} = getder(exprs{ei},m.xSyms);
            dexprvars{3*(ei-1)+1} = [exprvars{ei} 'x'];
            dexprs{3*(ei-1)+2} = getder(exprs{ei},m.kSyms);
            dexprvars{3*(ei-1)+2} = [exprvars{ei} 'k'];
            dexprs{3*(ei-1)+3} = getder(exprs{ei},m.uSyms);
            dexprvars{3*(ei-1)+3} = [exprvars{ei} 'u'];
        end
    end

    function exprname = getDerivativeName(exprvar)
        % Converts short representation of derivative name to full name,
        % i.e., fxu -> d2fdudx
        dnum = length(exprvar)-1;
        if dnum == 0
            exprname = exprvar(1);
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
        exprname = [numstr denstr];
    end

end

