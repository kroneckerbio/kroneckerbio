function m = symbolic2PseudoKronecker(SymModel, uNames, yNames, yMembers, yValues, opts)
%symbolic2PseudoKronecker converts a symbolic model into a pseudo-kronecker
%   model, which interacts with the Kronecker Bio toolbox exactly like a
%   Kronecker model, but with substantial performance reductions and loss
%   of abilities to modify the topology.
% 
%   m = symbolic2PseudoKronecker(SymModel, uNames, yNames, yMembers, yValues, opts)
% 
%   Inputs
%   SymModel: [ symbolic model scalar ]
%       A symbolic model
%   uNames: [ string | cell vector of strings | symbolic vector {[]} ]
%       Additional species that will be designated as inputs
%   yNames: [ string | cell vector of strings | symbolic vector  | 
%             positive integer vector {1:nx} ]
%       Names for the designated outptus. If yMembers is not provided,
%       the names are also interpreted as the species that will be
%       represented by the outputs. Naturally, if yMembers is not
%       provided, then each output can only reflect the concentration
%       of a single species.
%   yMembers: [ cell vector of cell vectors of strings ]
%       Regular expressions corresponding the species that will contribute
%       to the value of the output
%   yValues: [ cell vector of nonegative vectors ]
%       Each numeric entry in this vector is associated with one of the
%       expressions. This value tells how much a species will contribute
%       when it matches the corresponding expressions. If a species is
%       matched by multiple expressions, the last expression to match
%       overrides all others.
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

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Clean up inputs
if nargin < 6
    opts = [];
    if nargin < 5
        yValues = [];
        if nargin < 4
            yMembers = [];
            if nargin < 3
                yNames = [];
            end
        end
    end
end

% Default options
defaultOpts.Order             = 2;
defaultOpts.VolumeToParameter = false;
defaultOpts.Verbose           = 0;

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

%% Extract symbolic values
Name    = SymModel.Name;

nv      = SymModel.nv;
nxu     = SymModel.nxu;
nk      = SymModel.nk;
nr      = SymModel.nr;

vSyms   = SymModel.vSyms;
vNames  = SymModel.vNames;
dv      = SymModel.dv;
vValues = SymModel.vValues;

xuSyms   = SymModel.xuSyms;
xuNames  = SymModel.xuNames;
vxuNames = SymModel.vxuNames;
xu0      = SymModel.xu0;
isu      = SymModel.isu;

kSyms   = SymModel.kSyms;
kNames  = SymModel.kNames;
k       = SymModel.k;

rNames = SymModel.rNames;
r      = SymModel.r;
S      = SymModel.S;

f      = SymModel.f;

% Convert compartment volumes to parameters or constants
if opts.VolumeToParameter
    nk = nk + nv;
    kNames = [kNames; vNames];
    k = [k; vValues];
    kSyms = [kSyms; vSyms];
else% ~VolumeToParameter
    r = subs(r, vSyms, vValues, 0);
    f = subs(f, vSyms, vValues, 0);
end

% Generate derivatives of desired order
if opts.Order >= 1
%Gradient of r with respect to x
drdx = jacobian(r, xuSyms);

%Gradient of r with respect to k
drdk = jacobian(r, kSyms);

%Gradient of f with respect to x
dfdx = S*drdx;

%Gradient of f with respect to k
dfdk = S*drdk;
else
drdx = '';
drdk = '';
dfdx = '';
dfdk = '';
end

if opts.Order >= 2
%Gradient of drdx with respect to x
d2rdx2 = jacobian(vec(drdx), xuSyms);

%Gradient of drdk with respect to k
d2rdk2 = jacobian(vec(drdk), kSyms);

%Gradient of drdx with respect to k
d2rdkdx = jacobian(vec(drdx), kSyms);

%Gradient of drdk with respect to x
d2rdxdk = jacobian(vec(drdk), xuSyms);

%Gradient of dfdx with respect to x
d2fdx2 = S*reshape(d2rdx2, nr,nxu*nxu);
d2fdx2 = reshape(d2fdx2, nxu*nxu,nxu);

%Gradient of dfdk with respect to k
d2fdk2 = S*reshape(d2rdk2, nr,nk*nk);
d2fdk2 = reshape(d2fdk2, nxu*nk,nk);

%Gradient of dfdx with respect to k
d2fdkdx = S*reshape(d2rdkdx, nr,nxu*nk);
d2fdkdx = reshape(d2fdkdx, nxu*nxu,nk);

%Gradient of dfdk with respect to x
d2fdxdk = S*reshape(d2rdxdk, nr,nk*nxu);
d2fdxdk = reshape(d2fdxdk, nxu*nk,nxu);
else
d2rdx2  = '';
d2rdk2  = '';
d2rdkdx = '';
d2rdxdk = '';
d2fdx2  = '';
d2fdk2  = '';
d2fdkdx = '';
d2fdxdk = '';
end

if opts.Order >= 3
%Gradient of d2rdx2 with respect to x
d3rdx3 = jacobian(vec(d2rdx2), xuSyms);

%Gradient of d2rdx2 with respect to k
d3rdkdx2 = jacobian(vec(d2rdx2), kSyms);

%Gradient of d2fdx2 with respect to x
d3fdx3 = S*reshape(d3rdx3, nr,nxu*nxu*nxu);
d3fdx3 = reshape(d3fdx3, nxu*nxu*nxu,nxu);

%Gradient of d2fdx2 with respect to k
d3fdkdx2 = S*reshape(d3rdkdx2, nr,nxu*nxu*nk);
d3fdkdx2 = reshape(d3fdkdx2, nxu*nxu*nxu,nk);
else
d3rdx3   = '';
d3rdkdx2 = '';
d3fdx3   = '';
d3fdkdx2 = '';
end

xuStrs = cell(nxu,1);
for i = 1:nxu
    xuStrs{i} = char(xuSyms(i));
end

kStrs = cell(nk,1);
for i = 1:nk
    kStrs{i} = char(kSyms(i));
end

%% Process additional inputs
% Standardize uNames as a cell array of strings
if isa(uNames, 'char')
    uAddInd = find(strcmp(xuNames, uNames));
elseif isa(uNames, 'sym')
    temp = uNames;
    uNames = cell(length(uNames),1);
    for i = 1:length(uNames)
        uNames{i} = char(temp(i));
        uAddInd = find(strcmp(xuNames, uNames{i}));
    end
elseif isa(uNames, 'numeric')
    uAddInd = uNames;
elseif isa(uNames, 'cell')
    uAddInd = zeros(length(uNames),1);
    for i = 1:length(uNames)
        uAddInd(i) = find(strcmp(xuNames, uNames{i}));
    end
end

% Add new inputs
isu(uAddInd) = true;

% Seperate species into states and inputs
xNames = xuNames(~isu);
uNames = xuNames(isu);

xStrs = xuStrs(~isu);
uStrs = xuStrs(isu);

x0 = xu0(~isu);
u  = xu0(isu);

vxNames = vxuNames(~isu);
vuNames = vxuNames(isu);

vxInd = lookup(vxuNames(~isu), vNames);
vuInd = lookup(vxuNames(isu), vNames);

xNamesFull = strcat(vxNames, '.', xNames);
uNamesFull = strcat(vuNames, '.', uNames);

nx = numel(xStrs);
nu = numel(uStrs);

%% Delete the rows and columns corresponding to the input
% f
f(isu,:) = [];

% dfdx, keeping dfdu
if ~isempty(dfdx)
    dfdx(isu,:) = [];
    dfdu        = dfdx(:,isu);
    dfdx(:,isu) = [];
end

% dfdk
if ~isempty(dfdk)
    dfdk(isu,:) = [];
end

if ~isempty(d2fdx2)
    d2fdxu2 = d2fdx2;
    
    % d2fdx2
    d2fdx2              = reshape(d2fdxu2, nxu,nxu,nxu);
    d2fdx2(isu,:,:)     = [];
    d2fdx2(:,isu,:)     = [];
    d2fdx2(:,:,isu)     = [];
    d2fdx2              = reshape(d2fdx2, nx*nx,nx);

    % d2fdu2
    d2fdu2              = reshape(d2fdxu2, nxu,nxu,nxu);
    d2fdu2(isu,:,:)     = [];
    d2fdu2(:,~isu,:)    = [];
    d2fdu2(:,:,~isu)    = [];
    d2fdu2              = reshape(d2fdu2, nx*nu,nu);

    % d2fdudx
    d2fdudx             = reshape(d2fdxu2, nxu,nxu,nxu);
    d2fdudx(isu,:,:)    = [];
    d2fdudx(:,isu,:)    = [];
    d2fdudx(:,:,~isu)   = [];
    d2fdudx             = reshape(d2fdudx, nx*nx,nu);

    % d2fdxdu
    d2fdxdu             = reshape(d2fdxu2, nxu,nxu,nxu);
    d2fdxdu(isu,:,:)    = [];
    d2fdxdu(:,~isu,:)   = [];
    d2fdxdu(:,:,isu)    = [];
    d2fdxdu             = reshape(d2fdxdu, nx*nu,nx);
else
    d2fdu2  = '';
    d2fdudx = '';
    d2fdxdu = '';
end

% d2fdk2
if ~isempty(d2fdk2)
    d2fdk2          = reshape(d2fdk2, nxu,nk,nk);
    d2fdk2(isu,:,:) = [];
    d2fdk2          = reshape(d2fdk2, nx*nk,nk);
end

% d2fdkdx
if ~isempty(d2fdkdx)
    d2fdkdx            = reshape(d2fdkdx, nxu,nxu,nk);
    d2fdkdx(isu,:,:)   = [];
    d2fdkdx(:,isu,:)   = [];
    d2fdkdx            = reshape(d2fdkdx, nx*nx,nk);
end

% d2fdxdk
if ~isempty(d2fdxdk)
    d2fdxdk          = reshape(d2fdxdk, nxu,nk,nxu);
    d2fdxdk(isu,:,:) = [];
    d2fdxdk(:,:,isu) = [];
    d2fdxdk          = reshape(d2fdxdk, nx*nk,nx);
end

% d3fdx3
if ~isempty(d3fdx3)
    d3fdx3            = reshape(d3fdx3, nxu,nxu,nxu,nxu);
    d3fdx3(isu,:,:,:) = [];
    d3fdx3(:,isu,:,:) = [];
    d3fdx3(:,:,isu,:) = [];
    d3fdx3(:,:,isu,:) = [];
    d3fdx3            = reshape(d3fdx3, nx*nx*nx,nx);
end

% d3fdkdx2
if ~isempty(d3fdkdx2)
    d3fdkdx2            = reshape(d3fdkdx2, nxu,nxu,nxu,nk);
    d3fdkdx2(isu,:,:,:) = [];
    d3fdkdx2(:,isu,:,:) = [];
    d3fdkdx2(:,:,isu,:) = [];
    d3fdkdx2            = reshape(d3fdkdx2, nx*nx*nx,nk);
end

% S
S(isu,:) = [];

% drdx
if ~isempty(drdx)
    drdu        = drdx(:,isu);
    drdx(:,isu) = [];
end

% d2rdx2
if ~isempty(d2rdx2)
    d2rdx2            = reshape(d2rdx2, nr,nxu,nxu);
    d2rdx2(:,isu,:)   = [];
    d2rdx2(:,:,isu)   = [];
    d2rdx2            = reshape(d2rdx2, nr*nx,nx);
end

% d2rdkdx
if ~isempty(d2rdkdx)
    d2rdkdx          = reshape(d2rdkdx, nr,nxu,nk);
    d2rdkdx(:,isu,:) = [];
    d2rdkdx          = reshape(d2rdkdx, nr*nx,nk);
end

% d2rdxdk
if ~isempty(d2rdxdk)
    d2rdxdk          = reshape(d2rdxdk, nr,nk,nxu);
    d2rdxdk(:,:,isu) = [];
    d2rdxdk          = reshape(d2rdxdk, nr*nk,nx);
end

%% Process the output
% Standardize yNames as a cell vector of strings
% and yMembers as a cell array of cell arrays of strings
if isempty(yMembers)
    % yNames are the yMembers
    if isempty(yNames)
        % All species are outputs
        yNames = 1:nx;
    end
    
    if isa(yNames, 'char')
        ny = 1;
        yNames = {yNames};
        yMembers = {yNames};
    elseif isa(yNames, 'sym')
        ny = length(yNames);
        temp = yNames;
        yNames = cell(ny,1);
        yMembers = cell(ny,1);
        for iy = 1:ny
            yNames{iy} = char(temp(iy));
            yMembers{iy} = yNames(iy);
        end
    elseif isa(yNames, 'numeric')
        ny = length(yNames);
        yMembers = cell(ny,1);
        for iy = 1:ny
            yMembers{iy} = strcat('^', xNamesFull(yNames(iy)), '$');
        end
        yNames = xNames(yNames);
    elseif isa(yNames, 'cell')
        ny = length(yNames);
        yNames = vec(yNames);
        yMembers = cell(ny,1);
        for iy = 1:ny
            yMembers{iy} = yNames(iy);
        end
    else
        error('KroneckerBio:symbolic2PseudoKronecker', 'Invalid type supplied for yNames')
    end
    
    yValues = repmat({1}, ny,1);
else
    % yMembers and yValues is provided
    if isa(yNames, 'char')
        ny = 1;
        yNames = {yNames};
    elseif isa(yNames, 'sym')
        ny = length(yNames);
        temp = yNames;
        yNames = cell(ny,1);
        for iy = 1:ny
            yNames{iy} = char(temp(iy));
        end
    elseif isa(yNames, 'numeric')
        ny = length(yNames);
        yNames = xNames(yNames);
    elseif isa(yNames, 'cell')
        ny = length(yNames);
        yNames = vec(yNames);
    else
        error('KroneckerBio:symbolic2PseudoKronecker', 'Invalid type supplied for yNames')
    end
    
    yMembers = vec(yMembers);
    yValues = vec(yValues);
end

%% Construct output matrices
% Entries in each sparse matrix for compartment conversion
nC1Entries = 0;
nC2Entries = 0;
ncEntries = 0;

C1Entries = zeros(0,2);
C1Values  = zeros(0,1);
C2Entries = zeros(0,2);
C2Values  = zeros(0,1);
cEntries  = zeros(0,2);
cValues   = zeros(0,1);

for iy = 1:ny
    nExpr = numel(yMembers{iy});
    for iExpr = 1:nExpr
        % Find states that match the expression
        match = find(~cellfun(@isempty, regexp(xNamesFull, yMembers{iy}{iExpr}, 'once')));
        nAdd = numel(match);
        nC1Entries = nC1Entries + nAdd;
        
        % Add more room in vector if necessary
        currentLength = size(C1Entries,1);
        if nC1Entries > currentLength
            addlength = max(currentLength, nAdd);
            C1Entries = [C1Entries; zeros(addlength,2)];
            C1Values  = [C1Values;  zeros(addlength,1)];
        end
        
        % Add entries
        C1Entries(nC1Entries-nAdd+1:nC1Entries,1) = iy;
        C1Entries(nC1Entries-nAdd+1:nC1Entries,2) = match;
        C1Values(nC1Entries-nAdd+1:nC1Entries) = yValues{iy}(iExpr);
        
        % Find inputs that match the expression
        match = find(~cellfun(@isempty, regexp(uNamesFull, yMembers{iy}{iExpr}, 'once')));
        nAdd = numel(match);
        nC2Entries = nC2Entries + nAdd;
        
        % Add more room in vector if necessary
        currentLength = size(C2Entries,1);
        if nC2Entries > currentLength
            addlength = max(currentLength, nAdd);
            C2Entries = [C2Entries; zeros(addlength,2)];
            C2Values  = [C2Values; zeros(addlength,1)];
        end
        
        % Add entries
        C2Entries(nC2Entries-nAdd+1:nC2Entries,1) = iy;
        C2Entries(nC2Entries-nAdd+1:nC2Entries,2) = match;
        C2Values(nC2Entries-nAdd+1:nC2Entries) = yValues{iy}(iExpr);
        
        % Find empty expressions, which are constants
        if isempty(yMembers{iy}{iExpr})
            ncEntries = ncEntries + 1;
            
            % Add more room in vector if necessary
            currentLength = size(cEntries,1);
            if nbEntries > currentLength
                addlength = max(currentLength, 1);
                cEntries = [cEntries; zeros(addlength,2)];
                cValues = [cValues; zeros(addlength,1)];
            end
            
            % Add entries
            cEntries(ncEntries,1) = iy;
            cEntries(ncEntries,2) = 1;
            cValues(ncEntries) = yValues{iy}(iExpr);
        end
    end
end

% Remove duplicate entries
[C1Entries ind] = unique(C1Entries(1:nC1Entries,:), 'rows');
C1Values = C1Values(ind);

[C2Entries ind] = unique(C2Entries(1:nC2Entries,:), 'rows');
C2Values = C2Values(ind);

[cEntries ind] = unique(cEntries(1:ncEntries,:), 'rows');
cValues = cValues(ind);

% Construct matrices
C1 = sparse(C1Entries(:,1), C1Entries(:,2), C1Values, ny, nx);
C2 = sparse(C2Entries(:,1), C2Entries(:,2), C2Values, ny, nu);
c  = sparse(cEntries(:,1),  cEntries(:,2),  cValues,  ny, 1);

%% Replace Verbose Names with Systematic names
if verbose; fprintf('Converting symbolics to strings...\n'); end
% Convert the symbolics into strings
f           = strtrim(evalc('disp(f)'));
if ~isempty(dfdx)
    dfdx    = strtrim(evalc('disp(dfdx)'));
    if nu > 0
        % Symbolic toolbox fails if uAddInd is empty
        dfdu    = strtrim(evalc('disp(dfdu)'));
    else
        % Use "zeros" if it will be empty
        dfdu    = sprintf('zeros(%d,0)',nx);
    end
else
    dfdu = '';
end
if ~isempty(dfdk)
    dfdk     = strtrim(evalc('disp(dfdk)'));
end
if ~isempty(d2fdx2)
    d2fdx2   = strtrim(evalc('disp(d2fdx2)'));
end
if ~isempty(d2fdu2)
    d2fdu2   = strtrim(evalc('disp(d2fdu2)'));
end
if ~isempty(d2fdk2)
    d2fdk2   = strtrim(evalc('disp(d2fdk2)'));
end
if ~isempty(d2fdudx)
    d2fdudx  = strtrim(evalc('disp(d2fdudx)'));
end
if ~isempty(d2fdxdu)
    d2fdxdu  = strtrim(evalc('disp(d2fdxdu)'));
end
if ~isempty(d2fdk2)
    d2fdk2   = strtrim(evalc('disp(d2fdk2)'));
end
if ~isempty(d2fdkdx)
    d2fdkdx  = strtrim(evalc('disp(d2fdkdx)'));
end
if ~isempty(d2fdxdk)
    d2fdxdk  = strtrim(evalc('disp(d2fdxdk)'));
end
if ~isempty(d3fdx3)
    d3fdx3   = strtrim(evalc('disp(d3fdx3)'));
end
if ~isempty(d3fdkdx2)
    d3fdkdx2 = strtrim(evalc('disp(d3fdkdx2)'));
end
r           = strtrim(evalc('disp(r)'));
if ~isempty(drdx)
    drdx    = strtrim(evalc('disp(drdx)'));
    if nu > 0
        % Symbolic toolbox fails if uAddInd is empty
        drdu    = strtrim(evalc('disp(drdu)'));
    else
        % Use "zeros" if it will be empty
        drdu    = sprintf('zeros(%d,0)',nx);
    end
else
    drdu = '';
end
if ~isempty(drdk)
    drdk     = strtrim(evalc('disp(drdk)'));
end
if ~isempty(d2rdx2)
    d2rdx2   = strtrim(evalc('disp(d2rdx2)'));
end
if ~isempty(d2rdk2)
    d2rdk2   = strtrim(evalc('disp(d2rdk2)'));
end
if ~isempty(d2rdkdx)
    d2rdkdx  = strtrim(evalc('disp(d2rdkdx)'));
end
if ~isempty(d2rdxdk)
    d2rdxdk  = strtrim(evalc('disp(d2rdxdk)'));
end

% Replace species names with vector index names
if verbose; fprintf('Replacing names with vectorized variables...\n'); end
if verbose; fprintf('   states...\n');end
for i = 1:nx
    name = sprintf('x(%d)', i);
    f        = regexprep(f, xStrs{i}, name, 0);
    dfdx     = regexprep(dfdx, xStrs{i}, name, 0);
    dfdu     = regexprep(dfdu, xStrs{i}, name, 0);
    dfdk     = regexprep(dfdk, xStrs{i}, name, 0);
    d2fdx2   = regexprep(d2fdx2, xStrs{i}, name, 0);
    d2fdu2   = regexprep(d2fdu2, xStrs{i}, name, 0);
    d2fdk2   = regexprep(d2fdk2, xStrs{i}, name, 0);
    d2fdudx  = regexprep(d2fdudx, xStrs{i}, name, 0);
    d2fdxdu  = regexprep(d2fdxdu, xStrs{i}, name, 0);
    d2fdkdx  = regexprep(d2fdkdx, xStrs{i}, name, 0);
    d2fdxdk  = regexprep(d2fdxdk, xStrs{i}, name, 0);
    d3fdx3   = regexprep(d3fdx3, xStrs{i}, name, 0);
    d3fdkdx2 = regexprep(d3fdkdx2, xStrs{i}, name, 0);
    r        = regexprep(r, xStrs{i}, name, 0);
    drdx     = regexprep(drdx, xStrs{i}, name, 0);
    drdu     = regexprep(drdu, xStrs{i}, name, 0);
    drdk     = regexprep(drdk, xStrs{i}, name, 0);
    d2rdx2   = regexprep(d2rdx2, xStrs{i}, name, 0);
    d2rdk2   = regexprep(d2rdk2, xStrs{i}, name, 0);
    d2rdkdx  = regexprep(d2rdkdx, xStrs{i}, name, 0);
    d2rdxdk  = regexprep(d2rdxdk, xStrs{i}, name, 0);
end

% Replace input namus with vector index names
if verbose; fprintf('   input...\n'); end
for i = 1:nu
    name = sprintf('u(%d)', i);
    f        = regexprep(f, uStrs{i}, name, 0);
    dfdx     = regexprep(dfdx, uStrs{i}, name, 0);
    dfdu     = regexprep(dfdu, uStrs{i}, name, 0);
    dfdk     = regexprep(dfdk, uStrs{i}, name, 0);
    d2fdx2   = regexprep(d2fdx2, uStrs{i}, name, 0);
    d2fdu2   = regexprep(d2fdu2, uStrs{i}, name, 0);
    d2fdk2   = regexprep(d2fdk2, uStrs{i}, name, 0);
    d2fdudx  = regexprep(d2fdudx, uStrs{i}, name, 0);
    d2fdxdu  = regexprep(d2fdxdu, uStrs{i}, name, 0);
    d2fdkdx  = regexprep(d2fdkdx, uStrs{i}, name, 0);
    d2fdxdk  = regexprep(d2fdxdk, uStrs{i}, name, 0);
    d3fdx3   = regexprep(d3fdx3, uStrs{i}, name, 0);
    d3fdkdx2 = regexprep(d3fdkdx2, uStrs{i}, name, 0);
    r        = regexprep(r, uStrs{i}, name, 0);
    drdx     = regexprep(drdx, uStrs{i}, name, 0);
    drdu     = regexprep(drdu, uStrs{i}, name, 0);
    drdk     = regexprep(drdk, uStrs{i}, name, 0);
    d2rdx2   = regexprep(d2rdx2, uStrs{i}, name, 0);
    d2rdk2   = regexprep(d2rdk2, uStrs{i}, name, 0);
    d2rdkdx  = regexprep(d2rdkdx, uStrs{i}, name, 0);
    d2rdxdk  = regexprep(d2rdxdk, uStrs{i}, name, 0);
end

% Replace parameters with vector index names
if verbose; fprintf('   parameters...\n');end
for i = 1:nk
    name = sprintf('k(%d)', i);
    f        = regexprep(f, kStrs{i}, name, 0);
    dfdx     = regexprep(dfdx, kStrs{i}, name, 0);
    dfdu     = regexprep(dfdu, kStrs{i}, name, 0);
    dfdk     = regexprep(dfdk, kStrs{i}, name, 0);
    d2fdx2   = regexprep(d2fdx2, kStrs{i}, name, 0);
    d2fdu2   = regexprep(d2fdu2, kStrs{i}, name, 0);
    d2fdk2   = regexprep(d2fdk2, kStrs{i}, name, 0);
    d2fdudx  = regexprep(d2fdudx, kStrs{i}, name, 0);
    d2fdxdu  = regexprep(d2fdxdu, kStrs{i}, name, 0);
    d2fdkdx  = regexprep(d2fdkdx, kStrs{i}, name, 0);
    d2fdxdk  = regexprep(d2fdxdk, kStrs{i}, name, 0);
    d3fdx3   = regexprep(d3fdx3, kStrs{i}, name, 0);
    d3fdkdx2 = regexprep(d3fdkdx2, kStrs{i}, name, 0);
    r        = regexprep(r, kStrs{i}, name, 0);
    drdx     = regexprep(drdx, kStrs{i}, name, 0);
    drdu     = regexprep(drdu, kStrs{i}, name, 0);
    drdk     = regexprep(drdk, kStrs{i}, name, 0);
    d2rdx2   = regexprep(d2rdx2, kStrs{i}, name, 0);
    d2rdk2   = regexprep(d2rdk2, kStrs{i}, name, 0);
    d2rdkdx  = regexprep(d2rdkdx, kStrs{i}, name, 0);
    d2rdxdk  = regexprep(d2rdxdk, kStrs{i}, name, 0);
end

if verbose; fprintf('   ...done.\n'); end

%% Convert strings into function handles
if verbose; fprintf('Converting expressions into handles...'); end

% Clear unnecessary variables from scope
clear SymModel vSyms xuSyms xSyms uSyms kSyms ...
    xuStrs xStrs uStrs kStrs xNamesFull uNamesFull vxNames vuNames ...
    xNames uNames ...
    C1Entries C1Values C2Entries C2Values cEntries cValues ...
    nC1Entries nC2Entries ncEntries defaultOpts opts ...
    addlength currentLength i iExpr ind iy match nAdd nExpr uAddInd name ...

u = @(t,q)repmat(u, 1,numel(t));

f = eval(['@(t,x,u,k) [' f ']']);

if ~isempty(dfdx)
    dfdx = regexprep(dfdx, '[', ''); %remove extra "[" from front of lines
    dfdx = regexprep(dfdx, ']', ''); %and "]"
    dfdx = strtrim(dfdx);            %trim excess new lines from the end
    dfdx = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' dfdx '])))']);
end

if ~isempty(dfdu)
    dfdu = regexprep(dfdu, '[', '');
    dfdu = regexprep(dfdu, ']', '');
    dfdu = strtrim(dfdu);
    dfdu = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' dfdu '])))']);
end

if ~isempty(dfdk)
    dfdk = regexprep(dfdk, '[', '');
    dfdk = regexprep(dfdk, ']', '');
    dfdk = strtrim(dfdk);
    dfdk = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' dfdk '])))']);
end

if ~isempty(d2fdx2)
    d2fdx2 = regexprep(d2fdx2, '[', '');
    d2fdx2 = regexprep(d2fdx2, ']', '');
    d2fdx2 = strtrim(d2fdx2);
    d2fdx2 = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2fdx2 '])))']);
end

if ~isempty(d2fdu2)
    d2fdu2 = regexprep(d2fdu2, '[', '');
    d2fdu2 = regexprep(d2fdu2, ']', '');
    d2fdu2 = strtrim(d2fdu2);
    d2fdu2 = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2fdu2 '])))']);
end

if ~isempty(d2fdk2)
    d2fdk2 = regexprep(d2fdk2, '[', '');
    d2fdk2 = regexprep(d2fdk2, ']', '');
    d2fdk2 = strtrim(d2fdk2);
    d2fdk2 = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2fdk2 '])))']);
end

if ~isempty(d2fdudx)
    d2fdudx = regexprep(d2fdudx, '[', '');
    d2fdudx = regexprep(d2fdudx, ']', '');
    d2fdudx = strtrim(d2fdudx);
    d2fdudx = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2fdudx '])))']);
end

if ~isempty(d2fdxdu)
    d2fdxdu = regexprep(d2fdxdu, '[', '');
    d2fdxdu = regexprep(d2fdxdu, ']', '');
    d2fdxdu = strtrim(d2fdxdu);
    d2fdxdu = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2fdxdu '])))']);
end

if ~isempty(d2fdkdx)
    d2fdkdx = regexprep(d2fdkdx, '[', '');
    d2fdkdx = regexprep(d2fdkdx, ']', '');
    d2fdkdx = strtrim(d2fdkdx);
    d2fdkdx = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2fdkdx '])))']);
end

if ~isempty(d2fdxdk)
    d2fdxdk = regexprep(d2fdxdk, '[', '');
    d2fdxdk = regexprep(d2fdxdk, ']', '');
    d2fdxdk = strtrim(d2fdxdk);
    d2fdxdk = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2fdxdk '])))']);
end

if ~isempty(d3fdx3)
    d3fdx3 = regexprep(d3fdx3, '[', '');
    d3fdx3 = regexprep(d3fdx3, ']', '');
    d3fdx3 = strtrim(d3fdx3);
    d3fdx3 = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d3fdx3 '])))']);
end

if ~isempty(d3fdkdx2)
    d3fdkdx2 = regexprep(d3fdkdx2, '[', '');
    d3fdkdx2 = regexprep(d3fdkdx2, ']', '');
    d3fdkdx2 = strtrim(d3fdkdx2);
    d3fdkdx2 = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d3fdkdx2 '])))']);
end

r = eval(['@(t,x,u,k) [' r ']']);

if ~isempty(drdx)
    drdx = regexprep(drdx, '[', ''); %remove extra "[" from front of lines
    drdx = regexprep(drdx, ']', ''); %and "]"
    drdx = strtrim(drdx);            %trim excess new lines from the end
    drdx = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' drdx '])))']);
end

if ~isempty(drdu)
    drdu = regexprep(drdu, '[', '');
    drdu = regexprep(drdu, ']', '');
    drdu = strtrim(drdu);
    drdu = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' drdu '])))']);
end

if ~isempty(drdk)
    drdk = regexprep(drdk, '[', '');
    drdk = regexprep(drdk, ']', '');
    drdk = strtrim(drdk);
    drdk = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' drdk '])))']);
end

if ~isempty(d2rdx2)
    d2rdx2 = regexprep(d2rdx2, '[', '');
    d2rdx2 = regexprep(d2rdx2, ']', '');
    d2rdx2 = strtrim(d2rdx2);
    d2rdx2 = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2rdx2 '])))']);
end

if ~isempty(d2rdk2)
    d2rdk2 = regexprep(d2rdk2, '[', '');
    d2rdk2 = regexprep(d2rdk2, ']', '');
    d2rdk2 = strtrim(d2rdk2);
    d2rdk2 = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2rdk2 '])))']);
end

if ~isempty(d2rdkdx)
    d2rdkdx = regexprep(d2rdkdx, '[', '');
    d2rdkdx = regexprep(d2rdkdx, ']', '');
    d2rdkdx = strtrim(d2rdkdx);
    d2rdkdx = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2rdkdx '])))']);
end

if ~isempty(d2rdxdk)
    d2rdxdk = regexprep(d2rdxdk, '[', '');
    d2rdxdk = regexprep(d2rdxdk, ']', '');
    d2rdxdk = strtrim(d2rdxdk);
    d2rdxdk = eval(['@(t,x,u,k) inf2big(nan2zero(sparse([' d2rdxdk '])))']);
end

if verbose; fprintf('done.\n'); end

%% Replace rate constants with their values
if verbose; fprintf('Evaluating handles...'); end

m.Type = 'Model.AnalyticKronecker';
m.Name = Name;

m.Compartments = struct('Name', vNames, 'Dimension', num2cell(dv));
m.Species      = struct('Name', xuNames, 'Compartment', vxuNames, 'IsInput', num2cell(isu), 'Value', num2cell(xu0));
m.Outputs      = struct('Name', yNames, 'Expressions', yMembers, 'Values', yValues);
m.Parameters   = struct('Name', kNames, 'Value', num2cell(k));
m.Reactions    = struct('Name', rNames);

m.nv = nv;
m.nx = nx;
m.nu = nu;
m.nq = 0;
m.ny = ny;
m.nk = nk;
m.nr = nr;

m.dv = dv;
m.k  = k;
m.x0 = x0;

m.u    = @(t)u(t,zeros(0,1));
m.q    = zeros(0,1);
m.dudq = zeros(0,0);
m.nqu  = zeros(nu,1);

m.isu   = isu;
m.xInd = find(~isu);
m.uInd = find(isu);
m.vxInd = vxInd;
m.vuInd = vuInd;
%rOrder
%krInd

%As
%Bs
m.C1 = C1;
m.C2 = C2;
m.c  = c;

clear nv nx nk nu ny vNames xuNames yNames kNames rNames

m.f         = @(t,x,u)f(t,x,u,k);
if ~isempty(dfdx)
    m.dfdx      = @(t,x,u)dfdx(t,x,u,k);
end
if ~isempty(dfdk)
    m.dfdk      = @(t,x,u)dfdk(t,x,u,k);
end
if ~isempty(dfdu)
    m.dfdu      = @(t,x,u)dfdu(t,x,u,k);
end
if ~isempty(d2fdx2)
    m.d2fdx2    = @(t,x,u)d2fdx2(t,x,u,k);
end
if ~isempty(d2fdu2)
    m.d2fdu2    = @(t,x,u)d2fdu2(t,x,u,k);
end
if ~isempty(d2fdk2)
    m.d2fdk2    = @(t,x,u)d2fdk2(t,x,u,k);
end
if ~isempty(d2fdudx)
    m.d2fdudx   = @(t,x,u)d2fdudx(t,x,u,k);
end
if ~isempty(d2fdxdu)
    m.d2fdxdu   = @(t,x,u)d2fdxdu(t,x,u,k);
end
if ~isempty(d2fdkdx)
    m.d2fdkdx   = @(t,x,u)d2fdkdx(t,x,u,k);
end
if ~isempty(d2fdxdk)
    m.d2fdxdk   = @(t,x,u)d2fdxdk(t,x,u,k);
end
if ~isempty(d3fdx3)
    m.d3fdx3    = @(t,x,u)d3fdx3(t,x,u,k);
end
if ~isempty(d3fdkdx2)
    m.d3fdkdx2  = @(t,x,u)d3fdkdx2(t,x,u,k);
end

m.S = S;
m.r = @(t,x,u)r(t,x,u,k);

if ~isempty(drdx)
    m.drdx      = @(t,x,u)drdx(t,x,u,k);
end
if ~isempty(drdk)
    m.drdk      = @(t,x,u)drdk(t,x,u,k);
end
if ~isempty(drdu)
    m.drdu      = @(t,x,u)drdu(t,x,u,k);
end
if ~isempty(d2rdx2)
    m.d2rdx2    = @(t,x,u)d2rdx2(t,x,u,k);
end
if ~isempty(d2rdk2)
    m.d2rdk2    = @(t,x,u)d2rdk2(t,x,u,k);
end
if ~isempty(d2rdkdx)
    m.d2rdkdx   = @(t,x,u)d2rdkdx(t,x,u,k);
end
if ~isempty(d2rdxdk)
    m.d2rdxdk   = @(t,x,u)d2rdxdk(t,x,u,k);
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

function varargout = update(newk, newx0, newq)
        % Apply changes
        k = newk;
        x0 = newx0;
        q = newq;
        
        m.x0 = x0;
        m.k  = k;
        m.q  = q;
        
        % Update function handles
        m.u             = @(t)u(t,q);
        m.f             = @(t,x,u)f(t,x,u,k);
        if ~isempty(dfdx)
            m.dfdx      = @(t,x,u)dfdx(t,x,u,k);
        end
        if ~isempty(dfdk)
            m.dfdk      = @(t,x,u)dfdk(t,x,u,k);
        end
        if ~isempty(dfdu)
            m.dfdu      = @(t,x,u)dfdu(t,x,u,k);
        end
        if ~isempty(d2fdx2)
            m.d2fdx2    = @(t,x,u)d2fdx2(t,x,u,k);
        end
        if ~isempty(d2fdu2)
            m.d2fdu2    = @(t,x,u)d2fdu2(t,x,u,k);
        end
        if ~isempty(d2fdk2)
            m.d2fdk2    = @(t,x,u)d2fdk2(t,x,u,k);
        end
        if ~isempty(d2fdudx)
            m.d2fdudx   = @(t,x,u)d2fdudx(t,x,u,k);
        end
        if ~isempty(d2fdxdu)
            m.d2fdxdu   = @(t,x,u)d2fdxdu(t,x,u,k);
        end
        if ~isempty(d2fdkdx)
            m.d2fdkdx   = @(t,x,u)d2fdkdx(t,x,u,k);
        end
        if ~isempty(d2fdxdk)
            m.d2fdxdk   = @(t,x,u)d2fdxdk(t,x,u,k);
        end
        if ~isempty(d3fdx3)
            m.d3fdx3    = @(t,x,u)d3fdx3(t,x,u,k);
        end
        if ~isempty(d3fdkdx2)
            m.d3fdkdx2  = @(t,x,u)d3fdkdx2(t,x,u,k);
        end
        m.r             = @(t,x,u)r(t,x,u,k);
        if ~isempty(dfdx)
            m.drdx      = @(t,x,u)drdx(t,x,u,k);
        end
        if ~isempty(dfdk)
            m.drdk      = @(t,x,u)drdk(t,x,u,k);
        end
        if ~isempty(dfdu)
            m.drdu      = @(t,x,u)drdu(t,x,u,k);
        end
        if ~isempty(d2fdx2)
            m.d2rdx2    = @(t,x,u)d2rdx2(t,x,u,k);
        end
        if ~isempty(d2fdk2)
            m.d2rdk2    = @(t,x,u)d2rdk2(t,x,u,k);
        end
        if ~isempty(d2fdkdx)
            m.d2rdkdx   = @(t,x,u)d2rdkdx(t,x,u,k);
        end
        if ~isempty(d2fdxdk)
            m.d2rdxdk   = @(t,x,u)d2rdxdk(t,x,u,k);
        end
        m.Update = @update;
        
        varargout{1} = m;
    end

end