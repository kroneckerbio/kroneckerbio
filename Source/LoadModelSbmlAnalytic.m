function kronModel = LoadModelSbmlAnalytic(simbioModel, yNames, yMembers, yValues, opts)
%LoadModelSbmlAnalytic Convert an SBML model into a pseudo-kronecker model
%   which interacts with much of the Kronecker toolbox much like a
%   Kronecker model.
% 
%   m = LoadModelSbmlAnalytic(SimbioModel, yNames, yMembers, yValues, opts)
%
%   Inputs
%   simbioModel: [ path string | SimBiology Model scalar | 
%                  symbolic model scalar ]
%       This can be a path to an SBML file, a Simbiology Model, or a
%       symbolic model.
%   yNames: [ string | cell vector of strings {{m.Species.Name}} ]
%       SBML has no concept of outputs; this is used to add them to the
%       KroneckerBio model. The names of the outputs. If yMembers is not
%       provided, the names are also interpreted as the species that will
%       be represented by the outputs.
%   yMembers: [ cell vector ny of cell vectors of strings | {} ]
%       The names of the species to be included in each output
%   yValues: [ cell vector of nonegative vectors ]
%       Each numeric entry in this vector is associated with one of the
%       expressions in yMembers. This value tells how much a species will
%       contribute when it matches the corresponding expressions. If a
%       species is matched by multiple expressions, the last expression to
%       match overrides all others.
%   opts: [ options struct scalar {} ]
%       .Verbose [ logical scalar {false} ]
%       	Print progress to command window
%       .Order [ 0 | 1 | {2} | 3 ]
%       	Determines how deep the derivatives should be taken. Each level
%       	increases the cost exponentially, but increases the number of
%       	Kronecker functions that can be run on the model.
%
%   Outputs
%   m: [ Kronecker model struct scalar ]
%       An analytic pseudo-Kronecker model
%
%   Inputs are any species that have "constant" or "boundaryCondition" set.
%
%   Seeds are generated from any parameter that appears in an
%   InitialAssignment. The expression must be a linear combination of
%   parameters.
%
%   Limitations:
%   Not all Simbiology features are compatible with this converter. This
%   function ignores any events, most rules, and functions of the model.

% (c) 2013 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean up inputs
if nargin < 5
    opts = [];
    if nargin < 4
        yValues = [];
        if nargin < 3
            yMembers = [];
            if nargin < 2
                yNames = [];
            end
        end
    end
end

% Load sbml model directly to symbolic model
symModelDirect = sbml2Symbolic(simbioModel);

% Load sbml if file is provided
if ischar(simbioModel)
    simbioModel = sbmlimport(simbioModel);
end

% Use sbml2Symbolic to convert an SBML model to a symbolic model
symModel = simbio2Symbolic(simbioModel, opts);

%%%%%%%% Add specified outputs to the symbolic model %%%%%%%%%

%%% First, standardize yNames, yMembers, and yValues to their fully expressed
%%% forms, using indexes to refer to states and inputs.

% If only yNames were provided, then each yName is a single species to be
% added as an output. Copy the yNames over to the yMembers, encapsulating
% them in cells.
if ~isempty(yNames) && isempty(yMembers)
    yMembers = num2cell(yNames);
end

ny = length(yNames);
nmemberspery = cellfun(@length,yMembers);

% If no yValues were provided, default the yValues to 1 for each yMember.
if ~isempty(yNames) && isempty(yValues)
    yValues = arrayfun(@(len)ones(len,1),nmemberspery,'UniformOutput',false);
end

% Check output specifiers
assert(iscell(yMembers) && all(cellfun(@iscell,yMembers)), 'yMembers should be a cell vector of cell vectors of strings')

% Set up output function
% Ensure that symModel has y and yNames fields
if isfield(symModel,'y')
    oldny = length(symModel.y);
else
    symModel.y = sym([]);
    symModel.yNames = {};
    oldny = 0;
end
symModel.y = [symModel.y; sym(zeros(ny,1))];
xNames = symModel.xNames;
uNames = symModel.uNames;
xSyms = symModel.xSyms;
uSyms = symModel.uSyms;
nx = symModel.nx;
nu = symModel.nu;
C1 = zeros(ny,nx);
C2 = zeros(ny,nu);
c = zeros(ny,1);

for yi = 1:ny
    
    thisyMembers = yMembers{yi};
    thisnyMembers = length(thisyMembers);
    
    for ymi = 1:thisnyMembers

        thisyMember = thisyMembers{ymi};
        thisyValue = yValues{yi}(ymi);
        
        % Record constant value if member string is empty
        if isempty(thisyMember)
            c(oldny+yi) = thisyValue;
        % Otherwise...
        else
            % Look for name among states
            thisyMemberindex = find(strcmp(xNames,thisyMember));
            % If not found there...
            if isempty(thisyMemberindex)
                % Look for name among inputs
                thisyMemberindex = find(strcmp(uNames,thisyMember));
                % If not found there...
                if isempty(thisyMemberindex)
                    % Name is invalid
                    error(['Unrecognized state or input ' thisyMember])
                % If found among inputs...
                else
                    % Record yValue in C2
                    C2(oldny+yi,thisyMemberindex) = thisyValue;
                end
            % If found among states...
            else
                % Record yValue in C1
                C1(oldny+yi,thisyMemberindex) = thisyValue;
            end
        end

    end
        
    symModel.y(oldny+yi) = C1(oldny+yi,:)*xSyms + C2(oldny+yi,:)*uSyms + c(oldny+yi);
        
end

symModel.yNames = [symModel.yNames; yNames];


% Use symbolic2Kronecker to convert a symbolic model to a psuedo kronecker model
kronModel = symbolic2PseudoKronecker(symModel, opts);
