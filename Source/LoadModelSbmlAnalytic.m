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

% Load sbml if file is provided
if ischar(simbioModel)
    simbioModel = sbmlimport(simbioModel);
end

% Use sbml2Symbolic to convert an SBML model to a symbolic model
symModel = simbio2Symbolic(simbioModel, opts);

% Use symbolic2Kronecker to convert a symbolic model to a psuedo kronecker model
kronModel = symbolic2PseudoKronecker(symModel, yNames, yMembers, yValues, opts);
