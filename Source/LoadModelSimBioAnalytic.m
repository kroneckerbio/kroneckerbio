function m = LoadModelSimBioAnalytic(simbio)
%LoadModelSimBioAnalytic Import Matlab Simbiology model and covert to
%   Kronecker analytic model. The model may be further modified by the
%   user.
%
%   m = LoadModelSimBioAnalytic(simbio)
%
%   Inputs
%   simbio: [ simbio model object | string ]
%       Matlab SimBiology Model object or name of SBML file to import
%
%   Outputs
%   m: [ Model.Analytic struct ]
%       An analytic kroneckerbio model
%
%   Notes
%   - All species are tracked in amount.
%   - Species that are set to constant or boundaryCondition are
%   converted to inputs
%   - Initial assignment rules of species are copied to the initial
%   condition expression.
%   - Repeated assignment rules of compartments are copied to the
%   compartment size.
%   - Repeated assignment rules of parameters and species get converted
%   into Kronecker rules.
%   - Rate rules are converted to reactions with a single product.
%
%   Limitations
%   - Not all Simbiology features are compatible with this converter. This
%   function ignores events, algebraic rules, and all functions in the
%   model.
%   - Simbiology does not expose the dimension of a compartment. All
%   compartments are given a dimension of 3. This does not affect the
%   behavior of the model.

% (c) 2015 Kevin Shi, David R Hagen & Bruce Tidor
% This work is released under the MIT license.

if ischar(simbio)
    simbio = sbmlimport(simbio);
end

m = simbio2analytic(simbio);
