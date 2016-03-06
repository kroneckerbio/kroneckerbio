function m = LoadModelSbmlAnalytic(sbml, opts)
%LoadModelSbmlAnalytic Import SBML model and covert to Kronecker analytic
%   model. The model may be further modified by the user.
%
%   m = LoadModelSbmlAnalytic(sbml, opts)
%
%   Inputs
%   sbml: [ libSBML model object | string ]
%       LibSBML Model object or name of SBML file to import
%   opts: [ options struct scalar {} ]
%       .Verbose [ logical scalar {false} ]
%       	Print progress to command window
%       .Validate [ logical scalar {false} ]
%           Whether to use libSBML's model validation tool
%       .ICsAsSeeds [ {true} | false ]
%           Whether to make all state initial conditions seeds (that aren't
%           substituted rults) or hardcode initial conditions.
%
%   Outputs
%   m: [ Model.Analytic struct ]
%       An analytic kroneckerbio model
%
%   Notes
%   - All species are tracked in amount. Concentration species are
%   converted with best effort.
%   - Species that are set to constant or boundaryCondition are
%   converted to inputs
%   - Initial assignment rules of species are copied to the initial
%   condition expression.
%   - Repeated assignment rules of compartments are copied to the
%   compartment size.
%   - Repeated assignment rules of parameters and species get converted
%   into Kronecker rules.
%   - Rate rules are converted to reactions with a single product.
%   - The ID of a component is used if the component has no name.
%   - opts.Verbose > 0 will enable interactive validation when using
%   TranslateSBML. Disable interactive behavior in scripts by decreasing
%   the verbosity.
%
%   Limitations
%   - Not all Simbiology features are compatible with this converter. This
%   function ignores events, algebraic rules, and all functions in the
%   model.

% (c) 2015 Kevin Shi and David R Hagen
% This work is released under the MIT license.

%% Clean up inputs
if nargin < 2
    opts = [];
end

% Default options
opts_.Verbose = 0;
opts_.Validate = false;
opts_.ICsAsSeeds = false;

opts = mergestruct(opts_, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

%% Call libSBML to import SBML model
if verbose; fprintf('Load SBML model using libSBML...'); end

if ischar(sbml)
    sbml = TranslateSBML(sbml, double(opts.Validate), opts.Verbose);
end

if verbose; fprintf('done.\n'); end

m = sbml2analytic(sbml, opts);

end
