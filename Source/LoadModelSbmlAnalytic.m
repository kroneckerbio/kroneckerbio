function m = LoadModelSbmlAnalytic(filename, opts)
%LoadModelSbmlAnalytic Import SBML model and covert to kroneckerbio analytic
%model. Modify the model and add outputs after calling this.
%
%   m = LoadModelSbmlAnalytic(model, yNames, yMembers, yValues, opts)
%
%   Inputs
%   filename: [ string ]
%       Path to SBML file. Usually has .sbml or .xml extension
%   opts: [ options struct scalar {} ]
%       .Verbose [ logical scalar {false} ]
%       	Print progress to command window
%       .Validate [ logical scalar {false} ]
%           Whether to use libSBML's model validation tool
%       .UseNames [logical scalar {false} ]
%           Whether to convert SBML IDs to Names and autogenerate new IDs
%           Use this when the supplied SBML model uses "nice" names as IDs
%
%   Outputs
%   m: [ Model.Analytic struct ]
%       An analytic kroneckerbio model
%
%   Notes: 
%   - Inputs are any species that have "constant" or "boundaryCondition" set.
%   - Seeds are generated from any parameter that appears in an InitialAssignment.
%
%   Limitations:
%   Not all Simbiology features are compatible with this converter. This
%   function ignores events, some rules, and all functions in the model.

%% Clean up inputs
if nargin < 2
    opts = [];
end

% Default options
opts_.Verbose = 0;
opts_.Validate = false;
opts_.UseNames = false;

opts = mergestruct(opts_, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

%% Call libSBML to import SBML model
if verbose; fprintf('Convert SBML model using libSBML...'); end

sbml = TranslateSBML(filename, double(opts.Validate), opts.Verbose);

if verbose; fprintf('done.\n'); end

%% Convert model
symbolic = sbml2symbolic(sbml, opts);

assert(isValidSymbolicModel(symbolic), 'LoadModelSbmlAnalytic:InvalidSymbolicModel', 'Symbolic model intermediate failed validation check')
    
m = symbolic2analytic(symbolic, opts);

end
