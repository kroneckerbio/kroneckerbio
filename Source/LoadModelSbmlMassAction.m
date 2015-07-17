function m = LoadModelSbmlMassAction(filename, opts)
%LoadModelSbmlMassAction Load a mass action model from an SBML file.Modify the model and add outputs after calling this.
%
%   m = LoadModelSbmlMassAction(SimbioModel, opts)
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
%   m: [ Model.MassActionAmount struct ]
%       A mass action kroneckerbio model with all values in amounts (not
%       concentrations)
%
%   Only a subset of SBML models can be converted to mass action kroneckerbio
%   models. This function will throw an error if invalid rate forms and other
%   features are detected.

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

%% Call libSBML to import SBML model
if verbose; fprintf('Convert SBML model using libSBML...'); end

sbml = TranslateSBML(filename, double(opts.Validate), opts.Verbose);

if verbose; fprintf('done.\n'); end

%% Convert model
symbolic = sbml2symbolic(sbml, opts);

assert(isValidSymbolicModel(symbolic), 'LoadModelSbmlAnalytic:InvalidSymbolicModel', 'Symbolic model intermediate failed validation check')

m = symbolic2massaction(symbolic, opts);

end
