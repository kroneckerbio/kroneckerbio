function m = LoadModelSbmlMassAction(SimbioModel, opts)
%LoadModelSbmlMassAction Load a mass action model from an SBML file, a
%   SimBiology model, or a symbolic model
%
%   m = LoadModelSbmlMassAction(SimbioModel, opts)
%
%   Inputs
%   SimbioModel: [ path string | SimBiology Model scalar | 
%                  symbolic model scalar]
%       This can be a path to an SBML file, a Simbiology Model, or a
%       symbolic model.
%   opts: [ options struct ]
%       Optional
%       .Verbose [ nonnegative integer scalar {1} ]
%           Bigger number displays more progress information
%
%   Outputs
%   m: [ Kronecker model struct scalar ]
%       A complete mass action Kronecker model
%
%   Only a subset of SBML models can be converted to mass action Kronecker
%   models. This function will not return an analytic psuedo-Kronecker
%   model if the conversion to mass action fails; it will simply error out.

% (c) 2011 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Load sbml if file is provided
if ischar(SimbioModel)
    SimbioModel = sbmlimport(SimbioModel);
end

% Use simbio2Symbolic to convert an SBML model to a symbolic model
if strcmp(SimbioModel.Type, 'sbiomodel')
    SimbioModel = simbio2Symbolic(SimbioModel, opts);
end

% Use symbolic2MassAction to convert a symbolic model to a Kronecker mass action model
m = symbolic2MassAction(SimbioModel,opts);