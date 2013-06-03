function m = LoadModel(files)
%LoadModel Load a model from files using various methods
%
%   m = LoadModel(files)
%
%   Inputs
%   files: [ string | cell array of strings ]
%       The files that are to be loaded
%
%   Outputs
%   m: [ model struct scalar ]
%       The model was that loaded
%
%   Depending on the extension, the file will be interpretted as a
%   Kronecker mass action file (.txt) or an SBML file (.xml). This function
%   will first try to convert the SBML model to mass action, but if that
%   fails, it will be converted to a psuedo-Kronecker analytic model.
%
%   To force a particular type of interpretation use the corresponding
%   LoadModel* functions.

% (c) 2011 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Standardize files a cell vector of strings
if ischar(files)
    files = {files};
end

if strcmp(files{1}(end-3:end), '.txt')
    % Kronecker mass action model
    m = LoadModelMassAction(files);
elseif strcmp(files{1}(end-3:end), '.xml') || strcmp(files{1}(end-4:end), '.sbml')
    assert(numel(files) == 1, 'KroneckerBio:LoadModel:OneFileForSbml', 'Only one SBML model can be loaded at once')
    % SBML model
    try
        m = LoadModelSbmlMassAction(files{1});
    catch
        m = LoadModelSbmlAnalytic(files{1});
    end
else
    error('KroneckerBio:LoadModel:UnknownModelType', 'Unknown file extension while trying to load a model')
end