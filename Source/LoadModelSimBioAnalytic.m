function m = LoadModelSimBioAnalytic(simbio)
%LoadModelSimBioAnalytic Load analytic model from a Matlab SimBiology model
%   Modify the model and add outputs after calling this.
%
%   m = LoadModelSimBioAnalytic(SimbioModel)
%
%   Inputs
%   simbio: [ simbio model object | string ]
%       Matlab SimBiology Model object or name of SBML file ot import

% (c) 2015 Kevin Shi, David R Hagen & Bruce Tidor
% This work is released under the MIT license.

if ischar(simbio)
    simbio = sbmlimport(simbio);
end

m = simbio2analytic(simbio);
