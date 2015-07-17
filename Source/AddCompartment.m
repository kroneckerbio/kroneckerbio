function m = AddCompartment(m, varargin)
%AddCompartment Add a compartment to a generic model
%
%   m = AddCompartment(m, ...)
%
%   This is a generic function. Select the appropriate help file below
%   depending on the type of the model.
%
%   Model.MassActionAmount
%       help addCompartmentMassActionAmount
%
%   Model.MassActionConcentration
%       help addCompartmentMassActionConcentration
%
%   Model.Analytic
%       help addCompartmentAnalytic

% (c) 2015 David R Hagen
% This work is released under the MIT license.

if is(m, 'Model.MassActionAmount')
    m = addCompartmentMassActionAmount(m, varargin{:});
elseif is(m, 'Model.MassActionConcentration')
    m = addCompartmentMassActionConcentration(m, varargin{:});
elseif is(m, 'Model.Analytic')
    m = addCompartmentAnalytic(m, varargin{:});
else
    error('KroneckerBio:AddCompartment:m', 'm must be a model')
end
