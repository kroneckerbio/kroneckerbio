function m = AddInput(m, name, compartment, func, parameters)
%AddInput Add an input species to a KroneckerBio model
%
%   m = AddInput(m, name, compartment, value, parameters)
%
%   Inputs
%   m: [ model struct scalar ]
%       The model to which the input species will be added
%   name: [ string ]
%       A name for the input. This is the name by which reactions will
%       refer to it.
%   compartment: [ string ]
%       The name of the compartment to whcih it will be added
%   func: [ nonnegative scalar | string | 
%            function handle @(t,q) returns nonnegative scalar {0} ]
%       A scalar indicates that the input has a constant concentration, a
%       string is turned into a function handle of t and q, and the
%       function handle returns the value of the function at all time.
%   parameters: [ vector nq {zeros(0,1)} ]
%       This vector of parameters is passed to value(t,q) when determining
%       the value of the input at a particular time t
%
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new input added.

% (c) 2013 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean-up inputs
if nargin < 5
    parameters = [];
    if nargin < 4
        func = [];
    end
end

% Increment counter
nu = m.add.nu + 1;
m.add.nu = nu;
m.add.Inputs = growInputs(m.add.Inputs, nu);

% Add item
m.add.Inputs(nu).Name = fixSpeciesName(name);
m.add.Inputs(nu).Compartment = fixCompartmentName(compartment);
m.add.Inputs(nu).Function = fixInputValueFunction(func);
m.add.Inputs(nu).Parameters = fixInputParameters(parameters);

m.Ready = false;
