function m = AddSpecies(m, name, compartment, value, units, isInput, parameters)
%AddSpecies Add a species to a KroneckerBio model
%
%   m = AddSpecies(m, name, compartment, value, units, isInput, parameters)
%
%   Each species exists inside a single compartment. Transport from one
%   compartment to another requires a reaction. A species name only needs
%   to be unique within a compartment; species may have the same name if
%   they reside in different compartments.
%
%   Inputs
%   m: [ model struct scalar ]
%       The model to which the species will be added
%   name: [ string ]
%       A name for the species. This is the name by which reactions will
%       refer to it.
%   compartment: [ string ]
%       The name of the compartment it will be added
%   value: [ nonnegative scalar | string | 
%            function handle @(t,q) returns nonnegative scalar {0} ]
%       Optional
%       If the species is a state, this can only be a number and it
%       specifies the initial amount of the state. If the species is an
%       input, a number indicates its constant concentration, a string is
%       turned into a function handle of t and q, and the function handle
%       returns the value of the function at all time.
%   units: [ '' ]
%       This is currently not implemented. The canonical representaion of
%       species values is in amount.
%   isInput: [ logical scalar {false} ]
%       Specifies if the species is an input or a state.
%   parameters: [ vector nq ]
%       This vector of parameters is passed to value(t,q) when determining
%       the value of the input at a particular time t.
%
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new species added.

% (c) 2011 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean-up inputs
if nargin < 7
    parameters = [];
    if nargin < 6
        isInput = [];
        if nargin < 5
            units = [];
            if nargin < 4
                value = [];
            end
        end
    end
end

% Default inputs
if isempty(isInput)
    isInput = false;
end

% Simply defer calculation to the approriate functions
if isInput
    m = AddInput(m, name, compartment, value, parameters, units);
else
    m = AddState(m, name, compartment, value, units);
end
