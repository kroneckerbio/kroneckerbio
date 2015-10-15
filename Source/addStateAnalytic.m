function m = addStateAnalytic(m, name, compartment, ic)
%AddState Add a state species to a KroneckerBio model. The state's initial
%condition can be a scalar value or an arbitrary expression of seeds.
%
%   m = AddState(m, name, compartment, ic)
%
%   Inputs
%   m: [ model struct scalar ]
%       The model to which the state will be added
%   name: [ string ]
%       A name for the state. This is the name by which reactions will
%       refer to it.
%   compartment: [ string ]
%       The name of the compartment to which it will be added
%   ic: [ string | nonnegative scalar {0} ]
%       String expression for the initial condition in terms of seeds or
%       nonnegative scalar initial amount.
%
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new state added.
%

% Clean up inputs
if nargin < 4
    ic = [];
end

% Increment counter
nx = m.nx + 1;
m.nx = nx;
m.States = growStates(m.States, nx);

% Add item
m.States(nx).Name         = fixSpeciesName(name);
m.States(nx).Compartment  = fixCompartmentName(compartment);
m.States(nx).InitialValue = fixStateICAnalytic(ic);

m.Ready = false;
