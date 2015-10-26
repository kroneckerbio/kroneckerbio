function AddState(this, name, compartment, ic)
%AddState Add a state species to a massaction model
%
%   Inputs
%   this: [ scalar MassActionAmountModel ]
%       The model to which the state will be added
%   name: [ string ]
%       A name for the state. This is the name by which reactions will
%       refer to it.
%   compartment: [ string ]
%       The name of the compartment to which it will be added
%   ic: [ string | nonnegative scalar | cell vector of strings
%           | cell matrix {0} ]
%       This is the name of the seed that will control the initial value of
%       this state. If it is an empty string, then no seed will control it,
%       and it will always have an initial condition of 0. If a number is
%       given, this will be the initial condition of the state, which will
%       not have an associated parameter.

% Clean-up inputs
if nargin < 4
    ic = [];
end

% Increment counter
nx = this.nx + 1;
this.nx = nx;
this.growStates;

% Add item
this.States(nx).Name         = FieldValidator.SpeciesName(name);
this.States(nx).Compartment  = FieldValidator.CompartmentName(compartment);
this.States(nx).InitialValue = FieldValidator.StateInitialConditionMassAction(ic);

this.Ready = false;
